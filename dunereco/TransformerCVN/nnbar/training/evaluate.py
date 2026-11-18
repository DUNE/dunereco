"""Export event and prong predictions from a trusted Lightning checkpoint."""

from argparse import ArgumentParser
from pathlib import Path

import h5py
import torch

from transformercvn.options import Options

# Sample tags of the integer ``model`` column written by the nnbar preprocessing
# (nnbar-production cvn/preprocess_atmnu.py): the column is an index into this list.
DEFAULT_MODEL_NAMES = ["atm", "ha_lfg", "hn_lfg", "ha_esf", "hn_esf", "ha_br", "hn_br"]


def export_model_column(output_file, input_file, split, options):
    """Copy the per-event ``model`` code of the exported rows and its name list into the predictions file.

    The exported rows follow the input file order: the whole file for the testing split, the validation
    range ``[dataset_limit * train_validation_split, dataset_limit)`` of the training file otherwise.
    """
    with h5py.File(input_file, "r") as source:
        if "model" not in source:
            return None
        num_events = len(source["model"])
        if split == "testing":
            lower, upper = 0, num_events
        else:
            lower = int(round(options.dataset_limit * options.train_validation_split * num_events))
            upper = int(round(options.dataset_limit * num_events))
        model = source["model"][lower:upper]
        names = source.attrs.get("model_names", None)
    names = [n.decode() if isinstance(n, bytes) else str(n) for n in names] if names is not None else DEFAULT_MODEL_NAMES
    with h5py.File(output_file, "a") as output:
        if output.attrs["num_events"] != len(model):
            raise ValueError(f"model column has {len(model)} rows, predictions have {output.attrs['num_events']}")
        output.create_dataset("event_model", data=model)
        output.attrs["model_names"] = names
    return model, names


def export_predictions(model, dataloader, output_file, device="cpu"):
    """Stream predictions to HDF5, flattening valid prongs to support variable batches.

    ``prong_event_index`` maps each flattened prong to an exported event row.
    Existing output files are never overwritten.
    """
    model = model.to(device).eval()
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    event_offset = 0
    with h5py.File(output_file, "x") as output, torch.inference_mode():
        output.attrs["format_version"] = 1
        output.attrs["prong_layout"] = "flat; prong_event_index refers to event rows"
        for batch in dataloader:
            batch = [value.to(device) for value in batch]
            event_targets, prong_targets, event_logits, prong_logits = model.shared_step(batch)
            mask = prong_targets >= 0
            event_probs = event_logits.softmax(-1)
            prong_probs = prong_logits.softmax(-1)[mask]
            event_index = (
                torch.arange(len(event_targets), device=prong_targets.device)
                .unsqueeze(1)
                .expand_as(prong_targets)[mask]
                + event_offset
            )
            arrays = {
                "event_probabilities": event_probs,
                "event_predictions": event_probs.argmax(-1),
                "event_targets": event_targets,
                "prong_probabilities": prong_probs,
                "prong_predictions": prong_probs.argmax(-1),
                "prong_targets": prong_targets[mask],
                "prong_event_index": event_index,
            }
            for name, tensor in arrays.items():
                values = tensor.cpu().numpy()
                if name not in output:
                    output.create_dataset(
                        name, data=values, maxshape=(None, *values.shape[1:]), chunks=True
                    )
                else:
                    dataset = output[name]
                    start = len(dataset)
                    dataset.resize(start + len(values), axis=0)
                    dataset[start:] = values
            event_offset += len(event_targets)
        output.attrs["num_events"] = event_offset
    return output_file


def build_parser():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--options", required=True, help="Saved training options.json")
    parser.add_argument("--checkpoint", required=True, help="Checkpoint to evaluate")
    parser.add_argument("--training-file", help="Override the original training data path")
    parser.add_argument("--testing-file", help="Override testing_file for the test split")
    parser.add_argument("--split", choices=("validation", "testing"), default="validation")
    parser.add_argument("--network", choices=("dense", "sparse", "sdxl"), default="dense")
    parser.add_argument("--device", default="cpu", help="cpu or cuda:N")
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--output", required=True, help="New HDF5 prediction file")
    return parser


def main(args=None):
    args = build_parser().parse_args(args)
    if args.batch_size <= 0:
        raise ValueError("--batch-size must be positive")
    if Path(args.output).exists():
        raise FileExistsError(args.output)
    if args.network == "sparse":
        from transformercvn.network.trainers.neutrino_full_sparse_trainer import (
            NeutrinoFullSparseTrainer as Network,
        )
    elif args.network == "sdxl":
        from transformercvn.network.trainers.neutrino_full_sdxl_trainer import (
            NeutrinoFullSDXLTrainer as Network,
        )
    else:
        from transformercvn.network.trainers.neutrino_full_dense_trainer import (
            NeutrinoFullDenseTrainer as Network,
        )
    options = Options.load(args.options)
    if args.training_file:
        options.training_file = args.training_file
    if args.testing_file:
        options.testing_file = args.testing_file
    if args.split == "testing" and not options.testing_file:
        raise ValueError("The testing split requires --testing-file or options.testing_file")
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    state = checkpoint["state_dict"]

    class SizedNetwork(Network):
        """Size the event and prong classifiers from the checkpoint, not from the labels present in
        the input file (a background-only or signal-only file would otherwise build a smaller head)."""

        def create_datasets(self):
            training, validation, testing = super().create_datasets()
            training.num_event_classes = state["network.event_decoder.hidden_layer.weight"].shape[0]
            training.num_prong_classes = state["network.prong_decoder.output_layer.weight"].shape[0]
            return training, validation, testing

    model = SizedNetwork(options)
    model.load_state_dict(state)
    dataset = model.testing_dataset if args.split == "testing" else model.validation_dataset
    loader_options = dict(model.dataloader_options)
    loader_options.update(batch_size=args.batch_size, num_workers=0, drop_last=False)
    dataloader = model.dataloader(dataset, **loader_options)
    path = export_predictions(model, dataloader, args.output, args.device)
    print(f"Saved predictions to {path}")
    input_file = options.testing_file if args.split == "testing" else options.training_file
    exported = export_model_column(path, input_file, args.split, options)
    if exported is not None:
        model_column, names = exported
        counts = {names[code] if 0 <= code < len(names) else str(code): int((model_column == code).sum())
                  for code in sorted(set(model_column.tolist()))}
        print("Exported event_model (index into model_names):", counts)


if __name__ == "__main__":
    main()
