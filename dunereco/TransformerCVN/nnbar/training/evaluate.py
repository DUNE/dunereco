"""Export event and prong predictions from a trusted Lightning checkpoint."""

from argparse import ArgumentParser
from pathlib import Path

import h5py
import torch

from transformercvn.options import Options


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
    model = Network(options)
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    model.load_state_dict(checkpoint["state_dict"])
    dataset = model.testing_dataset if args.split == "testing" else model.validation_dataset
    loader_options = dict(model.dataloader_options)
    loader_options.update(batch_size=args.batch_size, num_workers=0, drop_last=False)
    dataloader = model.dataloader(dataset, **loader_options)
    path = export_predictions(model, dataloader, args.output, args.device)
    print(f"Saved predictions to {path}")


if __name__ == "__main__":
    main()
