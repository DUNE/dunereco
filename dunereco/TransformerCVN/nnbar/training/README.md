<!-- Copied from https://github.com/KaiwenYu2001/dune-nnbar-transformercvn_v2 at commit deb1014
     (author: Kaiwen Yu). Kept verbatim below this note so that the directory can be synchronised with
     the upstream repository. Not included from upstream: nnbar_best.ckpt (70 MB; fetched by
     ../setup/setup_nnbar_cvn.sh) and the files the author keeps local (pyproject.toml, docs/, tests/,
     THIRD_PARTY_NOTICES.md), so the "pip install -e" and pytest instructions below do not apply here:
     use ../setup/setup_nnbar_cvn.sh and PYTHONPATH (see ../README.md). Upstream has not chosen a
     license yet (see "Attribution and license status"); redistribution under the dunereco license
     is to be confirmed with the author. -->

# TransformerCVN for NNbar classification

Research code for joint event and prong classification using pixel-map CNNs and a
transformer encoder. The default training path uses a dense CNN; sparse CNN and
SDXL-derived encoder variants are also included.

## Setup

The baseline environment is Python 3.10, PyTorch 2.0.1, Lightning 1.9.5,
TorchMetrics 0.11.4, and NumPy 1.26.4. Lightning 2 is intentionally excluded:
the trainers still use the Lightning 1.x epoch hooks.

```bash
python3.10 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev,evaluation]'
```

For CUDA training, install a PyTorch build matching your machine before installing
this project. The dense path does not require MinkowskiEngine. Sparse variants
require a separately built MinkowskiEngine installation compatible with PyTorch
and CUDA; they have not been validated as part of this cleanup. `.[sdxl]` and
`.[legacy]` supply optional dependencies for experimental modules; those version
ranges are provisional and are not a reproducibility lockfile.

## Data and configuration

Training data and trained weights are not distributed. Supply preprocessed HDF5
files matching [the input schema](docs/data.md). `option_files/example.json` is the
starting configuration. Adjust its paths and hyperparameters for your data and
hardware; it is not a benchmark result.

```bash
python train.py --options_file option_files/example.json \
  --training_file data/training.h5 --gpus 1 --name nnbar_dense --log_dir outputs
```

Set `CUDA_VISIBLE_DEVICES` in your shell to select physical GPUs. `--gpus 0` selects
CPU execution. `--sparse` and `--sdxl` select alternative encoders. Resume with
`--checkpoint path/to/checkpoint.ckpt`. `--eval` is the validation interval in
training batches (default 500); it must fit the number of batches per epoch.
Gradient accumulation remains 8 batches, as in the original experiments.

Runs save TensorBoard logs, checkpoints and the resolved `options.json` under
`outputs/<name>/version_<n>/`. Checkpoints monitor `val_epoch_nnbar_ovr_tpr`:
one-vs-rest TPR for the last class at the first ROC point whose FPR is at least
0.00054 (0.054%). This preserves the original convention; the actual selected
FPR may exceed the target. The last class must represent NNbar for this metric
to have the intended meaning.

## Evaluation

Choose a checkpoint explicitly and reuse its saved options:

```bash
python evaluate.py --options outputs/nnbar_dense/version_0/options.json \
  --checkpoint outputs/nnbar_dense/version_0/checkpoints/last.ckpt \
  --training-file data/training.h5 --split validation \
  --output outputs/validation_predictions.h5 --device cpu
```

For a separate test sample, use `--split testing --testing-file data/testing.h5`.
The current trainer reconstructs its training dataset and normalization statistics
on initialization, so the original training data remains required. Only load
checkpoints you trust: Lightning checkpoints use Python pickle serialization.

The export streams event probabilities, predictions and **true labels** to HDF5.
Valid prongs are flattened, with `prong_event_index` mapping them back to event
rows; this handles batches containing different maximum numbers of prongs. Output
files must be new. `Evaluate.ipynb` plots these exports without repeating model
loading and inference for each sample.

## Repository layout

- `train.py`, `evaluate.py`: command-line entry points.
- `transformercvn/options.py`: shared defaults and JSON overrides.
- `transformercvn/dataset/`: HDF5 readers and collation.
- `transformercvn/network/{layers,networks,trainers}/`: model building blocks and training logic.
- `option_files/`: portable example configuration.
- `tests/`: CPU regression and smoke tests.

Historical trainers and experimental architectures are retained to preserve import
paths and research work. Only the full dense path is covered by the baseline
checks. See [the data contract](docs/data.md) for input requirements and historical
split behavior. The learning-rate schedule counts batches without accounting for
the entry point's gradient accumulation factor of 8; this historical behavior is
retained. The full network also retains its use of the event position parameter
for prong positions to preserve checkpoint behavior. These details require
validation when reproducing experiments.

## Development

```bash
ruff check .
ruff format --check .
python -m pytest -q
```

GitHub Actions runs these checks on Python 3.10. Tests use synthetic inputs; passing
them does not establish training convergence or reproduce physics results.

## Attribution and license status

An overall project license has not yet been selected. Some files contain upstream
copyright or attribution notices; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
The repository owner must resolve the missing upstream license material and choose
an appropriate license before release. No permission from upstream authors is
implied by this code cleanup.
