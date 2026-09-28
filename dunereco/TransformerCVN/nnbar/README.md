# TransformerCVN/nnbar -- the n-nbar classifier

The DUNE n-nbar search classifies events with a TransformerCVN trained on
n-nbar signal against atmospheric neutrinos. This directory contains the
complete chain from the RecoEnergyS art dump to per-event scores, plus the
training code, kept together as the record of the analysis. It is a standalone
tree of python and ROOT scripts: nothing here is compiled, installed or run
inside art, and the dunereco build does not touch it. Use it from a checkout. The simulation up to the art dump
is documented in [linyan-w/nnbar-production](https://github.com/linyan-w/nnbar-production);
the network code is the one of
[KaiwenYu2001/dune-nnbar-transformercvn_v2](https://github.com/KaiwenYu2001/dune-nnbar-transformercvn_v2)
(copied into `training/`, see there for the upstream commit).

## Layout: one directory per step

| Directory | Stage | What it does |
|---|---|---|
| `setup/` | -- | `setup_nnbar_cvn.sh` (ROOT, python venv, checkpoint), pinned `requirements*.txt` |
| `pixelmap/` | 5 | `make_text_file_to_root_trks_shws.C`: RecoEnergyS tree `recoEnergy/WC` -> `pixelmap` tree, one row per pixel of the event image and of each prong image (3 planes x 350 x 350 around the Pandora vertex); drops events with fewer than 100 hits or no vertex |
| `preprocess/` | 6 | `preprocess.py`: pixelmap files -> one HDF5 file (sparse images, prong features, truth, precut variables) |
| `sparsify/` | 7 | `sparsify.py`: stage-6 HDF5 -> the network input format (script version of the author's `CreateFullySparseDataset.ipynb`); applies the analysis precut; carries `event_id`, `file_name` and the `genie/` truth group through. `sparsify_streaming.py`: same output in blocks, for files larger than the memory |
| `training/` | -- | the network (`transformercvn/`), `train.py`, `evaluate.py`, `option_files/example.json` (the configuration of the published checkpoint) |
| `inference/` | 8 | `run_nnbar_inference.sh` (chains stages 5-8 or 7-8), `attach_truth.py` (copies `event_id` and the truth into `predictions.h5`), `summarize_scores.py`, `find_cut.py` (score cut at a target total background efficiency, and the signal efficiencies at it) |
| `plots/` | -- | `make_score_table.py` + `plot_nnbar_scores.C` (ROOT): score distributions per sample, CVN efficiency per GENIE decay mode with Clopper-Pearson intervals, scores per mode, mode key |

Stage 4, the `RecoEnergyS` art analyzer that writes `*_cvnpreprocess.root`,
is `dunereco/RecoEnergyStudies` (branch `feature/lwan_recoenergystudies`).

## Environment

```bash
source setup/setup_nnbar_cvn.sh            # inference: CPU torch, downloads the checkpoint
source setup/setup_nnbar_cvn.sh --no-eval  # stages 5-6 only (ROOT + h5 tools)
source setup/setup_nnbar_cvn.sh --train    # CUDA torch 2.0.1 (cu118) for training
```

The script finds this tree next to itself or through `$NNBAR_CVN_DIR`.
It creates the venv `$NNBAR_CVN_VENV` (default
`/exp/dune/app/users/$USER/nnbar-cvn-venv`, gpvm home areas are small) on
first use, puts `training/` on `PYTHONPATH`, and fetches the checkpoint into
`$NNBAR_CVN_WEIGHTS` (default next to the venv) from the author's repository
at the pinned commit. The pins (torch 2.0.1, Lightning 1.9.5, torchmetrics
0.11.4, numpy 1.26.4, rich 13.3.5, transformers 4.33.3) are the ones the
network was trained with; newer rich and transformers releases break
Lightning 1.9.5 and torch 2.0.1. ROOT comes from the current environment, or
6.28.12 from the larsoft spack area on cvmfs.

## Inference

```bash
# 4. art dump (dunereco with the RecoEnergyStudies package)
lar -c recoenergys.fcl -s reco2.root                 # -> <reco2>_cvnpreprocess.root

# 5-8. host environment (not inside a container)
source setup/setup_nnbar_cvn.sh
inference/run_nnbar_inference.sh -i /path/to/cvnpreprocess/files -o /exp/dune/data/users/$USER/nnbar_eval \
    -t nnbar -j 8 -m ha_br -d cuda:0
inference/run_nnbar_inference.sh -6 stage6.h5 -o OUT -d cuda:0   # start from an existing stage-6 file
inference/summarize_scores.py OUT/predictions.h5 --sparse OUT/sparse.h5 --stage6 stage6.h5 [--cut X]
```

`-i` accepts one file, a directory of files or a text file listing them
(`/pnfs` paths are read through xrootd when a bearer token is present,
otherwise through the NFS mount). `-s pixelmap|h5|sparse` stops the chain
early, `-P` skips the analysis precut in stage 7, `-d cuda:0` evaluates on a
GPU, `-w` and `-O` override the checkpoint and its options file. Stage 5 is
resumable: existing outputs above 10 kB are skipped. Stage 7 holds the whole
stage-6 file in memory, like the notebook it comes from, unless the file is
larger than `$NNBAR_CVN_STREAM_GB` (default 8 GB) or `-S` is given: then
`sparsify_streaming.py` produces the identical output in blocks of events. The `-t` type must be `nnbar` for every
sample (signal and atmospheric background).

Output layout under `-o`:

```
inputs.txt                          list of stage-4 files processed
pixelmap/<name>.root, <name>.log    stage 5
pixelmap.h5                         stage 6
sparse.h5                           stage 7  (network input; event_id and file_name carried along)
predictions.h5                      stage 8  (event_probabilities [N,4], event_predictions, event_targets,
                                              prong_* flattened with prong_event_index; plus event_id,
                                              file_name, truesE, prescut, models and genie/* copied from
                                              sparse.h5 by attach_truth.py)
```

The event classes are, in order, other (NC and nu_tau), nu_mu CC, nu_e CC and
n-nbar; the n-nbar score is the last column of `event_probabilities`. Rows of
`predictions.h5` correspond one to one to the events of `sparse.h5`, and
`event_id` (run, subrun, event) plus the `genie/` group (GENIE record of the
event: interaction mode, kinematics, n-nbar decay channel, final state) are
carried along from the stage-6 file so that scores can be matched to the
truth without the pixel maps. The score cut of the analysis is set
from the ROC on the atmospheric sample for the checkpoint in use; it is not
stored here.

Timing on one RTX 3090 with 20 cores: stage 7 about 30 s (12 s streaming) and
stage 8 about 7 min (batch 64, 13 GB of GPU memory) per 95 000 events.

### Reader quirk

The dataset reader (`training/transformercvn/dataset/minkowski_dataset.py`)
never reads the last event of a file: it computes an inclusive maximum index
and uses it as an exclusive slice bound. `sparsify.py` therefore appends one
dummy event, which is the one dropped, unless `--no-pad-last` is given. The
same happens to the training and validation ranges of a training file (one
event lost per range, harmless).

### Working point

The analysis cut on the n-nbar score is set on the atmospheric background at a
target total background efficiency, precut times CVN. `run_nnbar_inference.sh`
does not need the whole background: any stage-6 part of it gives, with
`find_cut.py --target 2e-4 part_predictions.h5`, the precut efficiency of that
part, the CVN false-positive rate that the target implies, the score cut, and
the signal efficiencies at that cut. `run_atm_part.sh`-style drivers (copy a
part, stage 7 streaming, stage 8, attach truth, add a `part` column, delete the
input) are machine specific and not kept here; the prediction files store
`events_before_precut` and `events_after_precut` as attributes for `find_cut.py`.

## Precuts

Stage 5 keeps an event only if it has at least 100 hits and a Pandora
neutrino vertex. Stage 7 applies the analysis precut
max(N_track, N_shower) > 2 and E_hit < 2 GeV, using the `precut` dataset of
the stage-6 file (track count, shower count, total hit energy), so the
counts before and after it are always available. Efficiencies quoted for the
network must include both.

## Training

The training file is a stage-7 file built from a stage-6 file that mixes the
signal models and the atmospheric background (the production uses
`nnbar-production/cvn/preprocess_atmnu.py`, of which `preprocess/preprocess.py`
is the inference variant with the same layout). With the environment from
`--train`:

```bash
python sparsify/sparsify.py mixed_stage6.h5 training.h5 --jobs 16 --no-pad-last
cd training
python train.py --options_file option_files/example.json --training_file /path/to/training.h5 \
    --gpus 1 --name nnbar_dense --log_dir outputs
```

`option_files/example.json` is the configuration of the published checkpoint
(dense CNN, DenseNet blocks [3,6,12,6,3] with growth rate 32, hidden
dimension 128, 6 encoder layers, 8 heads, GELU, batch 12 with gradient
accumulation over 8, AdamW at 1e-4 with 16 cosine cycles, focal loss with
gamma 1 and event/prong proportion 0.9, pixel values scaled by 1/255). The
training/validation split is 90/10 by position in the file, so shuffle the
stage-6 file before sparsifying (`preprocess_atmnu.py` does; `preprocess.py`
keeps the input order unless `--shuffle`). Checkpoints monitor the one-vs-rest
TPR of the last class (n-nbar) at the first ROC point with FPR >= 0.054%;
the five best are kept under `outputs/<name>/version_<n>/checkpoints/`. Each
checkpoint stores the feature normalisation, so inference needs only the
checkpoint and its `options.json`. Resume with `--checkpoint`; `--gpus 0`
runs on the CPU. `training/README.md` is the author's README (the
`pip install -e` and pytest instructions there refer to files the author
keeps local; use `setup/` instead).

During training the class counts of the network are inferred from the labels
present in `--training_file`. `evaluate.py` here differs from the upstream
copy in one respect: it sizes the event and prong classifiers from the
checkpoint instead, so that a background-only or signal-only file can be
scored (upstream would build a 3-class head for a file without n-nbar
events and fail to load the 4-class checkpoint).

## Stage-6 HDF5 layout

One row per event (N events, at most 20 prongs each). Images are stored
sparse: `cvnmap_index` rows are `(event, plane, wire, tick)` and
`png_cvnmap_index` rows `(event, prong, plane, wire, tick)`, with the pixel
charge in `cvnmap_value` and `png_cvnmap_value`; `cvnmap_shape` and
`png_cvnmap_shape` give the dense shapes `(N, 3, 350, 350)` and
`(N, 20, 3, 350, 350)`. Prong-level arrays: `input_png3d` (8 reconstructed
features x 20: energy, length, start x/y/z, direction x/y/z),
`input_png3d_pad_mask`, `png_ShwTrk` (0 shower, 1 track), and the truth
`mc.png_label`, `mc.png_label_split_photons`, `mc.png_mother`, `png_trueE`,
`png_trueP`. Event-level: `event_id` (run, subrun, event), `input_slice`
(reconstructed energy and vertex), `precut`, `tpc_id`, `model`, `file_name`,
the truth `mc.inter`, `trueE`, `trueVertex`, `trueP`, `mc.genie_*`,
`mc.prim_*`, the standard CVN scores `cvnnue`, `cvnnumu`, ..., and the group
`genie/` with the GENIE record of the event (MCNeutrino and GTruth fields,
`decay_channel` and `g_decay_mode` for n-nbar, the padded final-state list
`fs_*`, `matched`); see the README of `nnbar-production`, section "cvn".
Truth branches that are absent from the input are filled with -1, so the
chain runs unchanged on data or truth-less simulation.

## Validation status (2026-09-28)

* Stages 5 and 6 reproduce the production outputs (identical pixelmap trees
  and HDF5 datasets on the reference files); the run from art output through
  `run_nnbar_inference.sh` on the gpvm is still to be exercised.
* Stage 7 is bit-identical to the author's notebook on a real stage-6 file
  (96 848 events, all datasets compared), and `sparsify_streaming.py` is
  bit-identical to `sparsify.py` on the same file, truth group included.
* Stages 7-8 were run on the four reprocessed hA-BR detector-variation samples
  with the published checkpoint; the CUDA path and `summarize_scores.py` were
  exercised there.
