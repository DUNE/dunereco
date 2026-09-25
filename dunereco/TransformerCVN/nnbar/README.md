# TransformerCVN/nnbar -- n-nbar inference chain

Scripts to run the DUNE n-nbar TransformerCVN classifier on RecoEnergyS output
files, from the art dump to per-event scores. Nothing here is compiled; the
scripts are installed into the dunereco `bin` directory so they are on the
PATH after `setup dunereco`. Training is **not** done here (see
[KaiwenYu2001/dune-nnbar-transformercvn_v2](https://github.com/KaiwenYu2001/dune-nnbar-transformercvn_v2));
the production of the simulation up to the art dump is documented in
[linyan-w/nnbar-production](https://github.com/linyan-w/nnbar-production).

| Stage | Script | Input -> output |
|---|---|---|
| 4 | `RecoEnergyS` art analyzer (`dunereco/RecoEnergyStudies`, branch `feature/lwan_recoenergystudies`) | reco2 art file -> `*_cvnpreprocess.root` with the TTree `recoEnergy/WC` (hits, wires, prongs, truth) |
| 5 | `make_text_file_to_root_trks_shws.C` | `recoEnergy/WC` -> `pixelmap` TTree: one row per pixel of the event image and of each prong image (3 planes x 350 x 350 window around the Pandora vertex), plus an `events` tree with the GENIE final state; applies the precuts |
| 6 | `preprocess.py` | pixelmap files -> one HDF5 file: sparse event and prong images, per-prong reconstructed features, truth and CVN scores |
| 7 | `sparsify.py` (TransformerCVN training toolkit, **to be added**) | stage-6 HDF5 -> the network input schema |
| 8 | `evaluate.py` (TransformerCVN training toolkit, **to be added**) | network input + trained weights -> event and prong probabilities |

`run_nnbar_inference.sh` chains stages 5 to 8 and stops, with a message, at
the first stage whose script is not present.


## Running

```bash
# 4. art dump (needs dunereco with the RecoEnergyStudies package)
lar -c recoenergys.fcl -s reco2.root          # writes <reco2>_cvnpreprocess.root

# 5-8. host environment (not inside the SL7 container): ROOT + python venv
source setup_nnbar_cvn.sh                     # add --no-eval to skip the PyTorch packages
run_nnbar_inference.sh -i /path/to/cvnpreprocess/files -o /exp/dune/data/users/$USER/nnbar_eval \
    -t nnbar -j 8 -m ha_br -w weights.ckpt -O options.json
```

`-i` accepts one file, a directory of files or a text file listing them
(`/pnfs` paths are read through xrootd when a bearer token is present,
otherwise through the NFS mount). `-t atmnu` selects the geometry used for
the atmospheric-neutrino background sample. `-s pixelmap|h5|sparse` stops the
chain early. Stage 5 is resumable: existing outputs above 10 kB are skipped.

Output layout under `-o`:

```
inputs.txt          list of stage-4 files processed
pixelmap/<name>.root, <name>.log     stage 5
pixelmap.h5         stage 6
sparse.h5           stage 7
predictions.h5      stage 8
```

The python environment lives in `$NNBAR_CVN_VENV` (default
`/exp/dune/app/users/$USER/nnbar-cvn-venv`, since the gpvm home areas are
small) and is created on first use from `requirements.txt` and
`requirements-eval.txt`. The version pins in the latter are the ones the
network was trained with; newer `rich` and `transformers` releases break
Lightning 1.9.5 and torch 2.0.1.


## Precuts

Stage 5 keeps an event only if it has at least 100 hits and a Pandora
neutrino vertex. The reconstructed track and shower counts and the total hit
energy of every kept event are stored in the `precut` dataset of the HDF5
file, so any later selection on them can be applied or quoted downstream.
Efficiencies quoted for the network must include these cuts.


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
`mc.prim_*`, and the standard CVN scores `cvnnue`, `cvnnumu`, ...

Truth branches that are absent from the input are filled with -1, so the
chain runs unchanged on data or truth-less simulation.

`preprocess.py` is the inference variant of
`nnbar-production/cvn/preprocess_atmnu.py`: same datasets and dtypes, but the
files are kept in the given order (`--shuffle` restores the training
behaviour), the `model` column comes from `--model-tag` instead of the file
name, and no atmospheric background is mixed in.


## Trained weights

Weights and the matching `options.json` are not kept in git. Until they are
placed on StashCache next to the other TransformerCVN models, pass their
location with `-w` and `-O`.
