# TransformerCVN -- event and prong classifier

`TransformerCVN` evaluates a joint event/particle classifier: a convolutional
encoder turns an image of the whole event and one image per Pandora track and
shower into embeddings, and a transformer encoder combines them to classify
the event and each prong at once. The images are produced in the same way as
for the RegCNN numu energy estimator (3 views x 400 wires x 280 ticks), and
the package is modelled on `RegCNN` (see PR
[#128](https://github.com/DUNE/dunereco/pull/128)).

This README describes the package structure, the art modules and their
configuration, where the trained networks live, and how the n-nbar training
pipeline that feeds this package is organised. Training code is **not** kept
in dunereco; see "Training pipeline" for the repositories that hold it.


## Package structure

* `art/`  -- art modules and fcl:
  `TransformerCVNMapper` (pixel-map producer), `TransformerCVNEvaluator`
  (network evaluation), `TransformerPixelMapProducer` (image construction,
  adapted from `RegPixelMapProducer`).
* `func/` -- data products with ROOT dictionaries:
  `TransformerPixelMap`, `TransformerCVNResult`, `TransformerCVNPID`, and the
  `PType` prong-class enumeration.

The package is only built when `LIBTORCH_DIR` is defined (libtorch is needed
by the evaluator); the `func/` library depends on `RegCNNFunc`.


## Modules

### TransformerCVNMapper

Producer. Builds the event image and one image per track and per shower from
the hits of the event, cropped around the Pandora neutrino vertex. Produces
two collections of `cnn::TransformerPixelMap`:

| Instance label (fcl) | Default | Contents |
|---|---|---|
| `ClusterPMLabel` | `transformercvneventmap` | one map for the whole event |
| `ClusterPMProngLabel` | `transformercvnprongmap` | one map per Pandora track, then one per shower |

Main parameters (`TransformerCVNMapper.fcl`, `standard_transformercvnmapper`):

| Parameter | Default | Meaning |
|---|---|---|
| `HitsModuleLabel` | `linecluster` | hit collection |
| `TrackModuleLabel` / `ShowerModuleLabel` | `pandoraTrack` / `emshower` | prong definitions via hit associations |
| `PFParticleModuleLabel`, `VertexModuleLabel`, `PandoraNuVertexModuleLabel` | `pandora` | vertex used to crop the images |
| `WireLength` x `TdcWidth` | 400 x 280 | image size per view |
| `WireResolution`, `TimeResolution` | 7, 24 | pixel size in wires and ticks (RegCNN numu-energy geometry) |
| `GlobalWireMethod` | 2 | global wire numbering (1 = nue-energy convention) |
| `ProngTagMethod` | 0 | pixel prong tags in the event map: 0 = by shower, 1 = by track |
| `MinClusterHits` | 1 | skip events with fewer hits |

Each `TransformerPixelMap` holds per view the pixel charge (`fPEX`, `fPEY`,
`fPEZ`), the truth purity and label of each pixel, and the prong tag of each
pixel (`fProngTagX/Y/Z`), so the same product serves both training dumps and
evaluation.

### TransformerCVNEvaluator

Producer. Loads a TorchScript model once in `beginJob` and, for every event
with an event map, concatenates the flattened event image and up to
`MaxProngs` prong images into one float tensor of length
`3 * 400 * 280 * (n_prongs + 1)`, calls `forward`, and stores the outputs:

| Product (`std::vector<cnn::TransformerCVNResult>`) | Instance label | Content of `fOutput` |
|---|---|---|
| event result | `EventResultLabel` (`transformercvneventresult`) | 4 event-class scores |
| prong result | `ProngResultLabel` (`transformercvnprongresult`) | 8 scores per prong, prongs concatenated in mapper order |

Parameters (`TransformerCVNEvaluator.fcl`, `standard_transformercvnevaluator`):
`Network` (path of the TorchScript file), `MaxProngs` (20), and the
mapper instance labels `PixelMapInput`, `EventPixelMapInput`,
`ProngPixelMapInput`.

The network is expected to return a tuple `(event_scores, prong_scores)`. The
number of outputs (4 event classes, 8 prong classes) is currently fixed in the
module; a network with a different head needs a code change. Prong class
labels follow `cnn::PType` in `func/TransformerCVNPType.h`.

### Running

`art/transformercvnevaluatejob.fcl` runs mapper and evaluator on reco2 files
of the FD HD 1x2x6 geometry:

```
lar -c transformercvnevaluatejob.fcl -s reco2.root
```


## Trained networks

Model files are **not** stored in git. They are distributed through StashCache
under `/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/CVN/TransformerCVN/`
and selected with the `Network` parameter of the evaluator.

| Model | Path | Notes |
|---|---|---|
| FD HD beam, 2018 production | `FDHD/2018/dune_transformercvn_fd_hd_beam_2018prod.torchscript` | default in `TransformerCVNEvaluator.fcl` |
| FD HD n-nbar vs atmospheric | to be added | weights and the exported TorchScript are pending from the network author |


## Training pipeline (n-nbar search)

Training and evaluation of the network are done outside larsoft. The n-nbar
chain is documented in
[linyan-w/nnbar-production](https://github.com/linyan-w/nnbar-production),
whose README lists every stage, the repositories and feature branches that
hold the code, and the data locations. In short:

| Stage | What | Where |
|---|---|---|
| 1-3 | GENIE n-nbar generation (SK 2021 branching ratios), AddGENIE, G4/detsim/reco | `nnbar-production`, GENIE fork, `dunesim`/`dunesw` feature branches |
| 4 | Hit/wire dump in CVN pixel-map coordinates (`RecoEnergyS` analyzer, TTree) | `dunereco` branch `feature/lwan_recoenergystudies` (`dunereco/RecoEnergyStudies`) |
| 5 | TTree -> flat `pixelmap` TTree (event and prong images, 350x350 window, precuts) | `nnbar-production/cvn/make_text_file_to_root_trks_shws.C` |
| 6 | pixelmap files -> one HDF5 file (sparse `cvnmap_index`/`cvnmap_value`, prong arrays, truth, event ids) | `nnbar-production/cvn/preprocess_atmnu.py` |
| 7 | Sparsify: HDF5 of stage 6 -> the network input schema (per-event prong features and masks, event/prong sparse pixel coordinates and values, targets) | TransformerCVN training toolkit of the network author (not yet public) |
| 8 | Training and evaluation (`train.py`, `evaluate.py`; PyTorch 2.0.1, Lightning 1.9.5) | [KaiwenYu2001/dune-nnbar-transformercvn_v2](https://github.com/KaiwenYu2001/dune-nnbar-transformercvn_v2) |
| 9 | Export of the trained checkpoint to TorchScript for this evaluator | to be added to the repository of stage 8 |

The feature normalisation constants are stored inside the training checkpoint,
so the exported TorchScript must wrap the trainer's `forward` (which applies
them), not the bare network.


## Open items before the n-nbar network can run in art

* The evaluator feeds the network pixel maps only. The n-nbar network also
  takes a feature vector per prong and event-level variables, and its prong
  head has a different number of classes; the evaluator input construction
  and the fixed output sizes need to be generalised.
* The training images (stage 5, 350x350 window around the vertex) and the
  images this mapper produces (400x280) must be made identical, including the
  prong tagging, before a network trained on one is evaluated on the other.
* `evaluate.py`, the n-nbar weights and the TorchScript export are pending
  from the network author; they will be referenced here once available.
