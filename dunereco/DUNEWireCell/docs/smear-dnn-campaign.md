# Truth-label smearing + 6-channel DNN-ROI campaign (smear-dnn)

Branch: `smear-dnn` (dunereco, based on develop v10_21_01d00).
Validation work area: `/exp/dune/app/users/yuhw/cffm-if/smear-dnn` (dunebuild03).
Runtime stack: cvmfs dunesw `v10_21_01d00` e26:prof (wirecell `v0_37_0`,
larwirecell `v10_03_07`) + this dunereco source via `FHICL_FILE_PATH` /
`WIRECELL_PATH` (see `setup-smear-dnn.sh` in the work area).

Scope: 6 experiments — ProtoDUNEs `protodunevd`, `pdhd`; FDs
`dune10kt-1x2x6`, `dune-vd` (1x8x14), `dune10kt-hd`, `dune10kt-vd` (320 CRM).

## Task 1 — truth labeling: DepoFluxWriter + reco-resolution smearing

Every jsonnet whose emitted graph used `wclsSimChannelSink` (or
`wclsDepoSetSimChannelSink`) was switched to `wclsDepoFluxWriter`, following
the `dune10kt-1x2x6/wcls-sim-drift-simchannel-nf-sp.jsonnet` pattern
(`wclsSimDepoSetSource` with `id_is_track: false` + `DepoSetDrifter` +
writer; the bagger drops out).  Jsonnets already using the writer kept it.
FCL `inputers`/`outputers` updated accordingly
(`wclsSimDepoSetSource:…`, `wclsDepoFluxWriter:postdrift`).

Switched: `protodunevd/wcls-sim-drift-simchannel{,-splusn,-y-simchannel}`,
`pdhd/wcls-sim-drift-simchannel-priorSCE`, and the nf-sp jsonnets of
`dune-vd`, `dune10kt-hd`, `dune10kt-vd`.  (The deposplat jsonnets only
reference the sink from dead code — untouched.)

### energy: 0 (bug fix)

Verified in `larwirecell/Components/DepoFluxWriter.cxx`: a **zero** `energy`
makes the writer store the true depo energy (equivalent to SimChannelSink
`use_energy: true`); any **non-zero** value is a *constant* energy per
deposit.  `protodunevd/wcls-sim-drift-depoflux.jsonnet` and
`pdhd/wcls-sim-drift-simchannel-priorSCE-depoflux.jsonnet` had
`energy: 1 # equivalent to use_energy = true` — i.e. every SimChannel IDE
carried 1 MeV.  Fixed to `energy: 0` with a corrected comment.

### Smearing values

Extra smearing added to every DepoFluxWriter so the truth labels register
onto the deconvolved (SP-filtered) reco charge, added in quadrature to the
post-drift diffusion.  Recipe (per experiment `sp-filters.jsonnet`):

* `smear_long` [ticks] = sigma_t / tick with sigma_t = 1/(2π·f), f =
  `Gaus_wide` = 0.12 MHz for all six experiments → sigma_t = 1.326 µs;
  tick = 0.5 µs → **2.6526** everywhere.
* `smear_tran` [pitch] = 1/(2√π·k) per plane with k from
  `Wire_ind` (U,V) / `Wire_col` (W):

| experiment        | k U/V | k W  | smear_tran                      |
|-------------------|-------|------|---------------------------------|
| dune10kt-1x2x6    | 0.75  | 3.0  | [0.37612, 0.37612, 0.09403]     |
| dune10kt-hd       | 0.75  | 3.0  | [0.37612, 0.37612, 0.09403]     |
| dune-vd (1x8x14)  | 0.75  | 3.0  | [0.37612, 0.37612, 0.09403]     |
| dune10kt-vd       | 0.75  | 3.0  | [0.37612, 0.37612, 0.09403]     |
| protodunevd       | 5.0   | 10.0 | [0.056419, 0.056419, 0.028209]  |
| pdhd              | 0.75  | 10.0 | [0.37612, 0.37612, 0.028209]    |

Runtime support: released larwirecell v10_03_07 already implements
`smear_long`/`smear_tran` (verified in `libWireCellLarsoft.so`).

## Task 2 — 6-channel DNN-ROI for the 4 FD experiments

Following `protodunevd/wcls-nf-sp-dnnroi.jsonnet` (explicit `dnnroi_*`
parameters in the entry jsonnet) and the per-plane `dnnroi_pp.jsonnet`
subgraph:

* New `dnnroi_pp.jsonnet` per FD folder: `dune-vd` and `dune10kt-vd` ported
  from `protodunevd/dnnroi_pp.jsonnet`; `dune10kt-1x2x6` and `dune10kt-hd`
  ported from `pdhd/dnnroi_pp.jsonnet` **without** the PDHD APA0 V-plane
  passthrough quirk (a PDHD hardware issue, not applicable to FD sim).
* Entry nf-sp jsonnets define explicit locals: `dnnroi_nchan = 6`,
  `dnnroi_tick_per_slice = 4`, `dnnroi_output_scale = 1.0` (6-ch .ts bakes
  its normalization), `dnnroi_mask_thresh = 0.2`, `dnnroi_nchunks = 1`,
  `dnnroi_device = 'cpu'`, and the model path:
  * FD-VD → PDVD model
    `/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/WireCell/dune/dnn-roi/pdvd/20260615/pipe_distill_transformer_6ch.ts`
  * FD-HD → PDHD model
    `/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/WireCell/dune/dnn-roi/pdhd/20260615/pipe_distill_transformer_6ch.ts`
* No L1SP for FD (per request).  The 1x2x6 TritonService branch is left
  unchanged (model still from the `dnn_model` fcl param there).
* `dune10kt-hd` gained the `use_dnnroi` option from scratch (extVar,
  conditional `sp_override`, dnnroi stage in the pipelines, `dnnsp` merge in
  the retagger); its fcl blocks now pass `use_dnnroi` (default false).
* The 1x2x6 dnnroi-mode `sp_override` no longer blanks `tight_lf_tag` —
  tight_lf is a 6-channel model input.
* `wirecell_dune.fcl`: new variation blocks
  `dunefd_horizdrift_1x2x6_sim_nfsp_dnnroi`,
  `dune10kt_horizdrift_sim_nfsp_dnnroi`,
  `dune10kt_dunefd_vertdrift_1x8x14_3view_sim_nfsp_dnnroi`,
  `dune10kt_vertdrift_sim_nfsp_dnnroi` (inherit the base block, set
  `use_dnnroi: true`, add `WireCellImg`/`WireCellPytorch` plugins).

## Findings / fixes required beyond the plan

1. **WCT `DNNROIFinding` needs channel padding for FD-VD** *(fix outside
   dunereco — needs to go upstream to wire-cell-toolkit)*.
   FD-VD CRMs have 286 U/V strips; the 6-ch transformer model has two
   stride-2 levels on the channel axis, so the channel count must be a
   multiple of 4 (PDVD's own 476 is; 286 is not → TorchScript
   "Expected size 143 but got size 144" and job abort).  Released wirecell
   v0_37_0 only pads the tick axis.  A `chan_pad_multiple` option
   (zero-pad channel rows before inference, crop after) was added to
   `pytorch/{inc/WireCellPytorch,src}/DNNROIFinding.{h,cxx}` in the local
   checkout `/exp/dune/app/users/yuhw/wct` (branch 0.37.x, **uncommitted**)
   and used at runtime via `WCT_BUILD_DIR`.  The FD-VD `dnnroi_pp.jsonnet`
   set `chan_pad_multiple: 4` (silently ignored by older releases).
2. **`tick_pad_multiple` must be 128** for the 6-ch transformer models
   (5 stride-2 levels post-rebin → post-rebin ticks divisible by 32).
   The pdhd-derived default (`tick_per_slice*16` = 64) only worked for
   6000-tick detectors by luck (6016 is a multiple of 128); dune10kt-hd
   (4492 ticks) crashed.  Both FD-HD `dnnroi_pp.jsonnet` now default 128.
3. **`dune10kt-hd/sp.jsonnet` defined no numbered debug tags** (and
   `dune10kt-1x2x6/sp.jsonnet` lacked `tight_lf_tag`), so OmnibusSigProc
   emitted its unnumbered C++ defaults and all 6-ch DNN inputs
   (`loose_lf%d` etc.) arrived as zeros — DNN output empty for U/V while
   the job "succeeded".  Both sp.jsonnet now set the full per-anode tag set
   (matching dune-vd/dune10kt-vd/pdhd).
4. **dune10kt-vd has no full-VD LArSoft geometry preset in this release**
   (`dunevd10kt_geo` maps to a 1x6x6 config that even loads the HD ROOT
   geometry).  The validation jobs configure it explicitly — and note that
   `services.Geometry` alone is NOT enough for g4: Geant4 steps through
   `services.LArG4Detector.gdmlFileName_`, which must also be set to
   `dunevd10kt_3view_30deg_v6_refactored.gdml` (with
   `volumeNames: [volTPCActive, volCryostat]`), else zero SimEnergyDeposits.
5. **dune10kt-vd WCT-vs-LArSoft X-origin mismatch** *(known issue, needs an
   upstream decision)*: `dune10kt-vd/params-10kt.jsonnet` places the drift
   volumes at x ∈ [-4, 1300] cm (bottom anode ≈ -4 cm, cathode ≈ 646 cm,
   top anode ≈ 1300 cm) while the v6 gdml spans x ∈ [-652, 652] (cathode at
   0).  Depos with x < -4 cm (LArSoft bottom volume) are silently dropped by
   the WCT drifter.  The validation keeps the muon inside the overlap
   (x > 0); a proper fix must align `params-10kt.jsonnet` (or the wires
   file) with the v6 gdml coordinates.
6. PD detsim validation paths were restricted to
   `physics.simulate: [ rns, tpcrawdecoder ]` to avoid the unrelated
   `duneopdet` PDS issues (same class of problem as noted in the issue-485
   campaign doc).

## Validation

Work area `/exp/dune/app/users/yuhw/cffm-if/smear-dnn`: per experiment and
angle, single 10 GeV µ⁻ in the XZ plane at central Y (X = drift, Z = beam;
θ w.r.t. Z = 0°, 45°, 80°, 85°; gun placement from DumpGeometry envelopes in
`geom/`).  Chains (all inside SL7, `setup-smear-dnn.sh`):

* FDs: `gen_<ang>.fcl` → standard g4 fcl → standard detsim fcl with the
  DNN-enabled producer (leaf-flip `use_dnnroi: true` + plugins for
  1x2x6 / dune-vd / dune10kt-hd; block swap to
  `dune10kt_vertdrift_sim_nfsp_dnnroi` for dune10kt-vd) → `sp.root`.
* protodunevd: gen → g4 stage1+2 → `protodunevd_detsim.fcl` (smeared
  depoflux) → `protodunevd_nfsp_dnnroi_mc` reco (`wclsdatavd`) → `sp.root`.
* pdhd: gen → g4 → `standard_detsim_protodunehd.fcl` (smeared depoflux) →
  `standard_reco_protodunehd_MC.fcl` restricted to `wclsdatahd`
  (`protodunehd_nfsp_dnnroi_mc`) → `sp.root`.

Comparison plots (`plot-waveforms2.py`, gallery + matplotlib):

* The APA/CRM with the largest total SimChannel charge is selected using
  exact channel→anode maps built from each workflow's WCT wires file
  (`build-chanmaps.py` → `chanmap-<exp>.json.gz`; the VD and PD channel
  numbering is plane-major in paired-anode blocks, NOT arithmetic).
* 1D (per plane U/V/W of that APA, the max-truth-charge channel): the
  recob::Wire waveform is scaled back to electrons by **dividing** out the
  `wclsFrameSaver` `frame_scale` (0.005 for the FDs and protodunevd, 0.001
  for pdhd) and overlaid with the SimChannel electrons/tick on a single
  shared y-axis (no auto-scaling), zoomed to the signal.
* 2D (per plane of the same APA): two aligned channel-vs-tick panels —
  reco electrons and truth electrons — with a shared color scale; the 1D
  channel is marked.

Files: `<exp>/deg<NN>/plots/<exp>_deg<NN>_{U,V,W}.png` (1D) and
`..._2D_{U,V,W}.png`; job records (`gen_<NN>.fcl`, `gen.root`, `g4.root`,
`sp.root`, logs) in each combo dir.

All 24 experiment×angle combinations completed: 72 1D + 72 2D plots, with
the truth (smeared SimChannel) and reco waveform well registered in time,
comparable widths, and peak amplitudes agreeing at the tens-of-percent
level on the shared electron scale — including the near-drift-direction
80°/85° tracks.  Notes: tracks placed near the readout-window edge lose
the acausal early part of the reco waveform (expected); dune10kt-vd
overlays are self-consistent within WCT but sit in the [x > 0] overlap
region due to finding (5).

## Files touched (branch smear-dnn, not committed)

dunereco `dunereco/DUNEWireCell/`:

* `wirecell_dune.fcl` — DepoFluxWriter inputers/outputers for the switched
  blocks; `use_dnnroi` for dune10kt-hd; four new `_dnnroi` variation blocks.
* `protodunevd/wcls-sim-drift-{simchannel,simchannel-splusn,y-simchannel,depoflux}.jsonnet`
* `pdhd/wcls-sim-drift-simchannel-priorSCE{,-depoflux}.jsonnet`
* `dune10kt-1x2x6/{sp.jsonnet,wcls-sim-drift-simchannel-nf-sp.jsonnet,dnnroi_pp.jsonnet*}`
* `dune10kt-hd/{sp.jsonnet,wcls-sim-drift-simchannel-nf-sp.{jsonnet,fcl},dnnroi_pp.jsonnet*}`
* `dune-vd/{wcls-sim-drift-simchannel-nf-sp.{jsonnet,fcl},dnnroi_pp.jsonnet*}`
* `dune10kt-vd/{wcls-sim-drift-simchannel-nf-sp.{jsonnet,fcl},dnnroi_pp.jsonnet*}`

(`*` = new file; the folder CMakeLists glob `*.jsonnet`, so the new files
install automatically.)

wire-cell-toolkit (`/exp/dune/app/users/yuhw/wct`, 0.37.x, **uncommitted**):

* `pytorch/inc/WireCellPytorch/DNNROIFinding.h`,
  `pytorch/src/DNNROIFinding.cxx` — `chan_pad_multiple` option (finding 1).
