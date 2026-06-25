# DUNE FD trigger campaign — detsim notes

## Issue #485: `NamedFactory: Failed to find instance 'simdigits0' of class 'wclsFrameSaver'`

Upstream: <https://github.com/WireCell/wire-cell-toolkit/issues/485>

### Symptom

Running `standard_detsim_dunevd10kt_1x8x14_3view_30deg.fcl` (dunesw VD,
1x8x14, 3-view 30°) aborts at WireCell graph construction with:

```
NamedFactory: Failed to find instance "simdigits0" of class "wclsFrameSaver"
```

(art exit status 66). Reproduced on dunebuild03 (dunesw v10_20_03d00 stack).

### Root cause

The producer block used by that fcl,
`dune10kt_dunefd_vertdrift_1x8x14_3view_30deg_sim_nfsp`
(in `dunereco/dunereco/DUNEWireCell/wirecell_dune.fcl`), inherits its
`outputers` list from the parent `dune10kt_dunefd_vertdrift_1x8x14_3view_sim_nfsp`
via `@table::`, but overrides `process_mode` to `"full"`.

The two were out of sync:

* the inherited `outputers` listed `wclsFrameSaver:simdigits0 .. simdigits111`
  (added by commit `1da20f34` "per-APA rawdigits"), which only exist when the
  jsonnet runs in `process_mode = "full-sim-sp"`;
* `process_mode = "full"` makes the jsonnet
  (`dune-vd/wcls-sim-drift-simchannel-nf-sp.jsonnet`) add **no** digit saver to
  the graph (the `full_sim_pipes` `else []` branch).

So the fcl told larwirecell to fetch `wclsFrameSaver:simdigits0` as an outputer,
but no such node existed in the WCT graph → `NamedFactory` failure.

**The `outputers` list and `process_mode` must always agree.** This was not
enforced and was easy to break by inheritance.

A second, latent bug: the `1x8x6` block
(`dune10kt_dunefd_vertdrift_1x8x6_3view_sim_nfsp`) sets `process_mode = "1x8x6"`,
which the dune-vd jsonnet's `tools` selector did not handle (it errored
`unsupported process_mode`), and it also inherited the multi-simdigits outputers.

### Fix (all in `dunereco`, no other repo needed)

Make `outputers` reusable and explicit per detector/configuration, and pair each
nfsp block with a matching `process_mode`.

New files — reusable, per-detector WireCell fcl fragments holding the outputer
(wclsFrameSaver) lists:

* `dunereco/DUNEWireCell/dune-vd/wirecell_dune-vd.fcl`
  * `dunevd_nfsp_outputers_nodigits`     — `postdrift`, `spsignals`
  * `dunevd_nfsp_outputers_singledigits` — + `simdigits`
  * `dunevd_nfsp_outputers_multidigits`  — + `simdigits0 .. simdigits111` (1x8x14)
* `dunereco/DUNEWireCell/dune10kt-vd/wirecell_dune10kt-vd.fcl`
  * `dune10ktvd_nfsp_outputers_nodigits`     — `postdrift`, `spsignals`
  * `dune10ktvd_nfsp_outputers_singledigits` — + `simdigits`
  * (dune10kt-vd jsonnet defines no per-CRM saver, so no `multidigits`)

Mapping (outputers list ⇔ jsonnet `process_mode`):

| list          | process_mode    | products saved                                  |
|---------------|-----------------|-------------------------------------------------|
| `nodigits`    | `full`          | `spsignals` recob::Wire only                    |
| `singledigits`| `single-sim-sp` | + one merged `simdigits` raw::RawDigit (1 APA)  |
| `multidigits` | `full-sim-sp`   | + per-CRM `simdigits<TPC>` raw::RawDigit         |

`dunereco/DUNEWireCell/wirecell_dune.fcl`:

* `#include "dune-vd/wirecell_dune-vd.fcl"` and `#include "dune10kt-vd/wirecell_dune10kt-vd.fcl"`.
* Base `dune10kt_dunefd_vertdrift_1x8x14_3view_sim_nfsp` → `outputers =
  @local::dunevd_nfsp_outputers_nodigits`, `process_mode = "full"`. This makes
  the standard `_30deg_sim_nfsp` block (inherits) and the `1x8x6` block
  (inherits, keeps `process_mode "1x8x6"`) consistent — both now save only
  `spsignals`. **The standard 1x8x14 30° block needs no edit of its own.**
* Added variations:
  * `dune10kt_dunefd_vertdrift_1x8x14_3view_30deg_sim_nfsp_singledigits`
    (`single-sim-sp`, `singledigits`)
  * `dune10kt_dunefd_vertdrift_1x8x14_3view_30deg_sim_nfsp_multidigits`
    (`full-sim-sp`, `multidigits`)
  * `dune10kt_vertdrift_sim_nfsp_singledigits` (`single-sim-sp`, `singledigits`)
* dune10kt-vd base `dune10kt_vertdrift_sim_nfsp` outputers refactored to
  `@local::dune10ktvd_nfsp_outputers_nodigits` (behavior unchanged).
* Net effect: file shrank 2322 → 2251 lines (the 112-line inline list moved into
  the included fcl).

`dunereco/DUNEWireCell/dune-vd/wcls-sim-drift-simchannel-nf-sp.jsonnet`:

* Added `"1x8x6"` to the `tools` selector (`then tools_all`), so
  `process_mode "1x8x6"` is supported and produces only `spsignals`.

> Note on how edits take effect without a rebuild: `wire-cell-cfg/pgrapher/
> experiment/{dune-vd,dune10kt-vd}` are symlinks into `DUNEWireCell`, and
> `DUNEWireCell` is prepended to `FHICL_FILE_PATH` by the setup. So FCL and
> jsonnet edits in the source tree are live. (`nf.jsonnet`, imported but never
> referenced in the emitted graph, is dead code under jsonnet lazy eval — its
> absence is harmless.)

### Validation

Environment: SL7 apptainer, `setup-dom-v10_14_00d00.sh` (dunesw v10_14_00d00).
Method per variation (own folder, records kept): write a test fcl that
`#include`s the standard detsim entry fcl and overrides
`physics.producers.tpcrawdecoder` to the variation block; `fhicl-dump` →
`flat.fcl` (sanity); `lar -n 1 -c flat.fcl -s <input> -o detsim.root`;
`lar -n 1 -c eventdump.fcl -s detsim.root` to confirm products.

dune-vd: input `g4.root`, entry
`standard_detsim_dunevd10kt_1x8x14_3view_30deg.fcl`. Its `simulate` sequence is
`dunefd_vertdrift_detsim_all_systems`, which crashes in an **unrelated** module
(see below), so the path was restricted to `physics.simulate: [ tpcrawdecoder ]`
to validate the WCT module in isolation. Folders under `dune-vd/issue485-tests/`.

dune10kt-vd: input `input.root`, entry `standard_detsim_dunevd10kt.fcl`. Here the
producer is wired through `dunefd_vertdrift_producers` in
`workflow_detsim_dunevd10kt.fcl`, so the override is one level deeper
(`physics.producers.tpcrawdecoder: @local::...`). Its `simulate` sequence is
already `dunefd_vertdrift_detsim_tpc_only` (no opdet) and runs to completion.
Folders under `dune10kt-vd/issue485-tests/`.

| # | detector    | variation    | block / process_mode                                              | APA   | result                                                                 |
|---|-------------|--------------|-------------------------------------------------------------------|-------|------------------------------------------------------------------------|
| 1 | dune-vd     | nodigits     | `..._1x8x14_3view_30deg_sim_nfsp` / `full`                        | all   | PASS — recob::Wire gauss/wiener/dnnsp (96768), no RawDigit             |
| 2 | dune-vd     | singledigits | `..._30deg_sim_nfsp_singledigits` / `single-sim-sp`              | 20    | PASS — recob::Wire gausstpc20/wienertpc20 (864) + RawDigit daqtpc20 (864) |
| 3 | dune-vd     | multidigits  | `..._30deg_sim_nfsp_multidigits` / `full-sim-sp`                 | all   | PASS — merged recob::Wire (96768) + 112× RawDigit daqtpc0..111 (daqtpc20/21=864) |
| 4 | dune10kt-vd | nodigits     | `dune10kt_vertdrift_sim_nfsp` / `full`                           | all   | PASS — recob::Wire gauss/wiener/dnnsp (491520), no RawDigit            |
| 5 | dune10kt-vd | singledigits | `dune10kt_vertdrift_sim_nfsp_singledigits` / `single-sim-sp`     | 142   | PASS — recob::Wire gausstpc142/wienertpc142 (1536) + RawDigit daqtpc142 (1536) |

All five: `fhicl-dump`, `lar` detsim, and `eventdump` exit 0; the
`NamedFactory` error is gone and products are named/typed as expected.
(`single-sim-sp` runs one APA, so `process_tpc_index` was pointed at a populated
CRM: 20 for dune-vd `g4.root`, 142 for dune10kt-vd `input.root`.)

### Unrelated finding: segfault in the photon-detector digitizer

When the **full** `dunefd_vertdrift_detsim_all_systems` path is run (dune-vd
entry fcl, default sequence), the job segfaults — but **not** in WireCell. The
WCT TPC chain completes and saves correct products first; the crash is:

```
opdet::WaveformDigitizerSim::AddPEsToWaveform (... nChannelsPerOpDet=0 ...)
  duneopdet/OpticalDetector/WaveformDigitizerSim_module.cc:540
```

i.e. the optical/photon-detector digitizer in `duneopdet`, with
`nChannelsPerOpDet=0` (a PD geometry/input issue with these g4 inputs). This is
independent of issue #485 and of dunereco/WireCell. It is the reason the dune-vd
validation restricted the path to the WCT module. Worth a separate ticket
against `duneopdet` / the PD geometry for these inputs.

### Reproduction / environment

* Setup script for this issue (SL7):
  `/exp/dune/app/users/yuhw/wct-ci/dune/dune-vd/setup-issue485.sh`
  (mirrors `setup-dom.sh`; uses dunesw v10_20_03d00 + dbrailsf VD localProducts;
  prepends the local dunereco to `FHICL_FILE_PATH`/`WIRECELL_PATH`).
* Validation setup (as used above):
  `/exp/dune/app/users/yuhw/wct-ci/dune/setup-dom-v10_14_00d00.sh`.
* Container wrapper: `/exp/dune/app/users/yuhw/wct-ci/dune/in-sl7-dom.sh`
  (honors `SETUP_DOM=<setup script>`).
* Per-variation test folders (logs, `flat.fcl`, `detsim.root`) kept under
  `/exp/dune/app/users/yuhw/wct-ci/dune/{dune-vd,dune10kt-vd}/issue485-tests/`.

### Status

Changes are staged in the working tree only — **not committed or pushed**
(pending review). Files touched:

* `dunereco/DUNEWireCell/wirecell_dune.fcl` (modified)
* `dunereco/DUNEWireCell/dune-vd/wcls-sim-drift-simchannel-nf-sp.jsonnet` (modified)
* `dunereco/DUNEWireCell/dune-vd/wirecell_dune-vd.fcl` (new)
* `dunereco/DUNEWireCell/dune10kt-vd/wirecell_dune10kt-vd.fcl` (new)

Follow-up for production install (not needed for source-tree validation): ensure
the new per-detector `wirecell_dune-vd.fcl` / `wirecell_dune10kt-vd.fcl` files are
installed so `#include "dune-vd/wirecell_dune-vd.fcl"` resolves in an installed
build (e.g. via the DUNEWireCell `install_fhicl()` rules / subdir CMake).
