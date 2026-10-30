# VD wires-file ↔ LArSoft channel-map compatibility

**TL;DR.** A WCT wires file hard-codes the LArSoft *offline channel number* on
every wire. For a Wire-Cell/LArSoft (wcls) job to work, the LArSoft channel map
loaded at run time must number channels the **same way** the wires file was
dumped. DUNE VD has several LArSoft channel-map classes; which one serves a given
geometry is set by `DUNEWireReadout_service.cc` from the geometry `Name`.

> **Correction (supersedes an earlier version of this note).** Current
> **devel/dunecore** splits the `dunevd10kt` branch: a *workspace* Name
> (containing `1x8x`/`1x6x`) → non-wrapping `geo::WireReadoutCRUGeom`; the
> *genuine* full detector (e.g. `dunevd10kt_3view_30deg_v6_refactored`, used by
> `dune10kt-vd`) → **`geo::CRPWireReadoutGeom`** with **bridged (wrapped)
> induction**. So the shipped wrapped v6 file **is consistent** with current
> dunecore and is reproduced natively (491 520 channels) — verified at runtime.
> The earlier "incompatible" conclusion was based on an *older* dunecore
> (`reference/dunecore`) whose `dunevd10kt` branch mapped to `WireReadoutCRUGeom`
> **unconditionally**; that remains the hazard only for an out-of-date install.

This note gives a per-file compatibility table and a way to check any file
yourself. It is a companion to
[`vd-larsoft-wct-identifiers.md`](vd-larsoft-wct-identifiers.md).

---

## 1. Why compatibility is a hard constraint

In wcls, LArSoft objects (`sim::SimChannel`, `raw::RawDigit`, `recob::Wire`)
are keyed by the **LArSoft offline channel**, and the WCT wires file carries
that same integer on each `Wire`. The two numbering schemes are *the same
integers* or nothing lines up. So a wires file is only usable with the LArSoft
channel map that produced its channel column.

The schemes differ in two decisive ways:

- **Wire wrapping.** In a "wrapped" readout, an induction conductor spans the
  two drift volumes of a Charge-Readout Plane (CRP) and the two physical
  segments **share one channel**. In a "non-wrapping" readout every wire is its
  own channel. This changes the **total channel count** for the same geometry.
- **TPC/CRU ordering** (which physical strip gets which number), set by the
  LArSoft *sorter* tool.

A file dumped from a wrapping map simply has *fewer channels than wires*; a file
from a non-wrapping map has *channels == wires*. That single comparison is the
most reliable compatibility fingerprint (§3).

---

## 2. The LArSoft dispatch: geometry name → channel map

`dunecore/Geometry/DUNEWireReadout_service.cc` selects the channel-map class
purely by substring of the detector name (no version branch for VD):

| name substring | channel-map class | source |
|----------------|-------------------|--------|
| `dunevd10kt` + `1x8x`/`1x6x` (workspace) | **`geo::WireReadoutCRUGeom`** (non-wrapping) | `DUNEWireReadout_service.cc:108-110` |
| `dunevd10kt` otherwise (genuine FD) | **`geo::CRPWireReadoutGeom`** (induction wraps) | `DUNEWireReadout_service.cc:113-114` |
| `dunevdcb1` | **`geo::ColdBoxWireReadoutGeom`** (first-induction wraps) | `:109-110` |
| `dunecrpcb` | **`geo::CRPWireReadoutGeom`** (induction wraps) | `:113-114` |
| `protodunevd` | **`geo::CRPWireReadoutGeom`** (induction wraps) | `:138-139` |

Sorter tool (from `dunecore/Geometry/geometry_dune.fcl`):

| geometry | sorter | drift grouping |
|----------|--------|----------------|
| `dunevd10kt_*` (FD + workspaces) | `GeoObjectSorterCRU` | sorts by **Z then Y**, not drift |
| `dunevdcb1*` | `GeoObjectSorterCRU` | — |
| `dunecrpcb2` | `GeoObjectSorterCRU60D` | — |
| `protodunevd_*` | `GeoObjectSorterCRU60D` + `SortTPCPDVD:true` | sorts by **drift** first |

**The key fact:** the current `dunevd10kt` branch keys on the `1x8x`/`1x6x`
substring — workspaces get the non-wrapping `WireReadoutCRUGeom`, the genuine FD
gets the wrapping `CRPWireReadoutGeom`. (An *older* dunecore had no such split and
mapped every `dunevd10kt` name to `WireReadoutCRUGeom`.) So against an *old*
install a wrapped `dunevd10kt` file cannot be numbered the same as what that
dunecore produces.

---

## 3. Fingerprinting a wires file

Decompress and compare wires vs. channels:

```python
import bz2, json
d = json.load(bz2.open("FILE.json.bz2"))["Store"]
wires = [w["Wire"] for w in d["wires"]]
nwires = len(wires)
nchan  = len(set(w["channel"] for w in wires))
faces_per_anode = {len(a["Anode"]["faces"]) for a in d["anodes"]}
print(nwires, nchan, "WRAPPED" if nwires > nchan else "no-wrap", faces_per_anode)
```

| fingerprint | readout style | compatible map family |
|-------------|---------------|-----------------------|
| `nwires == nchan`, 1 face/anode | **independent-CRU, non-wrapping** | `WireReadoutCRUGeom` |
| `nwires > nchan`, 2 faces/anode, segments 0&1 | **wrapped CRP** (2-face + segment representation) | `CRPWireReadoutGeom` family |
| `nwires > nchan`, 1 face/anode, multi-segment channels | **wrapped** (1-face + segment representation) | `ColdBoxWireReadoutGeom` |
| `nwires > nchan`, 1 face/anode, channels shared *across anodes* | **wrapped "split-anode"** (each drift volume a separate anode) | `CRPWireReadoutGeom` / `ColdBoxWireReadoutGeom` |

---

## 4. Per-file compatibility table

All files under `dunereco/dunereco/DUNEWireCell/`. Measured directly.

| wires file | views | anodes | faces/anode | wires | **channels** | wrap? | compatible LArSoft map | sorter | geometry name |
|------------|:-----:|:------:|:-----------:|------:|-------------:|:-----:|------------------------|--------|---------------|
| **dune10kt-vd/** | | | | | | | | | |
| `dunevd10kt_3view_30deg_v6_refactored` | 3v 30° | 320 | 2 | 552 960 | 491 520 | **wrap** | `CRPWireReadoutGeom` ✅ | `GeoObjectSorterCRU60D` | `dunevd10kt` (genuine) |
| `dunevd10kt_3view_30deg_v6_refactored-2023` | 3v 30° | 320 | 2 | 553 600 | 491 520 | **wrap** | `CRPWireReadoutGeom` ✅ | `GeoObjectSorterCRU60D` | `dunevd10kt` (genuine) |
| **dune-vd/** | | | | | | | | | |
| `dunevd10kt_3view_30deg_v7_refactored_1x8x14` | 3v 30° | 112 | 1 | 96 768 | 96 768 | no-wrap | `WireReadoutCRUGeom` ✅ | `GeoObjectSorterCRU` | `dunevd10kt_*` |
| `dunevd10kt-1x8x14-3view-wires-v1` | 3view | 112 | 1 | 96 768 | 96 768 | no-wrap | `WireReadoutCRUGeom` ✅ | `GeoObjectSorterCRU` | `dunevd10kt_*` |
| `dunevd10kt_3view_v2_refactored_1x8x14ref` | 3view | 112 | 1 | 96 768 | 96 768 | no-wrap | `WireReadoutCRUGeom` ✅ | `GeoObjectSorterCRU` | `dunevd10kt_*` |
| `dunevd10kt_3view_30deg_v3_refactored_1x8x14ref` | 3v 30° | 112 | 1 | 96 768 | 96 768 | no-wrap | `WireReadoutCRUGeom` ✅ | `GeoObjectSorterCRU` | `dunevd10kt_*` |
| `dunevd10kt_3view_v2_refactored_1x8x6ref` | 3view | 48 | 1 | 41 472 | 41 472 | no-wrap | `WireReadoutCRUGeom` ✅ | `GeoObjectSorterCRU` | `dunevd10kt_*` |
| `dunevd10kt_3view_30deg_v2_refactored_1x8x6ref` | 3v 30° | 48 | 1 | 43 200 | 43 200 | no-wrap | `WireReadoutCRUGeom` ✅ | `GeoObjectSorterCRU` | `dunevd10kt_*` |
| `dunevd10kt_3view_30deg_v3_refactored_1x8x6ref` | 3v 30° | 48 | 1 | 41 472 | 41 472 | no-wrap | `WireReadoutCRUGeom` ✅ | `GeoObjectSorterCRU` | `dunevd10kt_*` |
| `dunevd10kt_2view_v2_refactored_1x8x6ref` | 2view | 36 | 1 | 33 408 | 33 408 | no-wrap | `WireReadoutCRUGeom` ✅ | `GeoObjectSorterCRU` | `dunevd10kt_*` |
| `dunevd10kt-1x6x6-2view-wires-v1` | 2view | 36 | 1 | 33 408 | 33 408 | no-wrap | `WireReadoutCRUGeom` ✅ | `GeoObjectSorterCRU` | `dunevd10kt_*` |
| `dunevd10kt-1x6x6-3view-wires-v1` | 3view | 36 | 1 | 31 104 | 31 104 | no-wrap | `WireReadoutCRUGeom` ✅ | `GeoObjectSorterCRU` | `dunevd10kt_*` |
| `dunevd10kt-1x6x6-3view30deg-wires-v1` | 3v 30° | 36 | 1 | 32 400 | 32 400 | no-wrap | `WireReadoutCRUGeom` ✅ | `GeoObjectSorterCRU` | `dunevd10kt_*` |
| `dunevd10kt_3view_30deg_v5_refactored_1x8x6ref` | 3v 30° | 24 | 2 | 41 520 | 36 864 | **wrap** | ⚠ file wraps, but the `1x8x` Name → non-wrapping `WireReadoutCRUGeom` | — | *mismatch, see §5* |
| `dunevdcb1-3view-wires-v1` | 3view | 2 | 1 | 3 456 | 3 200 | **wrap** (multi-seg) | `ColdBoxWireReadoutGeom` ✅ | `GeoObjectSorterCRU` | `dunevdcb1*` |
| **dune-vd-coldbox/** | | | | | | | | | |
| `dunevdcb1-3view-wires-v2-splitanode` | 3view | 4 | 1 | 5 120 | 3 200 | **wrap** (split) | `ColdBoxWireReadoutGeom` ✅ | `GeoObjectSorterCRU` | `dunevdcb1*` |
| **dunevd-crp2/** | | | | | | | | | |
| `dunevdcrp2-wires-larsoft-v1` | 3view | 4 | 1 | 4 976 | 3 072 | **wrap** (split) | `CRPWireReadoutGeom` ✅ | `GeoObjectSorterCRU60D` | `dunecrpcb*` |
| **protodunevd/** | | | | | | | | | |
| `protodunevd-wires-larsoft-v5` | 3v 30° | 8 | 2 | 13 840 | 12 288 | **wrap** | `CRPWireReadoutGeom` ✅ | `GeoObjectSorterCRU60D`+`SortTPCPDVD` | `protodunevd_*` |
| `protodunevd-wires-larsoft-v3` | 3v 30° | 8 | 2 | 13 840 | 12 288 | **wrap** | `CRPWireReadoutGeom` ✅ | `GeoObjectSorterCRU60D`+`SortTPCPDVD` | `protodunevd_*` |

✅ = structure matches the map that current dunecore selects for that geometry
name. ⚠ = structural mismatch with the current-dunecore selection; see §5.

(HD files — `dune10kt-hd/`, `pdhd/`, `dune10kt-1x2x6/`, `iceberg/` — use the
HD/APA channel maps, not these VD classes, and are out of scope here.)

---

## 5. The genuine full FDVD (v6) — consistent with current dunecore

`dune10kt-vd/wcls-sim-drift-simchannel-nf-sp.fcl:61` loads
`files_wires: "dunevd10kt_3view_30deg_v6_refactored.json.bz2"` and
`:1 #include "services_dune.fcl"`.

- That wires file is **wrapped**: 552 960 wires but only **491 520 channels**
  (induction conductors bridge across the two CRMs of each CRU; see
  `vd-larsoft-wct-identifiers.md` §5).
- Current dunecore loads this geometry under the Name
  `dunevd10kt_3view_30deg_v6_refactored` (no `1x8x`/`1x6x`), so `DUNEWireReadout`
  selects **`geo::CRPWireReadoutGeom`** — the **wrapping** map — which produces
  exactly **491 520 channels** with bridged induction.

So **the v6 wires file and current dunecore agree** — verified at runtime: a
fresh build of the current source logs `Initializing CRPWireReadoutGeom` for this
Name and yields 491 520 channels, matching the file (identical channel set,
identical wrapping topology, wire↔channel agreement to ~0.1 mm). It is reproduced
natively by `gen-wct-wires/fdvd-3view30deg-v6-wrapped.json`.

### The real hazards (when it still bites)

- **An out-of-date dunecore install.** An older `DUNEWireReadout` lacking the
  `1x8x`/`1x6x` split maps *every* `dunevd10kt` Name to the non-wrapping
  `WireReadoutCRUGeom` (552 960 channels) → disagrees with the wrapped v6. This
  repo's `reference/dunecore` is such an older tree; always check the **installed**
  dunecore, not just the source.
- **The `v5` workspace file.** `dunevd10kt_3view_30deg_v5_refactored_1x8x6ref` is
  a *wrapped* file, but its `1x8x` Name routes to the non-wrapping
  `WireReadoutCRUGeom` — so that shipped file is self-inconsistent with the
  mapper its own Name selects.
- A pure-WCT (non-wcls) job never consults LArSoft, so any wires file is
  self-consistent in isolation; disagreement only bites when LArSoft channel
  numbers must line up (wcls sim/sigproc, truth matching).

### Recommended actions

1. **Confirm before trusting the FD wcls chain:** dump channels from the actual
   geometry you run (`CheckCRPGeometry_module` with `DumpWires:true`) and check
   `max(channel)+1`. If it is **552 960** (non-wrapping) your geometry needs the
   **v7/1x8x14-style** wires file; if **491 520** (wrapping) the **v6** file is
   correct.
2. If on current dunecore, switch `files_wires` to a **non-wrapping** file
   consistent with `WireReadoutCRUGeom`, or regenerate wires straight from the
   geometry GDML with `wcpy util gdml-to-wires -d fdvd DUNEVD.gdml`
   (implemented by `gdml.py::_convert_vdfd`, which reproduces
   `WireReadoutCRUGeom` exactly: single-face anodes, no wrap).
3. Record, alongside every VD wires file, the **dunecore tag + geometry fcl**
   it was dumped from. The channel numbering is only meaningful relative to that.

---

## 6. Quick reference — which file for which geometry

| you are running (LArSoft, current dunecore) | map | use this wires-file family |
|---------------------------|-----|----------------------------|
| genuine `dunevd10kt` FD (Name has no `1x8x`/`1x6x`) | `CRPWireReadoutGeom` (wrap) | `v6_refactored` (wrapped, 491 520 ch) |
| `dunevd10kt` *workspace* (`1x8x6` / `1x8x14` / `1x6x6`) | `WireReadoutCRUGeom` (no wrap) | the matching `1x8x*` / `1x6x6` file (non-wrapping) |
| `dunevdcb1` cold box | `ColdBoxWireReadoutGeom` | `dunevdcb1-3view-wires-v1` or `-v2-splitanode` |
| `dunecrpcb2` CRP2 cold box | `CRPWireReadoutGeom` | `dunevdcrp2-wires-larsoft-v1` |
| `protodunevd` | `CRPWireReadoutGeom` | `protodunevd-wires-larsoft-v5` (or v3) |

(Always check the **installed** dunecore: an old one without the `1x8x`/`1x6x`
split maps genuine `dunevd10kt` to non-wrapping `WireReadoutCRUGeom` and needs a
non-wrapping full file instead — which is not shipped.)
