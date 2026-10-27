# FDVD wires-from-geometry specs

Input specs for `larwirecell/scripts/wct-wires-file.sh`, which runs the
`WireSchemaFileWriter` art module to build a WCT wires file directly from the
LArSoft geometry (so it is consistent-by-construction with the channel map).
See `DUNEWireCell/docs/vd-wires-consistency-tooling.md` for the full story.

Each `*.json` names a GDML, the geometry `Name` that drives the `DUNEWireReadout`
channel-map choice, the sorter tool, and the grouping/correction.  The canonical
output name is `<era><technology>-<variant>`, e.g. `fdvd-...`.

The `gdml` value may be absolute or **relative**: if not an existing path it is
looked up first in the current working directory, then by walking the
`:`-separated `$GDML_PATH`.  An optional `gdml_sha1sum` pins the exact GDML
(fatal error on mismatch); without it the first filename match is used.

## Specs here (all verified to reproduce)

Each was run through `wct-wires-file.sh generate` and its output compared to the shipped
`dune-vd/<detector_name>.json.bz2` (put through `wcwires -c 4`): structure is
identical and geometry agrees to ~10 µm (`wires-compare.py --atol 0.05 --rtol 1e-3`
→ `EQUIVALENT`).

| spec | detector_name | anodes / channels |
|------|---------------|-------------------|
| `fdvd-3view30deg-v2-1x8x6.json`  | dunevd10kt_3view_30deg_v2_refactored_1x8x6ref  | 48 / 43 200 |
| `fdvd-3view30deg-v3-1x8x6.json`  | dunevd10kt_3view_30deg_v3_refactored_1x8x6ref  | 48 / 41 472 |
| `fdvd-3view30deg-v3-1x8x14.json` | dunevd10kt_3view_30deg_v3_refactored_1x8x14ref | 112 / 96 768 |
| `fdvd-3view30deg-v7-1x8x14.json` | dunevd10kt_3view_30deg_v7_refactored_1x8x14    | 112 / 96 768 |

All use `WireReadoutCRUGeom` (non-wrapping) + `GeoObjectSorterCRU`.

## Full-detector spec (reproduces the shipped wrapped v6)

Terminology (per the readout, not the mechanics): a **CRM** is one WCT anode
*face*.  Induction strips are jumpered so that **2 CRM = 1 CRU = one WCT anode
with two faces**, and the induction channels wrap across the two faces.  ("CRP"
is a mechanical unit, not a readout unit, and is avoided here.)

`fdvd-3view30deg-v6-wrapped.json` uses the full
`dunevd10kt_3view_30deg_v6_refactored.gdml` (pinned via `gdml_sha1sum`) +
`GeoObjectSorterCRU60D` + grouping `pair-tpc`.

**Current dunecore serves the genuine full FDVD with the wrapping map.**
`DUNEWireReadout_service.cc` routes a detector `Name` containing `dunevd10kt`:
if it also contains `1x8x`/`1x6x` (a *workspace*) → the non-wrapping
`WireReadoutCRUGeom`; otherwise (the *genuine* detector) →
`CRPWireReadoutGeom` with bridged induction.  The production `Name`
`dunevd10kt_3view_30deg_v6_refactored` has no `1x8x`/`1x6x`, so it selects
`CRPWireReadoutGeom`.  (No name hack is needed; `detector_name` is the real
production name.)  `CRPWireReadoutGeom` additionally requires a `PDMapTool`
(`pdmap_tool`/`pdmap_file`, resolved against `$PDMAP_PATH`) — irrelevant to the
wires but required by the service.

This **reproduces the channel map** of the shipped
`dune10kt-vd/dunevd10kt_3view_30deg_v6_refactored.json.bz2`:

| | anodes | faces | wires | unique channels |
|---|---:|---:|---:|---:|
| generated (`CRPWireReadoutGeom`, `pair-tpc`) | 320 | 640 | 552 960 | **491 520** |
| shipped v6 | 320 | 640 | 552 960 | **491 520** |

- identical channel set (0 … 491 519) and identical wrapping topology (every
  channel has the same number of wire segments; 192 shared channels per anode =
  96 U + 96 V);
- wire→channel assignment agrees to ~0.1 mm (correction-vintage / raw-endpoint
  precision).

The only non-geometric difference is the **WCT anode-ident ordering** (which CRU
is anode 0,1,2,…) — a `pair-tpc` labeling permutation, not a channel or geometry
difference.

### Runtime note for the wrapped spec

```bash
export GDML_PATH=.../dunecore/dunecore/Geometry/gdml
export PDMAP_PATH=.../dunecore/dunecore/ChannelMap      # holds PDVD_PDS_Mapping_*.json
```

> The dunecore *installed* in a given env must be current: an older install whose
> `DUNEWireReadout` lacks the `1x8x`/`1x6x` split maps `dunevd10kt`
> unconditionally to `WireReadoutCRUGeom` (non-wrapping), and would NOT reproduce
> the wrapped v6.

## Run

```bash
export GDML_PATH=.../dunecore/dunecore/Geometry/gdml          # where the GDMLs live
# dunecore Geometry plugins + fcl must be on the runtime path:
export LD_LIBRARY_PATH=.../builds/envs/wcls/dunecore/lib:$LD_LIBRARY_PATH
export CET_PLUGIN_PATH=.../builds/envs/wcls/dunecore/lib:$CET_PLUGIN_PATH

larwirecell/scripts/wct-wires-file.sh generate -i fdvd-3view30deg-v7-1x8x14.json -o out/
```

Produces `out/<canonical>.json.bz2` (the wires file) and
`out/<canonical>.sidecar.json` (input spec + generated FHiCL + channel-map/sorter;
**no** wire summary).  Produce a summary separately:

```bash
larwirecell/scripts/wct-wires-file.sh summary out/<canonical>.json.bz2 -o out/<canonical>.summary.json
#   default = sha256 + totals; add --with-planes and/or --with-anodes for the arrays
```

## Not included (and why)

- The **`v5`** workspace file is wrapped, but its `1x8x` Name routes to the
  non-wrapping `WireReadoutCRUGeom`, so that shipped file is self-inconsistent
  with the map its own Name selects — no spec reproduces it.
- The older **`3view`/`2view` v2** GDMLs use a `gdml_simple_extension` tag the
  current Geant4 GDML parser rejects, so they cannot be loaded here.
