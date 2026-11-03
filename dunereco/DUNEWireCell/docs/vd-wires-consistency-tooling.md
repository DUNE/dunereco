# Keeping WCT wires consistent with the LArSoft channel map

This note documents tooling added to guarantee that the Wire-Cell Toolkit (WCT)
wire geometry a job uses is **consistent by construction** with the LArSoft
channel map (`WireReadoutCRUGeom` / `CRPWireReadoutGeom` / `ColdBoxWireReadoutGeom`)
that same job loads.  It is the implementation follow-up to
[`vd-wires-larsoft-compatibility.md`](vd-wires-larsoft-compatibility.md), which
documents the mismatch this prevents.

Two mechanisms are provided (the "Option 1" and "Option 3" from the design
discussion); both share one transform and both read the geometry from the live
LArSoft services, so neither can disagree with the installed channel map.

## The shared idea

The anode/face grouping, the channel numbers and the wire segmentation are all
read from the **same** `geo::WireReadoutGeom` the rest of the job uses, via its
generic `NTPCsets` / `TPCsetToTPCs` / `ROPtoWirePlanes` / `ChannelToWire` API.
That API expresses each mapper's structure (CRU→1 TPC/anode no-wrap; CRP→4
TPC/anode wrapped; ColdBox→2 TPC/anode), so **one transform covers every
mapper** — no per-detector code paths, no hand-coded `tpc//2` convention.

Core transform: `larwirecell/Components/WireSchemaBuilder.{h,cxx}`
```cpp
WireCell::WireSchema::Store
wcls::build_wire_schema(const geo::GeometryCore&, const geo::WireReadoutGeom&,
                        const wcls::GroupingPolicy&);
```
It emits a WCT `WireSchema::Store`, then makes it compliant with the wire-schema
conventions (`util/docs/wire-schema.org`) "by construction" by calling the newly
exposed `WireCell::WireSchema::correct()` (same per-plane fixers `wcwires`/`load()`
apply), and finally `WireSchema::validate()`.

### Grouping policy (`wcls::GroupingPolicy.grouping`)

| value | anode / face | reproduces |
|-------|--------------|-----------|
| `tpcset` (default) | anode = LArSoft TPCset (CRP); faces = its TPCs | **mapper-faithful** — recommended |
| `per-tpc` | anode = TPC, one face | wire-cell-python `dunevd` 3view (single-face) |
| `pair-tpc` | anode = TPC/2, face = (1+TPC)%2 | legacy `multitpc` wrapped two-face (eg the shipped `dunevd10kt_..._v6` file) |

The grouping only changes the anode/face *organisation*; the **channel numbers
are identical for every policy** (they come straight from `PlaneWireToChannel`),
so any policy is channel-consistent with the loaded map.

## WCT enhancement (prerequisite)

`WireSchema::correct()` was added so a store built in memory (not from a file)
can be run through the exact fixers `wcwires` uses:

```cpp
// util/inc/WireCellUtil/WireSchema.h
void  correct(StoreDB&, Correction target = Correction::pitch, Correction have = Correction::load);
Store correct(const Store& raw, Correction target = Correction::pitch, Correction have = Correction::load);
```
`load()` was refactored to share the same `fixers_for()` helper, so file-load and
build-in-memory take an identical correction path.  Verified: a hand-built raw
store → `correct()` → `validate()` passes and comes out pitch-ordered.

## Option 1 — write a wires file as a precursor art job

`larwirecell/Modules/WireSchemaFileWriter_module.cc` (art `EDAnalyzer`).  Build
the file once from the job's geometry, then point the "real" jobs at it.

```fcl
#include "services_dune.fcl"            # selects geometry + channel map + sorter
services: @local::dunefd_vd_services    # (whatever your real job uses)

physics.analyzers.wireswriter: {
  module_type: WireSchemaFileWriter
  output:     "dunevd10kt-wires-fromgeo.json.bz2"
  grouping:   "tpcset"      # tpcset | per-tpc | pair-tpc
  correction: 4             # 1=load 2=order 3=direction 4=pitch
  validate:   true
}
physics.ana: [ wireswriter ]
```
Run over a single (even empty) event.  **Record the dunecore tag + geometry fcl
alongside the emitted file** — the channel numbering is only meaningful relative
to them.  The output is written with WCT's own `WireSchema::dump()`, so it loads
and validates with `wcwires` by construction.

## Option 3 — provide the geometry in-process (no file)

Two pieces (the split is required because art and WCT each construct and own
their own objects):

- art service, data owner: `larwirecell/Services/WireSchemaService.{h,_service.cc}`
  — builds the `Store` lazily from the live geometry services.
- WCT forwarder component: `larwirecell/Components/WireSchemaSource.{h,cxx}`
  — a trivial default-constructible `IWireSchema`+`IConfigurable` that, in
  `configure()`, pulls the `Store` from the service.  Registered as the WCT type
  `wclsWireSchemaSource`.

It is a drop-in for `Gen::WireSchemaFile`.  `AnodePlane` is unchanged — it still
just calls `IWireSchema::wire_schema_store()`.

FCL:
```fcl
services.WireSchemaService: {
  grouping: "tpcset"  correction: 4  validate: true
}
# ensure the larwirecell WCT plugin is loaded (it already is in wcls jobs):
#   physics.producers.<wcls>.wcls_main.plugins: [ ..., "WireCellLarsoft" ]
```

Jsonnet (replace the `WireSchemaFile` provider; cf. `cfg/layers/low/anodes.jsonnet`):
```jsonnet
local wires = { type: "wclsWireSchemaSource", name: "", data: {} };
[ { type: "AnodePlane", name: "%d" % vol.wires,
    data: { ident: vol.wires, wire_schema: wc.tn(wires), faces: vol.faces },
    uses: [wires] } for vol in volumes ]
```
The `uses: [wires]` entry makes WCT configure the provider before `AnodePlane`
(the geometry services are already up at that point, so `configure()` can reach
them).

## Build status (xerosere `wcls` env, gcc 12 / LArSoft 10.26 / art 3.14)

- WCT branch `wcls-wires` (commit adding `correct()`) **builds & installs** to
  `installs/envs/wcls`; both `correct()` overloads are present in the installed
  `libWireCellUtil.so`.  The `correct()` change is also unit-verified (raw store
  → `correct()` → `validate()` passes, wires come out pitch-ordered).
- larwirecell (this branch) **builds** against that WCT + the wcls view:
  `libWireCellLarsoft.so` (with `WireSchemaBuilder`, `WireSchemaSource`,
  `WireSchemaService` impl), the `WireSchemaFileWriter` module, and the
  `WireSchemaService` service plugin all compile and link
  (`builds/envs/wcls/larwirecell`).
- `wires-summary.py` / `wires-compare.py` are exercised on the shipped FDVD
  files via `wcwires -c 4` (the reference side of the reproduction test).

Build notes (so a rebuild is reproducible):
- Use the system `g++` 12.2 (== the wcls toolchain) with
  `CMAKE_PREFIX_PATH="installs/envs/wcls;extern/envs/wcls/view"` and
  `PKG_CONFIG_PATH` pointing at the view's `lib/pkgconfig` (lardata pulls FFTW3
  via pkg-config).
- larwirecell needs `-DBOOST_STACKTRACE_USE_NOOP -Wno-error=maybe-uninitialized`
  (the env has no libbacktrace; the view's boost-stacktrace header otherwise
  trips `-Werror`).
- The view's WCT headers (a dir-symlink into the spack WCT 0.37.0 install)
  shadow the install's `WireSchema.h`; the added `correct()` declaration was
  overlaid there so it resolves.  Prefer rebuilding/replacing the env's WCT with
  the `wcls-wires` branch to make this permanent.

## Reproduction test — DONE (FDVD)

The full `lar`-based reproduction was run for the FDVD `1x8x6ref v2` geometry
and **passes**: `WireSchemaFileWriter` built the wires file from the live
`WireReadoutCRUGeom` channel map, and it matches the shipped
`dunevd10kt_3view_30deg_v2_refactored_1x8x6ref.json.bz2` (put through the same
`wcwires -c 4` corrections):

- **structure identical** — 48 anodes, 48 faces, 144 planes, 43 200 wires,
  43 200 channels; **0** per-plane wire/channel count mismatches;
- **geometry equivalent** — max Δ(wire angle) = 1.2e-4°, Δ(pitch) = 3.3e-5 mm,
  Δ(wire-0 center) = 1.0e-2 mm (10 µm).  (The generated geometry is in fact
  *cleaner* — exactly −60.0° vs the shipped file's −60.0001° real-world residual.)
  `wires-compare.py --atol 0.05 --rtol 1e-3` reports `EQUIVALENT`.

### dunecore build recipe used (Geometry only)

dunecore is not a spack root in the wcls env, but its deps are (Geant4, GENIE,
larpandora/reco/ana, nug4/nugen/nusimdata, …).  Only the **Geometry** subdir is
needed for the channel maps/sorters, and it does not use the DUNE-DAQ
dataformats.  So:

1. configure the dunecore clone against the wcls view with:
   - `CMAKE_PREFIX_PATH = installs/envs/wcls ; <stub dir> ; builds/envs/wcls/{nuevdb,artdaq_core,HighFive,TRACE,gallery} ; extern/envs/wcls/view`
   - `Torch_DIR` = the env's py-torch `.../torch/share/cmake/Torch` (larpandora pulls Torch);
   - `PKG_CONFIG_PATH` → view `lib/pkgconfig` (FFTW3);
   - `-DCMAKE_CXX_FLAGS="-DBOOST_STACKTRACE_USE_NOOP -Wno-error=maybe-uninitialized"`;
   - **stub Configs** for `dunedaqdataformats` / `dunedetdataformats` (absent from
     the env; only used by RawDecoding/HDF5Utils/DuneObj/DAQTriggerSim, which are
     not built) — a tiny INTERFACE-IMPORTED target + `set(<pkg>_FOUND TRUE)`.
2. build only: `dunecore_Geometry DUNEWireReadout_service
   dunecore_Geometry_GeoObjectSorterCRU_tool
   dunecore_Geometry_GeoObjectSorterCRU60D_tool
   dunecore_Geometry_GeoObjectSorterCRM_tool`.

### Running more geometries

`scripts/wct-wires-file.sh generate -i <spec>.json -o out/` with
`LD_LIBRARY_PATH`/`CET_PLUGIN_PATH` extended to
`builds/envs/wcls/dunecore/lib`, then
`scripts/wires-compare.py out/<canon>.json.bz2 <wcwires-corrected-existing>`.

Notes learned:
- Use the **with-wires** GDML (the channel map reads wire volumes from it; a
  `_nowires` GDML gives "non-existant wire 0").
- The sorter (`GeoObjectSorterCRU` for the 1x8xN workspaces, `GeoObjectSorterCRU60D`
  for the full `v6`/`v7`) goes on the **Geometry** service.
- Current dunecore picks the VD channel map from the geometry `Name`: a Name with
  `1x8x`/`1x6x` (workspace) → non-wrapping `WireReadoutCRUGeom`; the genuine full
  Name (e.g. `dunevd10kt_3view_30deg_v6_refactored`) → **`CRPWireReadoutGeom`**
  (wrapped).  `CRPWireReadoutGeom` additionally needs a `PDMapTool`
  (`pdmap_tool`/`pdmap_file`).
- The full **wrapped v6** IS reproduced natively (491 520 channels) via
  `gen-wct-wires/fdvd-3view30deg-v6-wrapped.json`.  The **`v5`** workspace file is
  the exception: it is wrapped but its `1x8x` Name selects the non-wrapping map,
  so it is self-inconsistent — see `vd-wires-larsoft-compatibility.md`.
- **Check the installed dunecore**, not just the source: an older install without
  the `1x8x`/`1x6x` split maps all `dunevd10kt` to non-wrapping and will not
  reproduce the wrapped v6.

## Recommended validation (do this regardless of option)

Add a comparison that, for each supported geometry, builds the Store from the
live map and asserts it matches the wires file a job would otherwise load:
`max(channel)+1`, per-`(TPC,plane)` channel ranges, and wire endpoints.  This is
the test that catches an installed-dunecore mismatch (eg an old `DUNEWireReadout`
giving `WireReadoutCRUGeom` where the file expects `CRPWireReadoutGeom`), and it
is how you confirm a given grouping policy reproduces an existing file.

## Relationship to other work

This complements `docs/wct-direct-wires-component.md` (beads `ddm-kdr`), which
evaluated a C++ `IWireSchema` component drop-in; Option 3 here is a concrete
realization of that seam (`IWireSchema::wire_schema_store()`), and Option 1 is
the file-producing variant.

## Note on "Option 2"

The earlier idea of having an art module query the live `WireReadout` service for
its class name and feed that to Jsonnet as a lookup key **does not work**: Jsonnet
is evaluated from FHiCL `params`/`structs` (`WCLS_tool.cc`), and there is no path
for a running art module to inject a value into those `extVar`s.  A human would
have to type the class name in FHiCL — a second place to get it wrong — so Option
2 is not pursued.  Options 1 and 3 avoid the question entirely by reading the
geometry directly.
