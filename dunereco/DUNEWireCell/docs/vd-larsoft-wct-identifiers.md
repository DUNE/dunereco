# LArSoft ↔ Wire-Cell identifiers for the DUNE Vertical Drift detector

This note explains how the identifier conventions of **LArSoft** and the
**Wire-Cell Toolkit (WCT)** relate for the DUNE far-detector Vertical Drift
(VD, `dunevd10kt` / `dune10kt-vd`), with emphasis on the quantity that causes
the most confusion: the WCT **anode "ident"** and the related **`WirePlaneId`**
packing.

The companion file [`vd-larsoft-wct-identifier-map.csv`](vd-larsoft-wct-identifier-map.csv)
tabulates the full correspondence for every plane of every TPC in the
production wires file.

---

## 1. Why "TPC" is a dangerous word

The same three letters mean different things on the two sides:

| Term | In LArSoft | In Wire-Cell |
|------|-----------|--------------|
| **TPC** | one **drift volume**: a `geo::TPCID = (cryostat, tpc)`. In VD this is one *Charge-Readout Module* (CRM / "volTPC"). | an **anode face** (the readout of one drift volume). WCT config/code *sometimes loosely* calls a whole anode a "TPC". |
| **anode** | — (no exact equivalent; closest is a CRP/CRU) | the top-level readout object, identified by an integer **`ident`**. Has 1 or more **faces**. |
| **face** | ≈ one LArSoft **TPC** | a face of an anode; identified by a small **index** (0 or 1), *not* an ident. |

So a LArSoft **TPC ≈ a WCT anode face**, while a WCT **anode** groups one or
more faces. A "TPC ID" is therefore ambiguous unless you say *whose* TPC.
**In the VD production geometry each WCT anode bundles two LArSoft TPCs**
(the two halves of a wrapped Charge-Readout Plane) — see §4.

The WCT VD configuration avoids the word "TPC" and instead calls each anode a
**CRM** (`dune10kt-vd/params-10kt.jsonnet`: `ncrm = 320`, `volumes[n].name =
"crm%d" % n`, `volumes[n].wires = n` ← *that `wires` field is the anode
number/ident*). There are **320 CRMs = 320 WCT anodes**.

---

## 2. The authoritative source of the anode ident: the wires file

A WCT "wires file" (here
`dune10kt-vd/dunevd10kt_3view_30deg_v6_refactored.json.bz2`) is the *only*
place an anode `ident` is defined. It is a serialized
`wirecell.util.wires.schema.Store` with this nesting:

```
Detector → Anode(ident) → Face(ident) → Plane(ident) → Wire(ident, channel, segment, tail, head) → Point
```

Each object carries its own integer `ident`. For the production VD file:

| Object | Count | `ident` values |
|--------|------:|----------------|
| Detector | 1 | 0 |
| **Anode** | **320** | **0 … 319** (dense) |
| Face | 640 | 0 or 1 within each anode (2 faces/anode) |
| Plane | 1920 | 0, 1, 2 within each face (= U, V, Z) |
| Wire | 552 960 | per-plane wire-in-plane index |
| Channel | 491 520 unique | 0 … 491 519 |

So the **anode ident is simply the integer `0 … 319`** stored in the file.
The interesting question is *how the converter chose it* from the LArSoft dump.

---

## 3. How the anode ident is formed (LArSoft dump → WCT wires file)

### 3.1 The LArSoft dump

The text dump fed to Wire-Cell is produced by the LArSoft analyzer
`dunecore/Geometry/CheckCRPGeometry_module.cc` with `DumpWires: true`. It
writes one line **per `geo::WireID`**:

```
# channel  tpc  plane  wire   sx sy sz   ex ey ez
```

- `channel` — LArSoft **offline channel** number (`PlaneWireToChannel(wid)`).
- `tpc` — LArSoft **TPC** number within the cryostat (= one CRM / drift volume).
- `plane` — **0 = U** (first induction), **1 = V** (second induction),
  **2 = Z** (collection).
- `wire` — **wire-in-plane** index.
- `sx…ez` — wire endpoints in cm.

Repeated rows with the *same channel* are successive **segments** of one
wrapped conductor (see §5).

### 3.2 The converter

The production file is built by the `multitpc` loader with the `dunevd` face
style (`wire-cell-python` `wirecell/util/wires/multitpc.py`). The only lines
that set the identifiers are:

```python
face = (1 + tpc) % 2    # "front" face is 0
apa  = tpc // 2         # the anode ident ("apa") -- marked "pure conjecture"
wid  = wire             # wire ident = LArSoft wire-in-plane index
# channel and segment copied straight through from the dump
```

Therefore, **for the VD wires file:**

> **anode ident = ⌊ LArSoft_TPC / 2 ⌋**
> **face index  = (1 + LArSoft_TPC) mod 2**  → even TPC → face **1**, odd TPC → face **0**

and the inverse:

> **LArSoft_TPC = 2 · (anode ident) + (1 − face index)**

| LArSoft TPC | anode ident | face index |
|---:|---:|---:|
| 0 | 0 | 1 |
| 1 | 0 | 0 |
| 2 | 1 | 1 |
| 3 | 1 | 0 |
| 2k | k | 1 |
| 2k+1 | k | 0 |

Plane ident, wire ident, channel, and segment are **passed through
unchanged** — these are LArSoft quantities preserved verbatim:

- **WCT plane ident == LArSoft plane** (0/1/2 = U/V/Z).
- **WCT wire ident == LArSoft wire-in-plane index.**
- **WCT channel == LArSoft offline channel** (identity; see §5).
- **WCT segment == dump segment order** (0, 1 for wrapped conductors).

> ⚠️ **Caveat — the `tpc // 2` pairing is an assumption.** The converter
> comment labels `apa = tpc // 2` as *"pure conjecture"*. It is correct for
> this geometry *because* consecutive LArSoft TPCs are the two wrapped halves
> of one Charge-Readout Plane, but it is **not** a universal LArSoft rule — see
> §7.

### 3.3 A different (newer) VD path exists

`wire-cell-python` also has two *other* VD converters that produce **different**
anode numbering — do not confuse them with the production file:

- `wirecell/util/wires/dunevd.py :: load(type='3view_30deg')` forces
  `face = 0`, giving **one face per anode** (used for some workspace geometries).
- `wirecell/util/wires/gdml.py :: _convert_vdfd()` is a native GDML→wires
  generator reproducing LArSoft `WireReadoutCRUGeom`: **single-face anodes, one
  anode per CRU, no wrapping**, anode idents assigned by sorting CRUs by center
  **Z ascending then Y ascending**. This matches *current* dunecore FD, which
  differs from the production file described here (see §7).

---

## 4. The VD anode/face geometry as it actually appears in the file

Measured directly from `dunevd10kt_3view_30deg_v6_refactored.json.bz2`:

- **320 anodes**, each with **2 faces** (face ident 0 and 1), each with **3
  planes** (U, V, Z).
- Each **face carries 864 wires** = **286 U + 286 V + 292 Z**. These per-CRM
  counts match the LArSoft generator exactly
  (`generate_dunevd10kt_3view_30deg_*.pl`: `Ind1=286, Ind2=286, Col=292`),
  confirming **one WCT face = one LArSoft CRM / TPC**.
- The two faces of an anode sit at the **same anode-plane X** (both faces point
  the same way) — unlike an HD APA whose two faces straddle it.

**Two drift volumes** (from `params-10kt.jsonnet` and confirmed by wire X):

| anode ident | drift volume | electronics | anode-plane X |
|:-----------:|:------------:|:-----------:|:-------------:|
| 0 … 159 | **bottom** | **BDE** (cold) | −651.9 cm |
| 160 … 319 | **top** | **TDE** | +652.15 cm |

(The induction wire angles mirror between the two volumes: bottom U/V ≈ +60°/−60°
from +Y, top U/V ≈ −60°/+60°.)

---

## 5. Channel numbering and wrapped induction

Channels are **copied from LArSoft unchanged** — a WCT channel *is* the LArSoft
offline channel. Properties observed in the file:

- **491 520 unique channels** = 320 anodes × **1 536 channels/anode**.
- Within each `(TPC, plane)` the channel set is **contiguous** and
  **linear in wire ident**: `channel = channel_min + wire_ident`
  (the CSV gives `channel_min` per row; note the plane's *wire list* is ordered
  by pitch, so list order ≠ ident order).
- **Induction wires wrap; collection does not.** Per anode (CRP):
  - U: the two faces **share 96 channels** (the "Ind1Bot" wrapped conductors);
    476 unique U channels total.
  - V: likewise **96 shared**, 476 unique.
  - Z: **0 shared**, 584 unique.
  - total unique per anode = 476 + 476 + 584 = **1 536**.
- A wrapped conductor appears as **two Wire rows with the same channel** but
  different `segment` (0 on one face, 1 on the other). In the file the wrapped
  conductors carry `segment = 0` on the face-1 side and `segment = 1` on the
  face-0 side.

---

## 6. `WirePlaneId`: packing anode + face + layer

WCT identifies a plane/face context with a single packed integer,
`WireCell::WirePlaneId` (`iface/.../WirePlaneId.{h,cxx}`):

```
bits:   [ apa (anode ident) | face | layer ]
ident = (anode << 4) | (face << 3) | layer

layer_mask = 0x7   face_shift = 3   apa_shift = 4
```

The **layer** field is a *bit mask*, not an index:

| layer | value | meaning |
|-------|------:|---------|
| `kUnknownLayer` | 0 | no layer — pure (anode, face) packing |
| `kUlayer` | 1 | U |
| `kVlayer` | 2 | V |
| `kWlayer` | 4 | W (= LArSoft Z/collection) |
| `kAllLayers` | 7 | "all" — the blessed **anode+face context** marker |

> Note the plane **index** (0,1,2 as stored in the wires file) maps to the
> layer **bit** via `iplane2layer[] = {1, 2, 4}` when WCT reads the file. So a
> stored plane ident of 2 (collection) becomes layer bit 4 (W).

### Identifying an anode *face*

Setting **layer = 0** gives a value that encodes *only* the anode ident and the
face index:

> **`WirePlaneId(layer=0).ident() = (anode_ident << 4) | (face_index << 3)
>  = 16·anode + 8·face`**

and to decode: `anode = ident >> 4`, `face = (ident >> 3) & 1`.
(`WirePlaneId::apa()` returns the anode ident; `face()` returns the face index.)
The code *also* defines `kAllLayers = 7` as the canonical "anode+face context";
both are in use — `layer=0` is the minimal packing, `layer=7` the labelled one.

### Worked examples

| anode | face | face WPID (layer 0) | U WPID (+1) | V WPID (+2) | W WPID (+4) | name (`aNfMpX`) |
|------:|-----:|--------------------:|------------:|------------:|------------:|:----------------|
| 0 | 0 | 0 | 1 | 2 | 4 | `a0f0pU` … |
| 0 | 1 | 8 | 9 | 10 | 12 | `a0f1pU` … |
| 1 | 0 | 16 | 17 | 18 | 20 | `a1f0pU` … |
| 1 | 1 | 24 | 25 | 26 | 28 | — |
| 159 | 1 | 2 552 | 2 553 | 2 554 | 2 556 | — |
| 319 | 1 | 5 112 | 5 113 | 5 114 | 5 116 | — |

The `WirePlaneId::name()` format is `a{apa}f{face}p{layer}` (e.g. `a2f1pU`).

---

## 7. Important caveats about the LArSoft side

LArSoft has **three** different VD channel-map/sorter implementations, chosen
by geometry name in `dunecore/Geometry/DUNEWireReadout_service.cc`. They differ
on exactly the points above, so **confirm which produced your sample**:

| geometry | mapper | TPC sort | TPCset | induction wrap |
|----------|--------|----------|--------|----------------|
| `dunevd10kt` (current FD) | `WireReadoutCRUGeom` | Z then Y (`GeoObjectSorterCRU`), **single drift** | **TPCset == TPC** | **none** — every wire its own channel, hierarchical TPC→plane→wire |
| `protodunevd` | `CRPWireReadoutGeom` | by drift coord first (`…CRU60D`, `SortTPCPDVD`) | TPCset == CRP (4 TPC) | induction wraps, Z does not |
| `dunevdcb1` (cold box) | `ColdBoxWireReadoutGeom` | `GeoObjectSorterCRU` | TPCset == CRP, pairs `{0,2},{1,3}` | first induction wraps |

**The production file documented here (`…v6_refactored`) exhibits two drift
volumes, wrapped U/V induction, and consecutive-TPC CRP pairing** — i.e. it
matches the **CRP/wrapping-style** conventions (consistent with its per-CRM
counts 286/286/292/96 and the `protodunevd`-family readout), **not** the
current single-drift, no-wrap `WireReadoutCRUGeom` FD path. In particular:

- "CRP = TPC // 2" is a property of *this* wrapped file, **not** a universal
  LArSoft rule (current FD groups 4 TPC per CRP; the cold box pairs
  `{0,2}`/`{1,3}`).
- "TPCset" means different things per mapper (== TPC in current FD, == CRP
  elsewhere) — do not assume TPCset = CRP = WCT anode without checking.
- If you regenerate wires from *current* dunecore FD geometry you will get the
  single-face, no-wrap numbering of `gdml.py::_convert_vdfd` (§3.3), which is a
  **different** anode-ident convention than this file.

---

## 8. The CSV

[`vd-larsoft-wct-identifier-map.csv`](vd-larsoft-wct-identifier-map.csv) has one
row per **(anode, face, plane) = (LArSoft TPC, plane)** — 1 920 rows. Columns:

| column | meaning |
|--------|---------|
| `ls_cryostat` | LArSoft cryostat (0; single cryostat assumed) |
| `ls_tpc` | LArSoft TPC = `2·anode + (1−face)` |
| `ls_plane`, `ls_view` | LArSoft plane index / view letter (0 U, 1 V, 2 Z) |
| `ls_drift_volume`, `ls_readout_elec` | bottom/top, BDE/TDE |
| `wct_anode_ident` | WCT anode ident (0…319) |
| `wct_crm_name` | config CRM name `crm{ident}` |
| `wct_face_index` | WCT face index (0/1) |
| `wct_plane_ident` | WCT plane ident (0/1/2) |
| `wct_layer_name`, `wct_layer_bit` | U/V/W and its mask bit 1/2/4 |
| `wpid_face_layer0` | `WirePlaneId` with layer 0 = `16·anode + 8·face` |
| `wpid_plane` | full-plane `WirePlaneId` = `wpid_face_layer0 + layer_bit` |
| `n_wires`, `wire_ident_min/max` | wires in this plane and ident range |
| `channel_min/max`, `n_channels` | LArSoft channel range (contiguous) |
| `segments` | segment values present (`0`, or `0;1` for wrapped) |
| `wire_angle_deg` | wire direction, degrees of (Δz, Δy) |

Reconstruct any single channel: for a row, LArSoft `channel = channel_min +
wire_ident`, where `wire_ident` is the LArSoft wire-in-plane index (0…`n_wires−1`).
