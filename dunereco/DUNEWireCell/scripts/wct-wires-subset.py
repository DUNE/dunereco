#!/usr/bin/env python3
"""Subset a Wire-Cell Toolkit "wires" file to a chosen set of anodes.

A WCT wires file (JSON, optionally bz2/gz) holds a single whole-detector
geometry under a top-level "Store":

    Store = { detectors[], anodes[], faces[], planes[], wires[], points[] }

The arrays are *flat* and cross-referenced by POSITION (array index), while each
object also carries its own `ident` metadata:

    anode.faces  -> positions into faces[]      (anode also has .ident)
    face.planes  -> positions into planes[]     (face  also has .ident)
    plane.wires  -> positions into wires[]      (plane also has .ident)
    wire.tail/head -> positions into points[]   (wire  also has .ident, channel)

WCT's AnodePlane selects its anode BY IDENT (ws_store.anode(ident)), so a subset
need only (a) keep the chosen anodes with their original idents and (b) re-index
faces/planes/wires/points compactly while remapping every positional reference.
That yields a valid, internally-consistent store a fraction of the size.

Why bother: WCT's WireSchemaFile loads and parses the ENTIRE wires file at job
start regardless of how many anodes a job processes.  For FD-VD the whole-detector
file is ~244 MB of JSON (552,960 wires) and its jsoncpp DOM peaks at ~2.8 GB --
the dominant term in the per-job memory floor.  A per-anode (or per-group) subset
file shrinks that cost roughly in proportion to the anodes kept.

The anode idents selected here are matched against each anode's `ident` (for the
dune10kt-vd files these equal the process_anode_indices used by the WCLS config).

Examples
--------
  # keep anodes 0 and 1
  wct-wires-subset.py -a 0,1 full.json.bz2 sub-0_1.json.bz2

  # keep a mixed spec, read plain json, write bz2 (by extension)
  wct-wires-subset.py -a 0-3,10,20-22 full.json out.json.bz2
"""

import argparse
import bz2
import gzip
import json
import os
import re
import sys


def parse_anode_spec(spec):
    """Parse 'N', 'A-B', comma/space separated into a sorted set of ints."""
    ids = set()
    for item in spec.replace(",", " ").split():
        m = re.fullmatch(r"(\d+)-(\d+)", item)
        if m:
            lo, hi = int(m.group(1)), int(m.group(2))
            if lo > hi:
                raise ValueError(f"empty range '{item}' (lo > hi)")
            ids.update(range(lo, hi + 1))
        elif re.fullmatch(r"\d+", item):
            ids.add(int(item))
        else:
            raise ValueError(f"bad anode item '{item}' (want N or A-B)")
    if not ids:
        raise ValueError("no anodes given")
    return ids


def _open_read(path):
    if path.endswith(".bz2"):
        return bz2.open(path, "rt")
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _open_write(path):
    if path.endswith(".bz2"):
        return bz2.open(path, "wt")
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def _reindex(old_positions):
    """Return (sorted list of old positions, {old_pos: new_pos})."""
    ordered = sorted(old_positions)
    return ordered, {old: new for new, old in enumerate(ordered)}


def subset_store(store, keep_idents):
    """Return a new Store holding only the anodes whose ident is in keep_idents,
    with faces/planes/wires/points pruned and all positional refs remapped.
    Reference-list ORDER (pitch order of plane.wires, face order, etc.) is
    preserved; only array positions are compacted."""
    anodes = store["anodes"]
    faces = store["faces"]
    planes = store["planes"]
    wires = store["wires"]
    points = store["points"]

    kept_anodes = [a for a in anodes if a["Anode"]["ident"] in keep_idents]
    found = {a["Anode"]["ident"] for a in kept_anodes}
    missing = keep_idents - found
    if missing:
        raise ValueError("wires file has no anode(s) with ident(s): %s"
                         % ",".join(str(i) for i in sorted(missing)))

    face_pos = {f for a in kept_anodes for f in a["Anode"]["faces"]}
    plane_pos = {p for fp in face_pos for p in faces[fp]["Face"]["planes"]}
    wire_pos = {w for pp in plane_pos for w in planes[pp]["Plane"]["wires"]}
    point_pos = set()
    for wp in wire_pos:
        wd = wires[wp]["Wire"]
        point_pos.add(wd["tail"])
        point_pos.add(wd["head"])

    _, fmap = _reindex(face_pos)
    _, pmap = _reindex(plane_pos)
    wlist, wmap = _reindex(wire_pos)
    plist, ptmap = _reindex(point_pos)

    new_points = [points[i] for i in plist]

    new_wires = []
    for i in wlist:
        wd = dict(wires[i]["Wire"])
        wd["tail"] = ptmap[wd["tail"]]
        wd["head"] = ptmap[wd["head"]]
        new_wires.append({"Wire": wd})

    new_planes = []
    for i in sorted(plane_pos):
        pd = dict(planes[i]["Plane"])
        pd["wires"] = [wmap[w] for w in pd["wires"]]   # keep pitch order
        new_planes.append({"Plane": pd})

    new_faces = []
    for i in sorted(face_pos):
        fd = dict(faces[i]["Face"])
        fd["planes"] = [pmap[p] for p in fd["planes"]]
        new_faces.append({"Face": fd})

    new_anodes = []
    for a in sorted(kept_anodes, key=lambda a: a["Anode"]["ident"]):
        ad = dict(a["Anode"])
        ad["faces"] = [fmap[f] for f in ad["faces"]]
        new_anodes.append({"Anode": ad})

    # detectors reference anodes positionally; the FD-VD file ships none.  If a
    # file has detectors, keep only those whose anode refs all survived, remapped
    # to the new anode positions.
    ident_to_newpos = {ad["Anode"]["ident"]: n for n, ad in enumerate(new_anodes)}
    oldpos_to_ident = {n: a["Anode"]["ident"] for n, a in enumerate(anodes)}
    new_dets = []
    for det in store.get("detectors", []):
        refs = det["Detector"]["anodes"]
        idents = [oldpos_to_ident[r] for r in refs]
        if all(i in ident_to_newpos for i in idents):
            dd = dict(det["Detector"])
            dd["anodes"] = [ident_to_newpos[i] for i in idents]
            new_dets.append({"Detector": dd})

    return {
        "detectors": new_dets,
        "anodes": new_anodes,
        "faces": new_faces,
        "planes": new_planes,
        "wires": new_wires,
        "points": new_points,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Subset a WCT wires file to a set of anodes.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("input", help="input WCT wires file (.json/.json.bz2/.gz)")
    ap.add_argument("output", help="output wires file (compression by extension)")
    ap.add_argument("-a", "--anodes", required=True, metavar="ANODESPEC",
                    help="anode idents to keep: N and/or A-B, comma/space "
                         "separated (e.g. '0-3,10,20-22')")
    args = ap.parse_args(argv)

    try:
        keep = parse_anode_spec(args.anodes)
    except ValueError as e:
        ap.error(str(e))

    with _open_read(args.input) as fp:
        doc = json.load(fp)
    store = doc["Store"]
    nw0 = len(store["wires"])

    new_store = subset_store(store, keep)
    nw1 = len(new_store["wires"])

    with _open_write(args.output) as fp:
        json.dump({"Store": new_store}, fp)

    sz = os.path.getsize(args.output)
    print("wrote %s  (anodes %s: %d/%d anodes, %d/%d wires, %d points, %.1f MB)"
          % (args.output, ",".join(str(i) for i in sorted(keep)),
             len(new_store["anodes"]), len(store["anodes"]),
             nw1, nw0, len(new_store["points"]), sz / 1e6),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
