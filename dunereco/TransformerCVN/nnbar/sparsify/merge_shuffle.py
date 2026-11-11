#!/usr/bin/env python3
"""Merge several stage-7 (sparsify.py) files into one training file with the events shuffled.

    merge_shuffle.py OUT.h5 IN1.h5:LABEL1 IN2.h5:LABEL2 ... [--seed 2026] [--limit N]

The trainer splits a training file by row (the first `train_validation_split` of the rows for training, the
rest for validation), so the events of the different samples must be mixed before training. Every input
must be written with `sparsify.py --no-pad-last` (no dummy last event). Per-event datasets (first dimension =
number of events: features, extra, truesE, prescut, models, event_target, prong_target, prong_mask, event_id,
file_name and the genie/* passthrough group) are permuted; the pixel lists (event_pixels_*, prong_pixels_*)
are gathered in the new event order and their compressed indices rebuilt; shape datasets are copied from the
first input and must agree. A per-event `sample` dataset (int8, the position of the input on the command line)
and the attribute `samples` (the LABELs) record where each event comes from.
--limit N keeps only the first N events of each input (for tests).
"""
import argparse

import h5py
import numpy as np

PIXEL_SETS = {"event": ("event_pixels_coordinates", "event_pixels_values", "event_compressed_index"),
              "prong": ("prong_pixels_coordinates", "prong_pixels_values", "prong_compressed_index")}
SHAPE_SETS = ("event_pixels_shape", "prong_pixels_shape", "full_pixels_shape")


def per_event_names(f, n):
    """Datasets (also inside groups) whose first dimension is the number of events."""
    skip = {x for v in PIXEL_SETS.values() for x in v} | set(SHAPE_SETS)
    names = []
    f.visititems(lambda k, o: names.append(k) if isinstance(o, h5py.Dataset) and k not in skip
                 and o.ndim > 0 and o.shape[0] == n else None)
    return names


def gather_pixels(coords, values, cidx, order):
    """Pixel rows of the events in `order`, and the compressed index of the result."""
    lo, hi = cidx[order, 0], cidx[order, 1]
    lengths = hi - lo
    starts = np.repeat(lo - np.concatenate(([0], np.cumsum(lengths)[:-1])), lengths)
    rows = starts + np.arange(lengths.sum())
    ends = np.cumsum(lengths)
    return coords[rows], values[rows], np.stack((ends - lengths, ends), axis=1).astype(cidx.dtype)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output")
    ap.add_argument("inputs", nargs="+", help="FILE:LABEL")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    files = [s.rsplit(":", 1) if ":" in s else (s, s) for s in a.inputs]

    per_event, pixels, shapes, counts = {}, {k: [[], [], []] for k in PIXEL_SETS}, None, []
    for si, (path, label) in enumerate(files):
        with h5py.File(path, "r") as f:
            n = f["event_target"].shape[0]
            if f["event_id"][-1][0] == -1 if "event_id" in f else False:
                raise SystemExit(f"{path}: ends with a padding event; write the inputs with --no-pad-last")
            m = min(n, a.limit) if a.limit else n
            names = per_event_names(f, n)
            if si == 0:
                keys = set(names)
                shapes = {k: f[k][:] for k in SHAPE_SETS}
            else:
                if set(names) != keys:
                    raise SystemExit(f"{path}: per-event datasets differ: {sorted(set(names) ^ keys)}")
                for k in SHAPE_SETS:
                    if not np.array_equal(f[k][:], shapes[k]):
                        raise SystemExit(f"{path}: {k} differs")
            for k in names:
                per_event.setdefault(k, []).append(f[k][:m])
            per_event.setdefault("sample", []).append(np.full(m, si, dtype=np.int8))
            for kind, (c, v, ci) in PIXEL_SETS.items():
                cidx = f[ci][:m]
                lo, hi = int(cidx[0, 0]), int(cidx[-1, 1])
                off = sum(len(x) for x in pixels[kind][0])
                pixels[kind][0].append(f[c][lo:hi]); pixels[kind][1].append(f[v][lo:hi])
                pixels[kind][2].append(cidx - lo + off)
            counts.append(m)
            print(f"{label}: {m} events from {path}", flush=True)

    total = sum(counts)
    order = np.random.default_rng(a.seed).permutation(total)
    with h5py.File(a.output, "w") as out:
        for k, parts in per_event.items():
            out.create_dataset(k, data=np.concatenate(parts)[order])
            parts.clear()
        for kind, (c, v, ci) in PIXEL_SETS.items():
            coords, values, cidx = (np.concatenate(x) for x in pixels[kind])
            pixels[kind] = None
            coords, values, cidx = gather_pixels(coords, values, cidx, order)
            out.create_dataset(c, data=coords); out.create_dataset(v, data=values); out.create_dataset(ci, data=cidx)
            del coords, values, cidx
        for k, x in shapes.items():
            out.create_dataset(k, data=x)
        out.attrs["samples"] = np.array([lab for _, lab in files], dtype="S")
        out.attrs["sample_counts"] = np.array(counts)
        out.attrs["shuffle_seed"] = a.seed
    print(f"wrote {total} events to {a.output}", flush=True)
    for si, (_, lab) in enumerate(files):
        print(f"  sample {si} {lab}: {counts[si]}", flush=True)


if __name__ == "__main__":
    main()
