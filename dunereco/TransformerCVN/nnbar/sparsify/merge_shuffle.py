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
Memory: the pixel lists of every input are held once (no concatenated copy); the output is written in blocks
of --block-events events, so the peak is about the total pixel size of the inputs plus one block.
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


def rows_of(lo, hi):
    """Row numbers of the pixel ranges [lo, hi) in order."""
    lengths = hi - lo
    return np.repeat(lo - np.concatenate(([0], np.cumsum(lengths)[:-1])), lengths) + np.arange(lengths.sum())


def write_pixels(out, names, parts, order, src, block):
    """Write the pixel lists of the events in `order` (global event numbers) in blocks; `parts` holds per input
    (coords, values, cidx relative to that input's arrays); `src` maps global event -> (input, local event)."""
    c, v, ci = names
    first_c, first_v, first_ci = parts[0]
    lens_all = [p[2][:, 1] - p[2][:, 0] for p in parts]
    sidx, eidx = src[order, 0], src[order, 1]
    lengths = np.empty(len(order), dtype=np.int64)
    for s in range(len(parts)):
        m = sidx == s
        lengths[m] = lens_all[s][eidx[m]]
    ends = np.cumsum(lengths)
    total = int(ends[-1]) if len(ends) else 0
    dc = out.create_dataset(c, shape=(total,) + first_c.shape[1:], dtype=first_c.dtype)
    dv = out.create_dataset(v, shape=(total,) + first_v.shape[1:], dtype=first_v.dtype)
    out.create_dataset(ci, data=np.stack((ends - lengths, ends), axis=1).astype(first_ci.dtype))
    for b0 in range(0, len(order), block):
        b1 = min(b0 + block, len(order))
        r0 = int(ends[b0 - 1]) if b0 else 0
        r1 = int(ends[b1 - 1])
        bc = np.empty((r1 - r0,) + first_c.shape[1:], dtype=first_c.dtype)
        bv = np.empty((r1 - r0,) + first_v.shape[1:], dtype=first_v.dtype)
        pos = (ends[b0:b1] - lengths[b0:b1]) - r0
        for s, (pc, pv, pci) in enumerate(parts):
            m = sidx[b0:b1] == s
            if not m.any():
                continue
            ev = eidx[b0:b1][m]
            src_rows = rows_of(pci[ev, 0], pci[ev, 1])
            dst_rows = rows_of(pos[m], pos[m] + lengths[b0:b1][m])
            bc[dst_rows] = pc[src_rows]
            bv[dst_rows] = pv[src_rows]
        dc[r0:r1] = bc
        dv[r0:r1] = bv


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output")
    ap.add_argument("inputs", nargs="+", help="FILE:LABEL")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--block-events", type=int, default=20000)
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
                pixels[kind][0].append(f[c][lo:hi]); pixels[kind][1].append(f[v][lo:hi])
                pixels[kind][2].append(cidx - lo)
            counts.append(m)
            print(f"{label}: {m} events from {path}", flush=True)

    total = sum(counts)
    order = np.random.default_rng(a.seed).permutation(total)
    with h5py.File(a.output, "w") as out:
        for k, parts in per_event.items():
            out.create_dataset(k, data=np.concatenate(parts)[order])
            parts.clear()
        src = np.concatenate([np.stack((np.full(m, si), np.arange(m)), axis=1) for si, m in enumerate(counts)])
        for kind, names in PIXEL_SETS.items():
            parts = list(zip(*pixels[kind]))
            pixels[kind] = None
            write_pixels(out, names, parts, order, src, a.block_events)
            del parts
            print(f"{kind} pixels written", flush=True)
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
