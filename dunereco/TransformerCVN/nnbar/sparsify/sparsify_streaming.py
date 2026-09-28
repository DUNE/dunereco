#!/usr/bin/env python
"""Stage 7 for stage-6 files that do not fit in memory: same output as sparsify.py, produced in
blocks of events (memory: a few GB whatever the input size).

    sparsify_streaming.py IN.h5 OUT.h5 [--no-precut] [--select-mask MASK.npy] [--no-pad-last]
                          [--block-events N] [--block-rows N]

The per-event conversion (sparse_to_sparse) is imported from sparsify.py, so the pixel content
is identical; the streaming version was checked bit for bit against it on a 97k-event file.
The output is first written with extendable (chunked) datasets and then rewritten contiguous,
because the dataset reader memory-maps the pixel arrays.
"""
import argparse
import os
import sys
import time

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sparsify import PASSTHROUGH_GROUPS, pad_row, sparse_to_sparse  # noqa: E402

FNAME_ITEMSIZE = 256


def compute_event_ranges(index_dset, n_events, chunk_rows=2_000_000):
    """(starts, ends) of each event's rows in a sparse index dataset whose column 0 is a
    non-decreasing event number, without loading the dataset."""
    nnz = int(index_dset.shape[0])
    starts = np.zeros(n_events, dtype=np.int64)
    ends = np.zeros(n_events, dtype=np.int64)
    if nnz == 0:
        return starts, ends
    offset = None
    current_event = None
    current_start = 0
    pos = 0
    while pos < nnz:
        take = min(chunk_rows, nnz - pos)
        col0 = index_dset[pos:pos + take, 0].astype(np.int64, copy=False)
        if offset is None:
            offset = int(col0[0])
        col0 = col0 - offset
        changes = np.nonzero(np.diff(col0) != 0)[0] + 1
        bounds = np.concatenate(([0], changes, [take]))
        for s in bounds[:-1]:
            e = int(col0[s])
            seg_start = pos + int(s)
            if current_event is None:
                current_event, current_start = e, seg_start
            if e > current_event:
                starts[current_event], ends[current_event] = current_start, seg_start
                if current_event + 1 < e:
                    starts[current_event + 1:e] = seg_start
                    ends[current_event + 1:e] = seg_start
                current_event, current_start = e, seg_start
        pos += take
    starts[current_event], ends[current_event] = current_start, nnz
    if current_event + 1 < n_events:
        starts[current_event + 1:] = nnz
        ends[current_event + 1:] = nnz
    return starts, ends


def choose_block_end(i, n_events, cvn_starts, cvn_ends, png_starts, png_ends, max_events, max_rows):
    j = min(i + max_events, n_events)
    while j > i + 1:
        if (cvn_ends[j - 1] - cvn_starts[i] <= max_rows) and (png_ends[j - 1] - png_starts[i] <= max_rows):
            break
        j = i + max(1, (j - i) // 2)
    return max(i + 1, j)


class Grower:
    """Append-only HDF5 dataset."""

    def __init__(self, group, name, shape_tail, dtype, chunk_rows):
        self.d = group.create_dataset(name, shape=(0, *shape_tail), maxshape=(None, *shape_tail), dtype=dtype,
                                      chunks=(chunk_rows, *shape_tail))

    def append(self, arr):
        n0 = self.d.shape[0]
        self.d.resize(n0 + len(arr), axis=0)
        if len(arr):
            self.d[n0:] = arr
        return n0


def rewrite_contiguous(path, block_rows=1_000_000):
    tmp = path + ".contiguous.tmp"
    with h5py.File(path, "r") as fin, h5py.File(tmp, "w") as fout:
        def copy(name, src):
            dst = fout.create_dataset(name, shape=src.shape, dtype=src.dtype)
            for i in range(0, src.shape[0], block_rows):
                dst[i:i + block_rows] = src[i:i + block_rows]
        for k, v in fin.items():
            if isinstance(v, h5py.Group):
                for kk, vv in v.items():
                    copy(f"{k}/{kk}", vv)
            else:
                copy(k, v)
        for k, v in fin.attrs.items():
            fout.attrs[k] = v
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--no-precut", action="store_true")
    ap.add_argument("--select-mask", help="optional .npy boolean mask selecting events after the precut")
    ap.add_argument("--rename-prong-target", action="store_true")
    ap.add_argument("--no-pad-last", action="store_true")
    ap.add_argument("--block-events", type=int, default=5000)
    ap.add_argument("--block-rows", type=int, default=1_500_000, help="max sparse rows read per block")
    a = ap.parse_args()
    t0 = time.time()

    with h5py.File(a.input, "r") as fin:
        png_shape = np.array(fin["png_cvnmap_shape"][:], dtype=np.int64)          # (N, max_prongs, F, H, W)
        n_events, max_prongs, png_features = int(png_shape[0]), int(png_shape[1]), int(png_shape[2])
        png_hw = tuple(map(int, png_shape[-2:]))
        cvn_shape = np.array(fin["cvnmap_shape"][:], dtype=np.int64)              # (N, F, H, W)
        if len(cvn_shape) != 4:
            raise ValueError(f"unexpected cvnmap_shape {cvn_shape}")
        cvn_features = int(cvn_shape[1]); cvn_hw = tuple(map(int, cvn_shape[-2:]))
        cvn_pixels_shape_out = np.array([1, cvn_features, *cvn_hw], dtype=np.int64)
        shape_scale = np.ascontiguousarray(np.array([png_hw[0] * png_hw[1], png_hw[1], 1], dtype=np.int64))

        # ---- event selection (same as sparsify.py) ----
        keep = np.ones(n_events, dtype=bool)
        if not a.no_precut and "precut" in fin:
            pc = np.asarray(fin["precut"][:], dtype=np.float32).reshape(n_events, -1)
            if pc.shape[1] != 3:
                raise ValueError(f"precut must have 3 columns, got {pc.shape[1]}")
            keep &= (np.maximum(pc[:, 0].astype(np.int32), pc[:, 1].astype(np.int32)) > 2) & (pc[:, 2] < 2.0)
            print(f"precut keeps {keep.sum()} / {n_events} events", flush=True)
        if a.select_mask:
            sel = np.load(a.select_mask) > 0.5
            if len(sel) != keep.sum() and len(sel) + 1 == keep.sum():
                sel = np.append(sel, [False])
            full = np.zeros(n_events, dtype=bool); full[np.where(keep)[0][sel]] = True; keep = full
            print(f"select mask keeps {keep.sum()} events", flush=True)
        out_n = int(keep.sum())

        use_pad_mask = "input_png3d_pad_mask" in fin and any(
            np.any(fin["input_png3d_pad_mask"][i:i + 8192]) for i in range(0, n_events, 8192))
        cvn_starts, cvn_ends = compute_event_ranges(fin["cvnmap_index"], n_events)
        png_starts, png_ends = compute_event_ranges(fin["png_cvnmap_index"], n_events)
        target_dtype = fin["mc.inter"].dtype
        prong_target_dtype = fin["mc.png_label_split_photons"].dtype
        mask_dtype = fin["input_png3d_pad_mask"].dtype if "input_png3d_pad_mask" in fin else np.int8
        has_event_id, has_file_name = "event_id" in fin, "file_name" in fin
        passthrough = []
        for g in PASSTHROUGH_GROUPS:
            if g in fin and isinstance(fin[g], h5py.Group):
                passthrough += [f"{g}/{k}" for k, d in fin[g].items()
                                if isinstance(d, h5py.Dataset) and d.ndim >= 1 and d.shape[0] == n_events]
        print(f"{n_events} events -> {out_n}; carrying {len(passthrough)} truth datasets; ranges computed [{time.time()-t0:.0f}s]", flush=True)

        os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
        with h5py.File(a.output, "w") as fout:
            cr = 4096
            G = {
                "features": Grower(fout, "features", (max_prongs, 8), np.float32, 512),
                "extra": Grower(fout, "extra", (4,), np.float32, cr),
                "truesE": Grower(fout, "truesE", (1,), np.float32, cr),
                "prescut": Grower(fout, "prescut", (3,), np.float32, cr),
                "models": Grower(fout, "models", (1,), np.int8, cr),
                "event_target": Grower(fout, "event_target", (), target_dtype, cr),
                "prong_target": Grower(fout, "prong_target", (max_prongs,), prong_target_dtype, cr),
                "prong_mask": Grower(fout, "prong_mask", (max_prongs,), mask_dtype, cr),
                "event_compressed_index": Grower(fout, "event_compressed_index", (2,), np.int64, cr),
                "prong_compressed_index": Grower(fout, "prong_compressed_index", (2,), np.int64, cr),
                "event_pixels_coordinates": Grower(fout, "event_pixels_coordinates", (3,), np.int32, 200_000),
                "event_pixels_values": Grower(fout, "event_pixels_values", (cvn_features,), np.float32, 200_000),
                "prong_pixels_coordinates": Grower(fout, "prong_pixels_coordinates", (3,), np.int32, 200_000),
                "prong_pixels_values": Grower(fout, "prong_pixels_values", (png_features,), np.float32, 200_000),
            }
            if has_event_id:
                G["event_id"] = Grower(fout, "event_id", fin["event_id"].shape[1:], fin["event_id"].dtype, cr)
            if has_file_name:
                G["file_name"] = Grower(fout, "file_name", (), fin["file_name"].dtype, cr)
            for k in passthrough:
                G[k] = Grower(fout, k, fin[k].shape[1:], fin[k].dtype, cr)
            fout.create_dataset("event_pixels_shape", data=cvn_pixels_shape_out)
            fout.create_dataset("prong_pixels_shape", data=np.array(list(png_shape[1:]), dtype=np.int64))
            fout.create_dataset("full_pixels_shape", data=np.array(list(png_shape[2:]), dtype=np.int64))

            evt_ptr = prg_ptr = 0
            i = 0
            while i < n_events:
                j = choose_block_end(i, n_events, cvn_starts, cvn_ends, png_starts, png_ends, a.block_events, a.block_rows)
                blk = keep[i:j]
                if not blk.any():
                    i = j; continue
                n_blk = j - i
                feat = np.ascontiguousarray(np.asarray(fin["input_png3d"][i:j], dtype=np.float32).transpose(0, 2, 1)) \
                    if "input_png3d" in fin else np.zeros((n_blk, max_prongs, 8), np.float32)
                extra = np.asarray(fin["input_slice"][i:j], dtype=np.float32).reshape(n_blk, -1)[:, :4] \
                    if "input_slice" in fin else np.zeros((n_blk, 4), np.float32)
                def col(name, dtype, default):
                    return np.asarray(fin[name][i:j], dtype=dtype).reshape(n_blk, -1)[:, :1] if name in fin \
                        else np.full((n_blk, 1), default, dtype=dtype)
                truesE, prescut_b, models = col("trueE", np.float32, 0.0), \
                    (np.asarray(fin["precut"][i:j], dtype=np.float32).reshape(n_blk, -1) if "precut" in fin else np.zeros((n_blk, 3), np.float32)), \
                    col("model", np.int8, 0)
                target = fin["mc.inter"][i:j].reshape(n_blk)
                prong_target = fin["mc.png_label_split_photons"][i:j]
                if a.rename_prong_target:
                    prong_target[prong_target == 9] = 8
                mask = fin["input_png3d_pad_mask"][i:j] if use_pad_mask else (prong_target >= 0).astype(mask_dtype)

                G["features"].append(feat[blk]); G["extra"].append(extra[blk]); G["truesE"].append(truesE[blk])
                G["prescut"].append(prescut_b[blk]); G["models"].append(models[blk]); G["event_target"].append(target[blk])
                G["prong_target"].append(prong_target[blk]); G["prong_mask"].append(mask[blk])
                if has_event_id:
                    G["event_id"].append(fin["event_id"][i:j][blk])
                if has_file_name:
                    G["file_name"].append(np.asarray(fin["file_name"][i:j]).reshape(n_blk)[blk])
                for k in passthrough:
                    G[k].append(fin[k][i:j][blk])

                lo, hi = int(cvn_starts[i]), int(cvn_ends[j - 1])
                idx_cvn = fin["cvnmap_index"][lo:hi] if hi > lo else np.empty((0, 4), np.int64)
                val_cvn = fin["cvnmap_value"][lo:hi] if hi > lo else np.empty((0,), np.float32)
                plo, phi = int(png_starts[i]), int(png_ends[j - 1])
                idx_png = fin["png_cvnmap_index"][plo:phi] if phi > plo else np.empty((0, 5), np.int64)
                val_png = fin["png_cvnmap_value"][plo:phi] if phi > plo else np.empty((0,), np.float32)

                ev_ci, pr_ci = [], []
                ev_co, ev_va, pr_co, pr_va = [], [], [], []
                for local_e in np.nonzero(blk)[0]:
                    e = i + int(local_e)
                    s_, t_ = int(cvn_starts[e] - lo), int(cvn_ends[e] - lo)
                    bi = idx_cvn[s_:t_]
                    coords = np.ascontiguousarray(np.concatenate([np.zeros((len(bi), 1), dtype=bi.dtype), bi[:, 1:4]], axis=1))
                    c, v = sparse_to_sparse(coords, np.asarray(val_cvn[s_:t_], dtype=np.float32), 1, cvn_features, shape_scale)
                    ev_co.append(c.astype(np.int32)); ev_va.append(v); ev_ci.append((evt_ptr, evt_ptr + len(c))); evt_ptr += len(c)
                    s_, t_ = int(png_starts[e] - plo), int(png_ends[e] - plo)
                    coords = np.ascontiguousarray(idx_png[s_:t_][:, 1:5])
                    n_pr = int(max(mask[local_e].sum(), 1))
                    c, v = sparse_to_sparse(coords, np.asarray(val_png[s_:t_], dtype=np.float32), n_pr, png_features, shape_scale)
                    pr_co.append(c.astype(np.int32)); pr_va.append(v); pr_ci.append((prg_ptr, prg_ptr + len(c))); prg_ptr += len(c)
                G["event_pixels_coordinates"].append(np.concatenate(ev_co)); G["event_pixels_values"].append(np.concatenate(ev_va))
                G["prong_pixels_coordinates"].append(np.concatenate(pr_co)); G["prong_pixels_values"].append(np.concatenate(pr_va))
                G["event_compressed_index"].append(np.array(ev_ci, dtype=np.int64)); G["prong_compressed_index"].append(np.array(pr_ci, dtype=np.int64))
                if (i // a.block_events) % 20 == 0:
                    print(f"  {j}/{n_events} events [{time.time()-t0:.0f}s]", flush=True)
                i = j

            n = int(fout["event_target"].shape[0])
            assert n == out_n, (n, out_n)
            print(f"event targets present: {np.unique(fout['event_target'][:]).tolist()}", flush=True)
            if not a.no_pad_last:
                # same dummy event as sparsify.py: one zero pixel, one zero prong, event_id -1
                G["features"].append(np.zeros((1, max_prongs, 8), np.float32)); G["extra"].append(np.zeros((1, 4), np.float32))
                G["truesE"].append(-np.ones((1, 1), np.float32)); G["prescut"].append(np.zeros((1, 3), np.float32))
                G["models"].append(-np.ones((1, 1), np.int8)); G["event_target"].append(np.array([fout["event_target"][0]], dtype=target_dtype))
                G["prong_target"].append(-np.ones((1, max_prongs), prong_target_dtype))
                pm = np.zeros((1, max_prongs), mask_dtype); pm[0, 0] = 1; G["prong_mask"].append(pm)
                e0 = G["event_pixels_coordinates"].append(np.zeros((1, 3), np.int32)); G["event_pixels_values"].append(np.zeros((1, cvn_features), np.float32))
                p0 = G["prong_pixels_coordinates"].append(np.zeros((1, 3), np.int32)); G["prong_pixels_values"].append(np.zeros((1, png_features), np.float32))
                G["event_compressed_index"].append(np.array([[e0, e0 + 1]], np.int64)); G["prong_compressed_index"].append(np.array([[p0, p0 + 1]], np.int64))
                if has_event_id:
                    G["event_id"].append(-np.ones((1, *fin["event_id"].shape[1:]), fin["event_id"].dtype))
                if has_file_name:
                    G["file_name"].append(np.array([b"padding"], dtype=fin["file_name"].dtype))
                for k in passthrough:
                    G[k].append(pad_row(fout[k][:1]))
                n += 1
                print("appended one dummy event (event_id -1) for the reader's off-by-one", flush=True)
    print(f"rewriting {a.output} contiguous [{time.time()-t0:.0f}s]", flush=True)
    rewrite_contiguous(a.output)
    print(f"wrote {a.output}: {n} events [{time.time()-t0:.0f}s]", flush=True)
    return 0 if n > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
