#!/usr/bin/env python3
"""Stage 7 of the n-nbar TransformerCVN chain: convert the stage-6 HDF5 file (preprocess.py)
into the input format of the network (the "fully sparse" dataset read by
transformercvn/dataset/minkowski_dataset.py).

    sparsify.py IN.h5 OUT.h5 [--no-precut] [--jobs N] [--select-mask MASK.npy]

Script version of CreateFullySparseDataset.ipynb
(https://github.com/linyan-w/dune-nnbar-transformercvn/blob/nnbar/), with the same logic:
  * the analysis precut max(N_track, N_shower) > 2 and E_hit < 2 GeV is applied here
    (use --no-precut to keep every event);
  * the three planes of a pixel are merged into one 3-channel entry, and the rows of the
    second plane are flipped (y -> 349 - y) as in the training data;
  * prong images are indexed by prong slot; a slot that is valid but has no pixel gets one
    zero entry so that images and prong features stay aligned;
  * prong labels are mc.png_label_split_photons, event labels the raw mc.inter codes
    (the dataset reader maps them to the four event classes).
  * event_id, file_name and the per-event truth group genie/ are carried through unchanged
    (selected with the same mask; the padding event gets -999 / empty values).
For stage-6 files larger than the memory, use sparsify_streaming.py (same output).
Everything is held in memory, as in the notebook: budget roughly the size of the input file.
"""
import argparse
import os
import sys

import h5py
import numba
import numpy as np
import sparse
from joblib import Parallel, delayed


@numba.njit()
def sparse_to_sparse(coords, values, num_prongs, num_features, shape_scale):
    DIMS = 4
    max_items = values.shape[0] + num_prongs

    output_coordinates = np.zeros((max_items, DIMS - 1), dtype=np.int64)
    output_values = np.zeros((max_items, num_features), dtype=np.float32)

    index_mapping = dict()
    num_unique_indices = 0

    for i in range(values.shape[0]):
        v = values[i]
        if v == 0:
            continue

        b = coords[i, 0]
        c = coords[i, 1]
        y = coords[i, 2]
        x = coords[i, 3]

        if num_features == 3 and c == 1:
            y = 349 - y  # second plane is flipped, as in the training data

        flat_index = shape_scale[0] * b + shape_scale[1] * y + shape_scale[2] * x

        if flat_index not in index_mapping:
            output_index = num_unique_indices
            index_mapping[flat_index] = num_unique_indices
            num_unique_indices += 1
        else:
            output_index = index_mapping[flat_index]

        output_coordinates[output_index][0] = b
        output_coordinates[output_index][1] = y
        output_coordinates[output_index][2] = x

        output_values[output_index][c] = v

    non_zero_planes = set(output_coordinates[:num_unique_indices, 0])
    all_planes = set(np.arange(num_prongs, dtype=np.int64))
    missing_planes = all_planes.difference(non_zero_planes)

    for i in missing_planes:
        output_coordinates[num_unique_indices, 0] = i
        num_unique_indices += 1

    output_coordinates = output_coordinates[:num_unique_indices]
    output_values = output_values[:num_unique_indices]

    sorting_indices = np.argsort((output_coordinates * shape_scale).sum(1))
    output_coordinates = np.ascontiguousarray(output_coordinates[sorting_indices])
    output_values = np.ascontiguousarray(output_values[sorting_indices])

    return output_coordinates, output_values


def compress_first_index(indices, shape):
    assert np.all(np.diff(indices[:, 0]) >= 0)
    starts = np.searchsorted(indices[:, 0], np.arange(shape[0]))
    ends = np.concatenate((starts[1:], [indices.shape[0]]))
    return np.stack([starts, ends]).T


PASSTHROUGH_GROUPS = ("genie",)


def load_passthrough(file, num_events):
    """Per-event truth carried through unchanged: every dataset of the groups in PASSTHROUGH_GROUPS
    whose first axis is the event axis (h5py Groups and datasets of other lengths are skipped)."""
    out = {}
    for g in PASSTHROUGH_GROUPS:
        if g not in file or not isinstance(file[g], h5py.Group):
            continue
        for k, d in file[g].items():
            if isinstance(d, h5py.Dataset) and d.ndim >= 1 and d.shape[0] == num_events:
                out[f"{g}/{k}"] = d[:]
    return out


def pad_row(v):
    """One dummy row for a carried dataset: -999 for numbers, False for booleans, empty bytes for strings."""
    row = np.zeros((1, *v.shape[1:]), dtype=v.dtype)
    if v.dtype.kind in "iuf":
        row[...] = -999 if v.dtype.kind != "u" else 0
    return row


def coo_select(index, value, shape, keep):
    """Apply a boolean event mask to a (nnz, ndim) index array via sparse.COO."""
    coo = sparse.COO(np.ascontiguousarray(index.T), value, tuple(map(int, shape)))
    coo = coo[keep]
    return np.ascontiguousarray(coo.coords.T), coo.data, np.array(coo.shape)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="stage-6 HDF5 file (preprocess.py output)")
    ap.add_argument("output", help="network input file (overwritten)")
    ap.add_argument("--no-precut", action="store_true", help="do not apply the analysis precut")
    ap.add_argument("--jobs", type=int, default=max(1, min(8, os.cpu_count() or 1)), help="parallel workers")
    ap.add_argument("--select-mask", help="optional .npy boolean mask selecting events after the precut")
    ap.add_argument("--rename-prong-target", action="store_true", help="map prong label 9 to 8 (historical option)")
    ap.add_argument("--no-pad-last", action="store_true",
                    help="do not append the dummy event that compensates the reader dropping the last event of a file")
    a = ap.parse_args()

    file = h5py.File(a.input, "r")
    num_events, num_prongs, *_ = file["png_cvnmap_shape"][:]

    data = file["input_png3d"][:] if "input_png3d" in file else np.zeros((num_events, 24, num_prongs), dtype=np.float32)
    mask = file["input_png3d_pad_mask"][:] if "input_png3d_pad_mask" in file else np.zeros((num_events, num_prongs), dtype=bool)
    extra = file["input_slice"][:] if "input_slice" in file else np.zeros(num_events, dtype=np.float32)
    truesE = file["trueE"][:] if "trueE" in file else np.zeros(num_events, dtype=np.float32)
    models = file["model"][:] if "model" in file else np.zeros(num_events, dtype=np.int8)
    prescut = file["precut"][:] if "precut" in file else np.zeros(num_events, dtype=np.float32)
    event_id = file["event_id"][:] if "event_id" in file else None
    file_name = file["file_name"][:] if "file_name" in file else None
    passthrough = load_passthrough(file, num_events)   # genie/* and other per-event truth, carried unchanged

    target = file["mc.inter"][:]
    prong_targets = file["mc.png_label_split_photons"][:]

    data = np.ascontiguousarray(data.astype(np.float32).transpose((0, 2, 1)))
    extra = extra.reshape(data.shape[0], -1).astype(np.float32)
    truesE = truesE.reshape(data.shape[0], -1).astype(np.float32)
    models = models.reshape(data.shape[0], -1).astype(np.int8)
    prescut = prescut.reshape(data.shape[0], -1).astype(np.float32)
    if a.rename_prong_target:
        prong_targets[prong_targets == 9] = 8
    if not mask.any():
        mask = prong_targets >= 0

    cvnmap_index = file["cvnmap_index"][:]
    cvnmap_value = file["cvnmap_value"][:]
    cvnmap_shape = file["cvnmap_shape"][:]
    # Event numbers in both sparse indices count from the same origin. The notebook subtracted the
    # first row of each index separately, which shifts every prong image by one event when the
    # first event of the file has no prong (1% of atmospheric events have none); the event-image
    # index always has the first event, so its first row defines the origin for both.
    event_origin = int(cvnmap_index[0][0])
    cvnmap_index[:, 0] -= event_origin

    png_cvnmap_index = file["png_cvnmap_index"][:]
    png_cvnmap_value = file["png_cvnmap_value"][:]
    png_cvnmap_shape = file["png_cvnmap_shape"][:]
    png_cvnmap_index[:, 0] -= event_origin
    file.close()

    # event images get a prong-slot column (always 0) so both maps share (event, prong, plane, y, x)
    cvnmap_index = np.insert(cvnmap_index, 1, np.zeros(cvnmap_index.shape[0], dtype=cvnmap_index.dtype), axis=1)
    cvnmap_shape = np.insert(cvnmap_shape, 1, 1)
    print(f"{num_events} events, {cvnmap_index.shape[0]} event pixels, {png_cvnmap_index.shape[0]} prong pixels", flush=True)

    # canonical (sorted, deduplicated) coordinate order
    cvnmap_index, cvnmap_value, cvnmap_shape = coo_select(cvnmap_index, cvnmap_value, cvnmap_shape, slice(None))

    keep = np.ones(num_events, dtype=bool)
    if not a.no_precut:
        if prescut.shape[1] != 3:
            raise ValueError(f"precut must have 3 columns, got {prescut.shape[1]}")
        n_tracks, n_showers, hit_tot_energy = prescut[:, 0].astype(np.int32), prescut[:, 1].astype(np.int32), prescut[:, 2]
        keep &= (np.maximum(n_tracks, n_showers) > 2) & (hit_tot_energy < 2.0)
        print(f"precut keeps {keep.sum()} / {num_events} events", flush=True)
    if a.select_mask:
        sel = np.load(a.select_mask) > 0.5
        if len(sel) != keep.sum() and len(sel) + 1 == keep.sum():
            sel = np.append(sel, [False])
        full = np.zeros(num_events, dtype=bool)
        full[np.where(keep)[0][sel]] = True
        keep = full
        print(f"select mask keeps {keep.sum()} events", flush=True)

    if not keep.all():
        extra, data, mask, truesE, models, prescut = extra[keep], data[keep], mask[keep], truesE[keep], models[keep], prescut[keep]
        target, prong_targets = target[keep], prong_targets[keep]
        if event_id is not None:
            event_id = event_id[keep]
        if file_name is not None:
            file_name = file_name[keep]
        passthrough = {k: v[keep] for k, v in passthrough.items()}
        cvnmap_index, cvnmap_value, cvnmap_shape = coo_select(cvnmap_index, cvnmap_value, cvnmap_shape, keep)
        png_cvnmap_index, png_cvnmap_value, png_cvnmap_shape = coo_select(png_cvnmap_index, png_cvnmap_value, png_cvnmap_shape, keep)

    for lab in np.unique(models):
        print(f"model {int(lab)}: {int((models == lab).sum())} events", flush=True)

    event_compressed_index = compress_first_index(cvnmap_index, cvnmap_shape)
    prong_compressed_index = compress_first_index(png_cvnmap_index, png_cvnmap_shape)
    shape_scale = np.ascontiguousarray(np.cumprod((*png_cvnmap_shape[-2:], 1)[::-1])[::-1])
    num_features_event, num_features_prong = int(cvnmap_shape[2]), int(png_cvnmap_shape[2])

    def extract_event(i):
        lower, upper = event_compressed_index[i]
        return sparse_to_sparse(cvnmap_index[lower:upper, 1:], cvnmap_value[lower:upper], 1, num_features_event, shape_scale)

    def extract_prong(i):
        lower, upper = prong_compressed_index[i]
        return sparse_to_sparse(png_cvnmap_index[lower:upper, 1:], png_cvnmap_value[lower:upper],
                                max(int(mask[i].sum()), 1), num_features_prong, shape_scale)

    n = len(event_compressed_index)
    print(f"converting {n} events with {a.jobs} workers", flush=True)
    results = Parallel(n_jobs=a.jobs, verbose=0)(delayed(extract_event)(i) for i in range(n))
    lengths = np.array([0] + [len(r[0]) for r in results])
    event_compressed_index = np.stack((np.cumsum(lengths)[:-1], np.cumsum(lengths)[1:]), axis=1)
    all_event_coordinates = np.concatenate([r[0] for r in results]).astype(np.int32)
    all_event_values = np.concatenate([r[1] for r in results])
    del results, cvnmap_index, cvnmap_value

    results = Parallel(n_jobs=a.jobs, verbose=0)(delayed(extract_prong)(i) for i in range(n))
    lengths = np.array([0] + [len(r[0]) for r in results])
    prong_compressed_index = np.stack((np.cumsum(lengths)[:-1], np.cumsum(lengths)[1:]), axis=1)
    all_prong_coordinates = np.concatenate([r[0] for r in results]).astype(np.int32)
    all_prong_values = np.concatenate([r[1] for r in results])
    del results, png_cvnmap_index, png_cvnmap_value

    _, _, *full_pixels_shape = png_cvnmap_shape
    print(f"event targets present: {np.unique(target).tolist()}", flush=True)

    if not a.no_pad_last:
        # minkowski_dataset.py slices [min_index:max_index] and so never reads the last event of a
        # file; append one dummy event (one zero pixel, one zero prong) so that every real event is
        # evaluated. The dummy (event_id -1) is the one the reader drops, so the prediction rows
        # correspond one to one to the real events.
        n += 1
        data = np.concatenate([data, np.zeros((1, *data.shape[1:]), dtype=data.dtype)])
        extra = np.concatenate([extra, np.zeros((1, extra.shape[1]), dtype=extra.dtype)])
        truesE = np.concatenate([truesE, -np.ones((1, truesE.shape[1]), dtype=truesE.dtype)])
        prescut = np.concatenate([prescut, np.zeros((1, prescut.shape[1]), dtype=prescut.dtype)])
        models = np.concatenate([models, -np.ones((1, models.shape[1]), dtype=models.dtype)])
        target = np.concatenate([target, np.array([target[0]], dtype=target.dtype)])
        prong_targets = np.concatenate([prong_targets, -np.ones((1, prong_targets.shape[1]), dtype=prong_targets.dtype)])
        pad_mask = np.zeros((1, mask.shape[1]), dtype=mask.dtype); pad_mask[0, 0] = 1
        mask = np.concatenate([mask, pad_mask])
        e0 = len(all_event_coordinates); p0 = len(all_prong_coordinates)
        all_event_coordinates = np.concatenate([all_event_coordinates, np.zeros((1, 3), dtype=all_event_coordinates.dtype)])
        all_event_values = np.concatenate([all_event_values, np.zeros((1, all_event_values.shape[1]), dtype=all_event_values.dtype)])
        all_prong_coordinates = np.concatenate([all_prong_coordinates, np.zeros((1, 3), dtype=all_prong_coordinates.dtype)])
        all_prong_values = np.concatenate([all_prong_values, np.zeros((1, all_prong_values.shape[1]), dtype=all_prong_values.dtype)])
        event_compressed_index = np.concatenate([event_compressed_index, np.array([[e0, e0 + 1]], dtype=event_compressed_index.dtype)])
        prong_compressed_index = np.concatenate([prong_compressed_index, np.array([[p0, p0 + 1]], dtype=prong_compressed_index.dtype)])
        if event_id is not None:
            event_id = np.concatenate([event_id, -np.ones((1, event_id.shape[1]), dtype=event_id.dtype)])
        if file_name is not None:
            file_name = np.concatenate([file_name, np.array([b"padding"], dtype=file_name.dtype)])
        passthrough = {k: np.concatenate([v, pad_row(v)]) for k, v in passthrough.items()}
        print("appended one dummy event (event_id -1) for the reader's off-by-one", flush=True)
    with h5py.File(a.output, "w") as out:   # contiguous datasets: the reader memory-maps them
        out.create_dataset("features", data=data)
        out.create_dataset("extra", data=extra)
        out.create_dataset("truesE", data=truesE)
        out.create_dataset("prescut", data=prescut)
        out.create_dataset("models", data=models)
        out.create_dataset("event_target", data=target)
        out.create_dataset("event_pixels_coordinates", data=all_event_coordinates)
        out.create_dataset("event_pixels_values", data=all_event_values)
        out.create_dataset("event_pixels_shape", data=cvnmap_shape[1:])
        out.create_dataset("event_compressed_index", data=event_compressed_index)
        out.create_dataset("prong_target", data=prong_targets)
        out.create_dataset("prong_mask", data=mask)
        out.create_dataset("prong_pixels_coordinates", data=all_prong_coordinates)
        out.create_dataset("prong_pixels_values", data=all_prong_values)
        out.create_dataset("prong_pixels_shape", data=png_cvnmap_shape[1:])
        out.create_dataset("prong_compressed_index", data=prong_compressed_index)
        out.create_dataset("full_pixels_shape", data=np.array(full_pixels_shape))
        if event_id is not None:
            out.create_dataset("event_id", data=event_id)
        if file_name is not None:
            out.create_dataset("file_name", data=file_name)
        for k, v in passthrough.items():
            out.create_dataset(k, data=v)
    print(f"wrote {a.output}: {n} events", flush=True)
    return 0 if n > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
