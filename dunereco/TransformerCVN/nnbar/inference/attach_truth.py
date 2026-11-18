#!/usr/bin/env python
"""Copy the per-event identifiers and truth of a stage-7 file (event_id, file_name, truesE, prescut,
models and the genie/ group) into the predictions.h5 of stage 8, row by row.

    attach_truth.py sparse.h5 predictions.h5

evaluate.py writes one row per evaluated event, in file order; the reader never evaluates the
last event of the file (the padding event when sparsify was run with padding). Existing
datasets of the same name in predictions.h5 are replaced.
"""
import sys
import h5py
import numpy as np

COLUMNS = ["event_id", "file_name", "truesE", "prescut", "models"]


def main(sparse_path, pred_path):
    with h5py.File(sparse_path, "r") as s, h5py.File(pred_path, "a") as p:
        n = int(p.attrs.get("num_events", p["event_predictions"].shape[0]))
        n_sparse = int(s["event_target"].shape[0])
        padded = "event_id" in s and bool((s["event_id"][-1] < 0).all())
        if n > n_sparse:
            sys.exit(f"predictions has {n} rows but {sparse_path} only {n_sparse} events")
        if n != n_sparse - 1:
            print(f"WARNING: {n} predictions for {n_sparse} events in {sparse_path} (expected {n_sparse - 1})")
        if not padded:
            print("WARNING: no padding event in the stage-7 file: its last real event was not evaluated")
        names = [c for c in COLUMNS if c in s] + [f"genie/{k}" for k in s["genie"]] if "genie" in s else [c for c in COLUMNS if c in s]
        for name in names:
            if name in p:
                del p[name]
            p.create_dataset(name, data=s[name][:n])
        p.attrs["truth_source"] = sparse_path
        print(f"attached {len(names)} datasets ({n} rows) from {sparse_path} to {pred_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
