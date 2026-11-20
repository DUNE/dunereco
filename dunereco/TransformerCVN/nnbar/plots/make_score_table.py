#!/usr/bin/env python
"""Write the tables read by plot_nnbar_scores.C from prediction files with truth attached.

    make_score_table.py OUTDIR NAME=LABEL:predictions.h5 [NAME=LABEL:predictions.h5 ...]

Each argument is a sample: NAME (used in file names), a ROOT-latex LABEL for legends, and the
predictions.h5 of stage 8 (with event_id and genie/ attached). Writes OUTDIR/scores.csv
(sample index, GENIE decay mode, n-nbar score; comma separated, no header),
OUTDIR/mode_labels.txt and OUTDIR/detvar_labels.txt.
"""
import os
import sys
import h5py
import numpy as np


def main(outdir, specs):
    os.makedirs(outdir, exist_ok=True)
    labels = {}
    with open(os.path.join(outdir, "scores.csv"), "w") as out, open(os.path.join(outdir, "detvar_labels.txt"), "w") as dl:
        for i, spec in enumerate(specs):
            name, rest = spec.split("=", 1)
            label, path = rest.rsplit(":", 1)
            with h5py.File(path, "r") as f:
                score = f["event_probabilities"][:, -1]
                mode = f["genie/g_decay_mode"][:]
                for m, c in zip(*np.unique(np.stack([mode.astype(str), f["genie/decay_channel"][:].astype(str)], 1), axis=0).T):
                    labels.setdefault(int(m), c)
            np.savetxt(out, np.column_stack([np.full(len(score), i), mode, score]), fmt=["%d", "%d", "%.6f"], delimiter=",")
            dl.write(f"{i}\t{name}\t{label}\n")
            print(f"{name}: {len(score)} events, modes {mode.min()}..{mode.max()}")
    with open(os.path.join(outdir, "mode_labels.txt"), "w") as ml:
        for k in sorted(labels):
            ml.write(f"{k}\t{labels[k]}\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2:])
