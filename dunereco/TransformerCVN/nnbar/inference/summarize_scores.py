#!/usr/bin/env python
"""Summarise a predictions.h5 file of stage 8: class counts, n-nbar score quantiles and, with
--cut, the fraction of events above a score cut (with the precut counted from the stage-6 file).

    summarize_scores.py predictions.h5 [--sparse sparse.h5] [--stage6 pixelmap.h5] [--cut 0.95]
"""
import argparse
import h5py
import numpy as np

CLASSES = ["other (NC, nu_tau)", "nu_mu CC", "nu_e CC", "n-nbar"]

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("predictions")
p.add_argument("--sparse", help="stage-7 file (to count the events after the precut)")
p.add_argument("--stage6", help="stage-6 file (to count the events before the precut)")
p.add_argument("--cut", type=float, help="n-nbar score cut; must come from the ROC of the checkpoint in use")
a = p.parse_args()

with h5py.File(a.predictions, "r") as f:
    probs = f["event_probabilities"][:]
    targets = f["event_targets"][:]
    n_eval = int(f.attrs.get("num_events", len(probs)))
score = probs[:, -1]
print(f"predictions        : {a.predictions}")
if a.stage6:
    with h5py.File(a.stage6, "r") as f:
        print(f"events before precut: {f['precut'].shape[0]}")
if a.sparse:
    with h5py.File(a.sparse, "r") as f:
        print(f"events after precut : {f['event_target'].shape[0]}")
print(f"events evaluated    : {n_eval}")
print("true classes        :", {CLASSES[t]: int(n) for t, n in zip(*np.unique(targets, return_counts=True))})
print("predicted classes   :", {CLASSES[i]: int(n) for i, n in enumerate(np.bincount(probs.argmax(1), minlength=4))})
q = np.quantile(score, [0.1, 0.25, 0.5, 0.75, 0.9])
print("n-nbar score quantiles 10/25/50/75/90%:", " ".join(f"{x:.4f}" for x in q))
if a.cut is not None:
    n_pass = int((score > a.cut).sum())
    err = np.sqrt(n_pass * (1 - n_pass / n_eval)) / n_eval
    print(f"score > {a.cut:g}         : {n_pass}  ({n_pass / n_eval:.2%} +- {err:.2%})")
