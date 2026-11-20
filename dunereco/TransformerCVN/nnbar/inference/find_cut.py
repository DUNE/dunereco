#!/usr/bin/env python
"""Score cut at a target TOTAL background efficiency (precut x CVN), from atmospheric prediction files,
and the signal efficiencies at that cut.

    find_cut.py --target 2e-4 data/atm/part*_predictions.h5 [--signal data/predictions/*_predictions.h5]

The precut efficiency of the background is taken from the attributes events_before_precut /
events_after_precut written by run_atm_part.sh; the CVN false-positive rate is then
target / precut efficiency, and the cut is the score above which that fraction of the precut-passing
background lies.
"""
import argparse
import glob
import h5py
import numpy as np

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("background", nargs="+")
ap.add_argument("--target", type=float, default=2e-4, help="total background efficiency incl. precut (default 0.02%%)")
ap.add_argument("--signal", nargs="*", default=sorted(glob.glob("data/predictions/*_predictions.h5")))
a = ap.parse_args()

scores, before, after = [], 0, 0
for path in a.background:
    with h5py.File(path, "r") as f:
        scores.append(f["event_probabilities"][:, -1])
        before += int(f.attrs["events_before_precut"]); after += int(f.attrs["events_after_precut"])
        print(f"background {path}: {f.attrs['events_before_precut']} events, {f.attrs['events_after_precut']} after precut, {len(scores[-1])} scored")
s = np.sort(np.concatenate(scores))[::-1]
eps_pre = after / before
fpr = a.target / eps_pre
k = int(round(fpr * len(s)))                       # events allowed above the cut
cut = s[k - 1] if k > 0 else 1.0                   # k-th highest score: exactly k events have score >= cut
n_above = int((s > cut).sum()); n_atleast = int((s >= cut).sum())
print(f"\nbackground precut efficiency : {after}/{before} = {eps_pre:.4%}")
print(f"target total efficiency      : {a.target:.4%}  ->  CVN FPR target {fpr:.5%}  ({fpr*len(s):.1f} of {len(s)} scored events)")
print(f"score cut                    : {cut:.6f}   ({n_atleast} events >= cut, {n_above} events > cut)")
print(f"achieved: CVN FPR {n_above/len(s):.5%} (>cut), total background efficiency {eps_pre*n_above/len(s):.5%}")
err = np.sqrt(max(n_above, 1)) / len(s)
print(f"statistical uncertainty on the FPR at this cut: +- {err:.5%} ({n_above} events)")
if a.signal:
    print("\nsignal at this cut (events after precut; the samples overlap with the training set):")
    for path in a.signal:
        with h5py.File(path, "r") as f:
            ss = f["event_probabilities"][:, -1]
        n = len(ss); p = int((ss > cut).sum())
        print(f"  {path.split('/')[-1]:32s} {p:6d}/{n}  CVN eff {p/n:.2%} +- {np.sqrt(p*(1-p/n))/n:.2%}")
print(f"\nCUT {cut:.6f}")
