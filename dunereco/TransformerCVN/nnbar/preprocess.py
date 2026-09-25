#!/usr/bin/env python3
"""Stage 6 of the n-nbar TransformerCVN chain: merge "pixelmap" ROOT files (output of
make_text_file_to_root_trks_shws.C) into one HDF5 file in the layout consumed by the
sparsify step of the TransformerCVN training toolkit.

    preprocess.py OUT.h5 [--filelist LIST | FILE ...] [--model-tag TAG] [--shuffle]

Inference variant of nnbar-production/cvn/preprocess_atmnu.py: same output datasets and
dtypes, but
  * files are taken in the given order unless --shuffle is set (training shuffles);
  * the "model" label comes from --model-tag (or -1), not from the file name;
  * truth branches that are absent from the input (data, or truth-less MC) are filled
    with -1 instead of aborting;
  * arrays are created on the first file that is actually processed, so a skipped
    first file no longer breaks the run.
Event images are stored sparse as cvnmap_index (event, plane, wire, tick) / cvnmap_value,
prong images as png_cvnmap_index (event, prong, plane, wire, tick) / png_cvnmap_value;
cvnmap_shape and png_cvnmap_shape give the dense shapes (N, 3, 350, 350) and
(N, 20, 3, 350, 350). At most 20 prongs per event are kept.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import tables as tb
import uproot

MODEL_TAGS = ['atm', 'ha_lfg', 'hn_lfg', 'ha_esf', 'hn_esf', 'ha_br', 'hn_br']
FNAME_ITEMSIZE = 256
MAX_GENIE = 32          # padding of the GENIE final-state particle arrays
MAX_PRONGS = 20
IMAGE = (3, 350, 350)   # planes, wires, ticks of one image
PDG_LABELS = {321: 7, 211: 4, 13: 1, 11: 0, 22: 6, 2212: 2, 1000010020: 7, 1000020040: 7,
              1000130270: 7, 1000010030: 7, 3112: 7}
CVN_VARS = ['cvnnue', 'cvnnumu', 'cvnnutau', 'cvnnc', 'cvn0protons', 'cvn1protons', 'cvn2protons',
            'cvnNprotons', 'cvn0pions', 'cvn1pions', 'cvn2pions', 'cvnNpions', 'cvn0pizeros',
            'cvn1pizeros', 'cvn2pizeros', 'cvnNpizeros', 'cvn0neutrons', 'cvn1neutrons',
            'cvn2neutrons', 'cvnNneutrons']


class Tree:
    """Lazy, cached access to the branches of the pixelmap tree; missing branches -> fill."""

    def __init__(self, tree):
        self.tree = tree
        self.cache = {}
        self.missing = set()

    def __getitem__(self, name):
        return self.get(name)

    def get(self, name, fill=-1, dtype=None):
        if name not in self.cache:
            if name in self.tree:
                self.cache[name] = self.tree[name].array(library='np')
            else:
                self.missing.add(name)
                n = self.tree.num_entries
                self.cache[name] = np.full(n, fill, dtype=dtype or (np.int32 if isinstance(fill, int) else np.float32))
        return self.cache[name]


def genie_truth(finh, unique_ids):
    """Per-event GENIE record from the 'events' tree written by the ROOT macro, aligned to
    the event order of the pixelmap tree via unique_id. Returns None if the tree is absent."""
    if 'events' not in finh:
        return None
    ev = finh['events'].arrays(['unique_id', 'prim_pdg', 'prim_E', 'ndaughters', 'daughter_pdg',
                                'daughter_E', 'daughter_depoE', 'daughter_px', 'daughter_py',
                                'daughter_pz'], library='np')
    row = {int(u): i for i, u in enumerate(ev['unique_id'])}
    n = len(unique_ids)
    out = {
        'mc.genie_npart': np.zeros(n, dtype=np.int16),
        'mc.genie_pdg':   np.zeros((n, MAX_GENIE), dtype=np.int32),
        'mc.genie_E':     -np.ones((n, MAX_GENIE), dtype=np.float32),
        'mc.genie_depoE': -np.ones((n, MAX_GENIE), dtype=np.float32),
        'mc.genie_p':     np.zeros((n, 3, MAX_GENIE), dtype=np.float32),
        'mc.prim_pdg':    np.zeros(n, dtype=np.int32),
        'mc.prim_E':      -np.ones(n, dtype=np.float32),
    }
    for j, u in enumerate(unique_ids):
        i = row.get(int(u))
        if i is None:
            print(f"[warn] unique_id {u} missing from events tree", flush=True)
            continue
        k = min(int(ev['ndaughters'][i]), MAX_GENIE)
        out['mc.genie_npart'][j] = ev['ndaughters'][i]
        out['mc.genie_pdg'][j, :k] = ev['daughter_pdg'][i][:k]
        out['mc.genie_E'][j, :k] = ev['daughter_E'][i][:k]
        out['mc.genie_depoE'][j, :k] = ev['daughter_depoE'][i][:k]
        out['mc.genie_p'][j, 0, :k] = ev['daughter_px'][i][:k]
        out['mc.genie_p'][j, 1, :k] = ev['daughter_py'][i][:k]
        out['mc.genie_p'][j, 2, :k] = ev['daughter_pz'][i][:k]
        out['mc.prim_pdg'][j] = ev['prim_pdg'][i]
        out['mc.prim_E'][j] = ev['prim_E'][i]
    return out


def copy_local(fin, timeout_s=90):
    """Stage a (pnfs) file locally; returns the local path or None."""
    tmp = f"stg_{os.path.basename(fin)}"
    try:
        subprocess.run(["timeout", f"{timeout_s}s", "cp", "--", fin, tmp], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if os.path.getsize(tmp) > 0:
            return tmp
    except Exception:
        pass
    try:
        os.remove(tmp)
    except OSError:
        pass
    return None


class Output:
    """Extendable arrays in the HDF5 file, created on first use."""

    def __init__(self, fout):
        self.fout = fout

    def append(self, name, obj):
        obj = np.asarray(obj)
        if name in self.fout.root:
            getattr(self.fout.root, name).append(obj)
        else:
            self.fout.create_earray(self.fout.root, name, obj=obj)

    def append_strings(self, name, values):
        if name not in self.fout.root:
            self.fout.create_earray(self.fout.root, name, atom=tb.StringAtom(itemsize=FNAME_ITEMSIZE), shape=(0,))
        getattr(self.fout.root, name).append(np.array(values, dtype=f"S{FNAME_ITEMSIZE}"))


def new_event_prongs():
    return dict(pad_mask=np.zeros(MAX_PRONGS, dtype=np.int8),
                label=-np.ones(MAX_PRONGS, dtype=np.int8),
                label_split_photons=-np.ones(MAX_PRONGS, dtype=np.int8),
                mother=-np.ones(MAX_PRONGS, dtype=np.int16),
                ShwTrk=-np.ones(MAX_PRONGS, dtype=np.int8),
                E=-np.ones(MAX_PRONGS, dtype=np.float32),
                P=-np.ones((3, MAX_PRONGS), dtype=np.float32),
                input=np.zeros((8, MAX_PRONGS), dtype=np.float32))


def process_file(df, finh, base, model, total_events, out):
    """Convert one pixelmap tree; returns the number of events written (0 if skipped)."""
    evt_ids = df['unique_id']
    events = np.unique(evt_ids)
    prong_tag_orig = df['prong_tag']
    png_tags = prong_tag_orig.copy()
    iplane, rel_wire, rel_tick, wire_charge = df['iplane'], df['rel_wire'], df['rel_tick'], df['wire_charge']

    # --- locate the first pixel of every event and of every prong image -----------------
    # prong_tag: -1 = event image, 0.. = track i, 1000+i = shower i (showers come first in
    # the tree); tracks are renumbered to follow the showers: tag -> tag + last_shower + 1.
    evt_start_indices = [0]
    png_start_indices = []
    last_shw = -1
    for ipix, tag in enumerate(png_tags):
        if tag == -1:
            last_shw = -1
            if ipix != 0:
                if evt_ids[ipix] != evt_ids[ipix - 1]:
                    evt_start_indices.append(ipix)
                if (png_tags[ipix - 1] == -1 or png_tags[ipix - 1] == 1000) and evt_ids[ipix] != evt_ids[ipix - 1]:
                    png_start_indices.append(ipix)
        if tag > -1:
            if tag >= 1000:
                continue
            png_tags[ipix] = tag + last_shw + 1
            if ipix != 0 and png_tags[ipix] != png_tags[ipix - 1]:
                png_start_indices.append(ipix)
    evt_start_indices = np.array(evt_start_indices)
    png_start_indices = np.array(png_start_indices, dtype=np.int64)
    if len(events) != len(evt_start_indices):
        print(f"[skip] {base}: {len(events)} events but {len(evt_start_indices)} event starts", flush=True)
        return 0

    # --- per-prong summary arrays (one entry per event, padded to MAX_PRONGS) ------------
    png_start_tags = png_tags[png_start_indices]
    pr = {k: df.get(k)[png_start_indices] for k in
          ['prong_true_pdg', 'prong_true_pdg_mom', 'prong_tag', 'prong_eng', 'prong_length',
           'prong_startx', 'prong_starty', 'prong_startz', 'prong_px', 'prong_py', 'prong_pz',
           'prong_true_eng', 'prong_true_px', 'prong_true_py', 'prong_true_pz']}
    per_event = []      # list of dicts, one per event
    cur = new_event_prongs()
    for ipng, png in enumerate(png_start_tags):
        if png >= MAX_PRONGS:
            continue
        if ipng != 0:
            if png == -1 and png_start_tags[ipng - 1] == -1:
                per_event.append(cur)
                continue
            if png == -1:
                per_event.append(cur)
                cur = new_event_prongs()
                per_event.append(cur)
                continue
            if png <= png_start_tags[ipng - 1]:
                per_event.append(cur)
                cur = new_event_prongs()
        cur['pad_mask'][png] = 1
        label = PDG_LABELS.get(int(abs(pr['prong_true_pdg'][ipng])))
        if label is None:
            label = 7
            if pr['prong_true_pdg'][ipng] != -1:
                print("Unknown pdg: %i" % pr['prong_true_pdg'][ipng], flush=True)
        cur['label'][png] = label
        cur['mother'][png] = pr['prong_true_pdg_mom'][ipng]
        cur['E'][png] = pr['prong_true_eng'][ipng]
        cur['P'][:, png] = [pr['prong_true_px'][ipng], pr['prong_true_py'][ipng], pr['prong_true_pz'][ipng]]
        length = pr['prong_length'][ipng]
        if length == -1:
            length = 0
        cur['input'][:, png] = [pr['prong_eng'][ipng], length, pr['prong_startx'][ipng], pr['prong_starty'][ipng],
                                pr['prong_startz'][ipng], pr['prong_px'][ipng], pr['prong_py'][ipng], pr['prong_pz'][ipng]]
        if pr['prong_true_pdg'][ipng] == 22 and pr['prong_true_pdg_mom'][ipng] == 111:
            cur['label_split_photons'][png] = 5
        elif pr['prong_true_pdg'][ipng] == 22 and pr['prong_true_pdg_mom'][ipng] == 2112:
            cur['label_split_photons'][png] = 3
        else:
            cur['label_split_photons'][png] = cur['label'][png]
        if pr['prong_tag'][ipng] >= 1000:
            cur['ShwTrk'][png] = 0
        elif pr['prong_tag'][ipng] >= 0:
            cur['ShwTrk'][png] = 1
    per_event.append(cur)
    if len(per_event) != len(evt_start_indices):
        print(f"[skip] {base}: {len(per_event)} prong records but {len(evt_start_indices)} events", flush=True)
        return 0

    # --- event images -------------------------------------------------------------------
    evt_image_indices = np.where(prong_tag_orig == -1)[0]
    out.append("cvnmap_index", np.column_stack((evt_ids[evt_image_indices] - 1 + total_events,
                                                iplane[evt_image_indices], rel_wire[evt_image_indices],
                                                rel_tick[evt_image_indices])))
    out.append("cvnmap_value", wire_charge[evt_image_indices])
    out.append("input_png3d_pad_mask", [e['pad_mask'] for e in per_event])
    out.append("mc.png_label", [e['label'] for e in per_event])
    out.append("mc.png_label_split_photons", [e['label_split_photons'] for e in per_event])
    out.append("mc.png_mother", [e['mother'] for e in per_event])
    out.append("png_ShwTrk", [e['ShwTrk'] for e in per_event])
    out.append("png_trueE", [e['E'] for e in per_event])
    out.append("png_trueP", [e['P'] for e in per_event])
    out.append("input_png3d", [e['input'] for e in per_event])
    out.append_strings("file_name", [base] * len(evt_start_indices))

    # --- prong images (renumbered tags, at most MAX_PRONGS per event) ---------------------
    png_image_indices = np.where(prong_tag_orig != -1)[0]
    tags = png_tags[png_image_indices]
    keep = tags < MAX_PRONGS
    tags, png_image_indices = tags[keep], png_image_indices[keep]
    out.append("png_cvnmap_index", np.column_stack((evt_ids[png_image_indices] - 1 + total_events, tags,
                                                    iplane[png_image_indices], rel_wire[png_image_indices],
                                                    rel_tick[png_image_indices])))
    out.append("png_cvnmap_value", wire_charge[png_image_indices])

    # --- per-event scalars ------------------------------------------------------------------
    s = evt_start_indices
    ccnc = df.get('trueCCNC')[s]
    pdg = df.get('trueNuPDG')[s]
    inter = np.where(ccnc == 1, 13, pdg)
    for raw, code in ((14, 0), (-14, 1), (12, 4), (-12, 5), (16, 8), (-16, 9)):
        inter = np.where(inter == raw, code, inter)
    out.append("mc.inter", inter)
    out.append("model", np.full(len(s), model, dtype=np.int8))
    out.append("event_id", np.column_stack((df['run'][s], df['subrun'][s], df['ievent'][s])).astype(np.int32))
    genie = genie_truth(finh, evt_ids[s])
    if genie is not None:
        for k, v in genie.items():
            out.append(k, v)
    out.append("tpc_id", df['tpc_id'][s])
    out.append("trueE", df.get('trueEnergy', -1.0)[s])
    out.append("precut", np.column_stack((df['n_tracks'][s], df['NShw'][s], df['hit_tot_energy'][s])))
    out.append("trueVertex", np.column_stack((df.get('trueVertex_X', -1.0)[s], df.get('trueVertex_Y', -1.0)[s],
                                              df.get('trueVertex_Z', -1.0)[s])))
    out.append("trueP", np.column_stack((df.get('truePx', -1.0)[s], df.get('truePy', -1.0)[s], df.get('truePz', -1.0)[s])))
    out.append("input_slice", np.column_stack((df['recoNu_E'][s], df['recoVertex_X'][s], df['recoVertex_Y'][s],
                                               df['recoVertex_Z'][s])))
    out.append("isNonFlux", np.full(len(events), 1 if 'nue' in base else 0, dtype=np.int8))
    for var in CVN_VARS:
        out.append(var, df.get(var, -1.0)[s])
    return len(events)


def main(fout_path, fins, model=-1, shuffle=False, stage=True):
    if shuffle:
        import random
        random.shuffle(fins)
    fout = tb.open_file(fout_path, mode='w')
    out = Output(fout)
    total_events, nfiles, missing = 0, 0, set()
    try:
        for fin in fins:
            local = copy_local(fin) if stage else fin
            if local is None:
                print(f"[skip] copy timed out/bad: {fin}", flush=True)
                continue
            try:
                with uproot.open(local) as finh:
                    if 'pixelmap' not in finh:
                        print(f"[skip] no pixelmap tree: {fin}", flush=True)
                        continue
                    df = Tree(finh['pixelmap'])
                    if df.tree.num_entries == 0:
                        print(f"[skip] empty pixelmap tree: {fin}", flush=True)
                        continue
                    n = process_file(df, finh, os.path.basename(fin).lower(), model, total_events, out)
                    missing |= df.missing
                total_events += n
                nfiles += 1 if n else 0
                print(f"{nfiles} files, {total_events} events", flush=True)
            finally:
                if stage and local != fin:
                    try:
                        os.remove(local)
                    except OSError:
                        pass
        out.append("cvnmap_shape", np.array([total_events, *IMAGE]))
        out.append("png_cvnmap_shape", np.array([total_events, MAX_PRONGS, *IMAGE]))
    finally:
        fout.close()
    if missing:
        print("[info] branches absent from the input, filled with -1: " + ", ".join(sorted(missing)), flush=True)
    return total_events


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", help="output HDF5 file (overwritten)")
    ap.add_argument("files", nargs="*", help="pixelmap ROOT files")
    ap.add_argument("--filelist", help="text file with one pixelmap ROOT file per line")
    ap.add_argument("--model-tag", default=None,
                    help="sample tag stored in 'model' (%s -> 0..%d); default -1" % (", ".join(MODEL_TAGS), len(MODEL_TAGS) - 1))
    ap.add_argument("--shuffle", action="store_true", help="shuffle the input files (as done for training)")
    ap.add_argument("--no-stage", action="store_true", help="read inputs in place instead of copying them locally first")
    a = ap.parse_args()
    fins = list(a.files)
    if a.filelist:
        with open(a.filelist) as f:
            fins += [ln.strip() for ln in f if ln.strip()]
    if not fins:
        ap.error("no input files")
    model = -1
    if a.model_tag is not None:
        if a.model_tag not in MODEL_TAGS:
            ap.error(f"unknown --model-tag {a.model_tag}; choose from {MODEL_TAGS}")
        model = MODEL_TAGS.index(a.model_tag)
    print(f"{len(fins)} files -> {a.out}", flush=True)
    n = main(a.out, fins, model=model, shuffle=a.shuffle, stage=not a.no_stage)
    print(f"wrote {a.out}: {n} events", flush=True)
    sys.exit(0 if n > 0 else 1)
