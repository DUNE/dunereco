#!/usr/bin/env python3
"""Export XGBoost model JSON to compact binary format (.bin) for FlashMatchBDT in DUNE.

Binary format specification (BDT1):
  Header (24 bytes):
    char magic[4]      = "BDT1"
    uint32_t version   = 1
    uint32_t num_trees
    uint32_t num_features
    uint32_t total_nodes
    float    base_margin
  Tree Table (num_trees entries, 6 bytes each):
    uint32_t node_offset
    uint16_t num_nodes
  Node Table (total_nodes entries, 10 bytes each):
    int16_t  feat       (feature index 0..12, or -1 if leaf)
    float    threshold  (split value if internal, weight if leaf)
    int16_t  left       (left child relative index, or -1 if leaf)
    int16_t  right      (right child relative index, or -1 if leaf)

Author: Pierre Granger (pierre.granger@cern.ch)
Date: 2026
"""

import os
import sys
import json
import struct
import argparse
import numpy as np


def convert_json_to_bin(json_path: str, bin_path: str, num_features: int = 13):
    print(f"Loading {json_path} ...")
    with open(json_path, "r") as f:
        data = json.load(f)

    trees = data["learner"]["gradient_booster"]["model"]["trees"]
    base_score = float(data["learner"]["learner_model_param"]["base_score"])
    base_margin = float(np.log(base_score / (1.0 - base_score)))

    num_trees = len(trees)
    total_nodes = sum(len(t["split_indices"]) for t in trees)

    print(f"  Trees: {num_trees}, Total nodes: {total_nodes}")
    print(f"  Base score: {base_score:.8f}, Base margin: {base_margin:.8f}")

    os.makedirs(os.path.dirname(os.path.abspath(bin_path)), exist_ok=True)
    with open(bin_path, "wb") as f_out:
        # 1. Header (24 bytes)
        f_out.write(b"BDT1")
        f_out.write(struct.pack("<IIIIf", 1, num_trees, num_features, total_nodes, base_margin))

        # 2. Tree table
        curr_offset = 0
        for t in trees:
            n_nodes = len(t["split_indices"])
            f_out.write(struct.pack("<IH", curr_offset, n_nodes))
            curr_offset += n_nodes

        # 3. Node table
        for t in trees:
            n_nodes = len(t["split_indices"])
            splits = t["split_indices"]
            conds = t["split_conditions"]
            lefts = t["left_children"]
            rights = t["right_children"]
            for i in range(n_nodes):
                feat = int(splits[i])
                thresh = float(conds[i])
                left = int(lefts[i])
                right = int(rights[i])
                f_out.write(struct.pack("<hfh h", feat, thresh, left, right))

    file_size_kb = os.path.getsize(bin_path) / 1024.0
    print(f"Successfully generated {bin_path} ({file_size_kb:.1f} KB)\n")


def main():
    parser = argparse.ArgumentParser(description="Convert XGBoost JSON model to FlashMatchBDT binary (.bin) format.")
    parser.add_argument("input_json", nargs="?", default=None, help="Path to input XGBoost JSON model")
    parser.add_argument("output_bin", nargs="?", default=None, help="Path to output binary file")
    parser.add_argument("--features", type=int, default=13, help="Number of model features (default: 13)")
    args = parser.parse_args()

    if args.input_json and args.output_bin:
        convert_json_to_bin(args.input_json, args.output_bin, args.features)
    else:
        # Default conversion for standard project paths
        script_dir = os.path.dirname(os.path.abspath(__file__))
        flashmatch_dir = os.path.dirname(script_dir)
        models = [
            ("/home/pigrange/projects/analysis_lowE/out/test_retrained_bdt/bdt_augmented.json",
             os.path.join(flashmatch_dir, "bdt_augmented.bin")),
            ("/home/pigrange/projects/analysis_lowE/out/test_retrained_bdt/bdt_baseline.json",
             os.path.join(flashmatch_dir, "bdt_baseline.bin")),
        ]
        for json_p, bin_p in models:
            if os.path.exists(json_p):
                convert_json_to_bin(json_p, bin_p, args.features)
            else:
                print(f"Notice: Default model {json_p} not found. Specify input and output paths.")


if __name__ == "__main__":
    main()
