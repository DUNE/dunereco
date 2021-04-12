"""
Prepare training data in numpy format from TH2s.
"""

import os, argparse
import uproot
import numpy as np
from matplotlib import pyplot as plt


def root_to_numpy(filepath, induct_saveloc, collect_saveloc):
    file = uproot.open(filepath)

    for idx, key in enumerate(file['infilldata'].keys()):
        z, _ = file['infilldata'][key].numpy()

        # print(z.dtype)
        # plt.imshow(z.T, aspect='auto', vmin=-20, vmax=20, cmap='coolwarm')
        # plt.show()

        if ('ROP0' in key) or ('ROP1' in key): 
            with open(os.path.join(induct_saveloc, "{}.npy".format(key[:-2])), "w") as f:
                np.save(f, z)

        elif ('ROP2' in key) or ('ROP3' in key): 
            with open(os.path.join(collect_saveloc, "{}.npy".format(key[:-2])), "w") as f:
                np.save(f, z)

        if (idx + 1) % 50 == 0:
            print("{}/{}".format(idx + 1, len(file["dec"].keys())))


def main(input_file, output_dir):
    
    induct_saveloc = os.path.join(output_dir, "induction")
    if not os.path.exists(induct_saveloc):
        os.makedirs(induct_saveloc)

    collect_saveloc = os.path.join(output_dir, "collection")
    if not os.path.exists(collect_saveloc):
        os.makedirs(collect_saveloc)

    root_to_numpy(input_file, induct_saveloc, collect_saveloc)


def parse_arguments():
    
    parser = argparse.ArgumentParser()

    parser.add_argument("input_file")
    parser.add_argument("output_dir")

    args = parser.parse_args()

    return (args.input_file, args.output_dir)


if __name__ == "__main__":
    arguments = parse_arguments()
    main(*arguments)
