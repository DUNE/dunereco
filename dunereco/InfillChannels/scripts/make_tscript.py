"""
Trace trained infill model to produce a TorchScript. Resulting .pt can then be loaded direclty into
C++.
"""

import os, argparse

import torch
import numpy as np

from model import UnetInduction, UnetCollection

def main(input_file, out_name, example_dir, collection):

    if collection:
        model = UnetCollection()

    else:
        model = UnetInduction()

    DEVICE = torch.device("cpu")
    model = model.to(DEVICE)
    model.eval()
    pretrained_dict = torch.load(input_file, map_location="cpu")
    pretrained_dict = {key.replace("module.", ""): value for key, value in pretrained_dict.items()}
    model.load_state_dict(pretrained_dict)

    for filename in os.listdir(example_dir):
        if filename.endswith(".npy"):
            arr = np.load(os.path.join(example_dir, filename)).T
            maskpattern = [114, 273, 401] # Example maskpattern
            arr[:, maskpattern] = 0
            example_img = torch.FloatTensor(arr.reshape(1, *arr.shape))
            example_img = torch.stack([example_img])
            break
    
    with torch.no_grad():  
        traced_model = torch.jit.trace(model, example_img)
    traced_model.save(out_name + ".pt")
    print("TorchScript summary:\n{}".format(traced_model))


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("input_file", help=".pth file to be serialised")
    parser.add_argument("output_name")
    parser.add_argument("example_dir", help="Directory containing input data for the model being serialised")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--collection",action='store_true')
    group.add_argument("--induction",action='store_true')

    args = parser.parse_args()

    return (args.input_file, args.output_name, args.example_dir, args.collection)


if __name__ == "__main__":
    arguments = parse_arguments()
    main(*arguments)