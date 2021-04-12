"""
Plot images that have been infilled by trained models.
"""


import os, sys, time, argparse, itertools
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D

import torch

from infill_loss import InfillLossInduction, InfillLossCollection
from model import UnetInduction, UnetCollection
# from loss_dense_infill_collection import DenseInfillLoss


def predict(model, test_dir, N, trace, info):
    test_masked_lst, test_true_lst, filenames = [], [], [],
    maskpattern_pool = itertools.cycle(info["maskpatterns"])

    for filename in os.listdir(test_dir):
        if filename.endswith(".npy"):
            arr = np.load(os.path.join(test_dir, filename)).T

            maskpattern = next(maskpattern_pool)
            arr_mask = np.copy(arr)
            arr_mask[:, maskpattern] = 0

            test_masked_lst.append(torch.FloatTensor(arr_mask.reshape(1, *arr_mask.shape)))
            test_true_lst.append(torch.FloatTensor(arr.reshape(1, *arr.shape)))
            filenames.append(filename)

        if len(test_masked_lst) >= N:
            break 

    model.eval()
    with torch.no_grad():
        for idx, masked in enumerate(test_masked_lst):
            print(filenames[idx])
            masked_tensor = torch.stack([masked])
            true_tensor = torch.stack([test_true_lst[idx]])

            masked_tensor.to(info["DEVICE"])
            true_tensor.to(info["DEVICE"])

            start = time.time()
            outputs = model(masked_tensor)
            end = time.time()
            print("Inference time:{:.1f}s".format(end - start))

            loss = info["criterion"](outputs, true_tensor, masked_tensor)[5]
            print("Loss: {}".format(loss))

            img_pred = outputs.detach().numpy()[0, 0, :, :].T
            img_true = true_tensor.detach().numpy()[0, 0, :, :].T
            img_masked = masked_tensor.detach().numpy()[0, 0, :, :].T

            dead_ch = [ idx for idx, col in enumerate(img_masked) if np.all(col == 0) ]
            print("Dead channels: {}".format(dead_ch))

            mask = np.zeros_like(img_masked)
            mask[np.array(dead_ch), :] = 1
            img_masked = np.ma.masked_array(img_masked, mask)
            img_infill = np.ma.masked_array(img_pred, np.logical_not(mask))

            # Plots truth image
            # fig, ax = plt.subplots()
            # fig.set_size_inches(16, 10)
            # im1 = ax.imshow(img_true.T, aspect='auto', cmap='coolwarm', vmin=-30, vmax=30, interpolation='None')
            # plt.show()

            plt.rc('font', family='serif')
            fig, ax = plt.subplots()
            fig.set_size_inches(16, 10)
            im1 = ax.imshow(img_masked.T, aspect='auto', cmap='coolwarm', vmin=-30, vmax=30, interpolation='None')
            im2 = ax.imshow(img_infill.T, aspect='auto', cmap='PRGn', vmin=-30, vmax=30, interpolation='None')
            plt.show()

            if trace:
                for ch in dead_ch:
                    print("Channel {}".format(ch))
                    tick_adc_true = img_true[ch, :]
                    tick_adc_pred = img_pred[ch, :]
                    tick = np.arange(1, 6001) 

                    plt.hist(tick, bins=len(tick), weights=tick_adc_true, histtype='step', label="True", linewidth=0.7)
                    plt.hist(tick, bins=len(tick), weights=tick_adc_pred, histtype='step', label="Network", linewidth=0.7)
                    plt.xlim(1,6001)
                    plt.xlabel("time tick", fontsize=20)
                    plt.ylabel("ADC", fontsize=20)
                    plt.title("Channel {}".format(ch))
                    ax = plt.gca()
                    ax.tick_params(axis='both', which='major', labelsize=16)
                    ax.tick_params(axis='both', which='minor', labelsize=16)
                    # What is this stuff doing?
                    tx = ax.yaxis.get_offset_text()
                    tx.set_fontsize(20)
                    handles, labels = ax.get_legend_handles_labels()
                    new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
                    plt.legend(handles=new_handles, labels=labels, prop={'size': 20})
                    plt.show()



def main(weights, test_dir, collection, induction, n, traces):
    DEVICE = torch.device("cpu")

    if collection:
        model = UnetCollection()
        criterion = InfillLossCollection().to(device=DEVICE)
        maskpatterns = [[55, 66, 78, 81, 89, 238, 370],
            [217, 219, 221, 223, 225, 227, 250, 251, 252, 361, 363, 365, 367, 369, 371],
            [20, 95, 134, 147, 196, 351],
            [2, 3, 25, 27, 29, 31, 33, 35, 289, 409, 411, 413, 415, 417, 419, 456],
            [4, 13, 424, 436]]

    elif induction:
        model = UnetInduction()
        criterion = InfillLossInduction().to(device=DEVICE)
        # Change these patters to the ROP ones from decodeDigits
        maskpatterns = [[1, 2, 4, 94, 200, 202, 204, 206, 208, 325, 400, 401, 442, 447, 453, 455,
             456, 472, 477, 571, 573],
            [0, 1, 76, 191, 193, 195, 197, 199, 400, 734, 739, 746],
            [114, 273, 401],
            [181, 183, 301, 303, 701, 703, 781, 783],
            [5, 151, 201, 241, 243, 257, 280, 303],
            [212],
            [0, 1, 238, 400, 648, 661],
            [0, 21, 23, 341, 343, 781, 783],
            [457, 560, 667, 784],
            [163, 230, 417, 419, 423, 429, 477, 629, 639],
            [1, 201, 281, 563]]
        
    model = model.to(DEVICE)

    # Needed to infer on cpu if model was trained using DataParallel
    pretrained_dict = torch.load(weights, map_location="cpu")
    pretrained_dict = {key.replace("module.", ""): value for key, value in pretrained_dict.items()}
    model.load_state_dict(pretrained_dict)

    info = {
            "DEVICE" : DEVICE,
            "criterion" : criterion,
            "maskpatterns" : maskpatterns
           }

    predict(model, test_dir, n, traces, info)


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("weights", help="Takes a .pth file")
    parser.add_argument("test_dir")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--collection",action="store_true")
    group.add_argument("--induction",action="store_true")

    parser.add_argument("-n", nargs="?", type=int, action="store", default=5, dest="N", 
        help="Number of files to do inference on")
    parser.add_argument("-t", "--traces", action="store_true",
    help="Plot true and predicted ADC by time tick for each dead channel that has been infilled")

    args = parser.parse_args()

    return (args.weights, args.test_dir, args.collection, args.induction, args.N, args.traces)


if __name__ == "__main__":
    arguments = parse_arguments()
    main(*arguments)








