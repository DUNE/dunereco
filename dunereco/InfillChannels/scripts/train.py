"""
Train the infill networks.
"""

import os, datetime, random, argparse, itertools
import numpy as np
import yaml

import torch
import torch.nn as nn

from infill_loss import InfillLossInduction, InfillLossCollection
from model import UnetInduction, UnetCollection


def train(model, train_files, valid_files, maskpatterns, epochs, batchsize, info):
    overtain_cntr = 0
    train_losses, valid_losses = [], []
    summary = {}
    now = datetime.datetime.now().strftime("%d%m%Y-%H%M%S")

    batchsizes_train = [batchsize]*(int((len(train_files)/batchsize)))
    batchsizes_train.append(len(train_files) % batchsize)
    batchsizes_valid = [batchsize]*(int((len(valid_files)/batchsize)))
    batchsizes_valid.append(len(valid_files) % batchsize)
    if batchsizes_train[-1] == 0:
        batchsizes_train.pop()
    if batchsizes_valid[-1] == 0:
        batchsizes_valid.pop()

    for epoch in range(epochs):
        model.train()
        
        epoch_running_train_loss = 0.0
        running_loss = 0.0

        random.shuffle(train_files)
        files_for_batches = np.split(np.array(train_files), [ sum(batchsizes_train[:i]) for i in range(1, len(batchsizes_train)) ])

        for idx, batch_files in enumerate(files_for_batches):
            masked_tensor_batch_lst, true_tensor_batch_lst = [], []

            for batch_idx, filepath in enumerate(batch_files):
                arr = np.load(filepath).T

                maskpattern = random.sample(maskpatterns, 1)[0]
                offset = random.randint(1,info["width"] - 1) # Exclude offset = 0 for validation set.
                offset_maskpattern = [ i + offset if (i + offset) < info["width"] else i - info["width"] + offset for i in maskpattern ]
                arr_mask = np.copy(arr)
                arr_mask[:, offset_maskpattern] = 0

                masked_tensor_batch_lst.append(torch.FloatTensor(arr_mask.reshape(1, *arr_mask.shape)))
                true_tensor_batch_lst.append(torch.FloatTensor(arr.reshape(1, *arr.shape)))

            masked_tensor_batch = torch.stack(masked_tensor_batch_lst)
            del masked_tensor_batch_lst
            true_tensor_batch = torch.stack(true_tensor_batch_lst)
            del true_tensor_batch_lst

            masked_tensor_batch = masked_tensor_batch.to(info["DEVICE"])
            true_tensor_batch = true_tensor_batch.to(info["DEVICE"])

            info["optimizer"].zero_grad()
            outputs = model(masked_tensor_batch)
            loss = info["criterion"](outputs, true_tensor_batch, masked_tensor_batch)[5]
            loss.backward()
            info["optimizer"].step()
   
            del masked_tensor_batch
            del true_tensor_batch

            running_loss += loss.item()
            epoch_running_train_loss += loss.item()
            if (idx + 1) % 5 == 0:
                print('[{}, {:2.2%}] loss: {:.2f}'.format(epoch + 1, (idx*batchsize)/float(len(train_files)), running_loss/5))
                running_loss = 0.0

#        adjust_learning_rate(info["optimizer"], epoch, info["lr"]) # lr decay
        
        train_losses.append(epoch_running_train_loss/len(files_for_batches))
        
        model.eval()

        running_loss = 0.0

        files_for_batches = np.split(np.array(valid_files), [ sum(batchsizes_valid[:i]) for i in range(1, len(batchsizes_valid)) ])
        maskpattern_pool = itertools.cycle(maskpatterns)

        with torch.no_grad():
            for batch_files in files_for_batches:
                masked_tensor_batch_lst, true_tensor_batch_lst = [], []

                for batch_idx, filepath in enumerate(batch_files):
                    arr = np.load(filepath).T

                    maskpattern = next(maskpattern_pool) # Use true mask patterns for validation
                    arr_mask = np.copy(arr)
                    arr_mask[:, maskpattern] = 0

                    masked_tensor_batch_lst.append(torch.FloatTensor(arr_mask.reshape(1, *arr_mask.shape)))
                    true_tensor_batch_lst.append(torch.FloatTensor(arr.reshape(1, *arr.shape)))

                masked_tensor_batch = torch.stack(masked_tensor_batch_lst)
                del masked_tensor_batch_lst
                true_tensor_batch = torch.stack(true_tensor_batch_lst)
                del true_tensor_batch_lst

                masked_tensor_batch = masked_tensor_batch.to(info["DEVICE"])
                true_tensor_batch = true_tensor_batch.to(info["DEVICE"])

                outputs = model(masked_tensor_batch)
                loss = info["criterion"](outputs, true_tensor_batch, masked_tensor_batch)[5]

                del masked_tensor_batch
                del true_tensor_batch

                running_loss += loss.item()

        valid_losses.append(running_loss/len(files_for_batches))
        print("Validation loss: {:.2f}".format(running_loss/len(files_for_batches)))

        summary['train losses'] = train_losses
        summary['valid losses'] = valid_losses

        if epoch == 0:
            torch.save(model.module.state_dict(), info["model_name"] + '.pth')
            old_valid_loss = valid_losses[0]

        else:
            if (valid_losses[-1] - old_valid_loss) < 0:
                torch.save(model.module.state_dict(), info["model_name"] + '.pth')
                old_valid_loss = valid_losses[-1]
                overtain_cntr = 0 
                summary['best epoch'] = epoch
                summary['best valid loss'] = valid_losses[-1]

            else:
                overtain_cntr += 1

        if overtain_cntr > 5:
            break

        with open('training_summary_{}.yaml'.format(now), 'w') as f:
             yaml.dump(summary, f)

    print("best valid loss: {} (at epoch {})".format(summary['best valid loss'], summary['best epoch']))
    print("train losees: {}\n".format(train_losses))
    print("valid losses: {}\n".format(valid_losses))


def adjust_learning_rate(optimizer, epoch, lr):
    lr = lr * (0.5 ** (epoch // 20)) 
    for param_group in optimizer.param_groups:
        param_group['lr'] = lr


def main(train_dir, valid_dir, collection, induction, epochs, batchsize, model_name):
    DEVICE = torch.device("cuda:0")

    if collection:
        model = UnetCollection()
        criterion = InfillLossCollection().to(device=DEVICE)
        maskpatterns = [[55, 66, 78, 81, 89, 238, 370],
            [217, 219, 221, 223, 225, 227, 250, 251, 252, 361, 363, 365, 367, 369, 371],
            [20, 95, 134, 147, 196, 351],
            [2, 3, 25, 27, 29, 31, 33, 35, 289, 409, 411, 413, 415, 417, 419, 456],
            [4, 13, 424, 436]]
        width = 480

    elif induction:
        model = UnetInduction()
        criterion = InfillLossInduction().to(device=DEVICE)
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
        width = 800

    model = nn.DataParallel(model)
    model.to(DEVICE)

    lr = 1.0e-4
    momentum = 0.9
    weight_decay = 1.0e-4
    optimizer = torch.optim.RMSprop(model.parameters(), lr=lr, weight_decay=weight_decay)

    train_files = [ os.path.join(train_dir, filename) for filename in os.listdir(train_dir) if filename.endswith(".npy") ]
    valid_files = [ os.path.join(valid_dir, filename) for filename in os.listdir(valid_dir) if filename.endswith(".npy") ]

    info = { 
             "DEVICE" : DEVICE,
             "criterion" : criterion,
             "optimizer" : optimizer,
             "lr" : lr,
             "model_name" : model_name,
             "width" : width
           }

    train(model, train_files, valid_files, maskpatterns, epochs, batchsize, info)


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("train_dir")
    parser.add_argument("valid_dir")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--collection",action='store_true')
    group.add_argument("--induction",action='store_true')

    parser.add_argument("-e", "--epochs", nargs='?', type=int, default=10, action='store', dest='EPOCHS')
    parser.add_argument("-b", "--batchsize", nargs='?', type=int, default=12, action='store', dest='BATCHSIZE')
    parser.add_argument("--model_name", nargs='?', type=str, action='store', dest='MODEL_NAME',
                        default="{}".format(datetime.datetime.now().strftime("_%d%m%Y-%H%M%S")))

    args = parser.parse_args()

    return (args.train_dir, args.valid_dir, args.collection, args.induction, args.EPOCHS,
            args.BATCHSIZE, args.MODEL_NAME)


if __name__ == "__main__":
    arguments = parse_arguments()
    main(*arguments)
