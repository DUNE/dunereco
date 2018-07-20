"""
This is the dataset generator module.
"""

__version__ = '1.0'
__author__ = 'Saul Alonso-Monsalve'
__email__ = "saul.alonso.monsalve@cern.ch"

import numpy as np
import glob
import ast
import ntpath
import pickle
import configparser
import logging
import sys
import time
import random

from sklearn.utils import class_weight
from collections import Counter

'''
****************************************
************** PARAMETERS **************
****************************************
'''

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

config = configparser.ConfigParser()
config.read('config/config.ini')

# random

SEED = int(config['random']['seed'])

if SEED == -1:
    SEED = int(time.time())

np.random.seed(SEED)

# images

IMAGES_PATH = config['images']['path']

# dataset

DATASET_PATH = config['dataset']['path']
PARTITION_PREFIX = config['dataset']['partition_prefix']
LABELS_PREFIX = config['dataset']['labels_prefix']
UNIFORM = ast.literal_eval(config['dataset']['uniform'])
 
# train

TRAIN_FRACTION = float(config['train']['fraction'])
WEIGHTED_LOSS_FUNCTION = ast.literal_eval(config['train']['weighted_loss_function'])
CLASS_WEIGHTS_PREFIX = config['train']['class_weights_prefix']

# validation

VALIDATION_FRACTION = float(config['validation']['fraction'])

# test

TEST_FRACTION = float(config['test']['fraction'])

if((TRAIN_FRACTION + VALIDATION_FRACTION + TEST_FRACTION) > 1):
    logging.error('(TRAIN_FRACTION + VALIDATION_FRACTION + TEST_FRACTION) must be <= 1')
    exit(-1)

def normalize(value):
    if value > 2:
        return 3
    return value 

def normalize2(value):
    if value < 0:
        return 1
    return 0 

'''
****************************************
*************** DATASETS ***************
****************************************
'''

partition = {'train' : [], 'validation' : [], 'test' : []} # Train, validation, and test IDs
labels = {}                                                # ID : label
y_train = []
y1_class_weights = []
y2_class_weights = []

if UNIFORM:

    pass

# Iterate through label folders

logging.info('Filling datasets...')

for images_path in glob.iglob(IMAGES_PATH + '/*'):

    count_train, count_val, count_test = (0, 0, 0)

    print images_path
    files = list(glob.iglob(images_path + "/images/*"))
    random.shuffle(files)

    for imagefile in files:
        #print imagefile
        ID = imagefile.split("/")[-1][:-3] 
        infofile = images_path + '/info/' + ID + '.info'

        #print infofile

        info = open(infofile, 'r').readlines()
        fInt = int(info[0].strip())
        flavour = fInt // 4
        interaction = fInt % 4

        fNuEnergy = float(info[1].strip())
        fLepEnergy = float(info[2].strip())
        fRecoNueEnergy = float(info[3].strip())
        fRecoNumuEnergy = float(info[4].strip())
        fEventWeight = float(info[5].strip())

        fNuPDG = normalize2(int(info[6].strip()))
        fNProton = normalize(int(info[7].strip()))
        fNPion = normalize(int(info[8].strip()))
        fNPizero = normalize(int(info[9].strip()))
        fNNeutron = normalize(int(info[10].strip()))

        if fInt == 13:
            fNuPDG = -1
            flavour = 3
            interaction = -1

        '''
        print "Info:", info
        print "fInt:", fInt
        print "flavour:", flavour
        print "interaction:", interaction
        print "fNuEnergy:", fNuEnergy
        print "fLepEnergy:", fLepEnergy
        print "fRecoNueEnergy:", fRecoNueEnergy
        print "fRecoNumuEnergy:", fRecoNumuEnergy
        print "fEventWeight:", fEventWeight
        print "fNuPDG:", fNuPDG
        print "fNProton:", fNProton
        print "fNPion:", fNPion
        print "fNPizero:", fNPizero
        print "fNNeutron:", fNNeutron
        '''

        random_value = np.random.uniform(0,1)

        # Fill training set

        if(random_value < TRAIN_FRACTION):
            random_value_uni = np.random.uniform(0,1)
            if not UNIFORM or random_value_uni < fraction:
                partition['train'].append(ID)
                labels[ID] = [fNuPDG, flavour, interaction, fNProton, fNPion, fNPizero, fNNeutron]

                count_train += 1

        # Fill validation set

        elif(random_value < (TRAIN_FRACTION + VALIDATION_FRACTION)):
            partition['validation'].append(ID)
            labels[ID] = [fNuPDG, flavour, interaction, fNProton, fNPion, fNPizero, fNNeutron]

            count_val += 1

        # Fill test set

        elif(random_value < (TRAIN_FRACTION + VALIDATION_FRACTION + TEST_FRACTION)):
            partition['test'].append(ID)
            labels[ID] = [fNuPDG, flavour, interaction, fNProton, fNPion, fNPizero, fNNeutron]

            count_test += 1

    logging.debug('%d train images', count_train)
    logging.debug('%d val images', count_val)
    logging.debug('%d test images', count_test)
    logging.debug('%d total images', count_train + count_val + count_test)

# Print dataset statistics

logging.info('Number of training examples: %d', len(partition['train']))
logging.info('Number of validation examples: %d', len(partition['validation']))
logging.info('Number of test examples: %d', len(partition['test']))

# Serialize partition and labels

logging.info('Serializing datasets...')

with open(DATASET_PATH + PARTITION_PREFIX + '.p', 'w') as partition_file:
    pickle.dump(partition, partition_file)

with open(DATASET_PATH + LABELS_PREFIX + '.p', 'w') as labels_file:
    pickle.dump(labels, labels_file)
