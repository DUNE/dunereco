import numpy as np
import glob
import ntpath
import pickle
import configparser
import logging, sys

'''
****************************************
************** PARAMETERS **************
****************************************
'''

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

config = configparser.ConfigParser()
config.read('config.ini')

# random

np.random.seed(int(config['random']['seed']))

# images

IMAGES_PATH = config['images']['path']

# dataset

DATASET_PATH = config['dataset']['path']
PARTITION_PREFIX = config['dataset']['partition_prefix']
LABELS_PREFIX = config['dataset']['labels_prefix']

# train

TRAIN_FRACTION = float(config['train']['fraction'])

# validation

VALIDATION_FRACTION = float(config['validation']['fraction'])

# test

TEST_FRACTION = float(config['test']['fraction'])

if((TRAIN_FRACTION + VALIDATION_FRACTION + TEST_FRACTION) > 1):
    logging.error('(TRAIN_FRACTION + VALIDATION_FRACTION + TEST_FRACTION) must be <= 1')
    exit(-1)

'''
****************************************
*************** DATASETS ***************
****************************************
'''

partition = {'train' : [], 'validation' : [], 'test' : []} # Train, validation, and test IDs
labels = {}                                                # ID : label

logging.info('Filling datasets...')

# Iterate through the label folders

for path in glob.iglob(IMAGES_PATH + '/*'):

    logging.debug('path: %s', path)

    # Iterate through the files under the directory

    for fil in glob.iglob(path + '/*'):

        ID = fil[:-7].split('/')[-1]  # it is basically the filename without '.txt.gz'
        label = ntpath.basename(path) # label

        random_value = np.random.uniform(0,1)

        # Fill training set

        if(random_value < TRAIN_FRACTION):
            partition['train'].append(ID)

        # Fill validation set

        elif(random_value < (TRAIN_FRACTION + VALIDATION_FRACTION)):
            partition['validation'].append(ID)

        # Fill test set

        elif(random_value < (TRAIN_FRACTION + VALIDATION_FRACTION + TEST_FRACTION)):
            partition['test'].append(ID)

        # Label

        labels[ID] = label

# Print dataset statistics

logging.info('Number of training examples: %d', len(partition['train']))
logging.info('Number of validation examples: %d', len(partition['validation']))
logging.info('Number of test examples: %d', len(partition['test']))

# Serialize partition and labels

logging.info('Serializing datasets...')

partition_file = open(DATASET_PATH + PARTITION_PREFIX + '.p', 'w')
pickle.dump(partition, partition_file)
partition_file.close()

labels_file = open(DATASET_PATH + LABELS_PREFIX + '.p', 'w')
pickle.dump(labels, labels_file)
labels_file.close()

