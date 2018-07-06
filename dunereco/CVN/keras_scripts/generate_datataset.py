import numpy as np
import glob
import ast
import ntpath
import pickle
import configparser
import logging, sys
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
config.read('config.ini')

# random

SEED = int(config['random']['seed'])

if SEED == -1:
    SEED = int(time.time())

np.random.seed(SEED)

# images

IMAGES_PATH = config['images']['path']
INTERACTION_LABELS = ast.literal_eval(config['images']['interaction_labels'])

# dataset

DATASET_PATH = config['dataset']['path']
PARTITION_PREFIX = config['dataset']['partition_prefix']
LABELS_PREFIX = config['dataset']['labels_prefix']
UNIFORM = ast.literal_eval(config['dataset']['uniform'])
INTERACTION_TYPES = ast.literal_eval(config['dataset']['interaction_types'])
NEUTRINO_LABELS = ast.literal_eval(config['images']['neutrino_labels'])
 
if(INTERACTION_TYPES):

    # Interaction types (from 0 to 13 (12))

    N_LABELS = len(Counter(INTERACTION_LABELS.values()))

else:

    # Neutrino types (from 0 to 3)

    N_LABELS = len(Counter(NEUTRINO_LABELS.values()))

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

    # calculate uniform constants

    logging.info('Calculating uniform constants...')

    COUNTS = np.zeros(N_LABELS, dtype=int)
    MIN_DIR = [float("inf"), '']

    for path in glob.iglob(IMAGES_PATH + '/*'):

        label = ntpath.basename(path) # label
        label_dir = glob.iglob(path + '/*')
        len_dir = float(len(list(label_dir)))

        if INTERACTION_TYPES:

            # Interaction types (from 0 to 12 (13))

            COUNTS[INTERACTION_LABELS[label]] += 1

            if len_dir < MIN_DIR[0]:

                MIN_DIR[0] = len_dir
                MIN_DIR[1] = label

        else:

            # Neutrino types: from 0 to 3

            COUNTS[NEUTRINO_LABELS[label]] += 1

            if len_dir < MIN_DIR[0]:

                MIN_DIR[0] = len_dir
                MIN_DIR[1] = label

# Iterate through label folders

logging.info('Filling datasets...')

for path in glob.iglob(IMAGES_PATH + '/*'):

    logging.debug('path: %s', path)

    label = ntpath.basename(path) # label
    label_dir = glob.iglob(path + '/*')

    if UNIFORM:

        # Calculate fraction

        if INTERACTION_TYPES:

            # Interaction types (from 0 to 13 (12))

            fraction = (MIN_DIR[0]) * (COUNTS[INTERACTION_LABELS[MIN_DIR[1]]] / COUNTS[INTERACTION_LABELS[label]]) / len(list(label_dir))

        else:

            # Neutrino types (from 0 to 3)

            fraction = (MIN_DIR[0]) * (COUNTS[NEUTRINO_LABELS[MIN_DIR[1]]] / COUNTS[NEUTRINO_LABELS[label]]) / len(list(label_dir))

        logging.debug('fraction: %f', fraction)  
 
    count_train, count_val, count_test = (0, 0, 0)
    label_dir = glob.iglob(path + '/*')
 
    # Iterate through the files under the directory

    count = 0

    fils = list(label_dir)
    random.shuffle(fils)

    for fil in fils:

        if(fil[-3:] != '.gz'):
            continue

        # REMOVE
        if label == "13" and count  >= 886894:
            break

        count+=1
        # END REMOVE        


        ID = fil[:-7].split('/')[-1]  # it is basically the filename without '.txt.gz'

        random_value = np.random.uniform(0,1)

        # Fill training set

        if(random_value < TRAIN_FRACTION):

            random_value_uni = np.random.uniform(0,1)

            if not UNIFORM or random_value_uni < fraction:

                partition['train'].append(ID)
                count_train += 1

                # Store y train

                if INTERACTION_TYPES:

                    y_train.append(INTERACTION_LABELS[label])
                    y1_class_weights.append(NEUTRINO_LABELS[label])

                    if label != '13':
                        y2_class_weights.append(int(INTERACTION_LABELS[label])%4)

                else:

                    y_train.append(NEUTRINO_LABELS[label])
 
                labels[ID] = label

        # Fill validation set

        elif(random_value < (TRAIN_FRACTION + VALIDATION_FRACTION)):

            partition['validation'].append(ID)
            count_val += 1

            labels[ID] = label

        # Fill test set

        elif(random_value < (TRAIN_FRACTION + VALIDATION_FRACTION + TEST_FRACTION)):

            partition['test'].append(ID)
            count_test += 1

            labels[ID] = label

    logging.debug('%d train images', count_train)
    logging.debug('%d val images', count_val)
    logging.debug('%d test images', count_test)
    logging.debug('%d total images', count_train + count_val + count_test)

# Print dataset statistics

logging.info('Number of training examples: %d', len(partition['train']))
logging.info('Number of validation examples: %d', len(partition['validation']))
logging.info('Number of test examples: %d', len(partition['test']))

# Calculate class weights (used for a weighted loss function)

logging.info('Calculating class weights...')

#class_weights = dict(enumerate(class_weight.compute_class_weight('balanced', np.unique(y_train), y_train)))
class_weights1 = dict(enumerate(class_weight.compute_class_weight('balanced', np.unique(y1_class_weights), y1_class_weights)))
class_weights2 = dict(enumerate(class_weight.compute_class_weight('balanced', np.unique(y2_class_weights), y2_class_weights)))

class_weights = class_weights1

for i in range(4):
    class_weights[i+4] = class_weights2[i]

'''
# Manually set class weights

class_weights[0]*=0.8
class_weights[1]*=0.8
class_weights[2]*=1.2
class_weights[3]*=0.8
'''

logging.info('Class weights: %s', class_weights)

# Serialize partition and labels

logging.info('Serializing datasets...')

with open(DATASET_PATH + PARTITION_PREFIX + '.p', 'w') as partition_file:
    pickle.dump(partition, partition_file)

with open(DATASET_PATH + LABELS_PREFIX + '.p', 'w') as labels_file:
    pickle.dump(labels, labels_file)

if WEIGHTED_LOSS_FUNCTION:

    with open(DATASET_PATH + CLASS_WEIGHTS_PREFIX + '.p', 'w') as class_weights_file:
        pickle.dump(class_weights, class_weights_file)
