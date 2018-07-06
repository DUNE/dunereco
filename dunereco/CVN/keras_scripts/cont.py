import numpy as np
import pickle
import configparser
import ast
import re
import logging, sys
import os
import time
import glob

from sklearn.metrics import classification_report, confusion_matrix
from os import listdir
from os.path import isfile, join
from keras.models import load_model
from data_generator import DataGenerator
from collections import Counter

sys.path.append("/home/salonsom/cvn_tensorflow/loss")

import multitask_loss

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
SHUFFLE = ast.literal_eval(config['random']['shuffle'])

# images

IMAGES_PATH = config['images']['path']
VIEWS = int(config['images']['views'])
PLANES = int(config['images']['planes'])
CELLS = int(config['images']['cells'])
STANDARDIZE = ast.literal_eval(config['images']['standardize'])
INTERACTION_LABELS = ast.literal_eval(config['images']['interaction_labels'])
FILTERED = ast.literal_eval(config['images']['filtered'])

INTERACTION_TYPES = ast.literal_eval(config['dataset']['interaction_types'])

if(INTERACTION_TYPES):

    # Interaction types (from 0 to 13 (12))

    NEUTRINO_LABELS = []
    N_LABELS = len(Counter(INTERACTION_LABELS.values()))

else:

    # Neutrino types (from 0 to 3)

    NEUTRINO_LABELS = ast.literal_eval(config['images']['neutrino_labels'])
    N_LABELS = len(Counter(NEUTRINO_LABELS.values()))

N_LABELS = 8

# dataset

DATASET_PATH = config['dataset']['path']
PARTITION_PREFIX = config['dataset']['partition_prefix']
LABELS_PREFIX = config['dataset']['labels_prefix']

# model

CHECKPOINT_PATH = config['model']['checkpoint_path']
CHECKPOINT_PREFIX = config['model']['checkpoint_prefix']
CHECKPOINT_SAVE_MANY = ast.literal_eval(config['model']['checkpoint_save_many'])
CHECKPOINT_SAVE_BEST_ONLY = ast.literal_eval(config['model']['checkpoint_save_best_only'])
PRINT_SUMMARY = ast.literal_eval(config['model']['print_summary'])
BRANCHES = ast.literal_eval(config['model']['branches'])

# test

OUTPUT_PATH = config['test']['output_path']
OUTPUT_PREFIX = config['test']['output_prefix']
CUT_NUE = float(config['test']['cut_nue'])
CUT_NUMU = float(config['test']['cut_numu'])
CUT_NUTAU = float(config['test']['cut_nutau'])
CUT_NC = float(config['test']['cut_nc'])
TEST_BATCH_SIZE = int(config['test']['batch_size'])

# test params

test_values = []

TEST_PARAMS = {'planes': PLANES,
               'cells': CELLS,
               'views': VIEWS,
               'batch_size': TEST_BATCH_SIZE,
               'n_labels': N_LABELS,
               'interaction_labels': INTERACTION_LABELS,
               'interaction_types': INTERACTION_TYPES,
               'branches': BRANCHES,
               'filtered': FILTERED,
               'neutrino_labels': NEUTRINO_LABELS,
               'images_path': IMAGES_PATH,
               'standardize': STANDARDIZE,
               'shuffle': SHUFFLE,
               'test_values': test_values}


'''
****************************************
*************** DATASETS ***************
****************************************
'''

# Load datasets

logging.info('Reading datasets from serialized files...')

with open(DATASET_PATH + PARTITION_PREFIX + '.p', 'r') as partition_file:
    partition = pickle.load(partition_file)

with open(DATASET_PATH + LABELS_PREFIX + '.p', 'r') as labels_file:
    labels = pickle.load(labels_file)

# Print some dataset statistics

logging.info('Number of training examples: %d', len(partition['train']))
logging.info('Number of validation examples: %d', len(partition['validation']))
logging.info('Number of test examples: %d', len(partition['test']))

prediction_generator_train = DataGenerator(**TEST_PARAMS).generate(labels, partition['train'],  False)
prediction_generator_test  = DataGenerator(**TEST_PARAMS).generate(labels, partition['test'],   False)

'''
****************************************
***************** TEST *****************
****************************************
'''

logging.info('PERFORMING TEST...\n')

prefix = '/scratch/cvn/seresnet34/checkpoint/multitask7-*'

for model_name in sorted(list(glob.iglob(prefix))):

    print model_name

    model = load_model(model_name, custom_objects={'masked_loss': multitask_loss.masked_loss, 'multitask_loss': multitask_loss.multitask_loss})

    # TRAIN

    test_values = []
    TEST_PARAMS['test_values'] = test_values
    TEST_PARAMS['shuffle'] = True
    prediction_generator_train  = DataGenerator(**TEST_PARAMS).generate(labels, partition['train'],   False)

    print "TRAIN"

    Y_pred = model.predict_generator(generator = prediction_generator_train,
                                     steps = len(partition['test'])//TEST_BATCH_SIZE,
                                     verbose = 0
                                    )

    test_values = np.array(test_values[0:Y_pred.shape[0]]) # array with energies and weights

    y_pred1 = np.argmax(Y_pred[:,:4], axis=1).reshape((Y_pred.shape[0], 1))                 # 1-DIM array of predicted values (neutrino)
    y_pred2 = np.argmax(Y_pred[:,4:], axis=1).reshape((Y_pred.shape[0], 1))                 # 1-DIM array of predicted values (interaction)
    y_pred3 = np.zeros(y_pred1.shape, dtype=int)                                            # 1-DIM array of predicted values (categories)

    for i in range(y_pred1.shape[0]):

        # inter

        y_pred3[i] = y_pred2[i]
    
        # flavour

        y_pred3[i] += (y_pred1[i]*4)

        if y_pred1[i] == 3:
            y_pred3[i] = 12
            y_pred2[i] = 4

    y_test1 = np.zeros(y_pred1.shape, dtype=int)
    y_test2 = np.zeros(y_pred1.shape, dtype=int)
    y_test3 = np.array([aux['y_value'] for aux in test_values]).reshape(y_pred1.shape)      # 1-DIM array of test values

    for i in range(y_test1.shape[0]):

        y_test1[i] = y_test3[i]//4 
        y_test2[i] = y_test3[i]%4

        if y_test1[i] == 3:
            y_test2[i] = 4

    #np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)}, threshold=np.nan)

    neutrino_target_names = ['CC Numu', 'CC Nue', 'CC Nutau', 'NC']
    inter_target_names    = ['CC QE', 'CC Res', 'CC DIS', 'CC Other', 'NULL']
    catego_target_names   = ['category 0', 'category 1', 'category 2', 'category 3', 'category 4', 'category 5', 'category 6', 
                         'category 7', 'category 8', 'category 9', 'category 10', 'category 11', 'category 13']

    # NEUTRINO FLAVOUR

    logging.info('Neutrino flavour report:\n')

    print(classification_report(y_test1, y_pred1, target_names=neutrino_target_names))

    '''
    # Confusion matrix - interaction types 

    logging.info('Neutrino flavour confusion matrix (rows = predicted classes, cols = actual classes):\n')

    inter_conf_matrix = confusion_matrix(y_pred1, y_test1)
    print inter_conf_matrix, '\n'
    '''

    # INTERACTION

    logging.info('Interaction type report:\n')

    print(classification_report(y_test2, y_pred2, target_names=inter_target_names))

    '''
    # Confusion matrix - interaction types 

    logging.info('Interaction type confusion matrix (rows = predicted classes, cols = actual classes):\n')

    inter_conf_matrix = confusion_matrix(y_pred2, y_test2)
    print inter_conf_matrix, '\n'

    # CATEGORIES

    logging.info('Categories report:\n')

    print(classification_report(y_test3, y_pred3, target_names=catego_target_names))

    # Confusion matrix - interaction types 

    logging.info('Categories confusion matrix (rows = predicted classes, cols = actual classes):\n')

    inter_conf_matrix = confusion_matrix(y_pred3, y_test3)
    print inter_conf_matrix, '\n'
    '''

    # TEST

    test_values = []
    TEST_PARAMS['test_values'] = test_values
    TEST_PARAMS['shuffle'] = False
    prediction_generator_test  = DataGenerator(**TEST_PARAMS).generate(labels, partition['test'],   False)

    print "TEST"

    Y_pred = model.predict_generator(generator = prediction_generator_test,
                                     steps = len(partition['test'])//TEST_BATCH_SIZE,
                                     verbose = 0
                                    )

    test_values = np.array(test_values[0:Y_pred.shape[0]]) # array with energies and weights

    y_pred1 = np.argmax(Y_pred[:,:4], axis=1).reshape((Y_pred.shape[0], 1))                 # 1-DIM array of predicted values (neutrino)
    y_pred2 = np.argmax(Y_pred[:,4:], axis=1).reshape((Y_pred.shape[0], 1))                 # 1-DIM array of predicted values (interaction)
    y_pred3 = np.zeros(y_pred1.shape, dtype=int)                                            # 1-DIM array of predicted values (categories)

    for i in range(y_pred1.shape[0]):

        # inter

        y_pred3[i] = y_pred2[i]
    
        # flavour

        y_pred3[i] += (y_pred1[i]*4)

        if y_pred1[i] == 3:
            y_pred3[i] = 12
            y_pred2[i] = 4

    y_test1 = np.zeros(y_pred1.shape, dtype=int)
    y_test2 = np.zeros(y_pred1.shape, dtype=int)
    y_test3 = np.array([aux['y_value'] for aux in test_values]).reshape(y_pred1.shape)      # 1-DIM array of test values

    for i in range(y_test1.shape[0]):

        y_test1[i] = y_test3[i]//4 
        y_test2[i] = y_test3[i]%4

        if y_test1[i] == 3:
            y_test2[i] = 4

    #np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)}, threshold=np.nan)

    # NEUTRINO FLAVOUR

    logging.info('Neutrino flavour report:\n')

    print(classification_report(y_test1, y_pred1, target_names=neutrino_target_names))

    # INTERACTION

    logging.info('Interaction type report:\n')

    print(classification_report(y_test2, y_pred2, target_names=inter_target_names))
