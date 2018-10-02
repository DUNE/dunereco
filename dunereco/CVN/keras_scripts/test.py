"""
This is the test module.
"""

__version__ = '1.0'
__author__ = 'Saul Alonso-Monsalve'
__email__ = "saul.alonso.monsalve@cern.ch"

import numpy as np
import pickle
import configparser
import ast
import re
import logging
import sys
import os
import time

sys.path.append(os.path.join(sys.path[0], 'modules'))

from random import shuffle
from keras import optimizers
from sklearn.metrics import classification_report, confusion_matrix
from os import listdir
from os.path import isfile, join
from keras.models import load_model
from data_generator import DataGenerator
from collections import Counter
import my_losses

# manually specify the GPUs to use
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"]="0"

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
SHUFFLE = ast.literal_eval(config['random']['shuffle'])
SHUFFLE = False

# images

IMAGES_PATH = config['images']['path']
VIEWS = int(config['images']['views'])
PLANES = int(config['images']['planes'])
CELLS = int(config['images']['cells'])
STANDARDIZE = ast.literal_eval(config['images']['standardize'])

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
PARALLELIZE = ast.literal_eval(config['model']['parallelize'])
OUTPUTS = int(config['model']['outputs'])

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

TEST_PARAMS = {'planes':PLANES,
               'cells':CELLS,
               'views':VIEWS,
               'batch_size':TEST_BATCH_SIZE,
               'branches':BRANCHES,
               'outputs': OUTPUTS,
               'images_path':IMAGES_PATH,
               'standardize':STANDARDIZE,
               'shuffle':SHUFFLE,
               'test_values':test_values}


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

'''
****************************************
************** GENERATORS **************
****************************************
'''

prediction_generator = DataGenerator(**TEST_PARAMS).generate(labels, partition['test'], False)


'''
****************************************
************** LOAD MODEL **************
****************************************
'''

# Load model

logging.info('Loading model from disk...')

if CHECKPOINT_SAVE_MANY:
    # Load the last generated model
    files = [f for f in os.listdir(CHECKPOINT_PATH) if os.path.isfile(os.path.join(CHECKPOINT_PATH, f))]
    files.sort(reverse=True)

    r = re.compile(CHECKPOINT_PREFIX[1:] + '-.*-.*.h5')

    for fil in files:
        if r.match(fil) is not None:
            model = load_model(CHECKPOINT_PATH + '/' + fil, 
                               custom_objects={'masked_loss':my_losses.masked_loss, 
                                               'multitask_loss': my_losses.multitask_loss,
                                               'masked_loss_binary': my_losses.masked_loss_binary, 
                                               'masked_loss_categorical': my_losses.masked_loss_categorical})
            logging.info('Loaded model: %s', CHECKPOINT_PATH + '/' + fil)
            break
else:
    # Load the model
    model = load_model(CHECKPOINT_PATH + CHECKPOINT_PREFIX + '.h5',
                       custom_objects={'masked_loss':my_losses.masked_loss,
                                       'multitask_loss': my_losses.multitask_loss,
                                       'masked_loss_binary': my_losses.masked_loss_binary,
                                       'masked_loss_categorical': my_losses.masked_loss_categorical})
  
    logging.info('Loaded model: %s', CHECKPOINT_PATH + CHECKPOINT_PREFIX + '.h5')

# Print model summary

if(PRINT_SUMMARY):
    model.summary()

'''
****************************************
***************** TEST *****************
****************************************
'''

logging.info('PERFORMING TEST...\n')

is_antineutrino_target_names = ['neutrino', 'antineutrino', 'NULL']
flavour_target_names = ['CC Numu', 'CC Nue', 'CC Nutau', 'NC']
interaction_target_names = ['CC QE', 'CC Res', 'CC DIS', 'CC Other', 'NULL']
categories_target_names = ['category 0', 'category 1', 'category 2', 'category 3', 'category 4', 'category 5', 'category 6', 
                           'category 7', 'category 8', 'category 9', 'category 10', 'category 11', 'category 13']
protons_target_names = ['0', '1', '2', '>2']
pions_target_names = ['0', '1', '2', '>2']
pizeros_target_names = ['0', '1', '2', '>2']
neutrons_target_names = ['0', '1', '2', '>2']

# Predict results

Y_pred = model.predict_generator(generator = prediction_generator,
                                 steps = len(partition['test'])//TEST_BATCH_SIZE,
                                 verbose = 1
                                )

#np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)}, threshold=np.nan)

if OUTPUTS == 1:
    # Single-output network
    test_values = np.array(test_values[0:Y_pred.shape[0]]) # array with y true value, energies and weights

    y_pred_categories = np.argmax(Y_pred, axis=1).reshape((Y_pred.shape[0], 1)) # 1-DIM array of predicted values (categories)
    Y_pred_flavour = np.zeros((Y_pred.shape[0], 4))                             # 2-DIM array of predicted values (flavour)

    y_test_categories = np.array([12 if aux['y_value'] == 13 else aux['y_value'] for aux in test_values]).reshape(y_pred_categories.shape) # 1-DIM array of test values (categories)
    y_test_flavour = np.zeros(y_test_categories.shape, dtype=int) # 1-DIM array of test values (flavour)

    # manually set y_pred_flavour and y_test_flavour

    for i in range(Y_pred_flavour.shape[0]):
        y_test_flavour[i] = y_test_categories[i]//4 # from  0-13 to 0-3

        # Add the interaction types for each neutrino type
        p = Y_pred[i]
        Y_pred_flavour[i][0] = p[0] + p[1] + p[2]  + p[3]  # NUMU  (0,1,2,3)
        Y_pred_flavour[i][1] = p[4] + p[5] + p[6]  + p[7]  # NUE   (4,5,6,7)
        Y_pred_flavour[i][2] = p[8] + p[9] + p[10] + p[11] # NUTAU (8,9,10,11)
        Y_pred_flavour[i][3] = p[12]                       # NC    (13)

    y_pred_flavour = np.argmax(Y_pred_flavour, axis=1) # 1-DIM array of predicted values (flavour)

    # flavour 

    logging.info('flavour report:\n')
    print(classification_report(y_test_flavour, y_pred_flavour, target_names=flavour_target_names))
    logging.info('flavour confusion matrix (rows = predicted classes, cols = actual classes):\n')
    flavour_conf_matrix = confusion_matrix(y_pred_flavour, y_test_flavour)
    print flavour_conf_matrix, '\n'

    # categories

    logging.info('categories report:\n')
    print(classification_report(y_test_categories, y_pred_categories, target_names=categories_target_names))
    logging.info('categories confusion matrix (rows = predicted classes, cols = actual classes):\n')
    categories_conf_matrix = confusion_matrix(y_pred_categories, y_test_categories)
    print categories_conf_matrix, '\n'

    # Apply cuts

    logging.info('Applying a numu cut of %.2f, a nue cut of %.2f, a nutau cut of %.2f, and a NC cut of %.2f...\n' % (CUT_NUMU, CUT_NUE, CUT_NUTAU, CUT_NC))

    weighted_conf_matrix     = np.zeros((4,4), dtype='float32')
    cut_weighted_conf_matrix = np.zeros((4,4), dtype='float32')

    for sample in range(len(Y_pred_flavour)):
        pred_flavour = int(y_pred_flavour[sample])  # get predicted class of sample
        test_flavour = int(y_test_flavour[sample])  # get actual class of sample
        weight = test_values[sample]['fEventWeight'] # event weight
        if Y_pred_flavour[sample][0] >= CUT_NUMU:
            cut_weighted_conf_matrix[0][test_flavour] += weight
        if Y_pred_flavour[sample][1] >= CUT_NUE:
            cut_weighted_conf_matrix[1][test_flavour] += weight
        if Y_pred_flavour[sample][2] >= CUT_NUTAU:
            cut_weighted_conf_matrix[2][test_flavour] += weight
        if Y_pred_flavour[sample][3] >= CUT_NC:
            cut_weighted_conf_matrix[3][test_flavour] += weight
        weighted_conf_matrix[pred_flavour][test_flavour] += weight

    # Confusion matrix - neutrino flavour (weighted)

    logging.info('Neutrino flavour weighted confusion matrix (rows = predicted classes, cols = actual classes):\n')
    print weighted_conf_matrix.astype(int), '\n'
    logging.info('Neutrino flavour weighted confusion matrix (rows = predicted classes, cols = actual classes) after applying the nue and numu cuts:\n')
    print cut_weighted_conf_matrix.astype(int), '\n'
    float_formatter = lambda x: "%.4f" % x
    np.set_printoptions(formatter={'float_kind':float_formatter})

    # Purity confusion matrix

    purity_conf_matrix = np.copy(cut_weighted_conf_matrix)
    for i in range(cut_weighted_conf_matrix.shape[0]):
        row_sum = np.sum(purity_conf_matrix[i])
        if(row_sum > 0):
            for j in range(cut_weighted_conf_matrix.shape[1]):
                purity_conf_matrix[i][j] /= row_sum
    logging.info('Purity confusion matrix (rows = predicted classes, cols = actual classes)\n')
    print purity_conf_matrix, '\n'

    # Efficiency confusion matrix

    logging.info('Efficiency confusion matrix (rows = predicted classes, cols = actual classes):\n')
    efficiency_conf_matrix = cut_weighted_conf_matrix.astype('float32')/ np.add.reduce(weighted_conf_matrix)
    print efficiency_conf_matrix, '\n'
else:
    # Multi-output network
    test_values = np.array(test_values[0:Y_pred[0].shape[0]]) # array with y true values, energies and weights

    y_pred_is_antineutrino = np.around(Y_pred[0]).reshape((Y_pred[0].shape[0], 1)).astype(int) # 1-DIM array of predicted values (is_antineutrino)
    y_pred_flavour = np.argmax(Y_pred[1], axis=1).reshape((Y_pred[1].shape[0], 1))             # 1-DIM array of predicted values (flavour)
    y_pred_interaction = np.argmax(Y_pred[2], axis=1).reshape((Y_pred[2].shape[0], 1))         # 1-DIM array of predicted values (interaction)
    y_pred_categories = np.zeros(y_pred_flavour.shape, dtype=int)                              # 1-DIM array of predicted values (categories)
    y_pred_protons = np.argmax(Y_pred[3], axis=1).reshape((Y_pred[3].shape[0], 1))             # 1-DIM array of predicted values (protons)
    y_pred_pions = np.argmax(Y_pred[4], axis=1).reshape((Y_pred[4].shape[0], 1))               # 1-DIM array of predicted values (pions)
    y_pred_pizeros = np.argmax(Y_pred[5], axis=1).reshape((Y_pred[5].shape[0], 1))             # 1-DIM array of predicted values (pizeros)
    y_pred_neutrons = np.argmax(Y_pred[6], axis=1).reshape((Y_pred[6].shape[0], 1))            # 1-DIM array of predicted values (neutrons)

    y_test_is_antineutrino = np.array([aux['y_value'][0] for aux in test_values]).reshape(y_pred_is_antineutrino.shape)
    y_test_flavour = np.array([aux['y_value'][1] for aux in test_values]).reshape(y_pred_flavour.shape)
    y_test_interaction = np.array([aux['y_value'][2] for aux in test_values]).reshape(y_pred_interaction.shape)
    y_test_categories = np.zeros(y_test_flavour.shape, dtype=int)
    y_test_protons = np.array([aux['y_value'][3] for aux in test_values]).reshape(y_pred_protons.shape)
    y_test_pions = np.array([aux['y_value'][4] for aux in test_values]).reshape(y_pred_pions.shape)
    y_test_pizeros = np.array([aux['y_value'][5] for aux in test_values]).reshape(y_pred_pizeros.shape)
    y_test_neutrons = np.array([aux['y_value'][6] for aux in test_values]).reshape(y_pred_neutrons.shape)
  
    # manually set y_pred_categories and y_test_categories

    for i in range(y_pred_categories.shape[0]):
        # inter
        y_pred_categories[i] = y_pred_interaction[i]
        y_test_categories[i] = y_test_interaction[i]    

        # flavour
        y_pred_categories[i] += (y_pred_flavour[i]*4)
        y_test_categories[i] += (y_test_flavour[i]*4)

        if y_pred_flavour[i] == 3:
            y_pred_is_antineutrino[i] = 2
            y_pred_interaction[i] = 4
            y_pred_categories[i] = 12

        if y_test_flavour[i] == 3:
            y_test_is_antineutrino[i] = 2
            y_test_interaction[i] = 4
            y_test_categories[i] = 12

    #np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)}, threshold=np.nan)

    # is_antineutrino

    logging.info('is_antineutrino report:\n')
    print(classification_report(y_test_is_antineutrino, y_pred_is_antineutrino, target_names=is_antineutrino_target_names))
    logging.info('is_antineutrino confusion matrix (rows = predicted classes, cols = actual classes):\n')
    is_antineutrino_conf_matrix = confusion_matrix(y_pred_is_antineutrino, y_test_is_antineutrino)
    print is_antineutrino_conf_matrix, '\n'

    # flavour 

    logging.info('flavour report:\n')
    print(classification_report(y_test_flavour, y_pred_flavour, target_names=flavour_target_names))
    logging.info('flavour confusion matrix (rows = predicted classes, cols = actual classes):\n')
    flavour_conf_matrix = confusion_matrix(y_pred_flavour, y_test_flavour)
    print flavour_conf_matrix, '\n'

    # interaction

    logging.info('interaction report:\n')
    print(classification_report(y_test_interaction, y_pred_interaction, target_names=interaction_target_names))
    logging.info('interaction confusion matrix (rows = predicted classes, cols = actual classes):\n')
    interaction_conf_matrix = confusion_matrix(y_pred_interaction, y_test_interaction)
    print interaction_conf_matrix, '\n'

    # categories

    logging.info('categories report:\n')
    print(classification_report(y_test_categories, y_pred_categories, target_names=categories_target_names))
    logging.info('categories confusion matrix (rows = predicted classes, cols = actual classes):\n')
    categories_conf_matrix = confusion_matrix(y_pred_categories, y_test_categories)
    print categories_conf_matrix, '\n'

    # protons

    logging.info('protons report:\n')
    print(classification_report(y_test_protons, y_pred_protons, target_names=protons_target_names))
    logging.info('protons confusion matrix (rows = predicted classes, cols = actual classes):\n')
    protons_conf_matrix = confusion_matrix(y_pred_protons, y_test_protons)
    print protons_conf_matrix, '\n'

    # pions

    logging.info('pions report:\n')
    print(classification_report(y_test_pions, y_pred_pions, target_names=pions_target_names))
    logging.info('pions confusion matrix (rows = predicted classes, cols = actual classes):\n')
    pions_conf_matrix = confusion_matrix(y_pred_pions, y_test_pions)
    print pions_conf_matrix, '\n'

    # pizeros

    logging.info('pizeros report:\n')
    print(classification_report(y_test_pizeros, y_pred_pizeros, target_names=pizeros_target_names))
    logging.info('pizeros confusion matrix (rows = predicted classes, cols = actual classes):\n')
    pizeros_conf_matrix = confusion_matrix(y_pred_pizeros, y_test_pizeros)
    print pizeros_conf_matrix, '\n'

    # neutrons

    logging.info('neutrons report:\n')
    print(classification_report(y_test_neutrons, y_pred_neutrons, target_names=neutrons_target_names))
    logging.info('neutrons confusion matrix (rows = predicted classes, cols = actual classes):\n')
    neutrons_conf_matrix = confusion_matrix(y_pred_neutrons, y_test_neutrons)
    print neutrons_conf_matrix, '\n'

    # Apply cuts

    logging.info('Applying a numu cut of %.2f, a nue cut of %.2f, a nutau cut of %.2f, and a NC cut of %.2f...\n' % (CUT_NUMU, CUT_NUE, CUT_NUTAU, CUT_NC))

    weighted_conf_matrix     = np.zeros((4,4), dtype='float32')
    cut_weighted_conf_matrix = np.zeros((4,4), dtype='float32')

    for sample in range(len(Y_pred[1])):
        pred_flavour = int(y_pred_flavour[sample])   # get predicted class of sample
        test_flavour = int(y_test_flavour[sample])   # get actual class of sample
        weight = test_values[sample]['fEventWeight'] # event weight
        if Y_pred[1][sample][0] >= CUT_NUMU:
            cut_weighted_conf_matrix[0][test_flavour] += weight
        if Y_pred[1][sample][1] >= CUT_NUE:
            cut_weighted_conf_matrix[1][test_flavour] += weight
        if Y_pred[1][sample][2] >= CUT_NUTAU:
            cut_weighted_conf_matrix[2][test_flavour] += weight
        if Y_pred[1][sample][3] >= CUT_NC:
            cut_weighted_conf_matrix[3][test_flavour] += weight
        weighted_conf_matrix[pred_flavour][test_flavour] += weight

    # Confusion matrix - neutrino flavour (weighted)

    logging.info('Neutrino flavour weighted confusion matrix (rows = predicted classes, cols = actual classes):\n')
    print weighted_conf_matrix.astype(int), '\n'
    logging.info('Neutrino flavour weighted confusion matrix (rows = predicted classes, cols = actual classes) after applying the nue and numu cuts:\n')
    print cut_weighted_conf_matrix.astype(int), '\n'
    float_formatter = lambda x: "%.4f" % x
    np.set_printoptions(formatter={'float_kind':float_formatter})

    # Purity confusion matrix

    purity_conf_matrix = np.copy(cut_weighted_conf_matrix)
    for i in range(cut_weighted_conf_matrix.shape[0]):
        row_sum = np.sum(purity_conf_matrix[i])
        if(row_sum > 0):
            for j in range(cut_weighted_conf_matrix.shape[1]):
                purity_conf_matrix[i][j] /= row_sum
    logging.info('Purity confusion matrix (rows = predicted classes, cols = actual classes)\n')
    print purity_conf_matrix, '\n'

    # Efficiency confusion matrix

    logging.info('Efficiency confusion matrix (rows = predicted classes, cols = actual classes):\n')
    efficiency_conf_matrix = cut_weighted_conf_matrix.astype('float32')/ np.add.reduce(weighted_conf_matrix)
    print efficiency_conf_matrix, '\n'

    # Dump test information

    logging.info('Dumping test information to \'%s\'...\n' % (OUTPUT_PATH + OUTPUT_PREFIX + '.np'))

    test_info = {'test_values':test_values,                        # Energy and weight values
                 'Y_pred':Y_pred,                                  # 3-DIM array of original probability predicted values
                 'y_pred_is_antineutrino':y_pred_is_antineutrino,  # 1-DIM array of is_antineutrino predicted values
                 'y_test_is_antineutrino':y_test_is_antineutrino,  # 1-DIM array of is_antineutrino test values
                 'y_pred_flavour':y_pred_flavour,                  # 1-DIM array of flavour predicted values
                 'y_test_flavour':y_test_flavour,                  # 1-DIM array of flavour test values
                 'y_pred_interaction':y_pred_interaction,          # 1-DIM array of interaction predicted values
                 'y_test_interaction':y_test_interaction,          # 1-DIM array of interaction test values
                 'y_pred_categories':y_pred_categories,            # 1-DIM array of categories predicted values
                 'y_test_categories':y_test_categories,            # 1-DIM array of categories test values
                 'y_pred_protons':y_pred_protons,                  # 1-DIM array of protons predicted values
                 'y_test_protons':y_test_protons,                  # 1-DIM array of protons test values
                 'y_pred_pions':y_pred_pions,                      # 1-DIM array of pions predicted values
                 'y_test_pions':y_test_pions,                      # 1-DIM array of pions test values
                 'y_pred_pizeros':y_pred_pizeros,                  # 1-DIM array of pizeros predicted values
                 'y_test_pizeros':y_test_pizeros,                  # 1-DIM array of pizeros test values
                 'y_pred_neutrons':y_pred_neutrons,                # 1-DIM array of neutrons predicted values
                 'y_test_neutrons':y_test_neutrons,                # 1-DIM array of neutrons test values
                 'is_antineutrino_cm':is_antineutrino_conf_matrix, # is_antineutrino confusion matrix
                 'flavour_cm':flavour_conf_matrix,                 # flavour confusion matrix
                 'interaction_cm':interaction_conf_matrix,         # interaction confusion matrix
                 'categories_cm':categories_conf_matrix,           # categories confusion matrix
                 'protons_cm':protons_conf_matrix,                 # protons confusion matrix
                 'pions_cm':pions_conf_matrix,                     # pions confusion matrix
                 'pizeros_cm':pizeros_conf_matrix,                 # pizeros confusion matrix
                 'neutrons_cm':neutrons_conf_matrix,               # neutrons confusion matrix
                 'cut_weighted_cm':cut_weighted_conf_matrix,       # Weighted neutrino types confusion matrix (after the cut)
                 'purity_cm':purity_conf_matrix,                   # Purity confusion matrix
                 'efficiency_cm':efficiency_conf_matrix            # Efficiency confusion matrix
                }

    test_info = np.array(test_info)

    with open(OUTPUT_PATH + OUTPUT_PREFIX + '.np', 'w') as test_info_file:
        test_info.dump(test_info_file)
