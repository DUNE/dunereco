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

sys.path.append("/raid/cvn_tensorflow/modules")

from keras import optimizers
from sklearn.metrics import classification_report, confusion_matrix
from os import listdir
from os.path import isfile, join
from keras.models import load_model
from data_generator import DataGenerator
from collections import Counter
import my_losses

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
    model = load_model(CHECKPOINT_PATH + CHECKPOINT_PREFIX + '.h5', custom_objects={'masked_loss':my_losses.masked_loss,'my_losses':my_losses.my_losses})
  
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

# Predict results

Y_pred = model.predict_generator(generator = prediction_generator,
                                 #steps = len(partition['test'])//TEST_BATCH_SIZE,
                                 steps = 20,
                                 verbose = 1
                                )

#np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)}, threshold=np.nan)

test_values = np.array(test_values[0:Y_pred[0].shape[0]]) # array with energies and weights

y_pred0 = np.around(Y_pred[0]).reshape((Y_pred[0].shape[0], 1)).astype(int) # 1-DIM array of predicted values (fNuPDG)
y_pred1 = np.argmax(Y_pred[1], axis=1).reshape((Y_pred[1].shape[0], 1))     # 1-DIM array of predicted values (flavour)
y_pred2 = np.argmax(Y_pred[2], axis=1).reshape((Y_pred[2].shape[0], 1))     # 1-DIM array of predicted values (interaction)
y_pred3 = np.zeros(y_pred1.shape, dtype=int)                                # 1-DIM array of predicted values (categories)
y_pred4 = np.argmax(Y_pred[3], axis=1).reshape((Y_pred[3].shape[0], 1))     # 1-DIM array of predicted values (fNProton)
y_pred5 = np.argmax(Y_pred[4], axis=1).reshape((Y_pred[4].shape[0], 1))     # 1-DIM array of predicted values (fNPion)
y_pred6 = np.argmax(Y_pred[5], axis=1).reshape((Y_pred[5].shape[0], 1))     # 1-DIM array of predicted values (fNPizero)
y_pred7 = np.argmax(Y_pred[6], axis=1).reshape((Y_pred[6].shape[0], 1))     # 1-DIM array of predicted values (fNNeutron)

y_test0 = np.array([aux['y_value'][0] for aux in test_values]).reshape(y_pred0.shape)
y_test1 = np.array([aux['y_value'][1] for aux in test_values]).reshape(y_pred1.shape)
y_test2 = np.array([aux['y_value'][2] for aux in test_values]).reshape(y_pred2.shape)
y_test3 = np.zeros(y_test1.shape, dtype=int)
y_test4 = np.array([aux['y_value'][3] for aux in test_values]).reshape(y_pred4.shape)
y_test5 = np.array([aux['y_value'][4] for aux in test_values]).reshape(y_pred5.shape)
y_test6 = np.array([aux['y_value'][5] for aux in test_values]).reshape(y_pred6.shape)
y_test7 = np.array([aux['y_value'][6] for aux in test_values]).reshape(y_pred7.shape)

# manually set y_pred3 and y_test3

for i in range(y_pred3.shape[0]):
    # inter
    y_pred3[i] = y_pred2[i]
    y_test3[i] = y_test2[i]    

    # flavour
    y_pred3[i] += (y_pred1[i]*4)
    y_test3[i] += (y_test1[i]*4)

    if y_pred1[i] == 3:
        y_pred0[i] = 2
        y_pred2[i] = 4
        y_pred3[i] = 12

    if y_test1[i] == 3:
        y_test0[i] = 2
        y_test2[i] = 4
        y_test3[i] = 12



#np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)}, threshold=np.nan)

fNuPDG_target_names = ['neutrino', 'antineutrino', 'NULL']
flavour_target_names = ['CC Numu', 'CC Nue', 'CC Nutau', 'NC']
interaction_target_names = ['CC QE', 'CC Res', 'CC DIS', 'CC Other', 'NULL']
categories_target_names = ['category 0', 'category 1', 'category 2', 'category 3', 'category 4', 'category 5', 'category 6', 
                           'category 7', 'category 8', 'category 9', 'category 10', 'category 11', 'category 13']
fNProton_target_names = ['0', '1', '2', '>2']
fNPion_target_names = ['0', '1', '2', '>2']
fNPizero_target_names = ['0', '1', '2', '>2']
fNNeutron_target_names = ['0', '1', '2', '>2']

# fNuPDG

logging.info('fNuPDG report:\n')
print(classification_report(y_test0, y_pred0, target_names=fNuPDG_target_names))
logging.info('fNuPDG confusion matrix (rows = predicted classes, cols = actual classes):\n')
fNuPDG_conf_matrix = confusion_matrix(y_pred0, y_test0)
print fNuPDG_conf_matrix, '\n'

# flavour 

logging.info('flavour report:\n')
print(classification_report(y_test1, y_pred1, target_names=flavour_target_names))
logging.info('flavour confusion matrix (rows = predicted classes, cols = actual classes):\n')
flavour_conf_matrix = confusion_matrix(y_pred1, y_test1)
print flavour_conf_matrix, '\n'

# interaction

logging.info('interaction report:\n')
print(classification_report(y_test2, y_pred2, target_names=interaction_target_names))
logging.info('interaction confusion matrix (rows = predicted classes, cols = actual classes):\n')
interaction_conf_matrix = confusion_matrix(y_pred2, y_test2)
print interaction_conf_matrix, '\n'

# categories

logging.info('categories report:\n')
print(classification_report(y_test3, y_pred3, target_names=categories_target_names))
logging.info('categories confusion matrix (rows = predicted classes, cols = actual classes):\n')
categories_conf_matrix = confusion_matrix(y_pred3, y_test3)
print categories_conf_matrix, '\n'

# fNProton

logging.info('fNProton report:\n')
print(classification_report(y_test4, y_pred4, target_names=fNProton_target_names))
logging.info('fNProton confusion matrix (rows = predicted classes, cols = actual classes):\n')
fNProton_conf_matrix = confusion_matrix(y_pred4, y_test4)
print fNProton_conf_matrix, '\n'

# fNPion

logging.info('fNPion report:\n')
print(classification_report(y_test5, y_pred5, target_names=fNPion_target_names))
logging.info('fNPion confusion matrix (rows = predicted classes, cols = actual classes):\n')
fNPion_conf_matrix = confusion_matrix(y_pred5, y_test5)
print fNPion_conf_matrix, '\n'

# fNPizero

logging.info('fNPizero report:\n')
print(classification_report(y_test6, y_pred6, target_names=fNPizero_target_names))
logging.info('fNPizero confusion matrix (rows = predicted classes, cols = actual classes):\n')
fNPizero_conf_matrix = confusion_matrix(y_pred6, y_test6)
print fNPizero_conf_matrix, '\n'

# fNNeutron

logging.info('fNNeutron report:\n')
print(classification_report(y_test7, y_pred7, target_names=fNNeutron_target_names))
logging.info('fNNeutron confusion matrix (rows = predicted classes, cols = actual classes):\n')
fNNeutron_conf_matrix = confusion_matrix(y_pred7, y_test7)
print fNNeutron_conf_matrix, '\n'

# Apply cuts

logging.info('Applying a nue cut of %.2f, a numu cut of %.2f, a nutau cut of %.2f, and a NC cut of %.2f...\n' % (CUT_NUE, CUT_NUMU, CUT_NUTAU, CUT_NC))

weighted_conf_matrix     = np.zeros((4,4), dtype='float32')
cut_weighted_conf_matrix = np.zeros((4,4), dtype='float32')

for sample in range(len(Y_pred[1])):
    pred_flavour = int(y_pred1[sample])  # get predicted class of sample
    test_flavour = int(y_test1[sample])  # get actual class of sample
    weight = test_values[sample]['fEventWeight'] # event weight

    if Y_pred[1][sample][0] >= CUT_NUE:
        # accumulate if the probability of that label of the sample is >= CUT_NUE
        cut_weighted_conf_matrix[0][test_flavour] += weight

    if Y_pred[1][sample][1] >= CUT_NUMU:
        # accumulate if the probability of that label of the sample is >= CUT_NUMU
        cut_weighted_conf_matrix[1][test_flavour] += weight

    if Y_pred[1][sample][2] >= CUT_NUTAU:
        # accumulate if the probability of that label of the sample is >= CUT_NUTAU
        cut_weighted_conf_matrix[2][test_flavour] += weight

    if Y_pred[1][sample][3] >= CUT_NC:
        # accumulate if the probability of that label of the sample is >= CUT_NC
        cut_weighted_conf_matrix[3][test_flavour] += weight

    weighted_conf_matrix[pred_flavour][test_flavour] += weight

# Confusion matrix - neutrino types (weighted)

logging.info('Neutrino types weighted confusion matrix (rows = predicted classes, cols = actual classes):\n')

print weighted_conf_matrix.astype(int), '\n'

logging.info('Neutrino types weighted confusion matrix (rows = predicted classes, cols = actual classes) after applying the nue and numu cuts:\n')

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

test_info = {'test_values':test_values,                  # Energy and weight values
             'Y_pred':Y_pred,                            # 3-DIM array of original probability predicted values
             'y_pred0':y_pred0,                          # 1-DIM array of fNuPDG predicted values
             'y_pred1':y_pred1,                          # 1-DIM array of flavour predicted values
             'y_pred2':y_pred2,                          # 1-DIM array of interaction predicted values
             'y_pred3':y_pred3,                          # 1-DIM array of categories predicted values
             'y_pred4':y_pred4,                          # 1-DIM array of fNProton predicted values
             'y_pred5':y_pred5,                          # 1-DIM array of fNPion predicted values
             'y_pred6':y_pred6,                          # 1-DIM array of fNPizero predicted values
             'y_pred7':y_pred7,                          # 1-DIM array of fNNeutron predicted values
             'y_test0':y_test0,                          # 1-DIM array of fNuPDG test values
             'y_test1':y_test1,                          # 1-DIM array of flavour test values
             'y_test2':y_test2,                          # 1-DIM array of interaction test values
             'y_test3':y_test3,                          # 1-DIM array of categories test values
             'y_test4':y_test4,                          # 1-DIM array of fNProton test values
             'y_test5':y_test5,                          # 1-DIM array of fNPion test values
             'y_test6':y_test6,                          # 1-DIM array of fNPizero test values
             'y_test7':y_test7,                          # 1-DIM array of fNNeutron test values
             'fNuPDG_cm':fNuPDG_conf_matrix,             # fNuPDG confusion matrix
             'flavour_cm':flavour_conf_matrix,           # flavour confusion matrix
             'interaction_cm':interaction_conf_matrix,   # interaction confusion matrix
             'categories_cm':categories_conf_matrix,     # categories confusion matrix
             'fNProton_cm':fNProton_conf_matrix,         # fNProton confusion matrix
             'fNPion_cm':fNPion_conf_matrix,             # fNPion confusion matrix
             'fNPizero_cm':fNPizero_conf_matrix,         # fNPizero confusion matrix
             'fNNeutron_cm':fNNeutron_conf_matrix,       # fNNeutron confusion matrix
             'cut_weighted_cm':cut_weighted_conf_matrix, # Weighted neutrino types confusion matrix (after the cut)
             'purity_cm':purity_conf_matrix,             # Purity confusion matrix
             'efficiency_cm':efficiency_conf_matrix      # Efficiency confusion matrix
            }

test_info = np.array(test_info)

with open(OUTPUT_PATH + OUTPUT_PREFIX + '.np', 'w') as test_info_file:
    test_info.dump(test_info_file)
