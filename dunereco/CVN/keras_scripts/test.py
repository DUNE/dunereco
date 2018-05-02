import numpy as np
import pickle
import configparser
import ast
import re
import logging, sys
import os

from sklearn.metrics import classification_report, confusion_matrix
from os import listdir
from os.path import isfile, join
from keras.models import load_model
from data_generator import DataGenerator
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

np.random.seed(int(config['random']['seed']))
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

# test

OUTPUT_PATH = config['test']['output_path']
OUTPUT_PREFIX = config['test']['output_prefix']
CUT_NUE = float(config['test']['cut_nue'])
CUT_NUMU = float(config['test']['cut_numu'])
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

partition_file = open(DATASET_PATH + PARTITION_PREFIX + '.p', 'r')
partition = pickle.load(partition_file)
partition_file.close()

labels_file = open(DATASET_PATH + LABELS_PREFIX + '.p', 'r')
labels = pickle.load(labels_file)
labels_file.close()

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
            model = load_model(CHECKPOINT_PATH + '/' + fil)

            logging.info('Loaded model: %s', CHECKPOINT_PATH + '/' + fil)

            break

else:

    # Load the model

    model = load_model(CHECKPOINT_PATH + CHECKPOINT_PREFIX + '.h5')
  
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
                                 steps = len(partition['test'])//TEST_BATCH_SIZE,
                                 verbose = 1
                                )

test_values = np.array(test_values[0:Y_pred.shape[0]]) # array with energies and weights

y_pred = np.argmax(Y_pred, axis=1).reshape((Y_pred.shape[0], 1))                 # 1-DIM array of predicted values
y_test = np.array([aux['y_value'] for aux in test_values]).reshape(y_pred.shape) # 1-DIM array of test values

'''
print '###############'

for i in Y_pred:
    print i

print '###############'
'''

if INTERACTION_TYPES:
    
    # Interaction types (from 0 to 12 (13))

    inter_target_names = ['interac. 0', 'interac. 1', 'interac. 2', 'interac. 3', 'interac. 4', 'interac. 5', 
                         'interac. 6', 'interac. 7', 'interac. 8', 'interac. 9', 'interac. 10', 'interac. 11', 'interac. 13']

    logging.info('Classification report (interaction types):\n')

    print(classification_report(y_test, y_pred, target_names=inter_target_names))

    # Confusion matrix - interaction types 

    logging.info('Interaction types confusion matrix (rows = predicted classes, cols = actual classes):\n')

    inter_conf_matrix = confusion_matrix(y_pred, y_test)
    print inter_conf_matrix, '\n'

    Y_pred_neutrino = np.zeros((Y_pred.shape[0], 4))
    y_test_neutrino = np.zeros(y_test.shape)

    # From interaction types to neutrino types

    logging.info('From interaction types to neutrino types')

    for i in range(len(Y_pred)):

        p = Y_pred[i]     # array of probabilities
        inter = y_test[i] # interaction type

        # labels from 0 to 13 (12) labels to 0 to 3

        if inter <= 3:

            y_test_neutrino[i] = 1 # NUMU

        elif inter <= 7:

            y_test_neutrino[i] = 0 # NUE

        elif inter <= 11:

            y_test_neutrino[i] = 2 # NUTAU

        else:
 
            y_test_neutrino[i] = 3 # NC

        # Add the interaction types for each neutrino type

        Y_pred_neutrino[i][0] = p[4] + p[5] + p[6]  + p[7]  # NUE   (4,5,6,7)
        Y_pred_neutrino[i][1] = p[0] + p[1] + p[2]  + p[3]  # NUMU  (0,1,2,3)
        Y_pred_neutrino[i][2] = p[8] + p[9] + p[10] + p[11] # NUTAU (8,9,10,11)
        Y_pred_neutrino[i][3] = p[12]                       # NC    (13)

else:

    # Y_pred and y_test already store neutrino types data

    Y_pred_neutrino = Y_pred
    y_test_neutrino = y_test 

# Report

y_pred_neutrino = np.argmax(Y_pred_neutrino, axis=1) # 1-DIM array of predicted values

# Neutrino types (from 0 to 3)

neutrino_target_names = ['nue', 'numu', 'nutau', 'nc']

logging.info('Classification report (neutrino flavours):\n')

print(classification_report(y_test_neutrino, y_pred_neutrino, target_names=neutrino_target_names))

# Confusion matrix - neutrino types 

logging.info('Neutrino types confusion matrix (rows = predicted classes, cols = actual classes):\n')

neutrino_conf_matrix = confusion_matrix(y_pred_neutrino, y_test_neutrino)
print neutrino_conf_matrix, '\n'

# Apply cuts

logging.info('Applying a nue cut of %.2f and a numu cut of %.2f...\n' % (CUT_NUE, CUT_NUMU))

weighted_conf_matrix     = np.zeros((4,4), dtype='float32')
cut_weighted_conf_matrix = np.zeros((4,4), dtype='float32')

for sample in range(len(Y_pred_neutrino)):
    
    pred_flavour = int(y_pred_neutrino[sample]) # get predicted class of sample
    test_flavour = int(y_test_neutrino[sample]) # get actual class of sample
    weight = test_values[sample]['fEventWeight']

    if Y_pred_neutrino[sample][0] >= CUT_NUE:

        cut_weighted_conf_matrix[0][test_flavour] += weight

        # accumulate if the probability of that label of the sample is >= CUT_NUE
 
    if Y_pred_neutrino[sample][1] >= CUT_NUMU:

        cut_weighted_conf_matrix[1][test_flavour] += weight

        # accumulate if the probability of that label of the sample is >= CUT_NUE

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

test_info = {'test_values':     test_values,              # Energies and weights values
             'Y_pred':          Y_pred,                   # 2-DIM array of original probability predicted values
             'y_pred':          y_pred,                   # 1-DIM array or original predicted values
             'y_test':          y_test,                   # 1-DIM array of original test values
             'Y_pred_neutrino': Y_pred_neutrino,          # 2-DIM array of neutrino probability predicted values
             'y_pred_neutrino': y_pred_neutrino,          # 1-DIM array of neutrino predicted values
             'y_test_neutrino': y_test_neutrino,          # 1-DIM array of neutrino test values
             'interaction_cm':  inter_conf_matrix,        # Interaction types confusion matrix
             'neutrino_cm':     neutrino_conf_matrix,     # Neutrino types confusion matrix
             'weighted_cm':     weighted_conf_matrix,     # Weighted neutrino types confusion matrix
             'cut_weighted_cm': cut_weighted_conf_matrix, # Weighted neutrino types confusion matrix (after the cut)
             'purity_cm':       purity_conf_matrix,       # Purity confusion matrix
             'efficiency_cm':   efficiency_conf_matrix    # Efficiency confusion matrix
            }
test_info = np.array(test_info)

test_info_file = open(OUTPUT_PATH + OUTPUT_PREFIX + '.np', 'w')
test_info.dump(test_info_file)
test_info_file.close()

logging.info('Test finished!')
