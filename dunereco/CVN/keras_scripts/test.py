import numpy as np
import pickle
import configparser
import ast
import re
import logging, sys

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

CUT = float(config['test']['cut'])
TEST_BATCH_SIZE = int(config['test']['batch_size'])

# test params

y_test = []

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
               'shuffle': SHUFFLE,
               'y_test': y_test}


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

    files = [f for f in listdir(CHECKPOINT_PATH) if isfile(join(CHECKPOINT_PATH, f))]
    files.sort(reverse=True)

    r = re.compile(CHECKPOINT_PREFIX[1:] + '-.*-.*.h5')

    for fil in files:
        if r.match(fil) is not None:
            model = load_model(CHECKPOINT_PATH + '/' + fil)
            break

else:

    # Load the model

    model = load_model(CHECKPOINT_PATH + CHECKPOINT_PREFIX + '.h5')
   
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
y_pred = np.argmax(Y_pred, axis=1)
y_test = y_test[0:len(y_pred)]

# Add the interaction types for each neutrino type

cut = 0.7
reca_acum = 0.0

if INTERACTION_TYPES:

    logging.info('Adding the interaction types for each neutrino type...')

    Y_pred_neutrino = np.zeros((len(Y_pred), 4))
    y_test_neutrino = np.copy(y_test)

    for i in range(len(Y_pred)):

        p = Y_pred[i]

        # from 0 to 13 (12) labels to 0 to 13

        if y_test_neutrino[i] <= 3:
            y_test_neutrino[i] = 1
        elif y_test_neutrino[i] <= 7:
            y_test_neutrino[i] = 0
        elif y_test_neutrino[i] <= 11:
            y_test_neutrino[i] = 2
        else:
            y_test_neutrino[i] = 3

        Y_pred_neutrino[i][0] = p[4] + p[5] + p[6]  + p[7]  # (4,5,6,7)
        Y_pred_neutrino[i][1] = p[0] + p[1] + p[2]  + p[3]  # (0,1,2,3)
        Y_pred_neutrino[i][2] = p[8] + p[9] + p[10] + p[11] # (8,9,10,11)
        Y_pred_neutrino[i][3] = p[12]                       # (13)

else:
    Y_pred_neutrino = Y_pred
    y_test_neutrino = y_test 

# Report

y_pred_neutrino = np.argmax(Y_pred_neutrino, axis=1)

if INTERACTION_TYPES:
    
    # Interaction types (from 0 to 12 (13))

    inter_target_names = ['interac. 0', 'interac. 1', 'interac. 2', 'interac. 3', 'interac. 4', 'interac. 5', 
                         'interac. 6', 'interac. 7', 'interac. 8', 'interac. 9', 'interac. 10', 'interac. 11', 'interac. 13']

    logging.info('Classification report (interaction types):\n')

    print(classification_report(y_test, y_pred, target_names=inter_target_names))

    logging.info('Interaction types confusion matrix (rows = predicted classes, cols = actual classes):\n')

    print(confusion_matrix(y_pred, y_test))
    print ''

# Neutrino types (from 0 to 3)

neutrino_target_names = ['nue', 'numu', 'nutau', 'nc']

logging.info('Classification report (neutrino flavours):\n')

print(classification_report(y_test_neutrino, y_pred_neutrino, target_names=neutrino_target_names))

logging.info('Neutrino flavours confusion matrix (rows = predicted classes, cols = actual classes):\n')


conf2 = confusion_matrix(y_pred_neutrino, y_test_neutrino)
print conf2
print ''

# Apply cut

logging.info('Applying a cut of %.2f...\n' % CUT)

cut_acum = np.zeros((len(Y_pred_neutrino[0]), len(Y_pred_neutrino[0])))

y_pred_neutrino_cut = []
y_test_neutrino_cut = []

for sample in range(len(Y_pred_neutrino)):
    
    test_flavour = y_test_neutrino[sample] # get actual class of sample
    pred_flavour = y_pred_neutrino[sample] # get predicted class of sample

    for label in range(len(Y_pred_neutrino[sample])):
 
        if Y_pred_neutrino[sample][label] >= CUT:

            # accumulate if the probability of that label of the sample is >= CUT
            
            y_test_neutrino_cut.append(test_flavour)
            y_pred_neutrino_cut.append(pred_flavour)

logging.info('Classification report (neutrino flavours) after applying the cut:\n')

print(classification_report(y_test_neutrino_cut, y_pred_neutrino_cut, target_names=neutrino_target_names))

logging.info('Neutrino flavours (after applying the cut) confusion matrix (rows = predicted classes, cols = actual classes):\n')

conf3 = confusion_matrix(y_pred_neutrino_cut, y_test_neutrino_cut)
print conf3
print ''

logging.info('Applying factors... confusion matrix obtained:\n')

n_nue = sum(conf2[:,0])
n_numu = sum(conf2[:,1])
n_nutau = sum(conf2[:,2])
n_nc = sum(conf2[:,3])

float_formatter = lambda x: "%.4f" % x
np.set_printoptions(formatter={'float_kind':float_formatter})

conf4 = conf3.astype('float32') * np.array([0.114*n_numu/n_nue,1,0.0297*n_numu/n_nutau,1.0135*n_numu/n_nc])
print conf4.astype(int)
print ''

for i in range(len(conf4[0])):
    print 'Purity', neutrino_target_names[i], ': ', conf4[0][i] / sum(conf4[0])

print ''

logging.info('Efficiency:\n')

print conf3.astype('float32') / np.add.reduce(conf2)
#print conf3.astype('float32')/conf2.astype('float32')

