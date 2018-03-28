import numpy as np
import pickle
import configparser
import ast
import re
import logging, sys

from os import listdir
from os.path import isfile, join
from keras.models import load_model
from data_generator import DataGenerator


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
SHUFFLE = ast.literal_eval(config['random']['shuffle'])

# images

IMAGES_PATH = config['images']['path']
VIEWS = int(config['images']['views'])
PLANES = int(config['images']['planes'])
CELLS = int(config['images']['cells'])
LABELS = ast.literal_eval(config['images']['labels'])
N_LABELS = len(LABELS)

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

TEST_BATCH_SIZE = int(config['test']['batch_size'])

# test params

TEST_PARAMS = {'planes': PLANES,
               'cells': CELLS,
               'views': VIEWS,
               'batch_size': TEST_BATCH_SIZE,
               'n_labels': N_LABELS,
               'labels': LABELS,
               'images_path': IMAGES_PATH,
               'shuffle': SHUFFLE}


'''
****************************************
*************** DATASETS ***************
****************************************
'''

logging.info('Reading datasets from serialized files...')

partition_file = open(DATASET_PATH + PARTITION_PREFIX + '.p', 'r')
partition = pickle.load(partition_file)
partition_file.close()

labels_file = open(DATASET_PATH + LABELS_PREFIX + '.p', 'r')
labels = pickle.load(labels_file)
labels_file.close()

logging.info('Number of training examples: %d', len(partition['train']))
logging.info('Number of validation examples: %d', len(partition['validation']))
logging.info('Number of test examples: %d', len(partition['test']))


'''
****************************************
************** GENERATORS **************
****************************************
'''

test_generator = DataGenerator(**TEST_PARAMS).generate(labels, partition['test'], True)


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

logging.info('STARTING TEST...')

# Test model

score = model.evaluate_generator(generator = test_generator,
                                 steps = len(partition['test'])//TEST_BATCH_SIZE
                                )

print('%s: %.2f%%' % (model.metrics_names[1], score[1]*100))
