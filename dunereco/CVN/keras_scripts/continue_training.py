import numpy as np
import pickle
import configparser
import ast
import logging, sys

from keras.models import Sequential, load_model
from keras.layers import Dense, Dropout, Flatten, BatchNormalization
from keras import optimizers
from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.callbacks import CSVLogger, ModelCheckpoint, EarlyStopping
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

# log

LOG_PATH = config['log']['path']
LOG_PREFIX = config['log']['prefix']

# model

CHECKPOINT_PATH = config['model']['checkpoint_path']
CHECKPOINT_PREFIX = config['model']['checkpoint_prefix']
CHECKPOINT_SAVE_MANY = ast.literal_eval(config['model']['checkpoint_save_many'])
CHECKPOINT_SAVE_BEST_ONLY = ast.literal_eval(config['model']['checkpoint_save_best_only'])
CHECKPOINT_PERIOD = int(config['model']['checkpoint_period'])
PRINT_SUMMARY = ast.literal_eval(config['model']['print_summary'])

# train

LEARNING_RATE = float(config['train']['lr'])
MOMENTUM = float(config['train']['momentum'])
DECAY = float(config['train']['decay'])
TRAIN_BATCH_SIZE = int(config['train']['batch_size'])
EPOCHS = int(config['train']['epochs'])
EARLY_STOPPING_PATIENCE = int(config['train']['early_stopping_patience'])

# validation

VALIDATION_FRACTION = float(config['validation']['fraction'])
VALIDATION_BATCH_SIZE = int(config['validation']['batch_size'])

# train params

TRAIN_PARAMS =      {'planes': PLANES,
                     'cells': CELLS,
                     'views': VIEWS,
                     'batch_size': TRAIN_BATCH_SIZE,
                     'n_labels': N_LABELS,
                     'labels': LABELS,
                     'images_path': IMAGES_PATH,
                     'shuffle': SHUFFLE}

# validation params

VALIDATION_PARAMS = {'planes': PLANES,
                     'cells': CELLS,
                     'views': VIEWS,
                     'batch_size': VALIDATION_BATCH_SIZE,
                     'n_labels': N_LABELS,
                     'labels': LABELS,
                     'images_path': IMAGES_PATH,
                     'shuffle': SHUFFLE}


'''
****************************************
*************** DATASETS ***************
****************************************
'''

partition = {'train' : [], 'validation' : [], 'test' : []} # Train, validation, and test IDs
labels = {}                                                # ID : label

# Load datasets

logging.info('Loading datasets from serialized files...')

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

training_generator = DataGenerator(**TRAIN_PARAMS).generate(labels, partition['train'], True)
validation_generator = DataGenerator(**VALIDATION_PARAMS).generate(labels, partition['validation'], True)

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
*************** CALLBACKS **************
****************************************
'''

# Checkpointing

logging.info('Configuring checkpointing...')

# Checkpoint one CVN model only

filepath = CHECKPOINT_PATH + CHECKPOINT_PREFIX + '.h5'

if VALIDATION_FRACTION > 0:

    # Validation accuracy

    if CHECKPOINT_SAVE_MANY:
        filepath = CHECKPOINT_PATH + CHECKPOINT_PREFIX + '-{epoch:02d}-{val_acc:.2f}.h5'

    monitor = 'val_acc'

else:

    # Training accuracy

    if CHECKPOINT_SAVE_MANY:
        filepath = CHECKPOINT_PATH + CHECKPOINT_PREFIX + '-{epoch:02d}-{acc:.2f}.h5'

    monitor = 'acc' 

checkpoint = ModelCheckpoint(filepath, monitor=monitor, verbose=1, save_best_only=CHECKPOINT_SAVE_BEST_ONLY, mode='max', period=CHECKPOINT_PERIOD)


# Early stopping

logging.info('Configuring early stopping...')

early_stopping = EarlyStopping(monitor=monitor, patience=EARLY_STOPPING_PATIENCE, mode='auto')
  
# Configuring log file

csv_logger = CSVLogger(LOG_PATH + LOG_PREFIX + '.log', append=True)
 
# Callbacks

logging.info('Setting callbacks...')

callbacks_list = [checkpoint, early_stopping, csv_logger]


'''
****************************************
*************** TRAINING ***************
****************************************
'''

logging.info('STARTING TRAINING...')

if VALIDATION_FRACTION > 0:

    # Training with validation

    model.fit_generator(generator = training_generator,
                        steps_per_epoch = len(partition['train'])//TRAIN_BATCH_SIZE,
                        validation_data = validation_generator,
                        validation_steps = len(partition['validation'])//VALIDATION_BATCH_SIZE,
                        epochs = EPOCHS,
                        callbacks = callbacks_list
                       )

else:

    # Training without validation

    model.fit_generator(generator = training_generator,
                        steps_per_epoch = len(partition['train'])//TRAIN_BATCH_SIZE,
                        epochs = EPOCHS,
                        callbacks = callbacks_list
                       )


