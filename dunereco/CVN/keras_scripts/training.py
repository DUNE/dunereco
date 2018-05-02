import numpy as np
import pickle
import configparser
import ast
import logging, sys
import re
import os

from sklearn.utils import class_weight
from keras.models import Model, Sequential, load_model
from keras.layers import Input, Dense, Activation, ZeroPadding2D, Dropout, Flatten, BatchNormalization, SeparableConv2D
from keras import regularizers, optimizers
from keras.layers.convolutional import Conv2D, MaxPooling2D, AveragePooling2D
from keras.callbacks import LearningRateScheduler, ReduceLROnPlateau, CSVLogger, ModelCheckpoint, EarlyStopping
from data_generator import DataGenerator
from collections import Counter

sys.path.append("/home/salonsom/cvn_tensorflow/networks")
sys.path.append("/home/salonsom/cvn_tensorflow/callbacks")

import se_resnet, resnet, resnetpa, googlenet, my_model
import my_callbacks

from keras import backend as K
K.set_image_data_format('channels_last')

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

RESUME = ast.literal_eval(config['train']['resume'])
LEARNING_RATE = float(config['train']['lr'])
MOMENTUM = float(config['train']['momentum'])
DECAY = float(config['train']['decay'])
TRAIN_BATCH_SIZE = int(config['train']['batch_size'])
EPOCHS = int(config['train']['epochs'])
EARLY_STOPPING_PATIENCE = int(config['train']['early_stopping_patience'])
WEIGHTED_LOSS_FUNCTION = ast.literal_eval(config['train']['weighted_loss_function'])
CLASS_WEIGHTS_PREFIX = config['train']['class_weights_prefix']

# validation

VALIDATION_FRACTION = float(config['validation']['fraction'])
VALIDATION_BATCH_SIZE = int(config['validation']['batch_size'])

# train params

TRAIN_PARAMS =      {'planes': PLANES,
                     'cells': CELLS,
                     'views': VIEWS,
                     'batch_size': TRAIN_BATCH_SIZE,
                     'n_labels': N_LABELS,
                     'interaction_labels': INTERACTION_LABELS,
                     'interaction_types': INTERACTION_TYPES,
                     'filtered': FILTERED,
                     'neutrino_labels': NEUTRINO_LABELS,
                     'images_path': IMAGES_PATH,
                     'standardize': STANDARDIZE,
                     'shuffle': SHUFFLE}

# validation params

VALIDATION_PARAMS = {'planes': PLANES,
                     'cells': CELLS,
                     'views': VIEWS,
                     'batch_size': VALIDATION_BATCH_SIZE,
                     'n_labels': N_LABELS,
                     'interaction_labels': INTERACTION_LABELS,
                     'interaction_types': INTERACTION_TYPES,
                     'filtered': FILTERED,
                     'neutrino_labels': NEUTRINO_LABELS,
                     'images_path': IMAGES_PATH,
                     'standardize': STANDARDIZE,
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

if WEIGHTED_LOSS_FUNCTION:

    class_weights_file = open(DATASET_PATH + CLASS_WEIGHTS_PREFIX + '.p', 'r')
    class_weights = pickle.load(class_weights_file)
    class_weights_file.close()

else:

    class_weights = None

# Print some dataset statistics

logging.info('Number of training examples: %d', len(partition['train']))
logging.info('Number of validation examples: %d', len(partition['validation']))
logging.info('Number of test examples: %d', len(partition['test']))
logging.info('Class weights: %s', class_weights)

'''
****************************************
************** GENERATORS **************
****************************************
'''

training_generator = DataGenerator(**TRAIN_PARAMS).generate(labels, partition['train'], True)
validation_generator = DataGenerator(**VALIDATION_PARAMS).generate(labels, partition['validation'], True)


'''
****************************************
*************** CVN MODEL **************
****************************************
'''

if RESUME:
 
    # Resume a previous training

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

else:

    # Start a new training

    logging.info('Creating model...')

    # Input shape: (PLANES x CELLS x VIEWS)

    input_shape = [PLANES, CELLS, VIEWS]

    #model = se_resnet.SEResNet(input_shape=input_shape, classes=N_LABELS)
    #model = resnetpa.ResNetPreAct()

    model = resnet.ResnetBuilder.build_resnet_18(input_shape, N_LABELS)
    #model = resnet.ResnetBuilder.build_resnet_34(input_shape, N_LABELS)
    #model = resnet.ResnetBuilder.build_resnet_50(input_shape, N_LABELS)
    #model = resnet.ResnetBuilder.build_resnet_101(input_shape, N_LABELS)
    #model = resnet.ResnetBuilder.build_resnet_152(input_shape, N_LABELS)
    #model = my_model.my_model(input_shape=input_shape, classes=N_LABELS)

    # Optimizer: Stochastic Gradient Descent

    logging.info('Setting optimizer...')

    opt = optimizers.SGD(lr=LEARNING_RATE, momentum=MOMENTUM, decay=DECAY, nesterov=True)

    # Compile model

    logging.info('Compiling model...')

    model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])

    '''
    model.compile(loss='categorical_crossentropy',
              optimizer='adam',
              metrics=['accuracy'])
    '''

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

    monitor_acc = 'val_acc'
    monitor_loss = 'val_loss'

else:

    # Training accuracy

    if CHECKPOINT_SAVE_MANY:
        filepath = CHECKPOINT_PATH + CHECKPOINT_PREFIX + '-{epoch:02d}-{acc:.2f}.h5'

    monitor_acc = 'acc'
    monitor_loss = 'loss'

checkpoint = ModelCheckpoint(filepath, monitor=monitor_acc, verbose=1, save_best_only=CHECKPOINT_SAVE_BEST_ONLY, mode='max', period=CHECKPOINT_PERIOD)

# Learning rate reducer

logging.info('Configuring learning rate reducer...')

#lr_reducer = LearningRateScheduler(schedule=lambda epoch,lr: (lr*0.01 if epoch % 2 == 0 else lr))
lr_reducer = ReduceLROnPlateau(monitor=monitor_loss, factor=0.1, cooldown=0, patience=3, min_lr=0.5e-6, verbose=1)

# Early stopping

logging.info('Configuring early stopping...')

early_stopping = EarlyStopping(monitor=monitor_acc, patience=EARLY_STOPPING_PATIENCE, mode='auto')
  
# Configuring log file

csv_logger = CSVLogger(LOG_PATH + LOG_PREFIX + '.log', append=RESUME)

# My callbacks

my_callback = my_callbacks.MyCallback()

# Callbacks

logging.info('Setting callbacks...')

callbacks_list = [lr_reducer, checkpoint, early_stopping, csv_logger]
#callbacks_list = [lr_reducer, checkpoint, early_stopping, csv_logger, my_callback]


'''
****************************************
*************** TRAINING ***************
****************************************
'''

if RESUME:

    # Resuming training...

    try:

        # Open previous log file in order to get the last epoch

        with open(LOG_PATH + LOG_PREFIX + '.log', 'r') as logfile:

            # initial_epoch = last_epoch + 1

            initial_epoch = int(re.search(r'\d+', logfile.read().split('\n')[-2]).group()) + 1

    except IOError:

        # Previous log file does not exist. Set initial epoch to 0    
    
        initial_epoch = 0

    logging.info('RESUMING TRAINING...')

else:

    # Starting a new training...
    # initial_epoch must be 0 when starting a training (not resuming it)

    initial_epoch = 0

    logging.info('STARTING TRAINING...')

if VALIDATION_FRACTION > 0:

    # Training with validation

    model.fit_generator(generator = training_generator,
                        steps_per_epoch = len(partition['train'])//TRAIN_BATCH_SIZE,
                        validation_data = validation_generator,
                        validation_steps = len(partition['validation'])//VALIDATION_BATCH_SIZE,
                        epochs = EPOCHS,
                        class_weight = class_weights,
                        callbacks = callbacks_list,
                        initial_epoch = initial_epoch,
                        verbose = 1
                       )

else:

    # Training without validation

    model.fit_generator(generator = training_generator,
                        steps_per_epoch = len(partition['train'])//TRAIN_BATCH_SIZE,
                        epochs = EPOCHS,
                        class_weight = class_weights,
                        callbacks = callbacks_list,
                        initial_epoch = initial_epoch,
                        verbose = 1
                       )


