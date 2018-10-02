"""
This is the train module.
"""

__version__ = '1.0'
__author__ = 'Saul Alonso-Monsalve'
__email__ = "saul.alonso.monsalve@cern.ch"

import tensorflow as tf
import numpy as np
import pickle
import configparser
import ast
import logging
import os
import sys
import re
import time

# manually specify the GPUs to use
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"]="0,1,2,3" 

sys.path.append(os.path.join(sys.path[0], 'modules'))
sys.path.append(os.path.join(sys.path[0], 'networks'))

from keras import backend as K, regularizers, optimizers
from keras.utils import multi_gpu_model
from keras.models import Model, Sequential, load_model
from keras.layers import Input, Dense, Activation, ZeroPadding2D, Dropout, Flatten, BatchNormalization, SeparableConv2D
from keras.callbacks import LearningRateScheduler, ReduceLROnPlateau, CSVLogger, ModelCheckpoint, EarlyStopping
from keras.layers.convolutional import Conv2D, MaxPooling2D, AveragePooling2D
from keras.regularizers import l2
from collections import Counter
from sklearn.utils import class_weight
from data_generator import DataGenerator
import networks
import my_losses
import my_callbacks

sess = tf.Session()
init = tf.global_variables_initializer()
sess.run(init)
K.set_session(sess)
K.set_image_data_format('channels_last')

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
    SEED = int(time.time()) # random seed   

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

# log

LOG_PATH = config['log']['path']
LOG_PREFIX = config['log']['prefix']

# model

ARCHITECTURE = config['model']['architecture']
CHECKPOINT_PATH = config['model']['checkpoint_path']
CHECKPOINT_PREFIX = config['model']['checkpoint_prefix']
CHECKPOINT_SAVE_MANY = ast.literal_eval(config['model']['checkpoint_save_many'])
CHECKPOINT_SAVE_BEST_ONLY = ast.literal_eval(config['model']['checkpoint_save_best_only'])
CHECKPOINT_PERIOD = int(config['model']['checkpoint_period'])
PARALLELIZE = ast.literal_eval(config['model']['parallelize'])
GPUS = int(config['model']['gpus'])
PRINT_SUMMARY = ast.literal_eval(config['model']['print_summary'])
BRANCHES = ast.literal_eval(config['model']['branches'])
OUTPUTS = int(config['model']['outputs'])

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
MAX_QUEUE_SIZE = int(config['train']['max_queue_size'])

# validation

VALIDATION_FRACTION = float(config['validation']['fraction'])
VALIDATION_BATCH_SIZE = int(config['validation']['batch_size'])

# train params

TRAIN_PARAMS = {'planes':PLANES,
                'cells':CELLS,
                'views':VIEWS,
                'batch_size':TRAIN_BATCH_SIZE,
                'branches':BRANCHES,
                'outputs': OUTPUTS,
                'images_path':IMAGES_PATH,
                'standardize':STANDARDIZE,
                'shuffle':SHUFFLE}

# validation params

VALIDATION_PARAMS = {'planes':PLANES,
                     'cells':CELLS,
                     'views':VIEWS,
                     'batch_size':VALIDATION_BATCH_SIZE,
                     'branches':BRANCHES,
                     'outputs': OUTPUTS,
                     'images_path':IMAGES_PATH,
                     'standardize':STANDARDIZE,
                     'shuffle':SHUFFLE}


'''
****************************************
*************** DATASETS ***************
****************************************
'''

partition = {'train' : [], 'validation' : [], 'test' : []} # Train, validation, and test IDs
labels = {}                                                # ID : label

# Load datasets

logging.info('Loading datasets from serialized files...')

with open(DATASET_PATH + PARTITION_PREFIX + '.p', 'r') as partition_file:
    partition = pickle.load(partition_file)

with open(DATASET_PATH + LABELS_PREFIX + '.p', 'r') as labels_file:
    labels = pickle.load(labels_file)

if WEIGHTED_LOSS_FUNCTION:
    with open(DATASET_PATH + CLASS_WEIGHTS_PREFIX + '.p', 'r') as class_weights_file:
        class_weights = pickle.load(class_weights_file)
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

# Optimizer: Stochastic Gradient Descent

logging.info('Setting optimizer...')

opt = optimizers.SGD(lr=LEARNING_RATE, momentum=MOMENTUM, decay=DECAY, nesterov=True) # SGD
#opt = optimizers.RMSprop(lr=0.045, rho=0.9, decay=0.94, epsilon=None, clipnorm=2.0)    # RMSprop (rho = Decay factor, decay = Learning rate decay over each update)
#opt = optimizers.Adam(lr=1e-3)							       # Adam

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
                filename = CHECKPOINT_PATH + '/' + fil
                sequential_model = load_model(filename, 
                                              custom_objects={'tf':tf, 
                                                              'masked_loss':my_losses.masked_loss, 
                                                              'multitask_loss':my_losses.multitask_loss,
                                                              'masked_loss_binary':my_losses.masked_loss_binary, 
                                                              'masked_loss_categorical':my_losses.masked_loss_categorical})
                logging.info('Loaded model: %s', CHECKPOINT_PATH + '/' + fil)                
                break
    else:
        # Load the model
        filename = CHECKPOINT_PATH + CHECKPOINT_PREFIX + '.h5'
        sequential_model = load_model(filename, custom_objects={'masked_loss':my_losses.masked_loss, 'my_losses':my_losses.my_losses})

        logging.info('Loaded model: %s', CHECKPOINT_PATH + CHECKPOINT_PREFIX + '.h5') 
else:
    # Start a new training
    logging.info('Creating model...')

    # Input shape
    if BRANCHES:
        input_shape = [PLANES, CELLS, 1]     # Input: VIEWS x [PLANES x CELLS x 1]
    else:
        input_shape = [PLANES, CELLS, VIEWS] # Input: [PLANES x CELLS x VIEWS]

    '''
    Create model. Options (argument 'network'):

                                                   Total params  Trainable params  Non-trainable params

    'xception'            -> Xception                20,888,117        20,833,589                54,528
    'vgg16'               -> VGG-16                 503,412,557       503,412,557                     0
    'vgg19'               -> VGG-19                 508,722,253       508,722,253                     0
    'resnet18'            -> ResNet-18               11,193,997        11,186,189                 7,808
    'resnet34'            -> ResNet-34               21,313,293        21,298,061                15,232
    'resnet50'            -> ResNet-50               23,694,221        23,641,101                53,120
    'resnet101'           -> ResNet-101              42,669,453        42,571,789                97,664
    'resnet152'           -> ResNet-152              58,382,221        58,238,477               143,744
    'inceptionv3'         -> Inception-v3            21,829,421        21,794,989                34,432
    'inceptionv4'         -> Inception-v4            41,194,381        41,131,213                63,168
    'inceptionresnetv2'   -> Inception-ResNet-v2     54,356,717        54,296,173                60,544
    'resnext'             -> ResNeXt                 89,712,576        89,612,096               100,480
    'seresnet18'          -> SE-ResNet-18            11,276,224        11,268,416                 7,808
    'seresnet34'          -> SE-ResNet-34            21,461,952        21,446,720                15,232
    'seresnet50'          -> SE-ResNet-50            26,087,360        26,041,920                45,440
    'seresnet101'         -> SE-ResNet-101           47,988,672        47,887,936               100,736
    'seresnet154'         -> SE-ResNet-154           64,884,672        64,740,928               143,744
    'seresnetsaul'        -> SE-ResNet-Saul          22,072,768        22,055,744                17,024
    'seinceptionv3'       -> SE-Inception-v3         23,480,365        23,445,933                34,432
    'seinceptionresnetv2' -> SE-Inception-ResNet-v2  64,094,445        64,033,901                60,544
    'seresnext'           -> SE-ResNeXt              97,869,632        97,869,632               100,480
    'mobilenet'           -> MobileNet                3,242,189         3,220,301                21,888
    'densenet121'         -> DenseNet-121             7,050,829         6,967,181                83,648
    'densenet169'         -> DenseNet-169            12,664,525        12,506,125               158,400
    'densenet201'         -> DenseNet-201            18,346,957        18,117,901               229,056
    other                 -> Custom model
    '''

    aux_model = networks.create_model(network=ARCHITECTURE, input_shape=input_shape)
    aux_model.layers.pop() # remove the last layer of the model

    weight_decay = 1e-4

    x = [None]*OUTPUTS

    if OUTPUTS == 1:
        x[0] = Dense(13, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='softmax', name='categories')(aux_model.layers[-1].output)
    elif OUTPUTS == 5:
        x[0] = Dense(4, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='softmax', name='flavour')(aux_model.layers[-1].output)
        x[1] = Dense(4, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='softmax', name='protons')(aux_model.layers[-1].output)
        x[2] = Dense(4, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='softmax', name='pions')(aux_model.layers[-1].output)
        x[3] = Dense(4, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='softmax', name='pizeros')(aux_model.layers[-1].output)
        x[4] = Dense(4, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='softmax', name='neutrons')(aux_model.layers[-1].output)
    else:
        x[0] = Dense(1, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='sigmoid', name='is_antineutrino')(aux_model.layers[-1].output)
        x[1] = Dense(4, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='softmax', name='flavour')(aux_model.layers[-1].output)
        x[2] = Dense(4, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='softmax', name='interaction')(aux_model.layers[-1].output)
        x[3] = Dense(4, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='softmax', name='protons')(aux_model.layers[-1].output)
        x[4] = Dense(4, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='softmax', name='pions')(aux_model.layers[-1].output)
        x[5] = Dense(4, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='softmax', name='pizeros')(aux_model.layers[-1].output)
        x[6] = Dense(4, use_bias=False, kernel_regularizer=l2(weight_decay),
                     activation='softmax', name='neutrons')(aux_model.layers[-1].output)

    sequential_model = Model(inputs=aux_model.inputs, outputs=x, name='resnext')

if PARALLELIZE:
    # Parallelize the model (use all the available GPUs)
    try:
        model = multi_gpu_model(sequential_model, gpus=GPUS, cpu_relocation=True)
        num_outputs = len(sequential_model.output_names)
        model.layers[-(num_outputs+1)].set_weights(sequential_model.get_weights()) # set the same weights as the sequential model
        logging.info('Training using %d GPUs...', GPUS)
    except:
        model = sequential_model
        logging.info('Training using single GPU or CPU...')
else:
    model = sequential_model

if not model._is_compiled:
    # Compile model
    logging.info('Compiling model...')

    if OUTPUTS == 1:
        model_loss = {'categories':my_losses.masked_loss_categorical}
    elif OUTPUTS == 5:
        model_loss = {'flavour':my_losses.masked_loss_categorical,
                      'protons':my_losses.masked_loss_categorical,
                      'pions':my_losses.masked_loss_categorical,
                      'pizeros':my_losses.masked_loss_categorical,
                      'neutrons':my_losses.masked_loss_categorical}
    else:
        model_loss = {'is_antineutrino':my_losses.masked_loss_binary,
                      'flavour':my_losses.masked_loss_categorical,
                      'interaction':my_losses.masked_loss_categorical,
                      'protons':my_losses.masked_loss_categorical,
                      'pions':my_losses.masked_loss_categorical,
                      'pizeros':my_losses.masked_loss_categorical,
                      'neutrons':my_losses.masked_loss_categorical}

    #model.compile(loss=my_losses.masked_loss, optimizer=opt, metrics=['accuracy']) 
    model.compile(#loss={'neutrino': my_losses.masked_loss_binary, 'flavour': my_losses.masked_loss_categorical, 'interaction': my_losses.masked_loss_categorical}, 
                  loss=model_loss,
                  #loss_weights={'neutrino':0.25, 'flavour':1.0, 'interaction': 0.5},
                  #loss_weights={'neutrino':0.33, 'flavour':0.33, 'interaction':0.33},
                  optimizer=opt, 
                  metrics=['accuracy']) 
    #model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])

if(PRINT_SUMMARY):
    # Print model summary
    if PARALLELIZE:
        sequential_model.summary()
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
    if OUTPUTS == 1:
        # Validation accuracy
        if CHECKPOINT_SAVE_MANY:
            filepath = CHECKPOINT_PATH + CHECKPOINT_PREFIX + '-{epoch:02d}-{val_acc:.2f}.h5'
        monitor_acc = 'val_acc'
        monitor_loss = 'val_loss'
    else:
        # Validation accuracy
        if CHECKPOINT_SAVE_MANY:
            filepath = CHECKPOINT_PATH + CHECKPOINT_PREFIX + '-{epoch:02d}-{val_flavour_acc:.2f}.h5'
        monitor_acc = 'val_flavour_acc'
        monitor_loss = 'val_flavour_loss'
else:
    # Training accuracy
    if CHECKPOINT_SAVE_MANY:
        filepath = CHECKPOINT_PATH + CHECKPOINT_PREFIX + '-{epoch:02d}-{acc:.2f}.h5'
    monitor_acc = 'acc'
    monitor_loss = 'loss'

if PARALLELIZE:
    #checkpoint = my_callbacks.MultiGPUCheckpointCallback(filepath, sequential_model, monitor=monitor_acc, verbose=1, save_best_only=CHECKPOINT_SAVE_BEST_ONLY, mode='max', period=CHECKPOINT_PERIOD)
    checkpoint = my_callbacks.ModelCheckpointDetached(filepath, monitor=monitor_acc, verbose=1, save_best_only=CHECKPOINT_SAVE_BEST_ONLY, save_weights_only=False, mode='max', period=CHECKPOINT_PERIOD)
else:
    checkpoint = ModelCheckpoint(filepath, monitor=monitor_acc, verbose=1, save_best_only=CHECKPOINT_SAVE_BEST_ONLY, mode='max', period=CHECKPOINT_PERIOD)

# Learning rate reducer

logging.info('Configuring learning rate reducer...')

#lr_reducer = LearningRateScheduler(schedule=lambda epoch,lr: (lr*0.01 if epoch % 2 == 0 else lr))
#lr_reducer = ReduceLROnPlateau(monitor=monitor_loss, factor=0.1, cooldown=0, patience=2, min_lr=0.5e-6, verbose=1)
lr_reducer = ReduceLROnPlateau(monitor=monitor_acc, mode='max', factor=0.1, cooldown=0, patience=10, min_lr=0.5e-6, verbose=1)

# Early stopping

logging.info('Configuring early stopping...')

early_stopping = EarlyStopping(monitor=monitor_acc, patience=EARLY_STOPPING_PATIENCE, mode='auto')
  
# Configuring log file

csv_logger = CSVLogger(LOG_PATH + LOG_PREFIX + '.log', append=RESUME)

# My callbacks

my_callback = my_callbacks.MyCallback()
#my_callback = my_callbacks.InceptionV4Callback()
#my_callback = my_callbacks.IterationsCallback(validation_generator=validation_generator, validation_steps=len(partition['validation'])//VALIDATION_BATCH_SIZE)

# Callbacks

logging.info('Setting callbacks...')

#callbacks_list = [checkpoint, csv_logger]
#callbacks_list = [checkpoint, csv_logger, my_callback]
#callbacks_list = [lr_reducer, checkpoint, early_stopping, csv_logger]
callbacks_list = [lr_reducer, checkpoint, early_stopping, csv_logger, my_callback]


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
    validation_data = validation_generator
    validation_steps = len(partition['validation'])//VALIDATION_BATCH_SIZE
else:
    # Training without validation
    validation_data = None
    validation_steps = None

# TRAINING

model.fit_generator(generator=training_generator,
                    steps_per_epoch=len(partition['train'])//TRAIN_BATCH_SIZE,
                    validation_data=validation_data,
                    validation_steps=validation_steps,
                    epochs=EPOCHS,
                    class_weight=class_weights,
                    callbacks=callbacks_list,
                    initial_epoch=initial_epoch,
                    max_queue_size=MAX_QUEUE_SIZE,
                    verbose=1,
                    use_multiprocessing=False,
                    workers=1
                   )

