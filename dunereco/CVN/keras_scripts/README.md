# Documentation

* [1. Configuration file](#1-configuration-file)
 * [[random]](#1-1-random)
 * [[images]](#1-2-images)
 * [[dataset]](#1-3-dataset)
 * [[log]](#1-4-log)
 * [[model]](#1-5-model)
 * [[train]](#1-6-train)
 * [[validation]](#1-7-validation)
 * [[test]](#1-8-test)
* [2. Dataset generator](#2-dataset-generator)
* [3. Training](#3-training)
* [4. Test](#4-test)
 
<br />

## 1. Configuration file

First, you need to specify all the configuration parameters in the *config.ini* file (you can always regenerate the configuration file by executing the Python script *generate_config.py*):

<br />

<a name="1-1-random"></a>
### [random]

#### - seed

> seed = RANDOM_SEED

**Description:** value used to initialize a pseudo-random number generator (used for reproducibility).

**Type:** Integer. 

**Example:** `seed = 7`

#### - shuffle

> shuffle = BOOLEAN

**Description:** value used (if True) to randomly shuffle the dataset before training/validating/testing.

**Type:** Boolean. 

**Example:** `shuffle = True`

<br />

<a name="1-2-images"></a>
### [images]

#### - filtered

> filtered = BOOLEAN

**Description:** value used (if True) to specify that the input dataset will consist of filtered images.

**Type:** Boolean.

**Example:** `filtered = True`

*NOTE:* it must only be True if the user already filtered the images.

#### - standardize

> standardize = BOOLEAN

**Description:** value used (if True) to normalize the pixel values (floating-point values from 0 to 1).

**Type:** Boolean.

**Example:** `standardize = True`

#### - views 

> views = NUMBER_OF_VIEWS

**Description:** number of views of each event (it should be 3).

**Type:** Integer.

**Example:** `views = 3`

#### - cells

> cells = NUMBER_OF_CELLS

**Description:** number of cells of each image (currently, the images shape is 500 planes x 500 cells).

**Type:** Integer.

**Example:** `cells = 500`

#### - interaction_labels

> interaction_labels = INTERACTION_TYPE_LABELS

**Description:** correspondence between each interaction type label and a number from 0 to N-1 (being N the total number of interaction type labels).

**Type:** Dictionary.

**Example:** `interaction_labels = {'0': 0, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, '10': 10, '11': 11, '13': 12}`

#### - neutrino_labels

> neutrino_labels = NEUTRINO_FLAVOUR_LABELS

**Description:** correspondence between each interaction type label and a number representing a neutrino flavour (e.g., 0 = numu, 1 = nue, 2 = nutau, 3 = other).

**Type:** Dictionary.

**Example:** `neutrino_labels = {'0': 1, '1': 1, '2': 1, '3': 1, '4': 0, '5': 0, '6': 0, '7': 0, '8': 2, '9': 2, '10': 2, '11': 2, '13': 3}`

#### - planes

> planes = NUMBER_OF_PLANES

**Description:** number of planes of each image (currently, the images shape is 500 planes x 500 cells).

**Type:** Integer.

**Example:** `planes = 500`

#### - path

> path = IMAGES_PATH

**Description:** main path of the compressed (using zlib) images. The structure of the directory should be the following:

```
IMAGES_PATH
│
└───LABEL_1
│        NEUTRINO.listNUMBER.txt.gz
│        NEUTRINO.listNUMBER.txt.gz
│        ...
│   
└───LABEL_2
│        NEUTRINO.listNUMBER.txt.gz
│        NEUTRINO.listNUMBER.txt.gz
│        ...
│ 
│   ... 
│ 
│
└───LABEL_N
         NEUTRINO.listNUMBER.txt.gz
         NEUTRINO.listNUMBER.txt.gz
         ...
```

Where IMAGES_PATH corresponds to the previously described path. That path consists of 13 folders, each one named LABEL_X, (that corresponds with an interaction type label), each of those storing the compressed images (where NEUTRINO corresponds with one of nue, numu, nutau, anue, anumu, and anutau; and NUMBER corresponds with an integer number). For example, the name of an image compressed file would be *anumu.list19746.txt.gz*.

**Type:** String.

**Example:** `path = /scratch/devFilesRaw`

<br />

<a name="1-3-dataset"></a>
### [dataset]

> path = DATASET_PATH

#### - path

**Description:** path where the generated training, validation, and test sequences will be stored after executing the file *generate_dataset.py*.

**Type:** String.

**Example:** `path = /scratch/cvn/dataset`

#### - labels_prefix

> labels_prefix = LABELS_PREFIX

**Description:** represents the name prefix the (pickled) generated labels file will have (the full name would be LABELS_PREFIX + '.p') after executing the file *generate_dataset.py*.

**Type:** String.

**Example:** `labels_prefix = /labels`

#### - interaction_types

> interaction_types = BOOLEAN

**Description:** value used (if True) to use either the interaction type labels or (if False) the neutrino flavour labels.

**Type:** Boolean.

**Example:** `interaction_types = False`

#### - partition_prefix

> partition_prefix = PARTITION_PREFIX

**Description:** represents the name prefix the (pickled) dictionary file containing the generated training, validation, and test sequences will have (the full name would be PARTITION_PREFIX + '.p') after executing the file *generate_dataset.py*.

**Type:** String.

**Example:** `partition_prefix = /partition`

#### - uniform

> uniform = BOOLEAN

**Description:** value used (if True) to perform rejection sampling in order to use the same number of samples per label.

**Type:** Boolean.

**Example:** `uniform = True`

<br />

<a name="1-4-log"></a>
### [log]

> path = LOG_PATH

#### - path

**Description:** path where the training and validation outputs will be logged while training.

**Type:** String.

**Example:** `path = /scratch/cvn/log`

#### - prefix

> prefix = LOG_PREFIX

**Description:** represents the name prefix the log file will have (the full name would be LOG_PREFIX + '.log').

**Type:** String.

**Example:** `prefix = /output`

<br />

<a name="1-5-model"></a>
### [model]

#### - checkpoint_path

> checkpoint_path = CHECKPOINT_PATH

**Description:** path where the model (the network architecture and its parameters) will be stored while training.

**Type:** String.

**Example:** `checkpoint_path = /scratch/cvn/checkpoint`

#### - checkpoint_save_many

> checkpoint_save_many = BOOLEAN

**Description:** value used (if True) to save several model files while training (e.g., save the model always in a different file at the end of each epoch).

**Type:** Boolean.

**Example:** `checkpoint_save_many = False`

#### - checkpoint_save_best_only

> checkpoint_save_best_only = BOOLEAN

**Description:** value used (if True) to only overwrite the latest stored model if the current one improves the quantity monitored. That means only one model will be stored after the training.

**Type:** Boolean.

**Example:** `checkpoint_save_best_only = True`

#### - checkpoint_prefix

> checkpoint_prefix = CHECKPOINT_PREFIX

**Description:** represents the name prefix the checkpoint file(s) will have (the full name would be CHECKPOINT_PREFIX [+ '-EPOCH-MONITORED_METRIC_VALUE'] + '.h5').

**Type:** String.

**Example:** `checkpoint_prefix = /model`

#### - print_summary

> print_summary = BOOLEAN

**Description:** value used (if True) to print a summary of the model before performing either the training or testing. An example of the summary of a model would be the following:

```
_________________________________________________________________
Layer (type)                 Output Shape              Param #
=================================================================
conv2d_1 (Conv2D)            (None, 64, 125, 125)      23296
_________________________________________________________________
max_pooling2d_1 (MaxPooling2 (None, 64, 31, 31)        0
_________________________________________________________________
batch_normalization_1 (Batch (None, 64, 31, 31)        256
_________________________________________________________________
conv2d_2 (Conv2D)            (None, 32, 16, 16)        100384
_________________________________________________________________
max_pooling2d_2 (MaxPooling2 (None, 32, 15, 15)        0
_________________________________________________________________
batch_normalization_2 (Batch (None, 32, 15, 15)        128
_________________________________________________________________
dropout_1 (Dropout)          (None, 32, 15, 15)        0
_________________________________________________________________
flatten_1 (Flatten)          (None, 7200)              0
_________________________________________________________________
dense_1 (Dense)              (None, 1000)              7201000
_________________________________________________________________
dropout_2 (Dropout)          (None, 1000)              0
_________________________________________________________________
dense_2 (Dense)              (None, 1000)              1001000
_________________________________________________________________
dropout_3 (Dropout)          (None, 1000)              0
_________________________________________________________________
dense_3 (Dense)              (None, 4)                 4004
=================================================================
Total params: 8,330,068
Trainable params: 8,329,876
Non-trainable params: 192

```

**Type:** Boolean.

**Example:** `print_summary = True`

#### - checkpoint_period

> checkpoint_period = CHECKPOINT_PERIOD

**Description:** interval (number of epochs) between checkpoints.

**Type:** Integer.

**Example:** `checkpoint_period = 1`

<br />

<a name="1-6-train"></a>
### [train]

#### - epochs

> epochs = TRAINING_EPOCHS

**Description:** number of training epochs.

**Type:** Integer.

**Example:** `epochs = 100`

*NOTE:* please, do not mix up the terms iteration and epoch: 

- An iteration one pass using [mini-batch size] number of examples. Thus, the number of iterations consists of the number of passes, each pass using [mini-batch size] number of examples. To be clear, one pass = one forward pass + one backward pass (we do not count the forward pass and backward pass as two different passes). 
- An epoch is an arbitrary cutoff, generally defined as "one pass over the entire dataset" (one forward pass and one backward pass of all the training examples), used to separate training into distinct phases, which is useful for logging and periodic evaluation.

Some Deep Learning frameworks (e.g., Caffe) do not use the term epoch; however, Tensorflow does.

#### - lr

> lr = LEARNING_RATE

**Description:** learning rate of the SGD (Stochastic Gradient Descent) algorithm.

**Type:** Floating-point.

**Example:** `lr = 0.001`

#### - batch_size

> batch_size = TRAINING_MINIBATCH_SIZE

**Description:** training mini-batch size.

**Type:** Integer.

**Example:** `batch_size = 8`

*NOTE:* a mini-batch consist of a set of N samples. The samples in a mini-batch are processed independently, in parallel. If training, a mini-batch results in only one update to the model (one forward/backward pass).

#### - fraction

> fraction = TRAINING_FRACTION

**Description:** fraction of the total amount of data used for training.

**Type:** Floating-point.

**Example:** `fraction = 0.6`

*NOTE:* TRAINING_FRACTION + VALIDATION_FRACTION + TEST_FRACTION should be less or equal to 1.

#### - resume

> resume = BOOLEAN

**Description:** value used (if True) to resume a previous training.

**Type:** Boolean.

**Example:** `resume = False`

#### - weighted_loss_function

> weighted_loss_function = BOOLEAN

**Description:** value used (if True) to use a (balanced) weighted loss function.

**Type:** Boolean.

**Example:** `weighted_loss_function = False`

#### - class_weights_prefix

> class_weights_prefix = CLASS_WEIGHTS_PREFIX

**Description:** represents the name prefix the (pickled) generated file that stores the balanced class weights will have (the full name would be CLASS_WEIGHTS_PREFIX + '.p').

**Type:** String.

**Example:** `class_weights_prefix = /class_weights_aux4`

#### - decay

> decay = DECAY

**Description:** learning rate decay over each update.

**Type:** Floating-point.

**Example:** `decay = 0.000005`

#### - early_stopping_patience

> early_stopping_patience = EARLY_STOPPING_PATIENCE

**Description:** number of epochs with no improvement after which training will be stopped.

**Type:** Integer.

**Example:** `early_stopping_patience = 5`
*NOTE:* if the user does not want to use Early Stopping, just set early_stopping_patience = epochs

#### - momentum

> momentum = MOMENTUM

**Description:** parameter that accelerates SGD (Stochastic Gradient Descent) in the relevant direction and dampens oscillations.

**Type:** Floating-point.

**Example:** `momentum = 0.9`

<br />

<a name="1-7-validation"></a>
### [validation]

#### - fraction

> fraction = VALIDATION_FRACTION

**Description:** fraction of the total amount of data used for validation.

**Type:** Floating-point.

**Example:** `fraction = 0.2`

*NOTE:* TRAINING_FRACTION + VALIDATION_FRACTION + TEST_FRACTION must be less or equal to 1. If VALIDATION_FRACTION = 0, the system does not compute validation.

#### - batch_size

> batch_size = VALIDATION_MINIBATCH_SIZE

**Description:** validation mini-batch size.

**Type:** Integer.

**Example:** `batch_size = 32`

*NOTE:* a mini-batch consist of a set of N samples. The samples in a mini-batch are processed independently, in parallel.

<br />

<a name="1-8-test"></a>
### [test]

#### - cut_nue

> cut_nue = CUT_NUE

**Description:** nue cut to apply during the analysis.

**Type:** Floating-point.

**Example:** `cut_nue = 0.7`


#### - cut_numu

> cut_numu = CUT_NUMU

**Description:** numu cut to apply during the analysis.

**Type:** Floating-point.

**Example:** `cut_numu = 0.7`

#### - output_prefix

> output_prefix = OUTPUT_PREFIX

**Description:** represents the name prefix the NumPy output file will have (the full name would be OUTPUT_PREFIX + '.np').

**Type:** String.

**Example:** `prefix = /test_info'`

#### - output_path

> output_path = OUTPUT_PATH

**Description:** path where the test output NumPy object will be stored.

**Type:** String.

**Example:** `output_path = /scratch/cvn/output`

#### - fraction

> fraction = TEST_FRACTION

**Description:** fraction of the total amount of data used for test.

**Type:** Floating-point.

**Example:** `fraction = 0.2`

*NOTE:* TRAINING_FRACTION + VALIDATION_FRACTION + TEST_FRACTION must be less or equal to 1.

#### - batch_size

> batch_size = TEST_MINIBATCH_SIZE

**Description:** test mini-batch size.

**Type:** Integer.

**Example:** `batch_size = 32`

*NOTE:* a mini-batch consist of a set of N samples. The samples in a mini-batch are processed independently, in parallel.

<br />

## 2. Dataset generator

Once you specified all the parameters in the *config.ini* file, you should generate the datasets by executing the Python script *dataset_generator.py* as follows:

> ./execute.sh generate_dataset.py

It will create the partition and labels pickled files under the dataset directory specified in the configuration file.

<br />

## 3. Training

Once the datasets are created, you should start the training by executing the Python script *training.py* as follows:

> ./execute.sh training.py

<br />

## 4. Test

Once the training is finished, you should start the testing by executing the Python script *test.py* as follows:

> ./execute.sh test.py
