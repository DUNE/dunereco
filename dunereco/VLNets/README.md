# VLNets -- Variable Length Array Neural Networks

This module is designed to evaluate neural networks that operate on variable
length array inputs.

At this moment, `VLNets` module supports evaluation of the RNN based neutrino
energy estimator. The RNN energy estimator uses information from PFParticles
identified by Pandora, and tries to predict neutrino energy with a help a of
Recurrent Neural Network.


## Module Overview

Each submodule of `VLNets` contains a README file with a short description of
the submodule structure and possibly some additional notes. This README file
contains description of the `VLNets` module structure, overview of the RNN
energy estimator training and evaluation procedures.

Please take a look at the Bugs section of this document if you intend to run
`VLNets` jobs on the grid.


## Module Structure

The `VLNets` module has the following structure:
* `art` -- _art_ specific code.
    Includes modules to evaluate RNN energy estimator and create training
    dataset.
* `data`   -- A collection of data structures and data format converters.
* `models` -- Wrappers to simplify evaluation of _TensorFlow_ networks.


## RNN Energy Estimator Training

Training procedure of the RNN energy estimator consists of two steps:
1. Generation of the training dataset.
2. Training of the neural network.

### 1. Generation of the Training Dataset

Training dataset can be extracted from the DUNE _art_ files by running
`VLNEnergyDataGen` module that is located under `art/data_generators`. This
module extracts relevant variables from the _art_ files and saves them in a
_csv_ format that is more convenient for training of neural networks.

### 2. Training of the Neural Network

The easiest way to train the RNN energy estimator is to use [lstm_ee][lstm_ee]
package. Please refer to the documentation of the `lstm_ee` package for the
details.

[lstm_ee]: https://github.com/usert5432/lstm_ee "LSTM EE package"


## RNN Energy Estimator Evaluation

Modules to evaluate RNN energy estimators are located under the
`art/evaluators` subdirectory.  One can use the `VLNEnergyProducer` module in
order to evaluate RNN energy estimator and save results in the `VLNEnergy` data
product.

The evaluation process consists of two three. At first step, input variables
are extracted from the DUNE _art_ dataset and saved to a `VarDict` structure.
Then, input tensors for the neural network are constructed from this structure.
Finally, neural network is evaluated on these input tensors.


## Bugs

There are two bugs that you might encounter while running this module on the
grid:
1. "Illegal Instruction..." crash on the grid
2. Jobs eating all available memory

### 1. Illegal Instruction

_TensorFlow_ is compiled with support of modern SIMD instructions. However,
some grid nodes have old CPUs that do not support modern instructions. If
_TensorFlow_ job lands on a node with the old CPU, then the job will crash with
"Illegal instruction..." exception.

At this moment, the only offending SIMD instruction sets are `ssse3` and
`sse4.1`. So, if you run the `VLNets` evaluation on the grid, make sure to
filter out nodes that do not support these SIMD instruction sets. The grid
nodes can be filtered with the help of a following jobsub argument:

```
--append_condor_requirements=((TARGET.has_ssse3==true)&&(TARGET.has_sse4.1==true))
```

Note, that new versions of _TensorFlow_ library may be compiled with newer SIMD
instruction sets (e.g. `avx2`) and the condor requirements would need to be
updated.

### 2. Jobs eating all memory

There is a bug in _art_/ROOT that results in jobs eating all available memory
and eventually being held/killed for such behavior. This but is not specific to
`VLNets` and apparently is about to get fixed, see
[Issue Link](https://cdcvs.fnal.gov/redmine/issues/25349).

Before this bug is actually fixed, you can run `MVASelect` producer in your
job, since running it somehow makes memory issue go away.


## Notes

`VLNets` module relies on a _TensorFlow_ library to evaluate neural networks.
Dependency on _TensorFlow_ is confined to the `models` submodule, and
abstracted away by various models. This allows _TensorFlow_ framework to be
easily replaced by other deep learning frameworks (e.g. _PyTorch_, _MXNet_,
etc).

At this moment `VLNets` relies on information from PFParticles, but it can be
easily modified to rely on information extracted showers or/and tracks found
without Pandora.


