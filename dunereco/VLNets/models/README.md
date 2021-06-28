# VLNets/models

This module contains _TensorFlow_ wrappers to simplify evaluation of the neural
networks. The most generic _TensorFlow_ wrapper is named `TFModel` and located
under `tf_model` subdirectory. Models specialized for a particular application
(e.g. neutrino energy estimation) are located under the `zoo` subdirectory.


## TFModel

`TFModel` wrapper is used to load _TensorFlow_ networks from a protobuf graph,
and perform network evaluation.

Constructor of `TFModel` receives path to a directory that contains trained
_TensorFlow_ network. This directory is expected to have the following
contents:
* `model.pb`    -- _TensorFlow_ network saved in a protobuf format.
* `config.json` -- Network configuration (c.f. "TFModel Config Format").

`TFModel` parses configuration file (`config.json`) and loads _TensorFlow_
graph (`model.pb`) lazily, only when the actual network evaluation is
requested.


## TFModel Config Format

Network configuration file `config.json` is a _json_ dictionary that specifies
input/output nodes of the _TensorFlow_ graph, names of the input variables and
their order. This information is used at the evaluation time to construct
proper input tensors for the network from a `VarDict` structure.

The `config.json` is specific for each model. For example, the `VLNEnergyModel`
defined in the `zoo` subdirectory, expects configuration file of the
following form:

```
{
    "input_png3d":    `Name of TF node that receives particle level inputs`,
    "input_slice":    `Name of TF node that receives event level inputs`,
    "target_primary": `Name of TF node that returns primary energy`,
    "target_total":   `Name of TF node that returns total energy`,
    "vars_png3d":     `List of particle level input variables`,
    "vars_slice":     `List of event level input variables`,
}
```

Cryptic names like `png3d` and `slice` are unfortunate result of porting RNN
energy estimator from NOvA and may be changed in the future to something more
meaningful.

