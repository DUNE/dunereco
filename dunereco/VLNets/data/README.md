# VLNets/data

This module contains definitions of structures that `VLNets` can produce
(in `data/structs`), and implementation of various data exporters (under
`data/exporters`).

Data exporters are objects that can serialize variables from a `VarDict`
structure into formats suitable for training of neural networks.  At this
moment only _csv_ data exporter is available.


## VarDict

The `VarDict` structure is defined in `data/structs/VarDict.h` and deserves
mentioning here. This structure is used as a temporary holder of variables
that have been extracted from the _art_ files. Variables, stored in this
structure, can later be used to construct input tensors for neural networks,
or can be exported to other formats with the help of data exporters.

