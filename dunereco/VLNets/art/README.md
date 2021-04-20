# VLNets/art

This subdirectory contains libraries that directly link to _art_.


## VLNets/art/data_generators

The `data_generators` subdirectory contains modules to construct training
datasets from the _art_ files and save them to simpler formats, more
appropriate for training of neural networks.

At present moment, `data_generators` contains only `VLNEnergyDataGen` module
that is capable of extracting training dataset for the RNN energy estimator.
The extracted dataset is saved to a _csv_ format.

If you wish to run `VLNEnergyDataGen` on a grid with _project.py_, then you
should specify extension (.csv) of the output files in the _project.py_ stage
configuration,
e.g.:
```
  <stage name="training_dataset_creation">
...
    <datafiletypes>csv</datafiletypes>
...
  </stage>
```


## VLNets/art/evaluators

The `evaluators` subdirectory contains _art_ modules to evaluate neural
networks.

The `VLNEnergyProducer` module can evaluate RNN  energy estimator and save
evaluation result in the _art_ files. It is intended to be a primary producer.

The `VLNEnergyAnalyzer` module extracts information necessary for RNN energy
estimator evaluation from the _art_ files, and saves it to the _csv_ file.
Additionally, it performs network evaluation itself and stores results of
network evaluation in the _csv_ file as well. This module can be used to
cross validate that _art_ level network evaluation produces correct results.


## VLNets/art/products

This subdirectory contains _art_ instructions to construct ROOT dictionaries
for products that `VLNets` can produce.


## VLNets/art/var_extractors

`var_extractors` are a collection of objects to extract information from
the _art_ files and store it in the temporary dictionary `VarDict`.

