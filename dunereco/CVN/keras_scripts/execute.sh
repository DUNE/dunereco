#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 file" >&2
  exit 1
fi

#singularity exec --nv -e -p /home/singularity/ML/ubuntu16-ML-keras-210.simg python $1
#singularity exec --nv -e -p /home/singularity/ML/ubuntu16-ML-keras-210-tf4.simg python $1
singularity exec --nv -e -p /home/singularity/ML/ubuntu1604-cuda-90-ML-tf1.8.simg python $1
