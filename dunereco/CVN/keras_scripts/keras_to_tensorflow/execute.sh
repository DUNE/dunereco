#!/bin/bash

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 arg1 arg2 arg3" >&2
  exit 1
fi

singularity exec -e -p /home/singularity/ML/ubuntu16-ML-keras-210-tf4.simg python $1 $2 $3
