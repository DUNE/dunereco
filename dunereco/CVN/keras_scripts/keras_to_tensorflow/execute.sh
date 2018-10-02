#!/bin/bash

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 arg1 arg2 arg3" >&2
  exit 1
fi

python $1 $2 $3
