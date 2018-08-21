#!/bin/bash

  for j in {0..19}
  do
     sbatch run.sh $1 j 1 1
 done
