#!/bin/bash

export LANG=C
export JOREK_HOST=jenkins
export compilethreads=1
export BATCHCOMMAND=""
export PRERUN="export OMP_NUM_THREADS=2"
export MPIRUN="mpirun.mpich2 -np "
