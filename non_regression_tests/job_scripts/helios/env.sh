#!/bin/bash

module purge
module load intel/14.0.4.211
module load bullxmpi/1.2.8.2
module load hdf5/1.8.14
export LANG=C
export JOREK_HOST=helios
export MPIRUN="mpirun -np"
export SUBCMD="sbatch -A JOREKIRQ"
