#!/bin/bash
#PBS -N JOREK
#PBS -V
#PBS -j oe
#PBS -l walltime=3:00:00
#PBS -l select=4:ncpus=8:ompthreads=8:mpiprocs=1

TRKDIR=~/sept12
MODEL=199
PREFIX=sp

source /data/JOREK/jorek_env/jorek_env.sh
cd $TRKDIR
make cleanall
util/compile_tests.sh ${MODEL}

export PRERUN="export OMP_NUM_THREADS=8"
export MPIRUN="mpirun -n 4"
export BASEDIR="/data/GL310811"
mpdboot -n 4 -f $PBS_NODEFILE
util/launch_tests.sh ${PREFIX}${MODEL}
mpdallexit
