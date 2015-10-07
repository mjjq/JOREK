#!/bin/bash

module purge  

module use /work/imas/opt/EasyBuild/modules/all
module load cURL

module use /work/imas/etc/modulefiles
module load hdf5

module load intel/12.0.2
module load mpich2/3.1.3-intel
module load scotch/5.1.12b
module load metis/5.1.0
module load mumps/4.10.0
module load pastix/5.2.2.16
module load ppplib/14.4.8
module load fftw/3.3.4

export LANG=C
export JOREK_HOST=iter-hpc
export compilethreads=8
export PRERUN="export OMP_NUM_THREADS=2"
export MPIRUN="mpirun -np "
export BATCHCOMMAND="qsub"

export http_proxy=10.10.240.146:8080
export https_proxy=10.10.240.146:8080
