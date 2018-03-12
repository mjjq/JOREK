#!/bin/bash
eval `tclsh /work/imas/opt/modules-tcl/modulecmd.tcl $(basename $SHELL) autoinit`

module purge

module use /work/imas/opt/EasyBuild/modules/all

module load cURL/7.40.0-GCC-4.8.3
module load zlib/1.2.8-gompi-1.5.16
module load intel/12.0.2
module load mpich2/3.1.3-intel
module load scotch/5.1.12b
module load metis/5.1.0
module load mumps/4.10.0
module load pastix/5.2.2.16
module load fftw/3.3.4
module load hdf5
module unload GCC

module list

export ZLIB_HOME=/work/imas/opt/EasyBuild/software/zlib/1.2.8-gompi-1.5.16

export LANG=C
export JOREK_HOST=iter-hpc
export compilethreads=4
export MAKEFLAGS="-j$compilethreads"
export PRERUN="export OMP_NUM_THREADS=4"
export MPIRUN="mpirun -np "
export BATCHCOMMAND="qsub"

export http_proxy=10.10.240.146:8080
export https_proxy=10.10.240.146:8080
