#!/bin/bash
eval `tclsh /work/imas/opt/modules-tcl/modulecmd.tcl $(basename $SHELL) autoinit`

module purge

module use /work/imas/opt/EasyBuild/modules/all
module use /work/imas/etc/modules/all

module load MUMPS/5.0.1-goolf-1.5.16-metis
module load PaStiX/5.2.2.22-goolf-1.5.16
module load HDF5/1.10.0-patch1-goolf-1.5.16
module load Szip/2.1-goolf-1.5.16
module load zlib/1.2.8-GCC-4.8.3
module load libibverbs/1.1.4
# patch libibverbs module on iter since it does not set LIBRARY_PATH correctly
export LIBRARY_PATH=$LIBRARY_PATH:/work/imas/opt/libibverbs/1.1.4/lib

export ZLIB_HOME=$EBROOTZLIB

export LANG=C
export JOREK_HOST=iter-hpc
export compilethreads=4
export MAKEFLAGS="-j$compilethreads"
export PRERUN="export OMP_NUM_THREADS=4"
export MPIRUN="mpirun -np "
export BATCHCOMMAND="qsub"
export CXXFLAGS=-O0 # problem with stdio library on ITER http://gcc.1065356.n8.nabble.com/g-4-8-fails-with-Ox-option-td953876.html

export http_proxy=10.10.240.146:8080
export https_proxy=10.10.240.146:8080
