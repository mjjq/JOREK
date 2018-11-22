#!/bin/bash
eval `tclsh /work/imas/opt/modules-tcl/modulecmd.tcl $(basename $SHELL) autoinit`

module purge

module use /work/imas/opt/EasyBuild/modules/all
module use /work/imas/etc/modules/all

module load MUMPS/5.1.2-foss-2018a-metis
module load PaStiX/5.2.3-foss-2018a
module load HDF5/1.10.1-foss-2018a

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

export MPIRUN='mpirun --allow-run-as-root -n'
