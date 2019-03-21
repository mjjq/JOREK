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
# OMP_NUM_THREADS>1 doesn't work with this GCC on the ITER platform...
export PRERUN="export OMP_NUM_THREADS=1"
export MPIRUN="mpirun -np "
export BATCHCOMMAND="qsub"

export http_proxy=10.10.240.146:8080
export https_proxy=10.10.240.146:8080

export MPIRUN='mpirun --allow-run-as-root -n'
