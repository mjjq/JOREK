#!/bin/bash
eval `tclsh /work/imas/opt/modules-tcl/modulecmd.tcl $(basename $SHELL) autoinit`

module purge

module use /work/imas/etc/modulefiles
module use /work/imas/opt/EasyBuild/modules/all

module load mpich2/3.1.3-gnu
module load SCOTCH/6.0.0_esmumps-goolf-1.5.16
module load METIS/4.0.3-goolf-1.5.16
module load MUMPS/4.10.0-goolf-1.5.16-metis
module load pastix/5.2.2.16 # TODO: get goolf version?
module load FFTW/3.3.4-gompi-1.5.16 
module load HDF5/1.8.9-goolf-1.5.16

export LANG=C
export JOREK_HOST=iter-hpc
export compilethreads=4
export MAKEFLAGS="-j$compilethreads"
export PRERUN="export OMP_NUM_THREADS=4"
export MPIRUN="mpirun -np "
export BATCHCOMMAND="qsub"

export http_proxy=10.10.240.146:8080
export https_proxy=10.10.240.146:8080
