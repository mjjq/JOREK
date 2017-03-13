#!/bin/bash

#module purge
#module load intel/12.1.0 mvapich2/1.8.1 mkl/10.3
#module load intel/12.1.0   intelmpi/5.0.1 mkl/11.2

export LANG=C
export JOREK_HOST=poincare
export compilethreads=8
export BATCHCOMMAND="llsubmit"
#export KMP_STACK_SIZE=16M
export KMP_AFFINITY="verbose,norespect," # if this is not set, the OpenMP threads are confined on one core
#export KMP_AFFINITY="verbose,norespect,compact"
#export OMP_PROC_BIND=true
#export HFI_NO_CPUAFFINITY=1
#export I_MPI_DEBUG=5
#export I_MPI_FABRICS="shm:tmi"
#export I_MPI_PIN_DOMAIN=16:compact
#export I_MPI_PIN_DOMAIN=omp:compact


export MKL_NUM_THREADS=1
export MKL_DYNAMIC=0
USE_INTELMPI=0
USE_MVAPICH2=1
#if [ ${USE_INTELMPI} -gt 0 ]; then
#HO=`basename $LOADL_HOSTFILE`
#uniq $LOADL_HOSTFILE > $PWD/$HO
which mpirun
# -nolocal -f $PWD/$HO
#--map-by ppr:1:node:PE=16 -f $LOADL_HOSTFILE -rmk llspawn.stdio
#export MPIRUN="mpirun -f $LOADL_HOSTFILE -ppn 1 -n"
#  echo $MPIRUN
#fi
export MV2_ENABLE_AFFINITY=0
export MPIRUN="mpiexec -envlist LANG,OMP_NUM_THREADS,KMP_STACK_SIZE,KMP_AFFINITY,MKL_NUM_THREADS,MKL_DYNAMIC,MV2_ENABLE_AFFINITY -launcher-exec /opt/ibmll/LoadL/scheduler/full/bin/llspawn.stdio -ppn 1 -binding -verbose -f $LOADL_HOSTFILE -n"
echo $MPIRUN

