#!/bin/sh

export LANG=C
export ARCH=marconi
export JOREK_HOST=marconi_knl
export compilethreads=32

export MKL_DOMAIN_NUM_THREADS="MKL_BLAS=1"
export MKL_NUM_THREADS="1"
export MKL_DYNAMIC="FALSE"
export MKL_ENABLE_INSTRUCTIONS=AVX512_MIC
#-genv I_MPI_HBW_POLICY=hbw_preferred,hbw_bind 
export PRERUN="uniq $PBS_NODEFILE > hosts.txt"
#-genv I_MPI_PIN=1 -genv I_MPI_PIN_DOMAIN=socket -genv I_MPI_PIN_ORDER=spread -genv OMP_PROC_BIND=1
if [ -n "$OMP_NUM_THREADS" ]; then
#  ((PPN=136/$OMP_NUM_THREADS))
  ((PPN=68/$OMP_NUM_THREADS))
else
  PPN=1
fi
export MPIRUN="mpiexec.hydra -genv I_MPI_HBW_POLICY=hbw_preferred,hbw_bind  -genv I_MPI_DEBUG=4 -genv KMP_AFFINITY=verbose,scatter -ppn $PPN -f hosts.txt -np"
export BATCHCOMMAND="qsub"
