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
export MPIRUN="mpiexe.hydra -genv I_MPI_HBW_POLICY=hbw_preferred,hbw_bind  -genv I_MPI_DEBUG=4 -genv KMP_AFFINITY=verbose,scatter -ppn 4 -f hosts.txt -np"
export BATCHCOMMAND="qsub"
