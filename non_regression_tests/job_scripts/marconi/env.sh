#!/bin/bash

export LANG=C
export JOREK_HOST=marconi
export compilethreads=32
export MPIRUN="mpirun -genv KMP_AFFINITY=verbose -genv I_MPI_DEBUG=4 -genv I_MPI_PIN=1 -genv I_MPI_PIN_DOMAIN=socket -genv I_MPI_PIN_ORDER=spread -genv OMP_PROC_BIND=1 --map-by ppr:1:socket:pe=$OMP_NUM_THREADS -np"
export BATCHCOMMAND="qsub"
export http_proxy=http://proxy:3128 
