#!/bin/bash

export LANG=C
export JOREK_HOST=marconi
export compilethreads=32
USE_INTELMPI=$(mpirun -V | grep -c -i "Intel.*MPI")
USE_OPENMPI=$(mpirun -V | grep -c -i "Open.*MPI")
MPIRUN="mpirun -np"
if [ ${USE_INTELMPI} -gt 0 ]; then
  MPIRUN="mpirun -genv KMP_AFFINITY=verbose -genv I_MPI_DEBUG=4 -genv I_MPI_PIN=1 -genv I_MPI_PIN_DOMAIN=socket -genv I_MPI_PIN_ORDER=spread -genv OMP_PROC_BIND=1 --map-by ppr:1:socket:pe=18 -np"
fi
if [ ${USE_OPENMPI} -gt 0 ]; then
  MPIRUN="mpirun -x KMP_AFFINITY=verbose -x OMP_PROC_BIND=1 -n"
fi
export MPIRUN
export BATCHCOMMAND="qsub"
