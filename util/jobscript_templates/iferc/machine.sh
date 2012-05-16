#!/bin/bash
#
# Used by prepare_run.sh. Do not execute manually.
#

# --- Check if variables are set correctly
if [ -z "$__JP__mpi" ] || [ -z "$__JP__omp" ] || [ -z "$__JP__nodes" ]; then
  echo "ERROR in machine.sh: Variable mpi and/or omp and/or nodes are not set." >&2
  exit 1
fi

# --- Compute derived values
add_param "ppn=$[$__JP__mpi/$__JP__nodes]"

# --- Perform some checks
if [ $__JP__nodes -gt 256 ]; then
  echo "WARNING: More than 256 compute nodes requested!" >&2
fi
if [ $[__JP__ppn*__JP__omp] -ne 16 ]; then
  echo "WARNING: PPN*OMP is not 16!" >&2
fi
