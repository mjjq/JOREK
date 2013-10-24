#!/bin/bash
#
# Used by prepare_run.sh. Do not execute manually.
#

# --- Check if variables are set correctly
if [ -z "$__JP__mpi" ] || [ -z "$__JP__omp" ] || [ -z "$__JP__nodes" ]; then
  echo "ERROR in machine.sh: Variable mpi and/or omp and/or nodes are not set." >&2
  echo "  Note that you cannot run this script manually." >&2
  echo "  It is sourced by the prepare_run.sh script." >&2
  exit 1
fi

# --- Compute derived values
add_param "procs=$[$__JP__nodes*16]"
add_param "tasks_per_node=$[$__JP__mpi/$__JP__nodes]"

# --- Perform some checks
if [ $__JP__nodes -gt 16 ]; then
  echo "WARNING: More than 16 compute nodes requested!" >&2
fi
if [ $__JP__omp -gt 16 ]; then
  echo "WARNING: Number of OpenMP threads larger than 16!" >&2
fi
