#!/bin/bash
# @ job_name = jorek
# @ class = clallmds+
# @ as_limit = 24gb
# @ error  = $(job_name).$(jobid).err
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL
# @ wall_clock_limit = 0:55:00
# @ job_type = mpich
# @ restart = no
#######################################################################
# Look to non-regression-tests/testcases/tearing_limiter_199
# to determine nb of nodes, MPI processes, and nb of threads
#######################################################################
# @ node = 1 
# @ total_tasks = 2
# @ resources = ConsumableCpus(8)
#####################################
# @ node_usage = not_shared
# @ queue

ulimit -d unlimited
ulimit -s unlimited
ulimit -Sa
ulimit -Ha
echo =0= $0
CASE=h5_tearing_limiter_199
JOREK_HOME=${HOME}/jorek

cd ${JOREK_HOME}/non-regression-tests
source examples/poincare/env.sh
#./run_test.sh -r $CASE
./run_test.sh -k -n $CASE
