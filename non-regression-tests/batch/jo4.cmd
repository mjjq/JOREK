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
# @ node = 1
# @ total_tasks = 2
# @ resources = ConsumableCpus(8)
# @ node_usage = not_shared
# @ queue

ulimit -d unlimited
ulimit -s unlimited
ulimit -Sa
ulimit -Ha
CASE=tearing_limiter_199
JOREK_HOME=/gpfshome/mds/staff/glatu/jorek
SCRATCHDIR=/gpfsdata/glatu
BASEDIR=${SCRATCHDIR}

cd ${JOREK_HOME}/non-regression-tests
export PRERUN='export KMP_AFFINTY=verbose'
export MPIRUN="mpiexec -launcher-exec /opt/ibmll/LoadL/scheduler/full/bin/llspawn.stdio -f $LOADL_HOSTFILE  -n "
./run_test.sh -k -n tearing_limiter_199
