. env.sh
SCRIPTDIR=`readlink -f $(dirname $0)/../..`
LIST="xpoint_equil xpoint_grid polar_grid tearing_limiter_199 tearing_limiter_303 ballooning_xpoint_303"
LIST="h5_tearing_limiter_199"
for case in ${LIST}; do
  echo $case
  if [ -f ${case}.cmd ]; then 
    ${SCRIPTDIR}/run_test.sh -r $case && llsubmit ${case}.cmd
  fi
done
