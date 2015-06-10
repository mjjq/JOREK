source env.sh
SCRIPTDIR=`readlink -f $(dirname $0)/../..`
LIST="h5_tearing_limiter_199 tearing_limiter_199 tearing_limiter_303 ballooning_xpoint_303 xpoint_equil xpoint_grid polar_grid"
${SCRIPTDIR}/testcases/download_ref.sh
for case in ${LIST}; do
  echo $case
  if [ -f ${case}.job ]; then 
    (${SCRIPTDIR}/run_test.sh -r $case 1>${case}_comp.out) && llsubmit ${case}.job
  fi
done
