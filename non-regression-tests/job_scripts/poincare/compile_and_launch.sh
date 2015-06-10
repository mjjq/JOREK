source env.sh
SCRIPTDIR=`readlink -f $(dirname $0)/../..`
LOCALDIR=$PWD
LIST="tearing_limiter_199"
#${SCRIPTDIR}/testcases/download_ref.sh
cd ${SCRIPTDIR}/..
for case in ${LIST}; do
  echo $case
  if [ -f ${case}.job ]; then 
    (${SCRIPTDIR}/run_test.sh -r $case 1>${case}_comp.out) && llsubmit ${LOCALDIR}/${case}.job
  fi
done
