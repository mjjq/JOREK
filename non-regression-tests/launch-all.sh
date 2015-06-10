#!/bin/bash

startdir=`readlink -f $(dirname $0)`
codedir=`readlink -f ${startdir}/..` # Assumption about source code location
if [ -z "$JOREK_HOST" ]; then
    echo "JOREK_HOST environment variable is not defined, can not launch test cases"
    exit 1
fi
cd  || exit 1
for dirname in ${startdir}/testcases/*; do
 if [ -d $dirname ]; then
   name=$(basename $dirname)
   echo "== Launch job $name"
   (cd ${startdir}/job_scripts/${JOREK_HOST}; llsubmit ${name}.job)
 fi
done
exit 0
