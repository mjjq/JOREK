#!/bin/bash

startdir=`readlink -f $(dirname $0)`
codedir=`readlink -f ${startdir}/..` # Assumption about source code location
cd  || exit 1
for dirname in ${startdir}/testcases/*; do
 if [ -d $dirname ]; then
   name=$(basename $dirname)
   echo "== Launch job $name"
   (cd ${startdir}/job_scripts/poincare; llsubmit ${name}.job)
 fi
done
exit 0
