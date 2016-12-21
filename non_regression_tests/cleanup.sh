#!/bin/bash

startdir=`pwd`


for i in `ls -1d non_regression_tests/testcases/*`; do
  if [ -d $i ]; then
    cd $i
    echo "Cleaning up '$i'"
    
    source settings.sh
    
    rm -f *.h5 *.rst jorek_model* rst_*bin* *.tgz ${extra_remote_files}
    
    cd $startdir
  fi
done
