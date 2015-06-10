#!/bin/bash

cd testcases || exit 1
for name in *; do
 if [ -d $name ]; then
   echo $name
   ./get-testcase-data.sh $name
 fi
done
exit 0
