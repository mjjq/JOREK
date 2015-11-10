#!/bin/bash
# See http://jorek.eu/wiki/doku.php?id=penning_test for details

# Make sure it is compiled first (in subshell so we stay in this dir)
JOREK_DIR=../../..
(cd $JOREK_DIR && make penning_test)

if [ $# -eq 0 ]
then
   cases="dt_cases/*"
else
   cases=$@
fi

mkdir -p dt_results

err=0
for testcase in $cases; do
   echo "-------------------------------------------------------------------------------"
   echo "- " $testcase
   echo "-------------------------------------------------------------------------------"
   outfile=dt_results/`basename $testcase`
   $JOREK_DIR/penning_test < $testcase | tee ${outfile}.log
   grep "RESULT:" ${outfile}.log > ${outfile}.dat
   if ! grep -q "Tests successfull" ${outfile}.log; then
      err=1
   fi
   #gnuplot -e "name='${outfile}'" dt_plot.gp
done

exit $err
