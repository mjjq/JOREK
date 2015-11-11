#!/bin/bash
# See http://jorek.eu/wiki/doku.php?id=penning_test for details

# Make sure it is compiled first (in subshell so we stay in this dir)
JOREK_DIR=../../..
(cd $JOREK_DIR && make penning_test)

if [ $# -eq 0 ]
then
   cases="n_radial_cases/*"
else
   cases=$@
fi

mkdir -p n_radial_results

err=0
for testcase in $cases; do
   echo "-------------------------------------------------------------------------------"
   echo "- " $testcase
   echo "-------------------------------------------------------------------------------"
   outfile=n_radial_results/`basename $testcase`
   $JOREK_DIR/penning_test < $testcase | tee ${outfile}.log
   grep "RESULT:" ${outfile}.log > ${outfile}.dat
   if ! grep -q "Tests successfull" ${outfile}.log; then
      err=1
   fi
done
# Gather all results into one file
cat n_radial_results/[0-9][0-9]*.dat > n_radial_results/n_radial.dat
gnuplot n_radial_plot.gp
