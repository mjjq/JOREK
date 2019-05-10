#!/bin/bash
# See http://jorek.eu/wiki/doku.php?id=penning_test for details
# Variables
OPTIND=1         # Reset in case getopts has been used previously in the shell.
verbose=0 # Default verbosity
cases="dt_cases/* n_radial_cases/*" # Default cases
stop_on_exit=0

# Help
show_help()
{
   echo "This file runs a testcase with the penning_test program"
   echo "See http://jorek.eu/wiki/doku.php?id=penning_test for details"
   echo ""
   echo "Optional arguments"
   echo "  -e  Quit on error in a case"
   echo "  -v  Enable verbose mode (show compilation and program output)"
   echo "  -h  Show this help"
}


# Exit on ctrl-C
control_c()
{
   exit $?
}
trap control_c SIGINT

# Colors
red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`

# Parse options
while getopts "h?ve" opt; do
   case "$opt" in
      h|\?)
      show_help
      exit 0
      ;;
      v)  verbose=1
      ;;
      e)  stop_on_exit=1
      ;;
   esac
done
# Remove options from $@
shift $((OPTIND-1))

# Make sure it is compiled first (in subshell so we stay in this dir)
echo -en "Compiling JOREK particle test..."
if [ $verbose -eq 1 ]; then
   echo "" # finish the line
   (j2p penning)
else
   (j2p penning 2>/dev/null >/dev/null)
fi
if [ $? -ne 0 ]
then
   echo -e "${red}FAIL${reset}"
   exit $?
else
   echo -e "${green}SUCCES${reset}"
fi

if [ $# -gt 0 ]
then
   cases=$@
fi

for testcase in $cases; do
   echo "-------------------------------------------------------------------------------"
   echo "- " $testcase
   echo "-------------------------------------------------------------------------------"
   outdir=`echo $testcase | cut -d/ -f1 | sed 's/cases/results/'`
   mkdir -p $outdir
   outfile=$outdir/`basename $testcase`
   if [ $verbose -eq 1 ]; then
      ./penning < $testcase | tee ${outfile}.log
      if [ ${PIPESTATUS[0]} -ne 0 ]
      then
	 echo -e "${red}FAIL${reset}"
	 if [ $stop_on_exit -eq 1 ]; then exit ${PIPESTATUS[0]}; fi
      else
	 echo -e "${green}SUCCES${reset}"
      fi
   else
      ./penning < $testcase > ${outfile}.log
      if [ $? -ne 0 ]
      then
	 echo -e "${red}FAIL${reset}"
	 if [ $stop_on_exit -eq 1 ]; then exit $?; fi
      else
	 echo -e "${green}SUCCES${reset}"
      fi
   fi
   grep "RESULT:" ${outfile}.log > ${outfile}.dat
done

# !!! extra processing step for n_radial results
cat n_radial_results/[0-9][0-9]*.dat > n_radial_results/n_radial.dat
