#!/bin/bash
#PBS -N energy_conservation
#PBS -V
#PBS -j oe
#PBS -l nodes=1
# See README.md for details
JOREK_DIR=~/jorek # absolute path required
OPTIND=1         # Reset in case getopts has been used previously in the shell.
verbose=0 # Default verbosity
analysis_only=0
no_recompile=0
no_analysis=0

# Help
show_help()
{
   echo "./run.sh [-vnah] directory"
   echo "This file runs a testcase with the jorek2_particles program"
   echo "See README.md for details"
   echo ""
   echo "Optional arguments"
   echo "  -v  Enable verbose mode (show compilation and program output)"
   echo "  -n  Do not recompile jorek2_particles"
   echo "  -a  Run analysis only (implies -n)"
   echo "  -h  Show this help"
}

control_c() {
   exit $?
}
trap control_c SIGINT

# Colors
red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`

# Do not parse options if $PBS_O_WORKDIR is set, and move there, and do not recompile
if [ ! -z "${PBS_O_WORKDIR}" ]; then
   echo "Running in batch mode"
   cd $PBS_O_WORKDIR
   verbose=1
   no_recompile=1

   if [ -d "$dir" ]; then
      echo "Switching to directory $dir"
      cd $dir
   else
      echo "Please give an existing directory as -v 'dir=DIRECTORY' parameter"
      exit 1
   fi
      
else
   # Parse options
   while getopts "h?vna" opt; do
      case "$opt" in
	 h|\?)
	 show_help
	 exit 0
	 ;;
	 v)  verbose=1
	 ;;
	 n)  no_recompile=1
	 ;;
	 a)  analysis_only=1; no_recompile=1
	 ;;
      esac
   done
   # Remove options from $@
   shift $((OPTIND-1))

   # Check if our target dir exists
   if [ "$#" -gt 0 ]; then
      if [ -d $1 ]; then
	 echo "Switching to directory $1"
	 cd $1
      fi
   else
      echo "Please give an existing directory as program parameter"
      show_help
      exit 1
   fi
fi

# Make sure it is compiled first (in subshell so we stay in this dir)
if [ "$no_recompile" -eq 0 ]; then
   echo -en "Compiling JOREK2_particles..."
   if [ $verbose -eq 1 ]; then
      echo "" # finish the line
      (cd $JOREK_DIR && make jorek2_particles)
   else
      (cd $JOREK_DIR && make jorek2_particles 2>/dev/null >/dev/null)
   fi
   if [ $? -ne 0 ]; then
      echo -e "${red}FAIL${reset}"
      exit $?
   else
      echo -e "${green}SUCCES${reset}"
   fi
fi


if [ ! -e jorek_in ]; then
   echo "${red}Please include jorek_in in the directory or set the desired directory${reset}"
   exit 2
fi
if [ ! -e jorek_restart.rst ]; then
   echo "${red}Please include jorek_restart.rst in the directory or give this command the parameter of the desired directory${reset}"
   exit 3
fi

if [ "$analysis_only" -eq 0 ]; then
   echo "Cleaning old files"
   rm *.dat *.png *.vtk
   echo "Running testcase jorek_restart.rst"
   # Run the testcase given as input parameter
   $JOREK_DIR/jorek2_particles < jorek_in | tee jorek_log
fi

if [ "$no_analysis" -eq 0 ]; then
   echo "Running analysis"
   # For now only 1-cpu support for the following part!
   grep "energy min" jorek_log | cut -d: -f2 > energy.dat
   grep "momentum min" jorek_log | cut -d: -f2 > momentum.dat

   # Awk magic to divide everything by the first line
   # See http://unix.stackexchange.com/questions/174371/calculate-and-divide-by-total-with-awk
   awkscr='FNR==NR{min=$1;mean=$2;max=$3;sd=$4;tot=$5;nextfile}{printf "%f\t%f\t%f\t%f\t%f\n", $1/min, $2/mean, $3/max, $4/sd, $5/tot}'
   awk "$awkscr" energy.dat energy.dat > energy_norm.dat
   awk "$awkscr" momentum.dat momentum.dat > momentum_norm.dat

   # Plot results TODO shading
   cat <<EOF > gnuplot_in
      set logscale y
      set ylabel "Energy"
      set xlabel "time"
      set terminal png
      set style data lines

      set out 'energy.png'
      plot 'energy.dat' u 1 t "Min", 'energy.dat' u 2 t 'Mean', 'energy.dat' u 3 t 'Max', 'energy.dat' u (\$2+2*\$4) t 'Upper 95% interval', 'energy.dat' u (\$2-2*\$4) t 'Lower 95% interval'

      unset logscale y
      set out 'momentum.png'
      set ylabel 'Toroidal momentum'
      plot 'momentum.dat' u 1 t "Min", 'momentum.dat' u 2 t 'Mean', 'momentum.dat' u 3 t 'Max', 'momentum.dat' u (\$2+2*\$4) t 'Upper 95% interval', 'momentum.dat' u (\$2-2*\$4) t 'Lower 95% interval'

      set out 'energy_norm.png'
      set ylabel 'Change in energy'
      plot 'energy_norm.dat' u 1 t "Min", 'energy_norm.dat' u 2 t 'Mean', 'energy_norm.dat' u 3 t 'Max', 'energy_norm.dat' u 4 t 'Sigma', 'energy_norm.dat' u 5 t 'Total'

      set out 'momentum_norm.png'
      set ylabel 'Change in momentum'
      plot 'momentum_norm.dat' u 1 t "Min", 'momentum_norm.dat' u 2 t 'Mean', 'momentum_norm.dat' u 3 t 'Max', 'momentum_norm.dat' u 4 t 'Sigma', 'momentum_norm.dat' u 5 t 'Total'
EOF

   # Create some simple graphs
   gnuplot gnuplot_in

   cp ../particle_stats.py .
   python particle_stats.py energy
   python particle_stats.py momentum
fi 
