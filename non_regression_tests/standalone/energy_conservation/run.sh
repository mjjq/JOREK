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
   echo "Running in batch mode, not running analysis"
   cd $PBS_O_WORKDIR
   verbose=1
   no_recompile=1
   no_analysis=1

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
   $JOREK_DIR/jorek2_particles < jorek_in
fi

if [ "$no_analysis" -eq 0 ]; then
   echo "Running analysis"
   # Find out how many output files exist now (assumes equal amount of momentum and energy)
   output_times=`ls -1 energy*_*dat | tr -d energydat. | cut -d_ -f1 | sort | uniq`
   output_cpus=` ls -1 energy*_*dat | tr -d energydat. | cut -d_ -f2 | sort | uniq`

   rm -f energy.dat momentum.dat
   # Join all output files together by cpu id
   for output_time in $output_times; do
      rm -f energy$output_time.dat momentum$output_time.dat
      for output_cpu in $output_cpus; do
	 cat energy${output_time}_${output_cpu}.dat >> energy$output_time.dat
	 cat momentum${output_time}_${output_cpu}.dat >> momentum$output_time.dat
      done

      echo -n $output_time " " >> energy.dat
      Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(min(d), max(d), mean(d), sd(d), sep=" ")' < energy$output_time.dat >> energy.dat
      echo "" >> energy.dat
      echo -n $output_time " " >> momentum.dat
      Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(min(d), max(d), mean(d), sd(d), sep=" ")' < momentum$output_time.dat >> momentum.dat
      echo "" >> momentum.dat
      echo $output_time
   done

   points="with points pointtype 6"
   # Plot results
   cat <<EOF > gnuplot_in
      set logscale y
      set ylabel "Energy"
      set xlabel "timestep"
      set terminal png

      set out 'energy.png'
      plot 'energy.dat' u 1:2 t "Min" w l, 'energy.dat' u 1:4 t 'Mean' w l, 'energy.dat' u 1:3 t 'Max' w l

      unset logscale y
      set out 'momentum.png'
      set ylabel 'Toroidal momentum'
      plot 'momentum.dat' u 1:2 t "Min" w l, 'momentum.dat' u 1:4 t 'Mean' w l, 'momentum.dat' u 1:3 t 'Max' w l
EOF

   # Create some simple graphs
   gnuplot gnuplot_in

   # Copy and run energy_momentum_plot.py for more complicated graphs
   cp ../energy_momentum_plot.py .
   echo "generating energy plot"
   python energy_momentum_plot.py energy
   echo "generating momentum plot"
   python energy_momentum_plot.py momentum
fi 
