#!/bin/bash
# Usage:                                                      
#   # Example For model 199
#    ./launch_test.sh sim199 # Compute the equilibrium
#    ./launch_test.sh sim199 # Compute the first 50 time steps
#    ./launch_test.sh sim199 # Compute the time steps 51 to 100

# ---------- Following variables has to be changed by the user ----------

# Trunk of jorek to be used
TRKDIR=/scratch/jorekgl/trunk2g
# Location of 'util' directory that contains 'setinput.sh'
UTILDIR=${TRKDIR}/util
# Location of namelist directory that contains input jorek files
NAMEDIR=${TRKDIR}/namelist
# Location of directory that contains executables of Jorek
EXEDIR=${TRKDIR}
# Location of target directory where simulation directory will be created
BASEDIR=/scratch/jorekgl
# Mpirun command
MPIRUN="mpiexec"

# List of executables used during following simulations
# The n-th executable has three parameters : 
#      model[n] ntor[n] nplane[n]
declare -a model=( 199 199 199 302 302)
declare -a ntor=(  1   3   7   1   3  )
declare -a nplane=(1   4   8   1   8  )

declare -a list_idx=(     0                  3       )
declare -a list_inputs=(  "in199"            "xx302" )
declare -a list_orig=(    "model199/intear"  "model300/inxflow" )

#-------------------------------------------------------------------------
SIMNAME="$1"
if [ ${#SIMNAME} -lt 4 ]; then
    echo "ERROR: the parameter (simulation directory name) is too short,"
    echo "try 'sim${model[0]}' for example."
    exit 0
fi
((l=${#SIMNAME}-3))
# Name of the prefix of the simulation directory
PREFIX=${SIMNAME:0:$l}
# Model number
MODNB=${SIMNAME:$l:3}
for (( j = 0, notfound=1 ; j < ${#list_idx[@]} ; j++ )); do
    i=${list_idx[$j]}
    if [ "${MODNB}" = "${model[$i]}" ]; then
	((notfound=0))
    fi
#    echo "/${model[$i]}/ /$MODNB/ $notfound"
done
if ((notfound)); then
    echo "ERROR: give as a parameter the simulation directory that you want,"
    echo "the last three digits must be a valid model number:"
    for (( j = 0, notfound=1 ; j < ${#list_idx[@]} ; j++ )); do
        i=${list_idx[$j]}
        printf "${model[$i]} "
    done
    echo ""
    exit 0
fi
echo "Launching simulation"
echo "Name of the simulation directory : $SIMNAME"
echo "Model number: $MODNB"


# Copy the executables to target simulation directory
# $WKDIR
for (( j = 0 ; j < ${#list_idx[@]} ; j++ )); do
   i=${list_idx[$j]}
   if [ ${model[$i]} = "$MODNB" ]; then
      WKDIR=${BASEDIR}/${PREFIX}${model[$i]}
      if [ ! -d $WKDIR ]; then
      	   mkdir $WKDIR
      	   if [ ! -f $NAMEDIR/${list_orig[$j]} ]; then
      	       echo "Pb $NAMEDIR/${list_orig[$j]} not found"
      	       exit 1
      	   fi
      	   sed "s/nstep.*=/nstep_n =/;s/tstep.*=/tstep_n =/;" ${NAMEDIR}/${list_orig[$j]} > ${WKDIR}/${list_inputs[$j]}
      fi
   fi
done
for (( i = 0 ; i < ${#model[@]} ; i++ )); do
   printf "model $i ${model[$i]}\n"
   if [ ${model[$i]} = "$MODNB" ]; then
     EXE=${EXEDIR}/j${model[$i]}_${ntor[$i]}
     WKDIR=${BASEDIR}/${PREFIX}${model[$i]}
     cp ${EXE} ${WKDIR}/
   fi
done

# First Case: model199
if [ "$MODNB" = "199" ]; then
  ((j=0))
  i=${list_idx[$j]}
  WKDIR=${BASEDIR}/${PREFIX}${model[$i]}
  cd $WKDIR
  cp macroscopic_vars.dat old_macros_vars.dat 2>/dev/null
  rm -f macroscopic_vars.dat 2>/dev/null
  INFILE=${list_inputs[$j]}
  
  if [ ! -f jorek_equil.rst ]; then
      ${UTILDIR}/setinput.sh ${INFILE} restart=.f. nstep_n=0 tstep_n=1
      EXE=j${model[$i]}_1
      ${MPIRUN} ./${EXE} < ${INFILE} | tee out_equil
      cp jorek_restart.rst jorek_equil.rst
  else
      ${UTILDIR}/setinput.sh ${INFILE} restart=.t. nstep_n=200 tstep_n=1000
      EXE=j${model[$i]}_3
      ${MPIRUN} ./${EXE} < ${INFILE} | tee out_loop
      
      ${UTILDIR}/extract_live_data.sh energies energies.dat 
  fi
  exit 0
fi 

# Second Case: model302 point X
if [ "$MODNB" = "302" ]; then
  ((j=1))
  i=${list_idx[$j]}
  WKDIR=${BASEDIR}/${PREFIX}${model[$i]}
  cd $WKDIR
  cp macroscopic_vars.dat old_macros_vars.dat
  rm -f macroscopic_vars.dat
  INFILE=${list_inputs[$j]}
  
  if [ ! -f jorek_equil.rst ]; then
      ${UTILDIR}/setinput.sh ${INFILE} restart=.f. nstep_n=0 
      EXE=j${model[$i]}_1
      ${MPIRUN} ./${EXE} < ${INFILE} 2>err0 | tee out_equil
      status=$?; if [ $status -eq 0 ]; then 
	  cp jorek_restart.rst jorek_equil.rst
      fi
  else
      if [ ! -f jorek_rst1.rst ]; then
	  cp jorek_equil.rst jorek_restart.rst 
	  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 10, 9, 9, 9, 4' 'tstep_n= 1e-3, 1e-2, 1e-1, 1, 2' nout=10
	  EXE=j${model[$i]}_1
	  ${MPIRUN} ./${EXE} < ${INFILE} 2>err1| tee out_loop1
	  status=$?; if [ $status -eq 0 ]; then 
	      cp jorek_restart.rst jorek_rst1.rst; 
	  fi
      else
	  if [ ! -f jorek_rst2.rst ]; then
	      cp jorek_rst1.rst jorek_restart.rst 
	      ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 5, 4' 'tstep_n= 2, 5' nout=10
	      EXE=j${model[$i]}_3
	      ${MPIRUN} ./${EXE} < ${INFILE} 2>err2 | tee out_loop2
	      status=$?; if [ $status -eq 0 ]; then 
		  cp jorek_restart.rst jorek_rst2.rst
	      fi
	  else
	      if [ ! -f jorek_rst3.rst ]; then
		  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 50' 'tstep_n= 5' nout=10
		  EXE=j${model[$i]}_3
		  ${MPIRUN} ./${EXE} < ${INFILE} 2>err3 | tee out_loop3
	      fi
	  fi
      fi
  fi
  ${UTILDIR}/extract_live_data.sh energies energies.dat 
  exit 0
fi

