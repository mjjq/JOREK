#!/bin/bash
# Location of 'util' directory that contains 'setinput.sh'
UTILDIR=/scratch/jorekgl/util
# Location of target directory where executables will be copied
EXEDIR=/scratch/jorekgl

# List of executables
# The n-th executable has three parameters : 
#      model[n] ntor[n] nplane[n]
declare -a model=( 199 199 199 302 302 302)
declare -a ntor=(  1   3   7   1   3   7  )
declare -a nplane=(1   4   8   1   4   8  )

declare -a list_idx=(     0        3      )
declare -a list_inputs=(  "in199"  "in302")

# Copy the executables to target directory.
# Reference executables have a _ref suffix,
# whereas new executables have a _new suffix.
for (( j = 0 ; j < ${#list_idx[@]} ; j++ )); do
   i=${list_idx[$j]}
   EXE=j${model[$i]}_${ntor[$i]}
   REFDIR=$EXEDIR/mod${model[$i]}ref
   NEWDIR=$EXEDIR/mod${model[$i]}new
   echo $i $EXE
   if [ ! -d $REFDIR ]; then
     mkdir $REFDIR
   fi
   if [ ! -f $EXEDIR/${list_inputs[$j]} ]; then
       echo "Pb $EXEDIR/${list_inputs[$j]} not found"
   fi
   cp $EXEDIR/${list_inputs[$j]} $REFDIR
   cp $EXEDIR/${EXE}_ref $REFDIR
   if [ ! -d $NEWDIR ]; then
     mkdir $NEWDIR
   fi
   cp $EXEDIR/${list_inputs[$j]} $NEWDIR
   cp $EXEDIR/${EXE}_new $NEWDIR
done

# First Case: model199
((j=0))
i=${list_idx[$j]}
for k in ref new; do
  cd $EXEDIR/mod${model[$i]}${k}
  rm -f macroscopic_vars.dat
  INFILE=${list_inputs[$j]}

  ${UTILDIR}/setinput.sh ${INFILE} restart=.f. nstep=0
  EXE=j${model[$i]}_1_${k}
  mpiexec omplace ../${EXE} < ${INFILE} | tee out_equil
  cp jorek_restart.rst jorek_equil.rst

#  cp jorek_equil.rst jorek_restart.rst
  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. nstep=180 tstep=1000
  EXE=j${model[$i]}_3_${k}
  mpiexec omplace ../${EXE} < ${INFILE} | tee out_loop
  ${UTILDIR}/extract_live_data.sh energies energies.dat 
done
 

# Second Case: model302
((j=1))
i=${list_idx[$j]}
for k in ref new; do
  cd $EXEDIR/mod${model[$i]}${k}
  rm -f macroscopic_vars.dat
  INFILE=${list_inputs[$j]}

  ${UTILDIR}/setinput.sh ${INFILE} restart=.f. nstep=0 
  EXE=j${model[$i]}_1_${k}
  mpiexec omplace ../${EXE} < ${INFILE} | tee out_equil
  cp jorek_restart.rst jorek_equil.rst

#  cp jorek_equil.rst jorek_restart.rst
  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 200' 'tstep_n= 2' nout=10
  EXE=j${model[$i]}_3_${k}
  mpiexec omplace ../${EXE} < ${INFILE} | tee out_loop
  ${UTILDIR}/extract_live_data.sh energies energies.dat 
done


