#!/bin/bash
# Location of new version of jorek to compile
NEWDIR=~/trunk2g
# Location of reference version of jorek to compile
REFDIR=~/trunkref
# Location of 'util' directory that contains 'setconfig.sh'
UTILDIR=~/util
# Location of target directory where executables will be copied
EXEDIR=/scratch/jorekgl

# Cleaning the old executables
if [ "$1" = "clean" ]; then
  for DIR in $REFDIR $NEWDIR; do
     echo "cleaning $DIR of old executables"
     cd $DIR
     rm -f j???_? j???_??
  done
fi

# List of executable to compile 
# The n-th executable has three parameters : 
#      model[n] ntor[n] nplane[n]
#declare -a model=( 199 199 199 302 302 302)
#declare -a ntor=(  1   3   7   1   3   7  )
#declare -a nplane=(1   4   8   1   4   8  )
declare -a model=( 199 199 199)
declare -a ntor=(  1   3   7  )
declare -a nplane=(1   4   8  )

# Launch compilation and generates all executables
for DIR in $REFDIR $NEWDIR; do
  for (( i = 0 ; i < ${#model[@]} ; i++ )); do
    cat <<EOF > go.sh
TG=jorek_model${model[$i]}
EXE=j${model[$i]}_${ntor[$i]}
if [ ! -f $DIR\$EXE ]; then
  echo ============= model=${model[$i]} n_tor=${ntor[$i]} n_period=1 n_plane=${nplane[$i]}
  cd $DIR
  $UTILDIR//setconfig.sh  model=${model[$i]} n_tor=${ntor[$i]} n_period=1 n_plane=${nplane[$i]}
  if [ -f \$TG ]; then rm -f \$TG; fi
  make -f NewMake.mk cleanall
  make -f NewMake.mk timing/trace.o
  make -f NewMake.mk -j 4 | tee log_${model[$i]}_${ntor[$i]}
  if [ -f \$TG ]; then mv \$TG \$EXE; fi
fi
EOF
    chmod +x go.sh
    ./go.sh
  done
done

# Copy the executables to target directory.
# Reference executables have a _ref suffix,
# whereas new executables have a _new suffix.
for (( i = 0 ; i < ${#model[@]} ; i++ )); do
   EXE=j${model[$i]}_${ntor[$i]}
   if [ -f $REFDIR/$EXE ]; then
      cp $REFDIR/$EXE $EXEDIR/${EXE}_ref
   else
      echo "ERROR: DO NOT FOUND " $REFDIR/$EXE
   fi
   if [ -f $NEWDIR/$EXE ]; then
      cp $NEWDIR/$EXE $EXEDIR/${EXE}_new
   else
      echo "ERROR: DO NOT FOUND " $NEWDIR/$EXE
   fi
done