#!/bin/bash
# Usage:                                                                
# ./compile_test.sh 199 # compile model 199 with several ntor/nplane settings
# ./compile_test.sh 302 # compile model 302 with several ntor/nplane settings

# ---------- Following variables has to be changed by the user ----------

# Location of trunk of jorek to compile
TRKDIR=/scratch/jorekgl/trunk2g
# Location of 'util' directory that contains 'setconfig.sh'
UTILDIR=${TRKDIR}/util
# Options when launching make
MAKEOPT="-j 4"

# List of executable to compile 
# The n-th executable has three parameters : 
#      model[n] ntor[n] nplane[n] nperiod[n]
if [ "$1" = "199" ]; then 
   declare -a model=(  199 199 )
   declare -a ntor=(    1   3  )
   declare -a nplane=(  1   4  )
   declare -a nperiod=( 1   1  ) 
fi
if [ "$1" = "302" ]; then 
   declare -a model=(  302 302 )
   declare -a ntor=(    3   1  )
   declare -a nplane=(  8   1  )
   declare -a nperiod=( 8   1  )
fi
if [ "$1" = "all" ]; then 
   declare -a model=(  199 199 199 302 302 )
   declare -a ntor=(   1   3   7   1   3   )
   declare -a nplane=( 1   4   8   1   8   )
   declare -a nperiod=(1   1   1   1   8   )
fi

#-------------------------------------------------------------------------

# Cleaning the old executables
if [ "$1" = "clean" ]; then
  echo "cleaning $DIR of old executables"
  cd $TRKDIR
  rm -f j???_? j???_??
fi

if [ ${#model[@]} -eq 0 ]; then
  echo "ERROR: no model number has be given in parameter."
  echo "Try for example '199' as parameter"
  exit 0
fi

if [ -z "$2" ]; then
    noselect=1
    targetid=0
else
    noselect=0
    targetid=$2
fi

for (( i = 0 ; i < ${#model[@]} ; i++ )); do
   # if a second parameter is given  in command line, it specifies
   # one subset of a model to be compiled.
  ((selectid=(i==$targetid)))
  if ((noselect || selectid)); then
    cat <<EOF > go.sh
TG=jorek_model${model[$i]}
EXE=j${model[$i]}_${ntor[$i]}
if [ ! -f $TRKDIR\$EXE ]; then
  echo ============= model=${model[$i]} n_tor=${ntor[$i]} n_period=1 n_plane=${nplane[$i]} n_period=${nperiod[$i]} 
  cd $TRKDIR
  $UTILDIR//setconfig.sh  model=${model[$i]} n_tor=${ntor[$i]} n_period=1 n_plane=${nplane[$i]} n_period=${nperiod[$i]}
  if [ -f \$TG ]; then rm -f \$TG; fi
  make -f NewMake.mk cleanall
  make -f NewMake.mk ${MAKEOPT} timing/trace.o
  make -f NewMake.mk ${MAKEOPT} | tee log_${model[$i]}_${ntor[$i]}
  if [ -f \$TG ]; then mv \$TG \$EXE; fi
fi
EOF
    chmod +x go.sh
    ./go.sh
  fi
done

