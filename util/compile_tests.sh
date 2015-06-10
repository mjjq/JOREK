#!/bin/bash
# Usage:                                                                
# ./compile_test.sh 199 # compile model 199 with several ntor/nplane settings
# ./compile_test.sh 302 # compile model 302 with several ntor/nplane settings
# ./compile_test.sh     # compile both models

# Options when launching make (use, e.g., -j 4 to compile in parallel)
MAKEOPT="-j 2"
# Location of 'util' directory that contains 'setconfig.sh'
UTILDIR=$(readlink -f `dirname $0`)
# Location of trunk of jorek to compile
TRKDIR="$(readlink -f $UTILDIR/..)"

# List of executable to compile 
# The n-th executable has three parameters : 
#      model[n] ntor[n] nplane[n] nperiod[n]
if [ "$1" = "199" ]; then 
   declare -a model=(  199 199 199)
   declare -a ntor=(    1   3  7)
   declare -a nplane=(  1   4  8)
   declare -a nperiod=( 1   1  1) 
fi
if [ "$1" = "302" ]; then 
   declare -a model=(  302 302 )
   declare -a ntor=(     1  3  )
   declare -a nplane=(   1  8  )
   declare -a nperiod=(  1  8  )
fi
if [ "$1" = "303" ]; then 
   declare -a model=(  303 303 303 )
   declare -a ntor=(    1   3   7  )
   declare -a nplane=(  1   8   16 )
   declare -a nperiod=( 1   2   2  )
fi
if [ "$1" = "all" ] || [ "$1" = "" ]; then 
   declare -a model=(  199 199 302 302 303 303 )
   declare -a ntor=(   1   3   1   3    1   3  )
   declare -a nplane=( 1   4   1   8    1   8  )
   declare -a nperiod=(1   1   1   8    1   2  )
fi

#-------------------------------------------------------------------------

# Cleaning the old executables
if [ "$1" = "clean" ]; then
  echo "cleaning $DIR of old executables"
  cd $TRKDIR
  rm -f j???_? j???_??
  exit 0
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
UTILDIR="$UTILDIR"
TRKDIR="$TRKDIR"
if [ ! -f \$TRKDIR/\$EXE ]; then
  echo ============= model=${model[$i]} n_tor=${ntor[$i]} n_period=1 n_plane=${nplane[$i]} n_period=${nperiod[$i]} 
  cd \$TRKDIR
  \$UTILDIR//setconfig.sh  model=${model[$i]} n_tor=${ntor[$i]} n_period=1 n_plane=${nplane[$i]} n_period=${nperiod[$i]}
  if [ -f \$TG ]; then rm -f \$TG; fi
  make cleanall
  make ${MAKEOPT} timing/trace.o  
  make ${MAKEOPT} | tee log_${model[$i]}_${ntor[$i]}
  if [ "\${PIPESTATUS[0]}" -ne 0 ]; then
    exit 1
  fi
  if [ -f \$TG ]; then mv \$TG \$EXE; fi
fi
EOF
    chmod +x go.sh
    ./go.sh || exit 1
  fi
done

