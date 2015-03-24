#!/bin/bash
if [ $# -ne 1 ]; then
cat <<EOF
Usage:                                                      
# Example for the model 199
# Compute the reference simulation in the directory sim199
export BASEDIR="/scratch/login/"
export PRERUN="export OMP_NUM_THREADS=8" 
export MPIRUN="mpiexec -machinefile $PBS_NODEFILE -n 4"

./launch_test.sh sim199 

EOF
exit 1
fi

# ---------- Set Suffix of restart file
function setsuf () {
    RST_HDF5=0
    RST_HDF5=`grep rst_hdf5 $INFILE | tr -s '=' ' ' |  awk '{print $2}'`
    if [ "$RST_HDF5" == "1" ]; then
	echo " Restart from HDF5 files "
	SUF="h5"
    else 
	echo " Restart from BINARY files "
	SUF="rst"
    fi
}
SUF="rst"
# ---------- Following variables has to be changed by the user ----------

# Trunk of jorek to be used
TRKDIR=$(readlink -f `dirname $0`)/..
# Location of 'util' directory that contains 'setinput.sh'
UTILDIR=${TRKDIR}/util
# Location of namelist directory that contains input jorek files
NAMEDIR=${TRKDIR}/namelist
# Location of directory that contains executables of Jorek
EXEDIR=${TRKDIR}
# Location of target directory where simulation directory will be created
if [ -z "$BASEDIR" ]; then
  BASEDIR=${TRKDIR}
fi
# Mpirun command
if [ -z "$PRERUN" ]; then
  PRERUN="sort -u $OAR_NODEFILE > mach"
fi
if [ -z "$MPIRUN" ]; then
  MPIRUN="mpiexec -launcher ssh -launcher-exec oarsh -f mach -iface ib0 -n 4"
fi


# List of executables used during following simulations
# The n-th executable has two parameters : 
#      model[n] ntor[n]
declare -a model=( 199 199 199 302 302 303 303 303)
declare -a ntor=(  1   3   7   1   3   1   3   7  )

declare -a list_idx=(     0                  3                  5               )
declare -a list_inputs=(  "in199"            "xx302"            "xx303"         )
declare -a list_orig=(    "model199/intear"  "model300/inxflow" "model303/jet_a")

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
# Detect if using murge
FIRSTLETTERS=${SIMNAME:0:2}
FIRSTLETTERS=$(echo "${FIRSTLETTERS}" | tr '[:lower:]' '[:upper:]')

if [ "${FIRSTLETTERS}" = "MU" ]; then
    SET_MURGE="use_murge=.t. use_murge_element=.t."
else
    SET_MURGE="use_murge=.f. use_murge_element=.f."
fi

UPPER_SIMNAME=$(echo "${SIMNAME}" | tr '[:lower:]' '[:upper:]')
echo ${UPPER_SIMNAME} | grep STARPU >/dev/null 2>&1
if [ $? -eq 0 ]
then
    SET_MURGE="${SET_MURGE} murge_with_starpu=.t."
fi
echo ${UPPER_SIMNAME} | grep NCUDA >/dev/null 2>&1
if [ $? -eq 0 ]
then
    NCUDA=$(echo ${UPPER_SIMNAME} | sed -e 's/.*NCUDA\([0-9]*\).*/\1/g')
    SET_MURGE="${SET_MURGE} murge_cuda_nbr=${NCUDA}"
fi

echo ${UPPER_SIMNAME} | grep NTOR >/dev/null 2>&1
if [ $? -eq 0 ]
then
    NTOR=$(echo ${UPPER_SIMNAME} | sed -e 's/.*NTOR\([0-9]*\).*/\1/g')
fi

# Model number
MODNB=${SIMNAME:$l:3}
echo "MODNB" $MODNB
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
      fi
      if [ ! -f  ${WKDIR}/${list_inputs[$j]} ]; then
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
if [ -f ${EXEDIR}/rst_bin2hdf5 ]; then
  cp  ${EXEDIR}/rst_bin2hdf5 $WKDIR
fi

# Jenkins small case: if simulation name contains jenkins199
if [ "$SIMNAME" != "${SIMNAME/JENKINS199/}" ]; then
  for c in $(seq 0 2); do
    ((j=0))
    i=${list_idx[$j]}
    WKDIR=${BASEDIR}/${PREFIX}${model[$i]}
    cd $WKDIR
    LASTNUM=$(ls  out_loop* | cut -c9- | sort -rn | head -1)
    if [ -z "$LASTNUM" ]; then
      ((LASTNUM=0))
    else
      ((LASTNUM++))
    fi
    cp macroscopic_vars.dat old_macros_vars.dat 2>/dev/null
    rm -f macroscopic_vars.dat 2>/dev/null
    INFILE=${list_inputs[$j]}
    setsuf
    ((RESTARTNB=89))
    ((RESTARTP1=RESTARTNB+1))
    TGRST=jorek000${RESTARTNB}.${SUF}
    EXPRST=jorek000${RESTARTNB}_export.${SUF}
    if [ -f $EXPRST  ]; then
      cp $EXPRST $TGRST
    fi
    if [ ! -f jorek_equil.${SUF} -a  ! -f ${TGRST} ]; then
        ${UTILDIR}/setinput.sh ${INFILE} restart=.f. nstep_n=0 tstep_n=1 n_flux=35 n_tht=14 ${SET_MURGE}
        EXE=j${model[$i]}_1
        eval ${PRERUN}
        ${MPIRUN} ./${EXE} < ${INFILE} | tee out_equil
        cp jorek_restart.${SUF} jorek_equil.${SUF}
    else if [ ! -f ${TGRST} ]; then
        cp jorek_equil.${SUF} jorek_restart.${SUF} 
        ${UTILDIR}/setinput.sh ${INFILE} restart=.t. nstep_n=89 tstep_n=1000  n_flux=35 n_tht=14 ${SET_MURGE}
        EXE=j${model[$i]}_3
        eval ${PRERUN}
        ${MPIRUN} ./${EXE} < ${INFILE} | tee out_loop${LASTNUM}        
        ${UTILDIR}/extract_live_data.sh energies energies.dat 
    else 
	cp ${TGRST} jorek_restart.${SUF}
        ${UTILDIR}/setinput.sh ${INFILE} restart=.t. nstep_n=1 tstep_n=1000  n_flux=35 n_tht=14 ${SET_MURGE}
        EXE=j${model[$i]}_3
        eval ${PRERUN}
        ${MPIRUN} ./${EXE} < ${INFILE} | tee out_loop${LASTNUM}
        cp ${TGRST} $EXPRST 2>/dev/null
        cp jorek000${RESTARTP1}.${SUF} jorek000${RESTARTP1}_export.${SUF} 2>/dev/null
        ${UTILDIR}/extract_live_data.sh energies energies.dat 
	exit 0
    fi
    fi
  done
  exit 0
fi 

# First Case: model199
if [ "$MODNB" = "199" ]; then
  for c in $(seq 0 2); do
    ((j=0))
    i=${list_idx[$j]}
    WKDIR=${BASEDIR}/${PREFIX}${model[$i]}
    cd $WKDIR
    LASTNUM=$(ls  out_loop* | cut -c9- | sort -rn | head -1)
    if [ -z "$LASTNUM" ]; then
      ((LASTNUM=0))
    else
      ((LASTNUM++))
    fi
    cp macroscopic_vars.dat old_macros_vars.dat 2>/dev/null
    rm -f macroscopic_vars.dat 2>/dev/null
    INFILE=${list_inputs[$j]}
    
    if [ ! -f jorek_equil.${SUF} ]; then
        ${UTILDIR}/setinput.sh ${INFILE} restart=.f. nstep_n=0 tstep_n=1 ${SET_MURGE}
        EXE=j${model[$i]}_1
        eval ${PRERUN}
        ${MPIRUN} ./${EXE} < ${INFILE} 2>&1 | tee out_equil
        cp jorek_restart.${SUF} jorek_equil.${SUF}
    else
        ${UTILDIR}/setinput.sh ${INFILE} restart=.t. nstep_n=200 tstep_n=1000 ${SET_MURGE}
        EXE=j${model[$i]}_3
        eval ${PRERUN}
        ${MPIRUN} ./${EXE} < ${INFILE} 2>&1 | tee out_loop${LASTNUM}
        
        ${UTILDIR}/extract_live_data.sh energies energies.dat 
    fi
  done
  exit 0
fi 

# Jenkins small case: if simulation name contains JENKINS302
if [ "$SIMNAME" != "${SIMNAME/JENKINS302/}" ]; then
  echo "JENKINS302  MINI TEST CASE"
  for c in $(seq 0 11); do
    ((j=1))
    i=${list_idx[$j]}
    WKDIR=${BASEDIR}/${PREFIX}${model[$i]}
    cd $WKDIR
    LASTNUM=$(ls  out_loop* | cut -c9- | sort -rn | head -1)
    if [ -z "$LASTNUM" ]; then
      ((LASTNUM=0))
    else
      ((LASTNUM++))
    fi

    cp macroscopic_vars.dat old_macros_vars.dat
    rm -f macroscopic_vars.dat
    INFILE=${list_inputs[$j]}
    setsuf
    COMMONOPT="n_flux=22 n_tht=30 n_open=7 n_leg=7 n_private=7"
    ((RESTARTNB=950))
    ((RESTARTP1=RESTARTNB+1))
    TGRST=jorek00${RESTARTNB}.${SUF}
    EXPRST=jorek00${RESTARTNB}_export.${SUF}
    if [ -f $EXPRST  ]; then
      cp $EXPRST $TGRST
    fi
    if [ ! -f jorek_equil.${SUF} ] && [ ! -f $TGRST ]; then
        ${UTILDIR}/setinput.sh ${INFILE} restart=.f. nstep_n=0 $COMMONOPT ${SET_MURGE}
        EXE=j${model[$i]}_1
        eval ${PRERUN}
        ${MPIRUN} ./${EXE} < ${INFILE} 2>err0 | tee out_equil
        status=$?; if [ $status -eq 0 ]; then
          cp jorek_restart.${SUF} jorek_equil.${SUF}
        fi
    else
        if [ ! -f jorek_rst1.${SUF} ]  && [ ! -f $TGRST ]; then
          cp jorek_equil.${SUF} jorek_restart.${SUF}
          ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 10, 9, 9, 9, 4' 'tstep_n= 1e-3, 1e-2, 1e-1, 1, 2' $COMMONOPT  nout=10 ${SET_MURGE}
          EXE=j${model[$i]}_3
          eval ${PRERUN}
          ${MPIRUN} ./${EXE} < ${INFILE} 2>err1| tee out_loop1
          status=$?; if [ $status -eq 0 ]; then
              cp jorek_restart.${SUF} jorek_rst1.${SUF};
          fi
        else
          if [ ! -f jorek_rst2.${SUF} ]  && [ ! -f $TGRST ]; then
              cp jorek_rst1.${SUF} jorek_restart.${SUF}
              ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 5, 4' 'tstep_n= 2, 5' $COMMONOPT  nout=10 ${SET_MURGE}
              EXE=j${model[$i]}_3
              eval ${PRERUN}
              ${MPIRUN} ./${EXE} < ${INFILE} 2>err2 | tee out_loop2
              status=$?; if [ $status -eq 0 ]; then
                  cp jorek_restart.${SUF} jorek_rst2.${SUF}
              fi
          else
              if [ ! -f $TGRST ]; then
                  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 100' 'tstep_n= 5' $COMMONOPT  nout=10 ${SET_MURGE}
                  EXE=j${model[$i]}_3
		  eval ${PRERUN}
                  ${MPIRUN} ./${EXE} < ${INFILE} 2>err3 | tee out_loop${LASTNUM}
	      else
		  cp $TGRST jorek_restart.${SUF}
                  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 1' 'tstep_n= 5' $COMMONOPT nout=1 ${SET_MURGE}
                  EXE=j${model[$i]}_3
		  eval ${PRERUN}
                  ${MPIRUN} ./${EXE} < ${INFILE} 2>err3 | tee out_loop${LASTNUM}
		  cp $TGRST $EXPRST 2>/dev/null
		  cp jorek00${RESTARTP1}.h5  jorek00${RESTARTP1}_export.h5 2>/dev/null	
		  exit 0
              fi
          fi
        fi
    fi
  done
  ${UTILDIR}/extract_live_data.sh energies energies.dat
  exit 0
fi

# Second Case: model302 point X
if [ "$MODNB" = "302" ]; then
  for c in $(seq 0 9); do
    ((j=1))
    i=${list_idx[$j]}
    WKDIR=${BASEDIR}/${PREFIX}${model[$i]}
    cd $WKDIR
    LASTNUM=$(ls  out_loop* | cut -c9- | sort -rn | head -1)
    if [ -z "$LASTNUM" ]; then
      ((LASTNUM=0))
    else
      ((LASTNUM++))
    fi
  
    cp macroscopic_vars.dat old_macros_vars.dat
    rm -f macroscopic_vars.dat
    INFILE=${list_inputs[$j]}
    
    if [ ! -f jorek_equil.${SUF} ]; then
        ${UTILDIR}/setinput.sh ${INFILE} restart=.f. nstep_n=0  ${SET_MURGE}
        EXE=j${model[$i]}_1
        eval ${PRERUN}
        ${MPIRUN} ./${EXE} < ${INFILE} 2>err0 | tee out_equil
        status=$?; if [ $status -eq 0 ]; then 
  	  cp jorek_restart.${SUF} jorek_equil.${SUF}
        fi
    else
        if [ ! -f jorek_rst1.${SUF} ]; then
  	  cp jorek_equil.${SUF} jorek_restart.${SUF} 
  	  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 10, 9, 9, 9, 4' 'tstep_n= 1e-3, 1e-2, 1e-1, 1, 2' nout=10 ${SET_MURGE}
  	  EXE=j${model[$i]}_1
        eval ${PRERUN}
  	  ${MPIRUN} ./${EXE} < ${INFILE} 2>err1| tee out_loop1
  	  status=$?; if [ $status -eq 0 ]; then 
  	      cp jorek_restart.${SUF} jorek_rst1.${SUF}; 
  	  fi
        else
  	  if [ ! -f jorek_rst2.${SUF} ]; then
  	      cp jorek_rst1.${SUF} jorek_restart.${SUF} 
  	      ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 5, 4' 'tstep_n= 2, 5' nout=10 ${SET_MURGE}
  	      EXE=j${model[$i]}_3
              eval ${PRERUN}
  	      ${MPIRUN} ./${EXE} < ${INFILE} 2>err2 | tee out_loop2
  	      status=$?; if [ $status -eq 0 ]; then 
  		  cp jorek_restart.${SUF} jorek_rst2.${SUF}
  	      fi
  	  else
  	      if [ ! -f jorek_rst3.${SUF} ]; then
  		  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 100' 'tstep_n= 5' nout=10 ${SET_MURGE}
  		  EXE=j${model[$i]}_3
                  eval ${PRERUN}
  		  ${MPIRUN} ./${EXE} < ${INFILE} 2>err3 | tee out_loop${LASTNUM}
  	      fi
  	  fi
        fi
    fi
  done
  ${UTILDIR}/extract_live_data.sh energies energies.dat 
  exit 0
fi

# Jenkins small case: if simulation name contains JENKINS303
if [ "$SIMNAME" != "${SIMNAME/JENKINS303/}" ]; then
  echo "JENKINS303 MINI TEST CASE"
  for c in $(seq 0 12); do
    ((j=2))
    i=${list_idx[$j]}
    WKDIR=${BASEDIR}/${PREFIX}${model[$i]}
    cd $WKDIR
    LASTNUM=$(ls  out_loop* | cut -c9- | sort -rn | head -1)
    if [ -z "$LASTNUM" ]; then
      ((LASTNUM=0))
    else
      ((LASTNUM++))
    fi
  
    cp macroscopic_vars.dat old_macros_vars.dat
    rm -f macroscopic_vars.dat
    INFILE=${list_inputs[$j]}
    setsuf
    sed -i "s/= *[^']none/= 'none'/" ${INFILE}
    COMMONOPT="n_flux=26 n_tht=28 n_open=8 n_leg=8 n_private=6"
    ((RESTARTNB=940))
    ((RESTARTP1=RESTARTNB+1))
    TGRST=jorek00${RESTARTNB}.${SUF}
    EXPRST=jorek00${RESTARTNB}_export.${SUF}
    if [ -f $EXPRST  ]; then
      cp $EXPRST $TGRST
    fi
    if [ ! -f jorek_equil.${SUF} ] && [ ! -f $TGRST ]; then
        ${UTILDIR}/setinput.sh ${INFILE} restart=.f. nstep_n=0 $COMMONOPT ${SET_MURGE}
        EXE=j${model[$i]}_1
        eval ${PRERUN}
        ${MPIRUN} ./${EXE} < ${INFILE} 2>err0 | tee out_equil
        status=$?; if [ $status -eq 0 ]; then 
  	  cp jorek_restart.${SUF} jorek_equil.${SUF}
        fi
    else
        if [ ! -f jorek_rst1.${SUF} ] && [ ! -f $TGRST ]; then
  	  cp jorek_equil.${SUF} jorek_restart.${SUF} 
  	  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 10, 10, 10, 10, 40, 40' 'tstep_n= 1e-3, 1e-2, 1e-1, 5e-1, 1, 2'  $COMMONOPT  nout=10 ${SET_MURGE}
  	  EXE=j${model[$i]}_3
          eval ${PRERUN}
  	  ${MPIRUN} ./${EXE} < ${INFILE} 2>err1| tee out_loop1
  	  status=$?; if [ $status -eq 0 ]; then 
  	      cp jorek_restart.${SUF} jorek_rst1.${SUF}; 
  	  fi
        else
  	  if [ ! -f jorek_rst2.${SUF} ] && [ ! -f $TGRST ]; then
  	      cp jorek_rst1.${SUF} jorek_restart.${SUF} 
  	      ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 20, 100' 'tstep_n= 2, 5'  $COMMONOPT  nout=10 ${SET_MURGE}
  	      EXE=j${model[$i]}_3
              eval ${PRERUN}
  	      ${MPIRUN} ./${EXE} < ${INFILE} 2>err2 | tee out_loop2
  	      status=$?; if [ $status -eq 0 ]; then 
  		  cp jorek_restart.${SUF} jorek_rst2.${SUF}
  	      fi
  	  else
              if [ ! -f $TGRST ]; then
                  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 100' 'tstep_n= 10' $COMMONOPT  nout=10 ${SET_MURGE}
                  EXE=j${model[$i]}_3
		  eval ${PRERUN}
                  ${MPIRUN} ./${EXE} < ${INFILE} 2>err3 | tee out_loop${LASTNUM}
	      else
		  cp $TGRST jorek_restart.${SUF}
                  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 1' 'tstep_n= 10' $COMMONOPT  nout=1 ${SET_MURGE}
                  EXE=j${model[$i]}_3
		  eval ${PRERUN}
                  ${MPIRUN} ./${EXE} < ${INFILE} 2>err3 | tee out_loop${LASTNUM}
		  cp $TGRST $EXPRST 2>/dev/null
		  cp jorek00${RESTARTP1}.h5  jorek00${RESTARTP1}_export.h5 2>/dev/null	
		  exit 0
             fi
  	  fi
        fi
    fi
  done
  ${UTILDIR}/extract_live_data.sh energies energies.dat 
  exit 0
fi


#Original serie
#nstep_n= 0.001 0.01 0.05 0.5 1   5   10
#tstep_n= 10    10   20   20  200 200 100
# Third Case: model303 point X
if [ "$MODNB" = "303" ]; then
  for c in $(seq 0 10); do
    ((j=2))
    i=${list_idx[$j]}
    WKDIR=${BASEDIR}/${PREFIX}${model[$i]}
    cd $WKDIR
    LASTNUM=$(ls  out_loop* | cut -c9- | sort -rn | head -1)
    if [ -z "$LASTNUM" ]; then
      ((LASTNUM=0))
    else
      ((LASTNUM++))
    fi
  
    cp macroscopic_vars.dat old_macros_vars.dat
    rm -f macroscopic_vars.dat
    INFILE=${list_inputs[$j]}
    
    if [ ! -f jorek_equil.${SUF} ]; then
        ${UTILDIR}/setinput.sh ${INFILE} restart=.f. nstep_n=0  ${SET_MURGE}
        EXE=j${model[$i]}_1
        eval ${PRERUN}
        ${MPIRUN} ./${EXE} < ${INFILE} 2>err0 | tee out_equil
        status=$?; if [ $status -eq 0 ]; then 
  	  cp jorek_restart.${SUF} jorek_equil.${SUF}
        fi
    else
        if [ ! -f jorek_rst1.${SUF} ]; then
  	  cp jorek_equil.${SUF} jorek_restart.${SUF} 
  	  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 10, 10, 10, 10, 40, 40' 'tstep_n= 1e-3, 1e-2, 1e-1, 5e-1, 1, 2' nout=10 ${SET_MURGE}
  	  EXE=j${model[$i]}_1
        eval ${PRERUN}
  	  ${MPIRUN} ./${EXE} < ${INFILE} 2>err1| tee out_loop1
  	  status=$?; if [ $status -eq 0 ]; then 
  	      cp jorek_restart.${SUF} jorek_rst1.${SUF}; 
  	  fi
        else
	    if [ -z $NTOR ]
	    then
		NTOR=3
	    fi
  	  if [ ! -f jorek_rst2.${SUF} ]; then
  	      cp jorek_rst1.${SUF} jorek_restart.${SUF} 
  	      ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 20, 100' 'tstep_n= 2, 5' nout=10 ${SET_MURGE}
  	      EXE=j${model[$i]}_$NTOR
        eval ${PRERUN}
  	      ${MPIRUN} ./${EXE} < ${INFILE} 2>err2 | tee out_loop2
  	      status=$?; if [ $status -eq 0 ]; then 
  		  cp jorek_restart.${SUF} jorek_rst2.${SUF}
  	      fi
  	  else
  	      if [ ! -f jorek_rst3.${SUF} ]; then
  		  ${UTILDIR}/setinput.sh ${INFILE} restart=.t. 'nstep_n= 100' 'tstep_n= 10' nout=10 ${SET_MURGE}
  		  EXE=j${model[$i]}_$NTOR
        eval ${PRERUN}
  		  ${MPIRUN} ./${EXE} < ${INFILE} 2>err3 | tee out_loop${LASTNUM}
  	      fi
  	  fi
        fi
    fi
  done
  ${UTILDIR}/extract_live_data.sh energies energies.dat 
  exit 0
fi

