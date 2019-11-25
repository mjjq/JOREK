#!/bin/bash

SCRIPTDIR=`dirname $0`; SCRIPTDIR=`readlink -f $SCRIPTDIR`

 itv=$1
 nvmax=$2
 idelta=$3

 rm -f Average_Terms_Voriticity_Gauss.dat
 rm -f Average_Terms_Voriticity_Gauss_IT.dat

 echo "It" > temp.dat

 while [ $itv -le $nvmax ]
 do

   if [ $itv -lt 10 ]
   then
     itv=0000$itv
   else
     if [ $itv -lt 100 ]
     then
       itv=000$itv
     else
       if [ $itv -lt 1000 ]
       then
         itv=00$itv
       else
         if [ $itv -lt 10000 ]
         then
           itv=0$itv
         fi
       fi
     fi
   fi

   vtkfile='jorekGauss'$itv'.vtk'
   
   # ''Intelligently'' determine whether to use .rst or .h5
   . ${SCRIPTDIR}/detect_rst_type.sh
   if [ "$RST_TYPE" != "h5" ] && [ "$RST_TYPE" != "rst" ]; then
     echo "ERROR: RST_TYPE not detected properly: $RST_TYPE"
     stop
   fi
   ext=$RST_TYPE
   
   echo $rstfile$ext

   cp $rstfile jorek_restart$ext

   ./jorek2vtk_GaussVortTerms < INPUT

   cp jorek_tmp.vtk $vtkfile
   
   echo $itv >> temp.dat

   itv=`echo "$itv + $idelta" | bc`

  done
  
  paste temp.dat Average_Terms_Voriticity_Gauss.dat > Average_Terms_Voriticity_Gauss_IT.dat
  rm temp.dat
