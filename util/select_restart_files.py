# -*- coding: utf-8 -*-

#
# Purpose: Auxiliary script for select_restart_files.sh to select restart files for given times,
#          using the data extracted from macroscopic_vars.dat
#
# Date: 2019-11-03
# Author: Fabian Wieschollek, IPP Garching
#

import numpy as np
import os
import sys

def getIndex(time):

  #Determine stepnumber from time
  step=np.abs(time-TimeS)
  step=np.where(step==min(step))[0][0]

  #Determine index of file from step number
  index=np.abs(step-FileS)
  index=np.where(index==min(index))[0][0]

  return index



# Evaluates the arguments
if len(sys.argv)<5 :
  print("Not enough arguments!")
  exit()

full=int(sys.argv[3])
name_times=sys.argv[1]
name_steps=sys.argv[2]
onlyT=sys.argv[4]
listfile=int(sys.argv[5])



# Loads the list of times and step numbers from the existing restart files
FileS=np.loadtxt(name_steps,dtype=int)
TimeS=np.loadtxt(name_times)



# Creates the list of selected times, by evaluating onlyT
SelectedTimeS=np.array([])
for onlyT_0 in onlyT.split(","):
  if "-" in onlyT_0:
    onlyT_1=onlyT_0.split("-")
    if float(onlyT_1[2]) > TimeS[-1] :
      SelectedTimeS=np.concatenate((SelectedTimeS,np.arange(float(onlyT_1[0]),TimeS[-1],float(onlyT_1[1]))))
      SelectedTimeS=np.append(SelectedTimeS,TimeS[-1])
    else:
      SelectedTimeS=np.concatenate((SelectedTimeS,np.arange(float(onlyT_1[0]),float(onlyT_1[2]),float(onlyT_1[1]))))
  else:
    SelectedTimeS=np.append(SelectedTimeS,float(onlyT_0))



# For each selected time, the correspoding restart file is being identified
SelectedFileS=[FileS[getIndex(a)] for a in SelectedTimeS] 
SelectedFileS=np.unique(SelectedFileS) #remove possible duplicates



#Prints list of absolute paths of the selected restart files or only their step numbers
if full:
  for SelectedFile in SelectedFileS:
    print(os.environ["sourceDir"]+"/jorek"+str(SelectedFile).zfill(5)+"."+os.environ["RST_TYPE"])
else:
  for SelectedFile in SelectedFileS:
    print(str(SelectedFile).zfill(5))



#Stores selected step numbers and corresponding times in a file
if listfile:
  with open("selected_files_"+onlyT+".dat","w+") as f:
    f.write('"step"   "time"\n')
    for SelectedFile in SelectedFileS:
      f.write(str(SelectedFile).zfill(5)+"    "+str(TimeS[SelectedFile])+"\n")
