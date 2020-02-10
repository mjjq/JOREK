# -*- coding: utf-8 -*-

#
# Purpose: Auxiliary script for select_restart_files.sh to select restart files for given times,
#          using the data extracted from macroscopic_vars.dat
#
# Date: 2019-11-03
# Author: Fabian Wieschollek, IPP Garching
#

try:
  import numpy as np
except:
  print("ERROR: Python module numpy is not available!")
  exit()
import os
import sys

def get_step(time):
  #Determine stepnumber from time
  temp   = np.abs(time-available_times)
  index  = np.where(temp==min(temp))[0][0]	
  return files[index]



# --- Process the arguments
if len(sys.argv)<7 :
  print("Not enough arguments!")
  exit()

full            = int(sys.argv[3])
name_times      = sys.argv[1]
name_steps      = sys.argv[2]
selected_times  = sys.argv[4]
listfile        = int(sys.argv[5])
unit            = float(sys.argv[6])*1000
ms              = int(sys.argv[7])



# --- Loads the list of times and step numbers from the existing restart files
#     Ignores restart files with indices greater than highest index in macroscopic_vars.dat
files           = np.loadtxt(name_steps,dtype=int)
times           = np.loadtxt(name_times)
available_times = times[files[files < len(times)]]



# --- Creates the list of selected times, by evaluating selected_times
times_array     = np.array([])
for selected_times_0 in selected_times.split(","):
  if "-" in selected_times_0:
    selected_times_1 = np.array(selected_times_0.split("-"),dtype=float)
    if ms:
      selected_times_1 = selected_times_1/unit
    if selected_times_1[2] > times[-1] :
      times_array = np.concatenate((times_array,np.arange(selected_times_1[0],times[-1],selected_times_1[1])))
      times_array = np.append(times_array,times[-1])
    else:
      times_array = np.concatenate((times_array,np.arange(selected_times_1[0],selected_times_1[2],selected_times_1[1])))
  else:
    times_array   = np.append(times_array,float(selected_times_0))



# --- For each selected time, the correspoding restart file is being identified
selected_files = [get_step(a) for a in times_array] 
selected_files = np.unique(selected_files) #remove possible duplicates



# --- Prints list of absolute paths of the selected restart files or only their step numbers
if full:
  for selected_file in selected_files:
    print(os.environ["sourceDir"]+"/jorek"+str(selected_file).zfill(5)+"."+os.environ["RST_TYPE"])
else:
  for selected_file in selected_files:
    print(str(selected_file).zfill(5))



# --- Stores selected step numbers and corresponding times in a file
if listfile:
  if ms:
    pre="_ms"
  else:
    pre=""
  with open("selected_files_"+selected_times+pre+".dat","w+") as f:
    f.write('"step"    "time/sqrt(mu_0rho_0)"    "time/ms"\n')
    for selected_file in selected_files:
      f.write(str(selected_file).zfill(5)+"    "+'%.7E' % times[selected_file]+"    "+'%.7E' % (times[selected_file]*unit)+"\n")
