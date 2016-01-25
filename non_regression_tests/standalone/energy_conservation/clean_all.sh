#!/bin/bash

read -r -p "Are you sure? [y/N] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then
   rm -fv */*.{png,dat,vtk,xml,dump} */jorek_log energy_conservation.o*
else
   exit 0
fi
