#!/bin/bash

read -r -p "Are you sure? [y/N] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then
   rm -rfv *_out energy_conservation.o*
else
   exit 0
fi
