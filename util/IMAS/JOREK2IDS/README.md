# Writing IMAS MHD IDS

PyQt5 is required for selecting input directory of HDF5 files

IMAS needs to be set up with

    module load imasenv
    imasdb jorek
    imasdb

Sample run can be saved with

     python3 jorekHDF5toIDS.py --shot=303 --run=1 --user=${USER} --database=jorek --occurrence=0


after selecting input directory such as inxflow_nper6


MHD IDS is written under ~/public/imasdb/jorek/3/0/ as 

     ids_3030001.characteristics  ids_3030001.datafile  ids_3030001.tree


MHD IDS can be visualised with 
https://git.iter.org/projects/VIS/repos/paraview-ggd-plugin under 
https://git.iter.org/projects/BND/repos/smiter or ParaView/5.8.0+

 
