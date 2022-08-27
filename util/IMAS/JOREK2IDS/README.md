# Writing IMAS MHD IDS

IMAS needs to be set up with

    module load IMAS # at least 3.37 required
    imasdb jorek
    imasdb

## Preparing python environment


	python3 -m venv local
	local/bin/python -m pip install --upgrade pip
	local/bin/python -m pip install h5py VTK pytz 
	local/bin/python jorekHDF5toIDS.py --help

Sample run can be saved with

    local/bin/python3 jorekHDF5toIDS.py --shot=303 --run=1 --user=${USER} --database=jorek --occurrence=0 jorek00000.h5


and created for visualisation with ParaView using

	local/bin/python3 IDS_to_VTK.py --shot=303 --run=1 --user=${USER} --database=jorek --occurrence=0
	paraview jorek..vtu


MHD IDS is written under `~/public/imasdb/jorek/3/0/` as 

    ids_3030001.characteristics  ids_3030001.datafile  ids_3030001.tree

## Generating VTK

    local/bin/python IDS_to_VTK.py
    module load ParaView
    paraview jorek..vtu

If you start typing sub in ParaView Advanced Properties pane, the 
`NonLinear Subdivision Level` can be increased from 1 to 3 
to get recomputted at finer FEM.


## Preparing IMAS and h5py with SMITER 1.6.4

    cd ~/smiter
    make env_launch.sh
    source env_launch.sh
    cd ~/jorek/util/IMAS/JOREK2IDS
    python -m venv mypy
    mypy/bin/pip install --upgrade pip
    mypy/bin/pip install h5py
    mypy/bin/python jorekHDF5toIDS.py -h
    cp jorek00000.h5 /tmp/jorek_restart.h5
    mypy/bin/python jorekHDF5toIDS.py
    mypy/bin/python IDS_to_VTK.py
    
     
 
