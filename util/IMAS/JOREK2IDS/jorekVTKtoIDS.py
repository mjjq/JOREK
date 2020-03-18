#!/usr/bin/env python

#  Name : jorekHDF5toIDS.py
#
#  Description :
#           A script which reads a single JOREK VTK output (ASCII. With
#           BINARY only a single data field could be read) file and writes its
#           contents (grid geometry, data fields) to mhd IDS.
#
#           Requirements:
#           - pip3 install --user vtk
#           - pip3 install --user pyvtk
#           - pip3 install --user pyqt5
#           - IMAS
#
#  Author :
#         Dejan Penko
#  E-mail :
#         dejan.penko@lecad.fs.uni-lj.si
#
#****************************************************
#     Copyright(c) 2020- D. Penko

import numpy as np
import math
import logging

from vtkUtilities import vtkReadUtilities, convert
from idsUtilities import basicIDS, writeIDS

def setVtkReaderData():
    """Set and return vtk readerData.
    """

    # Object for "VTK utilities" module. The constructor takes care of the basic
    # procedures such as file selection etc.
    vru = vtkReadUtilities()

    # Populate readerData object with representation of the VTK grid geometry
    vru.readVtkUnstructuredGrid()
    # Populate readerData object with representation of the VTK field data
    # (as a dictionary)
    vru.readVtkDataFields()

    # Add checks
    num_obj_2D = len(vru.readerData.obj_2D_list)
    if num_obj_2D == 0:
        logging.warning(" The list of 2D elements is empty!")
        return
    num_obj_3D = len(vru.readerData.obj_3D_list)
    if num_obj_3D > 0:
        logging.warning(" The grid is 3D!")
        return

    return vru.readerData

if  __name__ == "__main__":

    import imas

    # FIRST FILE:
    print("READING FIRST VTK FILE")
    readerData = setVtkReaderData()

    print("-------------------------")

    # # SECOND FILE:
    # print("READING SECOND VTK FILE")
    # readerData_2 = setVtkReaderData()

    # Quick comparison of the geometry/coordinated of the first point of each
    # grid slice
    print("First grid slice: geometry/coordinates of the first point: ",
          readerData.points_geo[0])

    shot = 1001
    run = 1
    username = 'penkod'
    database = 'jorek'

    # b_ids = basicIDS(shot, run, username, database)
    # b_ids.createIMASdatabase()

    w_ids = writeIDS(shot, run, username, database)
    comment = "Written JOREK output VTK file."
    # w_ids.createDatasetDescriptionIDS(user = 'penkod',
    #                                   source = 'JOREK',
    #                                   comment = '')
    w_ids.createGridGGD('mhd', 1)

    # Changing arrays notation from Python (starting with 0) to
    # Fortran notation (starting with 1)
    obj_0D_list_f90 = np.array(readerData.obj_0D_list)
    obj_0D_list_f90 = obj_0D_list_f90 + 1
    obj_1D_list_f90 = np.array(readerData.obj_1D_list)
    obj_1D_list_f90 = obj_1D_list_f90 + 1
    obj_2D_list_f90 = np.array(readerData.obj_2D_list)
    obj_2D_list_f90 = obj_2D_list_f90 + 1
    obj_3D_list_f90 = np.array(readerData.obj_3D_list)
    obj_3D_list_f90 = obj_3D_list_f90 + 1

    # Write grid geometry
    w_ids.writeMeshToSlice(points_geo = readerData.points_geo,
                           obj_0D_list = obj_0D_list_f90,
                           obj_1D_list = obj_1D_list_f90,
                           obj_2D_list = obj_2D_list_f90,
                           obj_3D_list = obj_3D_list_f90,
                           n = 0,
                           label = 'JOREK output VTK file grid')

    # Resize GGD
    w_ids.imas_obj.mhd.ggd.resize(1)
    # Set empty IDS path for quantity tree node
    IDSQuantityPath = None

    print("Dictionary of quantity arrays: ", readerData.field_dict)
    for key in readerData.field_dict:
        print("Current quantity array: ", key)
        if key == 'Flux':
            # Resize/Allocate
            w_ids.imas_obj.mhd.ggd[0].psi.resize(1)
            # Set IDS path
            IDSQuantityPath = w_ids.imas_obj.mhd.ggd[0].psi[0]
            # Write quantities
            w_ids.ggdWriteQuantityArray(IDSQuantityPath, readerData.field_dict[key], 1)
        elif key == 'Potential':
            # Resize/Allocate
            w_ids.imas_obj.mhd.ggd[0].phi_potential.resize(1)
            # Set IDS path
            IDSQuantityPath = w_ids.imas_obj.mhd.ggd[0].phi_potential[0]
            # Write quantities
            w_ids.ggdWriteQuantityArray(IDSQuantityPath, readerData.field_dict[key], 1)
        elif key == 'Current':
            print("Writing quantity array: ", key)
            # Resize/Allocate
            w_ids.imas_obj.mhd.ggd[i_slice].j_tor.resize(1)
            # Set IDS path
            IDSQuantityPath = w_ids.imas_obj.mhd.ggd[i_slice].j_tor[0]
            # Write quantities
            w_ids.ggdWriteQuantityArray(IDSQuantityPath, f.val[i,:], 1)
        elif key == 'Vorticity':
            # Resize/Allocate
            w_ids.imas_obj.mhd.ggd[0].vorticity.resize(1)
            # Set IDS path
            IDSQuantityPath = w_ids.imas_obj.mhd.ggd[0].vorticity[0]
            # Write quantities
            w_ids.ggdWriteQuantityArray(IDSQuantityPath, readerData.field_dict[key], 1)
        elif key == 'Density':
            # # Resize/Allocate
            # w_ids.imas_obj.mhd.ggd[0].electrons.density.resize(1)
            # # Set IDS path
            # IDSQuantityPath = w_ids.imas_obj.mhd.ggd[0].electrons.density[0]
            # # Write quantities
            # w_ids.ggdWriteQuantityArray(IDSQuantityPath, readerData.field_dict[key], 1)
            pass
        elif key == 'Temperature':
            # Resize/Allocate
            w_ids.imas_obj.mhd.ggd[0].electrons.temperature.resize(1)
            # Set IDS path
            IDSQuantityPath = w_ids.imas_obj.mhd.ggd[0].electrons.temperature[0]
            # Write quantities
            w_ids.ggdWriteQuantityArray(IDSQuantityPath, readerData.field_dict[key], 1)
        elif key == 'V_parallel':
            # Resize/Allocate
            w_ids.imas_obj.mhd.ggd[0].velocity_parallel.resize(1)
            # Set IDS path
            IDSQuantityPath = w_ids.imas_obj.mhd.ggd[0].velocity_parallel[0]
            # Write quantities
            w_ids.ggdWriteQuantityArray(IDSQuantityPath, readerData.field_dict[key], 1)

    w_ids.ids.put()

    w_ids.idsClose()



