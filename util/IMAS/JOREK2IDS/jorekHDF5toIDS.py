#!/usr/bin/env python

#  Name : jorekHDF5toIDS.py
#
#  Description :
#           A script which reads a single JOREK HDF5 output file and writes its
#           contents (grid geometry, data fields) to mhd IDS.
#
#           JOREK HDF5 file contents:
#           https://www.jorek.eu/wiki/doku.php?id=hdf5-tools&s[]=hdf5
#
#           Requirements:
#               - pip3 install --user h2py
#               - IMAS
#
#  Author :
#         Dejan Penko
#  E-mail :
#         dejan.penko@lecad.fs.uni-lj.si
#
#****************************************************
#     Copyright(c) 2020- D. Penko

from os import listdir
from os.path import isfile, join
import numpy as np
import math
import logging
import sys
import vtk
import numpy as np
from PyQt5 import QtWidgets

sys.path.append('../../../util')
import jorek_read_h5 as jorek
import imas

from idsUtilities import basicIDS, writeIDS

def getHDF5FileDialog():
    """Run as standalone application.
    """
    # Set QApplication
    app = QtWidgets.QApplication([])
    # Set options
    options = QtWidgets.QFileDialog.Options()
    options |= QtWidgets.QFileDialog.DontUseNativeDialog
    # Set file dialog
    # Note: arguments are:  - parent (None),
    #                       - window title
    #                       - default file name
    #                       - file types selectable options
    #                       - options
    filePath, _ = QtWidgets.QFileDialog \
        .getOpenFileName(None,
                         "Select VTK file",
                         "",
                         "All Files (*);;H5 Files (*.h5);;HDF5 Files (*.hdf5)",
                         options=options)
    # Close application
    app.exit()

    print("Selected VTK file (full path): ", filePath)

    return filePath

if  __name__ == "__main__":

    # FIRST FILE:
    print("READING FIRST HDF5 FILE")
    filePath = getHDF5FileDialog()

    print(filePath.split("/"))
    fileName = filePath.split("/")[-1]
    print(fileName)
    fileDirPath = filePath.split(fileName)[0]
    print(fileDirPath)

    # onlyfiles = [f for f in listdir(fileDirPath) if isfile(join(fileDirPath, f))]
    # print(onlyfiles)

    filePathList = []

    for file in listdir(fileDirPath):
        if file.startswith("jorek0") and file.endswith(".h5"):
            # print(join(fileDirPath, file))
            filePathList.append(join(fileDirPath, file))

    # Sort list of files alphabetically
    filePathList = sorted(filePathList)

    f = jorek.fields()

    print("PREPARING IDS")
    shot = 1000
    run = 4
    username = 'penkod'
    database = 'jorek'

    b_ids = basicIDS(shot, run, username, database)
    b_ids.createIMASdatabase()

    w_ids = writeIDS(shot, run, username, database)
    comment = "Written results of multiple JOREK output HDF5 files/timeslices."
    # w_ids.createDatasetDescriptionIDS(user = 'penkod',
    #                                   source = 'JOREK',
    #                                   comment = '')
    w_ids.createGridGGD('mhd', 1)

    # Resize GGD to the number of timeslices
    w_ids.imas_obj.mhd.ggd.resize(len(filePathList))

    print("WRITING TO IDS")

    # Loop through the list of HDF5 files
    for i_slice in range(len(filePathList)):

        print("Slice: ", i_slice)

        # Read file
        f.read(filePathList[i_slice], variables=[0,1,2,3,4,5,6])

        print("Time step: ", f.tstep)
        print("Time: ", f.t_now)

        # Set angle for VTK object to be interpolated (3D grid geometry
        # and data fields)
        # Set phi=0 for regular 2D slice

        # Pure 2D grid with quad elements
        # grid = f.to_vtk(phi=0, quadratic=False, n_sub=4)
        grid = f.to_vtk(phi=0, quadratic=False, n_sub=2, n_plane=1)
        # grid = f.to_vtk(phi=[0,np.pi/2], quadratic=False, n_plane=1)

        # Denser grid (n_sub)
        # grid = f.to_vtk(phi=0, quadratic=False, n_sub=8)

        # Grid geometry is taken only from the first time slice. All other
        # time slices share the same geometry. No need to write the same grid
        # multiple times as all of them have the same.
        if i_slice == 0:

            # x_coord = f.xyz[:,0]
            # y_coord = f.xyz[:,1]

            list_vertex = f.vertex
            vtk_quad_conn_array = f.ien
            # Remove the vtk cell type ID from the matrix
            quad_conn_array = vtk_quad_conn_array[:,1:]

            # print("num_coord: ", len(f.xyz))
            # print("x_coord: ", x_coord)
            # print("y_coord: ", y_coord)
            # print("list_vertex: ", list_vertex)
            # print("vtk_quad_conn_array: ", vtk_quad_conn_array)
            # print("quad_conn_array: ", quad_conn_array)

            # # Changing arrays notation from Python (starting with 0) to
            # # Fortran notation (starting with 1)
            obj_2D_list_f90 = np.array(quad_conn_array)
            obj_2D_list_f90 = obj_2D_list_f90 + 1

            # Write grid geometry
            w_ids.writeMeshToSlice(points_geo = f.xyz,
                                   obj_0D_list = [],
                                   obj_1D_list = [],
                                   obj_2D_list = obj_2D_list_f90,
                                   obj_3D_list = [],
                                   n = 0,
                                   label = 'JOREK output HDF5 file grid with quantities')

        quantity_names_list = f.var_names
        quantities_array = f.val
        # print("quantity_names_list: ", quantity_names_list)
        # print("quantities_array: ", quantities_array)


        # Set empty IDS path for quantity tree node
        IDSQuantityPath = None

        # Set time
        w_ids.imas_obj.mhd.ggd[i_slice].time = f.t_now

        print("Array of quantity labels: ", quantity_names_list)
        for i in range(len(quantity_names_list)):
            label = quantity_names_list[i]
            print("Current quantity array: ", label)
            if label == 'psi': # Flux / poloidal magnetic flux
                print("Writing quantity array: ", label)
                # Resize/Allocate
                w_ids.imas_obj.mhd.ggd[i_slice].psi.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.imas_obj.mhd.ggd[i_slice].psi[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, f.val[i,:], 1)
            elif label == 'u': # Potential / electric potential
                print("Writing quantity array: ", label)
                # Resize/Allocate
                w_ids.imas_obj.mhd.ggd[i_slice].phi_potential.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.imas_obj.mhd.ggd[i_slice].phi_potential[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, f.val[i,:], 1)
            elif label == 'j': # Current / toroidal current density
                print("Writing quantity array: ", label)
                # Resize/Allocate
                w_ids.imas_obj.mhd.ggd[i_slice].j_tor.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.imas_obj.mhd.ggd[i_slice].j_tor[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, f.val[i,:], 1)
            elif label == 'w': # Vorticity
                print("Writing quantity array: ", label)
                # Resize/Allocate
                w_ids.imas_obj.mhd.ggd[i_slice].vorticity.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.imas_obj.mhd.ggd[i_slice].vorticity[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, f.val[i,:], 1)
            elif label == 'rho': # Mass Density
                # print("Writing quantity array: ", label)
                # # Resize/Allocate
                # w_ids.imas_obj.mhd.ggd[i_slice].electrons.density.resize(1)
                # # Set IDS path
                # IDSQuantityPath = w_ids.imas_obj.mhd.ggd[i_slice].electrons.density[0]
                # # Write quantities
                # w_ids.ggdWriteQuantityArray(IDSQuantityPath, f.val[i,:], 1)
                pass
            elif label == 'T': # Temperature / total temperature
                print("Writing quantity array: ", label)
                # Resize/Allocate
                w_ids.imas_obj.mhd.ggd[i_slice].electrons.temperature.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.imas_obj.mhd.ggd[i_slice].electrons.temperature[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, f.val[i,:], 1)
            elif label == 'v_par': # V_parallel / parallel velocity
                print("Writing quantity array: ", label)
                # Resize/Allocate
                w_ids.imas_obj.mhd.ggd[i_slice].velocity_parallel.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.imas_obj.mhd.ggd[i_slice].velocity_parallel[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, f.val[i,:], 1)

                print("v_par max value: ", max(f.val[i,:]))

    w_ids.ids.put()

    w_ids.idsClose()
