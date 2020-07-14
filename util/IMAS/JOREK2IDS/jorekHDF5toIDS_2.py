#!/usr/bin/env python

#   Name : jorekHDF5toIDS_bezier.py
#
#   Description :
#       A script which reads a single JOREK HDF5 output file and writes its
#       contents (Bezier finite elements: 'x', 'values', 'vertex', 'size') to
#       MHD IDS.
#
#       JOREK HDF5 file contents:
#           https://www.jorek.eu/wiki/doku.php?id=hdf5-tools&s[]=hdf5
#
#   Requirements:
#       - pip3 install --user h2py
#       - IMAS
#
#
#   Author :
#       Dejan Penko
#   E-mail :
#       dejan.penko@lecad.fs.uni-lj.si
#
# *****************************************************************************
#     Copyright(c) 2020- D. Penko

from os import listdir, getenv
from os.path import isfile, join
import numpy as np
from time import time
import sys
import platform
import numpy as np
from PyQt5 import QtWidgets
import h5py
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import datetime
from xml.etree import ElementTree
from xml.dom import minidom

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
                         "H5 Files (*.h5);;HDF5 Files (*.hdf5)",
                         options=options)
    # Close application
    app.exit()

    print("Selected VTK file (full path): ", filePath)

    return filePath


def getHDF5File(filePath):
    hdf5_file = h5py.File(filePath, 'r')
    return hdf5_file


def checkArguments():
    """ Check arguments when running from the terminal.
    """

    if (len(sys.argv) > 1):
        import argparse
        from argparse import RawTextHelpFormatter
        description = """Utility for storing JOREK output stored in HDF5 files
to IMAS (IDSs). Example command:

>>> python3 jorekHDF5toIDS_bezier.py --shot=303 --run=1 --user=penkod --database=jorek --occurrence=0


"""

        parser = argparse.ArgumentParser(description=description,
                                         formatter_class=RawTextHelpFormatter)

        parser.add_argument("-s", "--shot", type=int, required=True,
                            help="Case parameter: shot")
        parser.add_argument("-r", "--run", type=int, required=True,
                            help="Case parameter: run")
        parser.add_argument("-u", "--user", type=str, required=False,
                            help="Case parameter: username",
                            default=getenv("USER"))
        parser.add_argument("-d", "--database", type=str, required=True,
                            help="Case parameter: database")
        parser.add_argument("-o", "--occurrence", type=int, required=False,
                            help="Case parameter: occurrence",
                            default=0)

        args = parser.parse_args()
        IDS_parameters = {"shot": args.shot,
                          "run": args.run,
                          "user": args.user,
                          "database": args.database,
                          "occurrence": args.occurrence}
    else:
        # Default parameters
        print("Using default parameters")
        IDS_parameters = {"shot": 303,
                          "run": 1,
                          "user": "penkod",
                          "database": "jorek",
                          "occurrence": 0}

    return IDS_parameters


def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def setIDSProperties(mhd, hdf5_file):

    # Data Dictionary entry for homogeneous_time:
    #  This node must be filled (with 0, 1, or 2) for the IDS to be valid.
    #  If 1, the time of this IDS is homogeneous, i.e. the time values for
    #  this IDS are stored in the time node just below the root of this IDS.
    #  If 0, the time values are stored in the various time fields at lower
    #  levels in the tree. In the case only constant or static nodes are filled
    #  within the IDS, homogeneous_time must be set to 2
    mhd.ids_properties.homogeneous_time = 1  # will fill in e.g. mhd.ggd[:].time

    mhd.time.resize(1)  # mandatory

    print(f"Model: {hdf5_file['jorek_model'][0]}")

    mhd.ids_properties.comment = f"model{hdf5_file['jorek_model'][0]}"
    mhd.ids_properties.source = "JOREK"
    mhd.ids_properties.provider = getenv("USER")
    mhd.ids_properties.creation_date = str(datetime.datetime.now())
    mhd.ids_properties.version_put.data_dictionary = getenv("IMAS_VERSION")
    mhd.ids_properties.version_put.access_layer = getenv("UAL_VERSION")
    mhd.ids_properties.version_put.access_layer_language = \
        "Python " + platform.python_version()


def setCodeContents_unfinished(mhd, file=None):
    """Write code (input) parameters e.g. contents of the "intear" file of
    model199.
    """

    mhd.code.name = "jorekHDF5toIDS_2"
    mhd.code.repository = "https://git.iter.org/projects/STAB/repos/jorek"

    # Set code parameters (XML -> string)
    root = Element('root')
    root.append(Comment('Test comment'))

    head = SubElement(root, 'head')
    title = SubElement(head, 'title')
    title.text = 'Code parameters'

    generated_on = str(datetime.datetime.now())
    dc = SubElement(head, 'dateCreated')
    dc.text = generated_on
    dm = SubElement(head, 'dateModified')
    dm.text = generated_on

    mhd.code.parameters = prettify(root)


def setBezierGrid(mhd, ti, hdf5_file):
    """Set Bezier grid based on 'x', 'vertex' and 'size' matrices
    (found in HDF5 file).

    ti ... time slice
    """

    x = hdf5_file['x']
    x_T = np.array(x).transpose()
    size = hdf5_file['size']
    size_T = np.array(size).transpose()
    vertex = hdf5_file['vertex']
    vertex_T = np.array(vertex).transpose()
    values = hdf5_file['values']

    t_now = hdf5_file['t_now'][0]
    tstep = hdf5_file['tstep'][0]
    jorek_model = hdf5_file['jorek_model'][0]
    n_elements = hdf5_file['n_elements'][0]
    n_nodes = hdf5_file['n_nodes'][0]
    n_degrees = hdf5_file['n_degrees'][0]
    n_tor = hdf5_file['n_tor'][0]

    print("x shape: ", hdf5_file['x'].shape)
    print("x_T shape: ", x_T.shape)
    print("size shape: ", hdf5_file['size'].shape)
    print("size_T shape: ", size_T.shape)
    print("vertex shape: ", hdf5_file['vertex'].shape)
    print("vertex_T shape: ", vertex_T.shape)
    print("values shape: ", hdf5_file['values'].shape)

    print(f"t_now: {t_now}")
    print(f"tstep: {tstep}")
    print(f"jorek_model: {jorek_model}")
    print(f"n_elements: {n_elements}")
    print(f"n_nodes: {n_nodes}")
    print(f"n_degrees: {n_degrees}")
    print(f"n_tor: {n_tor}")

    g = mhd.grid_ggd[ti]
    g.identifier.name = "JOREK grid"
    g.index = 1
    g.description = "JOREK grid (Bezier finite elements)"
    g.time = hdf5_file['t_now'][0]

    # R: mu = 0
    # Z: mu = 1
    # p_k: dof = 0
    # u_k: dof = 1
    # v_k: dof = 2
    # w_k: dof = 3
    # x[mu][dof][k]
    # g.space[0].objects_per_dimension[0].object[k].geometry =
    #                                              [p_k_R, u_k_R, v_k_R, w_k_R,
    #                                               p_k_Z, u_k_Z, v_k_Z, w_k_Z]
    g.space.resize(2)

    # The coordinate types are specified in Data Dictionary (DD)
    # utilities/coordinate_identifier.xml
    # Four dofs: # p_k, u_k, v_k, w_k
    g.space[0].coordinates_type.resize(2)
    g.space[0].coordinates_type = np.array([4, 3])  # R, Z

    g.space[0].identifier.name = "Space R-Z"
    g.space[0].identifier.index = 1  # Fortran notation required
    g.space[0].identifier.description = """Space R-Z, four degrees of freedom
per coordinate: p_k, u_k, v_k and w_k.

- Information relevant to nodes is being stored to
grid_ggd[ti].space[0].objects_per_dimension[0].

- Information relevant to 2D elements (quads) is being stored to
grid_ggd[ti].space[0].objects_per_dimension[2].

- The dofs for each node 'k' are stored as
grid_ggd[ti].space[0].objects_per_dimension[0].object[k].geometry =
[p_k_R, u_k_R, v_k_R, w_k_R, p_k_Z, u_k_Z, v_k_Z, w_k_Z].

- The connectivity array for 2D elements (quads) is being stored as
grid_ggd[ti].space[0].objects_per_dimension[2].object[k].nodes =
[k_0, k_1, k_2, k_3]
where k_0, k_1, k_2, k_3 are indices of nodes forming this 2D element.

- The element properties - measures for the distances of the control points from
the element nodes d_{u,k}, d_{v_k}, are being stored as (single array):

grid_ggd[ti].space[0].objects_per_dimension[2].object[el].geometry =
[1.0, d_{u_0}, d_{v_0}, d_{u_0, v_0},
 1.0, d_{u_1}, d_{v_1}, d_{u_1, v_1},
 1.0, d_{u_2}, d_{v_2}, d_{u_2, v_2},
 1.0, d_{u_3}, d_{v_3}, d_{u_3, v_3}].

Note: A bit inconvenient because of the GGD structure limitations and
'geometry' node can hold only 1D array of float values.

"""
    g.space[0].objects_per_dimension.resize(3)
    s0opd0 = g.space[0].objects_per_dimension[0]
    s0opd0.object.resize(hdf5_file['n_nodes'][0])
    for k in range(n_nodes):
        s0opd0.object[k].geometry.resize(n_degrees*2)
        for dof in range(n_degrees):
            s0opd0.object[k].geometry = np.ravel(x_T[k], order='F')

    # dummy object_per_dimension[1]
    # Note: leaving mid AOS (Arrays of Structures) empty might result in
    #       the put command skipping all the next/following AOS altogether
    g.space[0].objects_per_dimension[1].object.resize(1)

    # vertex (connectivity array for 2D elements)
    s0opd2 = g.space[0].objects_per_dimension[2]
    s0opd2.object.resize(n_elements)
    for i_element in range(n_elements):
        s0opd2.object[i_element].nodes.resize(vertex_T.shape[1])
        s0opd2.object[i_element].nodes = vertex_T[i_element]

        # Data Dictionary entry for geometry:
        #  Geometry data associated with the object. Its dimension depends
        #  on the type of object, geometry and coordinate considered.
        # s0opd2.object[i_element].geometry.resize(size.shape[0]*size.shape[1])
        s0opd2.object[i_element].geometry.resize(size_T[0].size)
        # Note: storing "size" is problematic, as there are, for example,
        #       4 values per 4 nodes of the element -> 16 values
        # IDEA:
        #   s0opd2.object[i_element].geometry = [1.0, d_u1, d_v1, (d_u1 d_v1),
        #                                       [1.0, d_u2, d_v2, (d_u2 d_v2),
        #                                       [1.0, d_u3, d_v3, (d_u3 d_v3),
        #                                       [1.0, d_u4, d_v4, (d_u4 d_v4),]
        s0opd2.object[i_element].geometry = np.ravel(size_T[0], order='C')

    print("Distances of the control points from the element node as set to "
          "space[0].objects_per_dimension[2].object[0]: \n", s0opd2.object[0])
    print("space[0].objects_per_dimension[0].object[0].geometry: \n",
          s0opd0.object[0].geometry)
    print("space[0].objects_per_dimension[2].object[0].geometry: \n",
          s0opd2.object[0].geometry)


if __name__ == "__main__":

    # Set mandatory arguments
    IDS_parameters = checkArguments()

    # FIRST FILE:
    print("READING FIRST HDF5 FILE")
    filePath = getHDF5FileDialog()

    print(filePath.split("/"))
    fileName = filePath.split("/")[-1]
    print(fileName)
    fileDirPath = filePath.split(fileName)[0]
    print(fileDirPath)

    filePathList = []
    for file in listdir(fileDirPath):
        if file.startswith("jorek0") and file.endswith(".h5"):
            filePathList.append(join(fileDirPath, file))

    # Sort list of files alphabetically
    filePathList = sorted(filePathList)

    print("PREPARING IDS")
    shot = IDS_parameters['shot']
    run = IDS_parameters['run']
    username = IDS_parameters['user']
    database = IDS_parameters['database']
    # occurrence = IDS_parameters['occurrence']  # Not yet implemented

    b_ids = basicIDS(shot, run, username, database)
    b_ids.createIMASdatabase()

    # MHD ids
    mhd = b_ids.imas_obj.mhd
    # Physics Data Dictionary
    pdd = b_ids.imas_obj

    # w_ids = writeIDS(shot, run, username, database) # TODO
    comment = "WORK IN PROGRESS; Written part of the output results of multiple JOREK output " \
              "HDF5 files/timeslices ('x', 'size', 'vertex', values'), to MHD IDS."

    print("WRITING TO IDS")

    # Loop through the list of HDF5 files
    for ti in range(len(filePathList)):
    # for ti in range(1):  # hardcodded for taking only first slice for testing purposes

        print("Slice: ", ti)

        hdf5_file = getHDF5File(filePathList[ti])

        print("keys: ", hdf5_file.keys())

        # variables = [0, 1, 2, 3, 4, 5, 6]
        var_names = ["psi", "u", "j", "w", "rho", "T", "v_par"]

        setIDSProperties(mhd, hdf5_file=hdf5_file)

        # set code parameters to IDS
        setCodeContents_unfinished(mhd)

        mhd.grid_ggd.resize(len(filePathList))
        t1 = time()
        setBezierGrid(mhd, ti=ti, hdf5_file=hdf5_file)
        mhd.time.resize(len(filePathList))
        mhd.time[ti] = hdf5_file['t_now'][0]
        # t2 = time()
        print(f"Time required to set grid_ggd/ggd for time slice {ti}: {time()-t1:.2f}s")

    t1 = time()
    mhd.put()
    print(f"Time required for put() command to finish: {time()-t1:.2f}s")

    pdd.close()
