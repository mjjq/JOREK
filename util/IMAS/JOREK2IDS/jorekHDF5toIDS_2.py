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
import sys
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

    # Data Dictionary entry for homogeneous time:
    #  This node must be filled (with 0, 1, or 2) for the IDS to be valid.
    #  If 1, the time of this IDS is homogeneous, i.e. the time values for
    #  this IDS are stored in the time node just below the root of this IDS.
    #  If 0, the time values are stored in the various time fields at lower
    #  levels in the tree. In the case only constant or static nodes are filled
    #  within the IDS, homogeneous_time must be set to 2
    mhd.ids_properties.homogeneous_time = 0  # will fill in e.g. mhd.ggd[:].time

    mhd.time.resize(1) # mandatory

    # mhd.ids_properties.comment = ...
    mhd.ids_properties.source = "JOREK"
    mhd.ids_properties.provider = getenv("USER")
    mhd.ids_properties.creation_date = str(datetime.datetime.now())
    mhd.ids_properties.version_put = getenv("UAL_VERSION")


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


def setBezierGrid(mhd, slice, hdf5_file):
    """Set Bezier grid based on 'x', 'vertex' and 'size' matrices
    (found in HDF5 file)
    """

    g = mhd.grid_ggd[i_slice]
    g.time = hdf5_file['t_now'][0]

    # R: mu = 0
    # Z: mu = 1
    # p_k: dof = 0
    # u_k: dof = 1
    # v_k: dof = 2
    # w_k: dof = 3
    # x[mu][dof][i_node]
    # g.space[mu].objects_per_dimension[0].object[i_node].geometry = [p_k, u_k, v_k, w_k]
    g.space.resize(2)

    # The coordinate types are specified in $IMAS_PREFIX/include/cpp, though
    # the degrees of freedom (dofs) are not foreseen here.
    # IDEA: could the dofs be suggested to be put here? Using dummy indices for
    # now
    # Four dofs: # p_k, u_k, v_k, w_k
    g.space[0].coordinates_type = [1000, 1001, 1002, 1003]  # dummy IDs

    g.space[0].objects_per_dimension.resize(3)

    opd0 = g.space[0].objects_per_dimension[0]
    opd0.object.resize(hdf5_file['n_nodes'][0])






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
    # for i_slice in range(len(filePathList)):
    for i_slice in range(1):

        print("Slice: ", i_slice)

        hdf5_file = getHDF5File(filePathList[i_slice])

        print("keys: ", hdf5_file.keys())

        x = hdf5_file['x'][0]
        size = hdf5_file['size'][0]
        vertex = hdf5_file['vertex'][0]
        values = hdf5_file['values'][0]

        t_now = hdf5_file['t_now'][0]
        tstep = hdf5_file['tstep'][0]
        jorek_model = hdf5_file['jorek_model'][0]
        n_elements = hdf5_file['n_elements'][0]
        n_nodes = hdf5_file['n_nodes'][0]

        print(f"t_now: {t_now}")
        print(f"tstep: {tstep}")
        print(f"jorek_model: {jorek_model}")
        print(f"n_elements: {n_elements}")
        print(f"n_nodes: {n_nodes}")

        # variables = [0, 1, 2, 3, 4, 5, 6]
        var_names = ["psi", "u", "j", "w", "rho", "T", "v_par"]

        setIDSProperties(mhd, hdf5_file=hdf5_file)

        # set code parameters to IDS
        setCodeContents_unfinished(mhd)

        mhd.grid_ggd.resize(len(filePathList))

        setBezierGrid(mhd, slice=i_slice, hdf5_file=hdf5_file)

    mhd.put()

    pdd.close()
