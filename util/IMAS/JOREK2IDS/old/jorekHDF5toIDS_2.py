#!/usr/bin/env python

#   Name : jorekHDF5toIDS_bezier.py
#
#   Description :
#       This script reads the contents of the JOREK case run directory,
#       extracts the data from HDF5 output files and writes their
#       contents (Bezier finite elements: 'x', 'values', 'vertex', 'size') to
#       MHD IDS. Contents of each HDF5 file is being written under its own
#       time slice.
#
#       JOREK HDF5 file contents:
#           https://www.jorek.eu/wiki/doku.php?id=hdf5-tools&s[]=hdf5
#
#   Requirements:
#       - pip3 install --user h5py PyQt5
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
import logging
import inspect

from idsUtilities import basicIDS, writeIDS

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


class jorekRun2IDS():
    def __init__(self):

        # Set log parser
        self.log = logging.getLogger(__name__)
        self.log.addHandler(logging.StreamHandler())
        self.setDebugMode(0)  # Set to INFO

        # Set initial variables
        # - Path to JOREK (case) run directory
        self.rundir = None
        # - A list to hold all HDF5 files found in
        #   JOREK run directory
        self.HDF5pathList = []
        # - reference to Physics Data Dictionary (root)
        self.pdd = None
        # - reference to MHD IDS member
        self.mhd = None

    def writeLogDebug(self, instance, currentframe, msg):
        """ Print to DEBUG log.
        Arguments:
            instance     (obj) : Class instance (e.g. self).
            currentframe (obj) : Frame object for the caller’s stack frame.
            msg          (str) : Additional message (usually "START" or "END")
        """

        self.log.debug(f"DEBUG | {type(instance).__name__} | "
                       f"{currentframe.f_code.co_name}  | {msg}.")

    def writeLogError(self, instance, currentframe, msg):
        """ Print to ERROR log.
        Arguments:
            instance     (obj) : Class instance (e.g. self).
            currentframe (obj) : Frame object for the caller’s stack frame.
            msg          (str) : Additional message (usually "START" or "END")
        """

        self.log.debug(f"ERROR | {type(instance).__name__} | "
                       f"{currentframe.f_code.co_name}  | {msg}.", exc_info=True)

    def getLogger(self):
        return self.log

    def setDebugMode(self, state):
        """Set debug mode on/off.
        """
        if state:
            self.log.setLevel(logging.DEBUG)
        else:
            self.log.setLevel(logging.INFO)

    def setRunDirectory(self):

        self.writeLogDebug(self, inspect.currentframe(), "START")
        # Set QApplication (required for dialog to show)
        app = QtWidgets.QApplication([])
        # Set options
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        # get rundir
        self.rundir = str(QtWidgets.QFileDialog.getExistingDirectory(
            None, "Select JOREK run directory", options=options))

        # Close application
        app.exit()

        self.writeLogDebug(self, inspect.currentframe(), "END")

    def setHDF5filesList(self):
        self.writeLogDebug(self, inspect.currentframe(), "START")

        if self.rundir is None:
            print("No rundir set. Returning.")
            return

        for file in listdir(self.rundir):
            if file.startswith("jorek") and file.endswith(".h5"):
                self.HDF5pathList.append(join(self.rundir, file))

        # Sort list of files alphabetically
        self.HDF5pathList = sorted(self.HDF5pathList)

        if len(self.HDF5pathList) > 0:
            self.log.info(f"List of HDF5 files: {self.HDF5pathList}")
        else:
            self.log.warning(f"The list of HDF5 is EMPTY!")

        self.writeLogDebug(self, inspect.currentframe(), "END")

    def getHDF5File(self, hdf5_file_id):
        hdf5_file = h5py.File(self.HDF5pathList[hdf5_file_id], 'r')
        return hdf5_file

    def setIDS(self, shot, run, username, database):

        # set IDS object
        basic_ids = basicIDS(shot, run, username, database)

        basic_ids.createNewIMASdatabase()

        # Physics Data Dictionary (root)
        self.pdd = basic_ids.data_entry
        # MHD ids
        self.mhd = self.pdd.mhd
        # Allocate grid_ggd AOS (grid geometry description)
        self.mhd.grid_ggd.resize(len(self.HDF5pathList))
        # Allocate ggd AOS (quantities)
        self.mhd.ggd.resize(len(self.HDF5pathList))
        # Allocate MHD time array
        self.mhd.time.resize(len(self.HDF5pathList))

    def setIDSProperties(self, hdf5_file_id):

        hdf5_file = self.getHDF5File(hdf5_file_id=hdf5_file_id)

        # Data Dictionary entry for homogeneous_time:
        #  This node must be filled (with 0, 1, or 2) for the IDS to be valid.
        #  If 1, the time of this IDS is homogeneous, i.e. the time values for
        #  this IDS are stored in the time node just below the root of this IDS.
        #  If 0, the time values are stored in the various time fields at lower
        #  levels in the tree. In the case only constant or static nodes are filled
        #  within the IDS, homogeneous_time must be set to 2
        self.mhd.ids_properties.homogeneous_time = 1  # will fill in e.g. mhd.ggd[:].time

        # Allocate MHD time array
        self.mhd.time.resize(len(self.HDF5pathList))  # mandatory
        print(f"Model: {hdf5_file['jorek_model'][0]}")

        self.mhd.ids_properties.comment = f"model{hdf5_file['jorek_model'][0]}"
        self.mhd.ids_properties.source = "JOREK"
        self.mhd.ids_properties.provider = getenv("USER")
        self.mhd.ids_properties.creation_date = str(datetime.datetime.now())
        self.mhd.ids_properties.version_put.data_dictionary = getenv("IMAS_VERSION")
        self.mhd.ids_properties.version_put.access_layer = getenv("UAL_VERSION")
        self.mhd.ids_properties.version_put.access_layer_language = \
            "Python " + platform.python_version()

    def setCodeContents_unfinished(self):
        """Write code (input) parameters e.g. contents of the "intear" file of
        model199.
        """

        self.mhd.code.name = "jorekHDF5toIDS_2"
        self.mhd.code.repository = "https://git.iter.org/projects/STAB/repos/jorek"

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

        self.mhd.code.parameters = prettify(root)

    def setBezierGrid2IDS(self, hdf5_file_id):
        """Set Bezier grid based on 'x', 'vertex' and 'size' matrices
        (found in HDF5 file).

        ti ... time slice
        """

        hdf5_file = self.getHDF5File(hdf5_file_id=hdf5_file_id)

        x = hdf5_file['x']
        x_T = np.array(x).transpose()
        size = hdf5_file['size']
        size_T = np.array(size).transpose()
        vertex = hdf5_file['vertex']
        vertex_T = np.array(vertex).transpose()

        t_now = hdf5_file['t_now'][0]
        tstep = hdf5_file['tstep'][0]
        jorek_model = hdf5_file['jorek_model'][0]
        n_elements = hdf5_file['n_elements'][0]
        n_nodes = hdf5_file['n_nodes'][0]
        n_degrees = hdf5_file['n_degrees'][0]
        n_tor = hdf5_file['n_tor'][0]

        print("x shape: ", x.shape)
        print("x_T shape: ", x_T.shape)
        print("size shape: ", size.shape)
        print("size_T shape: ", size_T.shape)
        print("vertex shape: ", vertex.shape)
        print("vertex_T shape: ", vertex_T.shape)

        print(f"t_now: {t_now}")
        print(f"tstep: {tstep}")
        print(f"jorek_model: {jorek_model}")
        print(f"n_elements: {n_elements}")
        print(f"n_nodes: {n_nodes}")
        print(f"n_degrees: {n_degrees}")
        print(f"n_tor: {n_tor}")

        grid = self.mhd.grid_ggd[ti]
        grid.identifier.name = "JOREK grid"
        grid.index = 1
        grid.description = "JOREK grid (Bezier finite elements)"
        grid.time = hdf5_file['t_now'][0]

        # R: mu = 0
        # Z: mu = 1
        # p_k: dof = 0
        # u_k: dof = 1
        # v_k: dof = 2
        # w_k: dof = 3
        # x[mu][dof][k]
        # Four dofs: # p_k, u_k, v_k, w_k
        # grid.space[0].objects_per_dimension[0].object[k].geometry =
        #                                              [p_k_R, u_k_R, v_k_R, w_k_R,
        #                                               p_k_Z, u_k_Z, v_k_Z, w_k_Z]
        grid.space.resize(2)

        # The coordinate types are specified in Data Dictionary (DD)
        # utilities/coordinate_identifier.xml
        grid.space[0].coordinates_type.resize(2)
        grid.space[0].coordinates_type = np.array([4, 3])  # R, Z

        grid.space[0].identifier.name = "Space R-Z"
        grid.space[0].identifier.index = 1  # Fortran notation required
        grid.space[0].identifier.description = """Space R-Z, four degrees of freedom
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
        grid.space[0].objects_per_dimension.resize(3)
        s0opd0 = grid.space[0].objects_per_dimension[0]
        s0opd0.object.resize(hdf5_file['n_nodes'][0])
        for k in range(n_nodes):
            s0opd0.object[k].geometry.resize(n_degrees*2)
            for dof in range(n_degrees):
                s0opd0.object[k].geometry = np.ravel(x_T[k], order='F')

        # dummy object_per_dimension[1]
        # Note: leaving mid AOS (Arrays of Structures) empty might result in
        #       the put command skipping all the next/following AOS altogether
        grid.space[0].objects_per_dimension[1].object.resize(1)

        # vertex (connectivity array for 2D elements)
        s0opd2 = grid.space[0].objects_per_dimension[2]
        s0opd2.object.resize(n_elements)
        for i_element in range(n_elements):
            s0opd2.object[i_element].nodes.resize(vertex_T.shape[1])
            s0opd2.object[i_element].nodes = vertex_T[i_element]

            # Data Dictionary entry for geometry:
            #  Geometry data associated with the object. Its dimension depends
            #  on the type of object, geometry and coordinate considered.
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

    def setBezierValues2IDS(self, hdf5_file_id):
        """Set Bezier grid based values from 'values' matrix (found in HDF5 file).

        ti ... time slice
        """

        hdf5_file = self.getHDF5File(hdf5_file_id=hdf5_file_id)

        values = hdf5_file['values']
        values_T = np.array(values).transpose()
        t_now = hdf5_file['t_now'][0]
        n_var = hdf5_file['n_var'][0]
        n_tor = hdf5_file['n_tor'][0]
        n_elements = hdf5_file['n_elements'][0]
        n_nodes = hdf5_file['n_nodes'][0]
        n_degrees = hdf5_file['n_degrees'][0]
        jorek_model = hdf5_file['jorek_model'][0]

        print("values.shape: ", values.shape)
        print("values_T.shape: ", values_T.shape)
        print(f"n_elements: {n_elements}")
        print(f"n_nodes: {n_nodes}")
        print(f"n_degrees: {n_degrees}")

        ggd = self.mhd.ggd[ti]
        ggd.time = t_now

        if (jorek_model == 303 or jorek_model == 307):
            variables = ["flux", "potential", "current", "vorticity", "density",
                         "temperature", "v_par"]
        else:
            print(f"No support for given model {jorek_model}. Returning.")
            return

        for v in range(n_var):
            if variables[v] == "flux":  # v = 0
                print("Writing quantity array: ", variables[v])
                ggd.psi.resize(n_tor)
                for l in range(n_tor):
                    ggd.psi[l].grid_index = 1
                    ggd.psi[l].grid_subset_index = -l
                    ggd.psi[l].coefficients.resize(n_nodes, n_degrees)
                    ggd.psi[l].coefficients = values[0][:][l]
            elif variables[v] == "potential":  # v = 1
                print("Writing quantity array: ", variables[v])
                ggd.phi_potential.resize(n_tor)
                for l in range(n_tor):
                    ggd.phi_potential[l].grid_index = 1
                    ggd.phi_potential[l].grid_subset_index = -l
                    ggd.phi_potential[l].coefficients.resize(n_nodes, n_degrees)
                    ggd.phi_potential[l].coefficients = values[0][:][l]
            elif variables[v] == "current":  # v = 2
                print("Writing quantity array: ", variables[v])
                ggd.j_tor.resize(n_tor)
                for l in range(n_tor):
                    ggd.j_tor[l].grid_index = 1
                    ggd.j_tor[l].grid_subset_index = -l
                    ggd.j_tor[l].coefficients.resize(n_nodes, n_degrees)
                    ggd.j_tor[l].coefficients = values[0][:][l]
            elif variables[v] == "vorticity":  # v = 3
                print("Writing quantity array: ", variables[v])
                ggd.vorticity.resize(n_tor)
                for l in range(n_tor):
                    ggd.vorticity[l].grid_index = 1
                    ggd.vorticity[l].grid_subset_index = -l
                    ggd.vorticity[l].coefficients.resize(n_nodes, n_degrees)
                    ggd.vorticity[l].coefficients = values[0][:][l]
            elif variables[v] == "density":  # v = 4
                print("Writing quantity array: ", variables[v])
                ggd.mass_density.resize(n_tor)
                for l in range(n_tor):
                    ggd.mass_density[l].grid_index = 1
                    ggd.mass_density[l].grid_subset_index = -l
                    ggd.mass_density[l].coefficients.resize(n_nodes, n_degrees)
                    ggd.mass_density[l].coefficients = values[0][:][l]
            elif variables[v] == "temperature":  # v = 5
                print("Writing quantity array: ", variables[v])
                ggd.electrons.temperature.resize(n_tor)
                for l in range(n_tor):
                    ggd.electrons.temperature[l].grid_index = 1
                    ggd.electrons.temperature[l].grid_subset_index = -l
                    ggd.electrons.temperature[l].coefficients.resize(n_nodes, n_degrees)
                    ggd.electrons.temperature[l].coefficients = values[0][:][l]
            elif variables[v] == "v_par":  # v = 6
                print("Writing quantity array: ", variables[v])
                ggd.velocity_parallel.resize(n_tor)
                for l in range(n_tor):
                    ggd.velocity_parallel[l].grid_index = 1
                    ggd.velocity_parallel[l].grid_subset_index = -l
                    ggd.velocity_parallel[l].coefficients.resize(n_nodes, n_degrees)
                    ggd.velocity_parallel[l].coefficients = values[0][:][l]
            else:
                print(f"WARNING! Unrecognized variable {variables[v]} found. Skipping.")

            print(ggd)

    def write2IDS(self):

        self.log.info("WRITING DATA TO IDS (put)")
        # Write the set data to IDS
        tt = time()
        self.pdd.put(self.mhd)
        print(f"Time required for put() command to finish: {time()-tt:.2f}s")
        # Close the database
        self.pdd.close()


if __name__ == "__main__":

    # Set mandatory arguments
    IDS_parameters = checkArguments()

    # FIRST FILE:
    print("READING RUN DIRECTORY")
    j2IDS = jorekRun2IDS()
    # j2IDS.setDebugMode(1)  # Set debug mode ON
    j2IDS.setRunDirectory()
    j2IDS.setHDF5filesList()

    print("PREPARING IDS")
    shot = IDS_parameters['shot']
    run = IDS_parameters['run']
    username = IDS_parameters['user']
    database = IDS_parameters['database']
    # occurrence = IDS_parameters['occurrence']  # Not yet implemented

    j2IDS.setIDS(shot, run, username, database)
    # j2IDS.mhd.time.resize(len(j2IDS.HDF5pathList))
    # j2IDS.mhd.time = [0]*len(j2IDS.HDF5pathList)
    print("WRITING TO IDS")
    j2IDS.log.info(f"Shot: {shot}")
    j2IDS.log.info(f"Run: {run}")
    j2IDS.log.info(f"Username: {username}")
    j2IDS.log.info(f"Database: {database}")

    t0 = time()

    # set code parameters to IDS
    j2IDS.setCodeContents_unfinished()

    # Loop through the list of HDF5 files
    for ti in range(len(j2IDS.HDF5pathList)):

        hdf5_file = h5py.File(j2IDS.HDF5pathList[ti], 'r')

        print("Slice: ", ti)

        # hdf5_file = getHDF5File(j2IDS.HDF5pathList[ti])

        print("keys: ", hdf5_file.keys())

        # variables = [0, 1, 2, 3, 4, 5, 6]
        # var_names = ["psi", "u", "j", "w", "rho", "T", "v_par"]

        j2IDS.setIDSProperties(hdf5_file_id=ti)

        # set grid_ggd
        # j2IDS.mhd.grid_ggd.resize(len(j2IDS.HDF5pathList))
        t1 = time()
        j2IDS.setBezierGrid2IDS(hdf5_file_id=ti)
        # j2IDS.mhd.time.resize(len(j2IDS.HDF5pathList))
        j2IDS.mhd.time[ti] = hdf5_file['t_now'][0]
        print(f"Time required to set grid_ggd for time slice {ti}: {time()-t1:.2f}s")

        # set ggd
        # j2IDS.mhd.ggd.resize(len(j2IDS.HDF5pathList))
        t2 = time()
        j2IDS.setBezierValues2IDS(hdf5_file_id=ti)
        print(f"Time required to set ggd for time slice {ti}: {time()-t2:.2f}s")

    j2IDS.write2IDS()
    # t1 = time()
    # mhd.put()
    # print(f"Time required for put() command to finish: {time()-t1:.2f}s")

    print(f"Total time: {time()-t0:.2f}s")

    # pdd.close()
