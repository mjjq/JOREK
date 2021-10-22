#!/usr/bin/env python
#   Name : jorekHDF5toIDS.py
#
#   Description :
#       A script which reads JOREK HDF5 output file(s) and writes its
#       contents (grid geometry as Bezier elements, quad connectivity,
#       data fields with Fourer harmonics) to IMAS MHD IDS.
#  Main features:
#      1. There are two grid_ggd spaces:
#         - Space 1 is two-dimensional (R,Z) space of unstructured grid
#           with geometry_2d associated to nodes and cells.
#         - Space 2 is one-dimensional space with $\phi$ and
#           geometry_type.index = 1 (Fourier) with geometry
#           objects (1, 2, 3, ..., number of harmonics)
#     2. Values (Te, n, w, ...), stored under ggd, with the above two
#        spaces form a "structured" (implicitly defined) grid of node
#        values, where explicitly RZ values of first harmonics are saved
#        first, then RZ values of the second harmonics follow up to the
#        last RZ harmonics. Similarly, coefficients on the nodes are saved.
#        This definition follows column major (FORTRAN) notation,
#        meaning that with varying first index (R) the values are close
#        together in the memory and that the last index
#        (phi in Fourier space) defines RZ block of values in memory.
#
#       JOREK HDF5 file contents:
#           https://www.jorek.eu/wiki/doku.php?id=hdf5-tools&s[]=hdf5
#
#   Requirements:
#       - h5py
#       - IMAS
#       - mkdir -p $HOME/public/imasdb/smiter/3/0
#
#   Author :
#       Leon Kos and Miha Radež
#   E-mail :
#       leon.kos@lecad.fs.uni-lj.si
#
# *****************************************************************************

from os import listdir, getenv
from os.path import isfile, join
import numpy as np
import sys
import getpass
import argparse
import h5py
import vtk

import imas
from imas import imasdef
import logging, logging.config
from idsUtilities import basicIDS, writeIDS

prec=np.float32
vtk_prec=vtk.VTK_FLOAT

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Convert JOREK HDF5 file(s) to IMAS (IDSs)",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--shot", type=int, default=1,
                        help="Shot number")
    parser.add_argument("-r", "--run", type=int, default=7,
                        help="Run number")
    parser.add_argument("-u", "--user", type=str, default=getpass.getuser(),
                        help="Location of ~$USER/public/imasdb")
    parser.add_argument("-d", "--database", type=str, default="smiter",
                        help="Database name under public/imasdb/")
    parser.add_argument("-o", "--occurrence", type=int, default=0,
                        help="Occurrence number")
    parser.add_argument("-f", "--backend", type=int, default=imasdef.MDSPLUS_BACKEND,
                        help="Database format: 12=MDSPLUS, 13=HDF5")
    parser.add_argument("hdf5files", metavar='jorek?????.h5', nargs='*',
                        help="JOREK HDF5 file(s)", default=["/tmp/jorek_restart.h5"])

    args = parser.parse_args()

    logging.config.fileConfig(fname='log.conf')
    logger = logging.getLogger('JOREK_HDF5_2_IDS')
    np.set_printoptions(
        formatter={'float': '{: 0.3f}'.format})

    logger.info("PREPARING IDS")
    w_ids = writeIDS(args.shot, args.run, args.user, args.database, args.backend)
    comment = "Written results of multiple JOREK output HDF5 files/timeslices."
    # w_ids.createDatasetDescriptionIDS(user = 'penkod',
    #                                   source = 'JOREK',
    #                                   comment = '')
    w_ids.createGridGGD('mhd', 1)

    # w_ids.data_entry.mhd.time.resize(len(filePathList))
    allTimeValues = np.array([0]*len(args.hdf5files))

    # Resize GGD to the number of timeslices
    w_ids.ids.ggd.resize(len(args.hdf5files))

    logger.info("WRITING TO IDS")

    # Loop through the list of HDF5 files
    for i_slice in range(len(args.hdf5files)):

        logger.info("Slice: %d", i_slice)

        # Read file
        #f.read(args.hdf5files[i_slice], variables=[0, 1, 2, 3, 4, 5, 6])
        with h5py.File(args.hdf5files[0], 'r') as hf:
            n_var        = hf.get('n_var')[0]
            n_period     = hf.get('n_period')[0]
            n_tor        = hf.get('n_tor')[0]
            n_vertex_max = hf.get('n_vertex_max')[0]
            n_elements   = hf.get('n_elements')[0]
            vertex       = np.array(hf.get('vertex'))
            x            = np.array(hf.get('x'))
            size         = np.array(hf.get('size'))
            values       = np.array(hf.get('values'))
            tstep        = hf['tstep'][0]
            t_now        = hf['t_now'][0]
            
        logger.info("Time step: %d", tstep)
        logger.info("Time: %d", t_now)
        logger.info("n_period = %d", n_period)

        # Grid geometry is taken only from the first time slice. All other
        # time slices share the same geometry. No need to write the same grid
        # multiple times as all of them have the same.
        if i_slice == 0:

            # x_coord = f.xyz[:,0]
            # y_coord = f.xyz[:,1]
            logger.info (vertex)
            #logger.info ("* list_vertex: %d", vertex)
            logger.info ("* len(list_vertex[0]):%d", len(vertex[0]))
            
            #logger.info ("* list_vertex[0]: %d", vertex[0])
            logger.info (vertex[0])
            # Remove the vtk cell type ID from the matrix

            #ien for quad
            ien = np.swapaxes(vertex, 1, 0) - 1
            vtk_quad_conn_array = np.insert(ien, 0, 4, axis=1)
            quad_conn_array = ien

            #xyz for quad
            x0 = x[0, 0]
            y0 = x[1, 0]
            xyz = np.zeros((np.shape(x0)[0], 3))
            xyz[:, 0] = x0
            xyz[:, 1] = y0
            

            #val for quad
            val0 = values[:, 0, :, :]
            tor = [1, 1, 0, 1, 0]
            val = np.einsum('ijk,j->ik', val0, tor)
            
            logger.info ("* vtk_quad_conn_array:")
            logger.info(vtk_quad_conn_array)
            logger.info ("* vtk_quad_conn_array.shape:")
            logger.info (vtk_quad_conn_array.shape)
            logger.info("* len(f.xyz): %d", len(xyz))
            logger.info ("* f.xyz:")
            logger.info (xyz)
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
            w_ids.writeMeshToSlice(points_geo=xyz[:, :],
                                   obj_0D_list=[],
                                   obj_1D_list=[],
                                   obj_2D_list=obj_2D_list_f90,
                                   obj_3D_list=[],
                                   n=0,
                                   s=2,
                                   label='JOREK output HDF5 file grid with quantities')

            # w_ids.grid_ggd.array[0].space.resize(2)  # Assure second space is allocated
            toroidal_space = w_ids.grid_ggd[0].space[1]  # TODO multiple time slices
            toroidal_space.coordinates_type.resize(1)  # One dimensional (toroidal) space
            toroidal_space.coordinates_type[0] = 5  # Toroidal angle (phi) in radians
            toroidal_space.geometry_type.index = n_period  # Fourier with periodicity
            toroidal_space.identifier.description = "Toroidal space"

            toroidal_space.objects_per_dimension.resize(1)  # We have only one dimension of
            toroidal_space.objects_per_dimension[0].object.resize(n_tor)  # toroidal harmonics

            for tor in range(n_tor):
                toroidal_space.objects_per_dimension[0].object[tor].geometry = np.array([tor+1])

            gr2d = w_ids.grid_ggd[0].space[0]
            gr2d.geometry_type.index = 0  # Standard geometry (non Fourier)
            x_shape = np.shape(x)[2]
            # coordinate and derivates (s, t, mixed)
            for j in range(x_shape):
                gr2d.objects_per_dimension[0].object[j].geometry_2d = x[:, :, j]

            # size 1, d_{uk}, d_{vk}, d{uv}d{vk} as in Daan Van Vugt thesis
            size_shape = np.shape(size)[2]
            for i in range(size_shape):
                gr2d.objects_per_dimension[2].object[i].geometry_2d = size[:, :, i]

            
        quantity_names_list = ["psi", "u", "j", "w", "rho", "T", "v_par"]
        quantities_array = val
        # print("quantity_names_list: ", quantity_names_list)
        # print("quantities_array: ", quantities_array)

        # Set empty IDS path for quantity tree node
        IDSQuantityPath = None

        # Set time
        w_ids.ids.ggd[i_slice].time = t_now
        # Add to array of all time values of all time slices
        allTimeValues[i_slice] = t_now

        #w_ids.ids.put()
        #w_ids.idsClose()
        #exit(0)

        logger.info("Array of quantity labels: ")
        logger.info(quantity_names_list)
        for i in range(len(quantity_names_list)):
            label = quantity_names_list[i]
            logger.info("Current quantity array: %s", label)
            if label == 'psi':  # Flux / poloidal magnetic flux
                logger.info("Writing quantity array: %s", label)
                # Resize/Allocate
                w_ids.ids.ggd[i_slice].psi.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.ids.ggd[i_slice].psi[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, val[i, :], 1)

                val_tor = w_ids.ids.ggd[i_slice].psi
                psi = values[0, :, :, :]
                psi = np.swapaxes(psi, 0, 1)
                a = np.shape(psi)
                psi1 = np.reshape(psi, (a[0], a[1] * a[2]))
                val_tor[0].coefficients = psi1
            elif label == 'u':  # Potential / electric potential
                logger.info("Writing quantity array: %s", label)
                # Resize/Allocate
                w_ids.ids.ggd[i_slice].phi_potential.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.ids.ggd[i_slice].phi_potential[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, val[i, :], 1)

                val_tor = w_ids.ids.ggd[i_slice].phi_potential
                u = values[1, :, :, :]
                u = np.swapaxes(u, 0, 1)
                a = np.shape(u)
                u1 = np.reshape(u, (a[0], a[1] * a[2]))
                val_tor[0].coefficients = u1
            elif label == 'j':  # Current / toroidal current density
                logger.info("Writing quantity array: %s", label)
                # Resize/Allocate
                w_ids.ids.ggd[i_slice].j_tor.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.ids.ggd[i_slice].j_tor[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, val[i, :], 1)

                val_tor = w_ids.ids.ggd[i_slice].j_tor
                j = values[2, :, :, :]
                j = np.swapaxes(j, 0, 1)
                a = np.shape(j)
                j1 = np.reshape(j, (a[0], a[1] * a[2]))
                val_tor[0].coefficients = j1
            elif label == 'w':  # Vorticity
                logger.info("Writing quantity array: %s", label)
                # Resize/Allocate
                w_ids.ids.ggd[i_slice].vorticity.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.ids.ggd[i_slice].vorticity[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, val[i, :], 1)

                val_tor = w_ids.ids.ggd[i_slice].vorticity
                w = values[3, :, :, :]
                w = np.swapaxes(w, 0, 1)
                a = np.shape(w)
                w1 = np.reshape(w, (a[0], a[1] * a[2]))
                val_tor[0].coefficients = w1
                
            elif label == 'rho':  # Mass Density
                logger.info("Writing quantity array: %s", label)
                # Resize/Allocate
                w_ids.ids.ggd[i_slice].mass_density.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.ids.ggd[i_slice].mass_density[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, val[i,:], 1)
                
                val_tor = w_ids.ids.ggd[i_slice].mass_density
                rho = values[4, :, :, :]
                rho = np.swapaxes(rho, 0, 1)
                a = np.shape(rho)
                rho1 = np.reshape(rho, (a[0], a[1] * a[2]))
                val_tor[0].coefficients = rho1
                pass
            elif label == 'T':  # Temperature / total temperature
                logger.info("Writing quantity array: %s", label)
                # Resize/Allocate
                w_ids.ids.ggd[i_slice].electrons.temperature.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.ids.ggd[i_slice].electrons.temperature[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, val[i, :], 1)

                #adding values and derivates for each Fourier harmonic
                #caution for adding toroidal coordinate we must resize ...temperature to n_tor
                val_tor = w_ids.ids.ggd[i_slice].electrons.temperature
                temp = values[5, :, :, :]
                temp = np.swapaxes(temp, 0, 1)
                a = np.shape(temp)
                temp1 = np.reshape(temp, (a[0], a[1] * a[2]))
                val_tor[0].coefficients = temp1

            elif label == 'v_par':  # V_parallel / parallel velocity
                logger.info("Writing quantity array: %s", label)
                # Resize/Allocate
                w_ids.ids.ggd[i_slice].velocity_parallel.resize(1)
                # Set IDS path
                IDSQuantityPath = w_ids.ids.ggd[i_slice].velocity_parallel[0]
                # Write quantities
                w_ids.ggdWriteQuantityArray(IDSQuantityPath, val[i, :], 1)

                logger.info("v_par max value:")
                logger.info(max(val[i, :]))
                val_tor = w_ids.ids.ggd[i_slice].velocity_parallel
                vp = values[6, :, :, :]
                vp = np.swapaxes(vp, 0, 1)
                a = np.shape(vp)
                vp1 = np.reshape(vp, (a[0], a[1] * a[2]))
                val_tor[0].coefficients = vp1
        w_ids.data_entry.put(w_ids.ids)
        w_ids.idsClose()

    # Set time array
    w_ids.ids.time = allTimeValues
