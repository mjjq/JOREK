#! /usr/bin/env python3

#
# Purpose: Plot 2D distribution of VxB_tor_n=0 (dynamo term) from .vtk files.
# VXB_tor = <uR * BZ - uZ * BR>, <...> indicates the n = 0 component
# Date: 2023-11-07
# Author: Haowei Zhang, IPP Garching
#

import vtk
import numpy as np
import matplotlib.pyplot as plt
import argparse
import copy
import plot_vtk as pltvtk
# import matplotlib.tri as mtri
# import scipy.interpolate as sci_intp
# import mpl_toolkits.mplot3d as plt3d
# import vtkmodules.util.numpy_support as vtk_nump

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot 2D distribution of VxB_tor from .vtk files.")
    parser.add_argument("folder_name",
                        help="The general name of folders contain the target vtk files, e.g., vtk_no0_iplane.")
    parser.add_argument("-np", "--nplane", default=8, 
                        help="Value of nplane.")
    parser.add_argument('-ns', "--nstep", default=4000,
                        help="Value of nstep.")
    parser.add_argument('-Rg', "--R_geo", default=1.65,
                        help="Value of R_geo, i.e., axis of the grids in JORKE code.")
    parser.add_argument('-Zg', "--Z_geo", default=0.0,
                        help="Value of Z_geo, i.e., axis of the grids in JORKE code.")
    parser.add_argument('-de', "--delta", default=0.005,
                        help="Value of delta, i.e., singular grids with distance away from (R_geo, Z_geo) smaller than delta will be set to zero.")
    parser.add_argument('-me', "--method", choices=["it", "sc"], default="it: interpolate",
                        help="Use the specified method to plot the data, including 1. it: interpolate; 2. sc: scatter.")

    args = parser.parse_args()

    reader = vtk.vtkDataSetReader()
    for itor in range(args.nplane):
        file_path = args.folder_name + str(int(itor + 1)) + '/jorek.' + str(int(args.nstep)).zfill(5) + '.vtk'
        
        data_list, coordniates, time = pltvtk.read_data(file_path, reader, ['UR', 'UZ', 'B_R', 'B_Z'], itor)
        UR = data_list[:, 0]
        UZ = data_list[:, 1]
        BR = data_list[:, 2]
        BZ = data_list[:, 3]

        if itor == 0:
            distance = np.sqrt((coordniates[:, 0] - args.R_geo) ** 2 + (coordniates[:, 1] - args.Z_geo) ** 2)
            # indices = np.where((coordniates[:, 0] > (args.R_geo - args.delta)) & (coordniates[:, 0] < (args.R_geo + args.delta)) 
            #                    & (coordniates[:, 1] > (args.Z_geo-args.delta)) & (coordniates[:, 1] < (args.Z_geo + args.delta)))
            indices = np.where(distance < args.delta)
            # print(np.shape(indices))
        BR[indices] = 0.0
        BZ[indices] = 0.0

        VXB_tor = (UR * BZ - UZ * BR) / args.nplane
        if itor == 0:
            VXB_tor_n0 = copy.deepcopy(VXB_tor)
        else:
            VXB_tor_n0 += VXB_tor # average to get the n = 0 components

        # pltvtk.plot2D_interpolated(
        #         coordniates[:, 0], coordniates[:, 1], VXB_tor, 'VxB_tor_' + str(int(itor + 1)), args.folder_name, time)
    if args.method == "sc":
        pltvtk.plot2D_scatter(
                coordniates[:, 0], coordniates[:, 1], VXB_tor_n0, r'$(V \times B)_{\varphi,\ n = 0}$', args.folder_name, time)
    else:
        pltvtk.plot2D_interpolated(
                coordniates[:, 0], coordniates[:, 1], VXB_tor_n0, r'$(V \times B)_{\varphi,\ n = 0}$', args.folder_name, time)

    plt.connect('motion_notify_event',
                lambda event: pltvtk.display_values(event, coordniates[:, 0], coordniates[:, 1], VXB_tor_n0))
    plt.show()
