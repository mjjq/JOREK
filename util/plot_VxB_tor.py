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
    parser.add_argument("-np", "--nplane", type=int, default=8, 
                        help="Value of nplane.")
    parser.add_argument('-ns', "--nstep", type=int, default=4000,
                        help="Value of nstep.")
    parser.add_argument('-Rg', "--R_geo", type=float, default=1.65,
                        help="Value of R_geo, i.e., axis of the grids in JORKE code.")
    parser.add_argument('-Zg', "--Z_geo", type=float, default=0.0,
                        help="Value of Z_geo, i.e., axis of the grids in JORKE code.")
    parser.add_argument('-de', "--delta", type=float, default=0.005,
                        help="Value of delta, i.e., singular grids with distance away from (R_geo, Z_geo) smaller than delta will be set to zero.")
    parser.add_argument('-me', "--method", choices=["it", "sc"], default="it: interpolate",
                        help="Use the specified method to plot the data, including 1. it: interpolate; 2. sc: scatter.")
    parser.add_argument("-fz", "--fontsize", type=int, default=15,
                        help="Define the fontsize of figure.")
    parser.add_argument('-xylim', nargs=4, type=float, metavar=('xmin', 'xmax', 'ymin', 'ymax'),
                        help='Specify xmin, xmax, ymin, ymax for the plot.')
    parser.add_argument("-qs", "--q_surfaces", nargs="+", type=float, default=0.0,
                        help="Plot the surface(s) for specified safety factor q (please note the pm sign).")
    parser.add_argument("-tm", "--multiplier_for_time", type=float, default=1.0,
                        help="Multiple the time-data with the specified value.")
    parser.add_argument('-clim', nargs=2, type=float, metavar=('cmin', 'cmax'),
                        help='Specify cmin, cmax for the contourf.')
    args = parser.parse_args()

    legend_name = r'$(\delta V \times \delta B)_{\varphi,\ n = 0}$'
    folder_name = args.folder_name

    # plot VxB from n = 0 components of (V, B), by subtract n /= 0 components from the total components.
    # e.g., plot_VxB_tor.py -np 32 -ns 40000 vtk_no0_iplane  -xylim 1.20 2.00 -0.20 0.30 -qs -1.1 -tm $time2si  -clim -0.00000018 0.00000018
    plot_eq_component = False
    more_outfile_digits = True

    zfill_num = 5
    if more_outfile_digits:
        zfill_num = 6

    reader = vtk.vtkDataSetReader()
    for itor in range(args.nplane):
        file_path = args.folder_name + str(int(itor + 1)) + '/jorek.' + str(int(args.nstep)).zfill(zfill_num) + '.vtk'
        
        data_list, coordniates, time = pltvtk.read_data(file_path, reader, ['UR', 'UZ', 'B_R', 'B_Z'], itor)
        time = time * args.multiplier_for_time
        UR = data_list[:, 0]
        UZ = data_list[:, 1]
        BR = data_list[:, 2]
        BZ = data_list[:, 3]

        if plot_eq_component:
            legend_name = r'$(V_0 \times B_0)_{\varphi}$'
            folder_name = 'vtk_iplane-' + args.folder_name
            file_path_all = 'vtk_iplane' + str(int(itor + 1)) + '/jorek.' + str(int(args.nstep)).zfill(zfill_num) + '.vtk'
            data_list_all, coordniates_all, time_all = pltvtk.read_data(file_path_all, reader, ['UR', 'UZ', 'B_R', 'B_Z'], itor)
            UR = data_list_all[:, 0] - data_list[:, 0]
            UZ = data_list_all[:, 1] - data_list[:, 1]
            BR = data_list_all[:, 2] - data_list[:, 2]
            BZ = data_list_all[:, 3] - data_list[:, 3]

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
        #         coordniates[:, 0], coordniates[:, 1], VXB_tor, 'VxB_tor_' + str(int(itor + 1)), folder_name, time)
    si_unit = False if args.multiplier_for_time == 1.0 else True
    if args.method == "sc":
        pltvtk.plot2D_scatter(
                coordniates[:, 0], coordniates[:, 1], VXB_tor_n0, legend_name, folder_name, time, args.fontsize, si_unit, args.clim)
    else:
        pltvtk.plot2D_interpolated(
                coordniates[:, 0], coordniates[:, 1], VXB_tor_n0, legend_name, folder_name, time, args.fontsize, si_unit, args.clim)
        
    if args.method != "sf" and not isinstance(args.q_surfaces, float):
        for q_surface in args.q_surfaces:
            pltvtk.read_and_plot_q(q_surface, file_path)

    # plt.connect('motion_notify_event',
    #             lambda event: pltvtk.display_values(event, coordniates[:, 0], coordniates[:, 1], VXB_tor_n0))
    if args.xylim:
        xmin, xmax, ymin, ymax = args.xylim
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
    plt.show()
