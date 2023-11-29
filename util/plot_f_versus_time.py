#! /usr/bin/env python3

#
# Purpose: Plot radial (y-axis) distribution of specified variable versus time (x-axis) from postproc files.
#
# Date: 2023-11-28
# Author: Haowei Zhang, IPP Garching
#
import os
import math
import argparse
import numpy as np
import copy
import matplotlib.pyplot as plt
import plot_vtk as pltvtk
import plot_multiple_files as pltmpf
import plot_q_versus_time as pltqvst

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot radial distribution of q (y-axis) versus time (x-axis) from postproc files.")
    parser.add_argument("-fp", "--file_path", default="postproc/",
                        help="Directory contains ???????_sxxxxx.dat.")
    parser.add_argument("-fn", "--file_name", default="exprs_midplane_s",
                        help="Prefix of the files to be plotted.")
    parser.add_argument("-fz", "--fontsize", type=int, default=15,
                        help="Define the fontsize of figure.")
    parser.add_argument('-xlim', nargs=2, type=float, metavar=('xmin', 'xmax'),
                        help='Specify xmin, xmax for the plot.')
    parser.add_argument('-ylim', nargs=2, type=float, metavar=('ymin', 'ymax'),
                        help='Specify ymin, ymax for the plot.')
    parser.add_argument('-clim', nargs=2, type=float, metavar=('cmin', 'cmax'),
                        help='Specify cmin, cmax for the contourf.')
    parser.add_argument("-qs", "--q_surfaces", nargs="+", type=float, default=[1.0],
                        help="Plot the lines where q equals the specified value.")
    parser.add_argument("-sk", "--skiprows", default=1,
                        help="Define how many rows will be skiped when reading the data.")
    parser.add_argument("-yc", "--column_plot", type=int, default=1,
                        help="The column to be plotted.")
    parser.add_argument("-rp", "--pow_for_r", type=float, default=1.0,
                        help="Set power for the radial-axis, e.g., r_n = pow(psi_N, 0.5).")
    parser.add_argument("-xl", "--xlabel", default="", help="Manually define the xlabel.")
    parser.add_argument("-yl", "--ylabel", default="", help="Manually define the ylabel.")
    parser.add_argument("-ti", "--title", default="", help="Manually define the title.")
    parser.add_argument("-cs", "--contour_surfaces", nargs="+", type=float, default=[],
                        help="Plot the lines where q equals the specified value.")
    
    args = parser.parse_args()
    
    font_size = args.fontsize
    directory = args.file_path

    # Get a list of file names in the directory
    file_names = [file for file in os.listdir(directory) if file.startswith(args.file_name) and file.endswith('.dat')]

    # Sort the file names based on the number in the filename
    sorted_file_names = sorted(file_names, key=pltqvst.extract_number)
    # print(sorted_file_names)

    # Initialize 1D and 2D arrays
    numbers_1d_array = []
    data_2d_array = []
    radial_1d_array = []

    # Read the second column of each file and store it in the 2D array
    for i, file_name in enumerate(sorted_file_names):
        file_path = os.path.join(directory, file_name)
        numbers_1d_array.append(pltqvst.extract_number(file_path))
        column_names, data = pltmpf.read_ascii_file(file_path, args.skiprows)
        data_2d_array.append(data[:, args.column_plot])
        if i == 0:
            radial_1d_array = data[:, 0]
            coord_name = column_names[0]
            variable_name = column_names[args.column_plot]

    timefile = 'times.dat'
    if os.path.exists(timefile):
        column_names, data = pltmpf.read_ascii_file(timefile, 1)
        step_all_1d = data[:,0]
        times_all_1d = data[:,1]
        times_1d_array = [0]
        for it in numbers_1d_array:
            for i, step in enumerate(step_all_1d):
                if it == step:
                    times_1d_array.append(times_all_1d[i])
    else:
        print("times.dat required, please firstly execute $ plot_live_data.sh -q time ")
        exit(-1)

    # numbers_1d = np.array(numbers_1d_array)
    times_1d = np.array(times_1d_array)
    radial_1d = np.array(radial_1d_array)
    data_2d = np.array(data_2d_array)
    data_2d_num = np.shape(data_2d)[0] * np.shape(data_2d)[1]
    if len(times_1d) != np.shape(data_2d)[0] or len(radial_1d) != np.shape(data_2d)[1]:
        print('mismatch of arrays, please check!')
        print('shape of times_1d = ', np.shape(times_1d))
        print('shape of radial_1d  = ', np.shape(radial_1d))
        print('shape of data_2d  = ', np.shape(data_2d))
        exit(-1)

    if args.pow_for_r != 1.0:
        radial_1d = radial_1d ** args.pow_for_r
        coord_name = coord_name + '^' + "{:.1f}".format(args.pow_for_r)

    radial_2d, times_2d = np.meshgrid(radial_1d, times_1d)
    data_2d_plt = copy.deepcopy(data_2d)

    if args.xlim:
        xmin, xmax = args.xlim
        data_2d_plt[np.where(times_2d < xmin)] = math.nan
        data_2d_plt[np.where(times_2d > xmax)] = math.nan
    if args.ylim:
        ymin, ymax = args.ylim
        data_2d_plt[np.where(radial_2d < ymin)] = math.nan
        data_2d_plt[np.where(radial_2d > ymax)] = math.nan

    plt.figure(figsize=(10, 6)) 
    n_level = 100
    if args.clim:
        plt.contourf(times_2d, radial_2d, data_2d_plt, n_level, levels = np.linspace(args.clim[0], args.clim[1], n_level), cmap='jet', extend = 'both')
    else:
        plt.contourf(times_2d, radial_2d, data_2d_plt, n_level, cmap='jet')

    plt.colorbar()
    plt.title(r'radial distribution of ' + variable_name + ' versus time' if args.title == "" else args.title, fontsize=font_size)
    plt.xlabel(r'normalized time' if args.xlabel == "" else args.xlabel, fontsize=font_size)        
    plt.ylabel(coord_name if args.ylabel == "" else args.ylabel, fontsize=font_size)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)

    if args.contour_surfaces:
        contour_line = plt.contour(times_2d, radial_2d, data_2d_plt, levels=args.contour_surfaces, colors='k', linestyles = 'dotted')
        plt.clabel(contour_line, inline=True, fontsize=8)
    plt.tight_layout()

    x = np.reshape(times_2d, data_2d_num)
    y = np.reshape(radial_2d, data_2d_num)
    z = np.reshape(data_2d, data_2d_num)

    if args.xlim:
        plt.xlim(xmin, xmax)
    if args.ylim:
        plt.ylim(ymin, ymax)

    plt.connect('motion_notify_event', 
                lambda event: pltvtk.display_values(event, x, y, z))

    plt.show()