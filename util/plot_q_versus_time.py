#! /usr/bin/env python3

#
# Purpose: Plot radial (y-axis) distribution of q versus time (x-axis) from postproc files.
#
# Date: 2023-11-28
# Author: Haowei Zhang, IPP Garching
#
import os
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import plot_vtk as pltvtk
import plot_multiple_files as pltmpf

# Function to extract the number from the filename
def extract_number(filename):
    return int(filename.split('_s')[1].split('.')[0])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot radial distribution of q (y-axis) versus time (x-axis) from postproc files.")
    parser.add_argument("-fp", "--file_path", default="postproc/",
                        help="Directory contains qprofile_sxxxxx.dat.")
    parser.add_argument("-fz", "--fontsize", type=int, default=15,
                        help="Define the fontsize of figure.")
    parser.add_argument("-qc", "--q_cut", type=float, default=1.3,
                        help="Define the cutoff of the q for plotting.")
    parser.add_argument('-xlim', nargs=2, type=float, metavar=('xmin', 'xmax'),
                        help='Specify xmin, xmax for the plot.')
    parser.add_argument('-ylim', nargs=2, type=float, metavar=('ymin', 'ymax'),
                        help='Specify ymin, ymax for the plot.')
    parser.add_argument("-cs", "--contour_surfaces", nargs="+", type=float, default=[1.0],
                        help="Plot the lines where q equals the specified value.")
    parser.add_argument("-xl", "--xlabel", default="", help="Manually define the xlabel.")
    parser.add_argument("-yl", "--ylabel", default="", help="Manually define the ylabel.")
    parser.add_argument("-ti", "--title", default="", help="Manually define the title.")
    parser.add_argument("-xm", "--multiplier_for_X", type=float, default=1.0,
                        help="Multiple the X-data with the specified value.")
    parser.add_argument('-clim', nargs=2, type=float, metavar=('cmin', 'cmax'),
                        help='Specify cmin, cmax for the contourf.')
    
    args = parser.parse_args()
    
    font_size = args.fontsize
    q_cut = args.q_cut
    directory = args.file_path

    # Get a list of file names in the directory
    file_names = [file for file in os.listdir(directory) if file.startswith('qprofile_s') and file.endswith('.dat')]

    # Sort the file names based on the number in the filename
    sorted_file_names = sorted(file_names, key=extract_number)
    # print(sorted_file_names)

    # Initialize 1D and 2D arrays
    numbers_1d_array = []
    data_2d_array = []
    psin_1d_array = []

    # Read the second column of each file and store it in the 2D array
    for i, file_name in enumerate(sorted_file_names):
        file_path = os.path.join(directory, file_name)
        numbers_1d_array.append(extract_number(file_path))
        column_names, data = pltmpf.read_ascii_file(file_path, 0)
        data_2d_array.append(data[:,1])
        if i == 0:
            psin_1d_array = data[:,0]

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
    psin_1d = np.array(psin_1d_array)
    data_2d = np.array(data_2d_array)
    data_2d_num = np.shape(data_2d)[0] * np.shape(data_2d)[1]
    if len(times_1d) != np.shape(data_2d)[0] or len(psin_1d) != np.shape(data_2d)[1]:
        print('mismatch of arrays, please check!')
        print('shape of times_1d = ', np.shape(times_1d))
        print('shape of psin_1d  = ', np.shape(psin_1d))
        print('shape of data_2d  = ', np.shape(data_2d))
        exit(-1)

    times_1d = times_1d * args.multiplier_for_X

    psin_2d, times_2d = np.meshgrid(psin_1d, times_1d)
    data_2d_plt = data_2d
    data_2d_plt[np.where(abs(data_2d) > q_cut)] = q_cut # math.nan

    plt.figure(figsize=(10, 5)) 
    if args.clim:
        plt.contourf(times_2d, np.sqrt(psin_2d), abs(data_2d_plt), 50, levels = np.linspace(args.clim[0], args.clim[1], 100), cmap='viridis_r')#, extend = 'both')
    else:
        plt.contourf(times_2d, np.sqrt(psin_2d), abs(data_2d_plt), 50, cmap='viridis_r')

    plt.colorbar()
    plt.title(r'radial distribution of $q$ versus time' if args.title == "" else args.title, fontsize=font_size)
    plt.xlabel(r'normalized time' if args.xlabel == "" else args.xlabel, fontsize=font_size)
    plt.ylabel(r'$\sqrt{\psi_n}$' if args.ylabel == "" else args.ylabel, fontsize=font_size)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)

    contour_line = plt.contour(times_2d, np.sqrt(psin_2d), abs(data_2d_plt), levels=args.contour_surfaces, colors='k', linestyles = 'dotted')
    plt.clabel(contour_line, inline=True, fontsize=8)
    plt.tight_layout()

    x = np.reshape(times_2d, data_2d_num)
    y = np.reshape(np.sqrt(psin_2d), data_2d_num)
    z = np.reshape(abs(data_2d), data_2d_num)

    if args.xlim:
        xmin, xmax = args.xlim
        plt.xlim(xmin, xmax)
    if args.ylim:
        ymin, ymax = args.ylim
        plt.ylim(ymin, ymax)

    plt.connect('motion_notify_event', 
                lambda event: pltvtk.display_values(event, x, y, z))

    plt.show()