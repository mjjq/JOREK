#! /usr/bin/env python3

#
# Purpose: Perform mode structures with given m, n, step, and variable name
#
# Date: 2023-04-16
# Author: Haowei Zhang, IPP Garching
#
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import plot_multiple_files as plmf
import plot_poincare_all as plpoinc
import re
from scipy import interpolate
from tqdm import tqdm
import random

def generate_unique_colors(length):
  """
  Generates a list of 'length' unique RGB color codes.

  Args:
      length: The desired number of unique colors.

  Returns:
      A list of length 'length' containing unique RGB color tuples.
  """
  # Set the random seed, to generate reproducible random
  random.seed(20)
  colors = set()
  while len(colors) < length:
    # Generate random RGB values
    r = random.random()
    g = random.random()
    b = random.random()
    # Add to set if unique, otherwise regenerate
    color = (r, g, b)
    if color not in colors:
      colors.add(color)
  return list(colors)

def extract_number(filename):
    # Extract the number string at the end
    match = re.search(r'\d+$', filename)

    if match:
        number_string = match.group()
        return int(number_string)
    else:
        print("No number string found at the end.")
        return -1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot lines from multiple ASCII files.")
    parser.add_argument("folder_path", default="four_results/", help="Path to the ASCII files.")
    parser.add_argument("-ti", "--title", default="", help="Manually define the title.")
    parser.add_argument("-xl", "--xlabel", default="$\sqrt{\psi_n}$", help="Manually define the xlabel.")
    parser.add_argument("-yl", "--ylabel", default="A.U.", help="Manually define the ylabel.")
    parser.add_argument("-dm", "--data_marker", default=".", help="Manually define the marker for the data.")
    parser.add_argument("-va", "--variable_list", nargs="+", default={
                        "rho"}, help="Plot the specified variable(s), such as [A3, AR, AZ, UR, UZ, Uphi, rho, T, B_R, B_Z, B_phi, etc.].")
    parser.add_argument("-st", "--step_list", nargs="+", default={
                        "00000"}, help="Plot the specified timesteps.")
    parser.add_argument("-ml", "--list_of_m", type=int, nargs="+", default=1,
                        help="List of m to be plotted (default is 1 for plotting m = 1)")
    parser.add_argument("-nl", "--list_of_n", type=int, nargs="+", default=1,
                        help="List of n to be plotted (default is 1 for plotting m = 1)")
    parser.add_argument("-me", "--method", choices=["am", "cm", "ph", "al"], default="al",
                        help="Use the specified method to plot the data, including 1. am: amplitude; 2. cm: complex amplitude; 3. ph: phase; 4. al: all (including 1 ~ 3)")
    parser.add_argument("-ym", "--multiplier_for_Y", type=float, default=1.0,
                        help="Multiple the Y-data with the specified value.")
    parser.add_argument("-fz", "--fontsize", default=15,
                        help="Define the fontsize of figure.")
    parser.add_argument("-lw", "--linewidth", default=2,
                        help="Define the linewidth of figure.")
    parser.add_argument('-xlim', nargs=2, type=float, metavar=('xmin', 'xmax'),
                        help='Specify xmin, xmax for the plot.')
    parser.add_argument('-ylim', nargs=2, type=float, metavar=('ymin', 'ymax'),
                        help='Specify ymin, ymax for the plot.')
    parser.add_argument("-tm", "--multiplier_for_time", type=float, default=1.0,
                        help="Multiple the time-data with the specified value.")
    parser.add_argument("-ex", "--extract_time_series", type=bool, default=False,
                        help="Extract the time series.")
    parser.add_argument("-er", "--radial_location_for_extract_time_series", type=float, default=0.5,
                        help="Radial location for extracting the time series.")
    args = parser.parse_args()

    # Define the number of colors you need
    num_colors = len(args.list_of_n) * len(args.list_of_m)
    # Get a list of predefined color names from Matplotlib
    # colors = plt.cm.get_cmap('Paired', num_colors)
    colors = generate_unique_colors(num_colors)

    if(args.extract_time_series):
        for variable in args.variable_list:
            for intor in args.list_of_n:
                for impol in args.list_of_m:
                    filename = variable + '_n' + "{:03d}".format(intor) + '_m' + "{:03d}".format(impol) + '_step'
                    # print("Looking for " + filename + " in " + args.folder_path)
                    # Get a list of file names in the directory
                    file_names = [file for file in os.listdir(args.folder_path) if file.startswith(filename)]
                    # Sort the file names based on the number in the filename
                    sorted_file_names = sorted(file_names, key=extract_number)

                    # Initialize 1D and 2D arrays
                    numbers_1d_array = []
                    data_2d_array = []
                    phase_2d_array = []
                    psin_1d_array = []

                    # Read the second column of each file and store it in the 2D array
                    for i, file_name in tqdm(enumerate(sorted_file_names), total = len(sorted_file_names), desc="Loading files"):
                        file_path = os.path.join(args.folder_path, file_name)
                        numbers_1d_array.append(extract_number(file_path))
                        column_names, data = plmf.read_ascii_file(file_path, 1)
                        data_2d_array.append(data[:,2] + data[:,3] * 1j) # read the complex value
                        phase_2d_array.append(data[:,4]) # read the phase
                        if i == 0:
                            psin_1d_array = data[:,0]                    

                    timefile = 'times.dat'
                    if os.path.exists(timefile):
                        column_names, data = plmf.read_ascii_file(timefile, 1)
                        new_row = np.zeros((1, data.shape[1]))
                        # Inserting the new row at the beginning
                        data = np.vstack((new_row, data))
                        step_all_1d = data[:,0]
                        times_all_1d = data[:,1]
                        times_1d_array = []
                        for it in tqdm(numbers_1d_array, desc="Assembling data"):
                            i = np.where(step_all_1d == it)[0]
                            if len(i) > 1:
                                print("Find multiple for step = ", it, " in " + timefile)
                            elif len(i) == 1:
                                times_1d_array.append(times_all_1d[i[0]])
                            else:
                                print("Do not find step = ", it, " in " + timefile)
                    else:
                        print("times.dat required, please firstly execute $ plot_live_data.sh -q time ")
                        exit(-1)

                    # numbers_1d = np.array(numbers_1d_array)
                    times_1d = np.array(times_1d_array)
                    psin_1d = np.array(psin_1d_array)
                    data_2d = np.array(data_2d_array)
                    phase_2d = np.array(phase_2d_array)
                    if len(times_1d) != np.shape(data_2d)[0] or len(psin_1d) != np.shape(data_2d)[1]:
                        print('mismatch of arrays, please check!')
                        print('shape of times_1d = ', np.shape(times_1d))
                        print('shape of psin_1d  = ', np.shape(psin_1d))
                        print('shape of data_2d  = ', np.shape(data_2d))
                        exit(-1)

                    times_1d = times_1d * args.multiplier_for_time

                    psin_2d, times_2d = np.meshgrid(psin_1d, times_1d)
                    r2d = np.sqrt(psin_2d)

                    id_target = np.abs(np.sqrt(psin_1d) - args.radial_location_for_extract_time_series).argmin()
                    data_slice_plt = data_2d[:,id_target]
                    phase_slice_plt = phase_2d[:,id_target]
                    # phase_slice_plt = np.angle(data_slice_plt)
                    filename = variable + '_n' + "{:03d}".format(intor) + '_m' + "{:03d}".format(impol)
                    file_save = args.folder_path + '/1D_' + filename + '_vs_time_at_rhop_eq_' + '{:.2f}'.format(np.sqrt(psin_1d[id_target])) + '.dat'
                    data_save=np.column_stack((times_1d, abs(data_slice_plt), data_slice_plt.real, data_slice_plt.imag, phase_slice_plt))
                    np.savetxt(file_save, data_save, delimiter=' ', comments='', header = 'time       amplitude       real       imag       phase')
                    print(filename + " vs time is saved to " + file_save)
                        
                    plt.figure(figsize=(10, 5))
                    plt.plot(times_1d, data_slice_plt.real, linestyle = '-' , linewidth = 2, label = 'real')
                    plt.plot(times_1d, data_slice_plt.imag, linestyle = '--', linewidth = 2, label = 'imag')
                    # plt.plot(times_1d, phase_slice_plt, linestyle = '-' , linewidth = 2, label = 'phase')
                    plt.xlabel(r'time', fontsize=args.fontsize)
                    plt.ylabel(filename + ' value', fontsize=args.fontsize)
                    plt.title(filename + r' vs time @ $\rho_p$ = ' + '{:.2f}'.format(np.sqrt(psin_1d[id_target])), fontsize=args.fontsize)
                    plt.xticks(fontsize=args.fontsize); plt.yticks(fontsize=args.fontsize)
                    plt.xlim(min(times_1d), max(times_1d))
                    plt.grid(); plt.legend(); plt.tight_layout(); plt.show()


    ivar_st = -1
    for variable in args.variable_list:
        for ist in args.step_list:
            ivar_st = ivar_st + 1
            for intor in args.list_of_n:
                inm = -1
                for impol in args.list_of_m:
                    inm = inm + 1
                    filename =  args.folder_path + "/" + variable + '_n' + "{:03d}".format(intor) + '_m' + "{:03d}".format(impol) + '_step' + ist
                    if os.path.exists(filename):
                        print("Plotting " + filename)
                    else:
                        print(filename + " does not exist.")
                        continue
                    
                    column_names, data = plmf.read_ascii_file(filename, 1)
                    x = np.sqrt(data[:, 0])  # Assuming the first column is the x-axis data
                    labelmn = 'm/n = ' + str(impol) + '/' + str(intor)
                    colormn = colors[inm]
                    if args.method == "am" or args.method == "al":
                        plt.figure(3*ivar_st + 1)
                        y = data[:, 1] * args.multiplier_for_Y
                        plt.plot(x, y, color = colormn, label= labelmn, linewidth=args.linewidth)

                    if args.method == "cm" or args.method == "al":
                        plt.figure(3*ivar_st + 2)
                        y = ( data[:, 2] + data[:, 3] * 1j ) * args.multiplier_for_Y
                        plt.plot(x, y.real, color = colormn, linewidth=args.linewidth, linestyle='-', label='Re(-), Im(--): ' + labelmn)
                        plt.plot(x, y.imag, color = colormn, linewidth=args.linewidth, linestyle='--')

                    if args.method == "ph" or args.method == "al":
                        plt.figure(3*ivar_st + 3)
                        y = data[:, 4]
                        plt.plot(x, y, color = colormn, label=labelmn, linewidth=args.linewidth)

                time = plpoinc.get_time(int(ist.lstrip('0')))
                if args.multiplier_for_time != 1.0:
                    time = time * args.multiplier_for_time
                if args.title:
                    title = args.title
                else:
                    title = "time = " + "{:.3f}".format(time) + ", mode structures: "
                if args.method == "am" or args.method == "al":
                    plt.figure(3*ivar_st + 1)
                    legend = plt.legend()
                    legend.set_draggable(True)
                    for label in legend.get_texts():
                        label.set_fontsize(10)
                    plt.grid(True)
                    plt.xlabel(args.xlabel, fontsize=args.fontsize)
                    plt.ylabel(args.ylabel, fontsize=args.fontsize)
                    plt.title(title + " amplitude of " + variable, fontsize=args.fontsize)
                    if args.xlim:
                        xmin, xmax = args.xlim
                        plt.xlim(xmin, xmax)
                    else:
                        plt.xlim(0, 1)
                    if args.ylim:
                        ymin, ymax = args.ylim
                        plt.ylim(ymin, ymax)
                    plt.tight_layout()

                if args.method == "cm" or args.method == "al":
                    plt.figure(3*ivar_st + 2)
                    legend = plt.legend()
                    legend.set_draggable(True)
                    for label in legend.get_texts():
                        label.set_fontsize(10)
                    plt.grid(True)
                    plt.xlabel(args.xlabel, fontsize=args.fontsize)
                    plt.ylabel(args.ylabel, fontsize=args.fontsize)
                    plt.title(title + " real/imag of " + variable, fontsize=args.fontsize)
                    if args.xlim:
                        xmin, xmax = args.xlim
                        plt.xlim(xmin, xmax)
                    else:
                        plt.xlim(0, 1)
                    if args.ylim:
                        ymin, ymax = args.ylim
                        plt.ylim(ymin, ymax)
                    plt.tight_layout()
            
                if args.method == "ph" or args.method == "al":
                    plt.figure(3*ivar_st + 3)
                    legend = plt.legend()
                    legend.set_draggable(True)
                    for label in legend.get_texts():
                        label.set_fontsize(10)
                    plt.grid(True)
                    plt.xlabel(args.xlabel, fontsize=args.fontsize)
                    plt.ylabel(args.ylabel, fontsize=args.fontsize)
                    plt.title(title + " phase of " + variable, fontsize=args.fontsize)
                    if args.xlim:
                        xmin, xmax = args.xlim
                        plt.xlim(xmin, xmax)
                    else:
                        plt.xlim(0, 1)
                    plt.ylim(-np.pi, np.pi)
                    plt.tight_layout()

                plt.show()





