#! /usr/bin/env python3

#
# Purpose: Plot 1D lines from multiple ASCII files.
#
# Date: 2023-09-18
# Author: Haowei Zhang, IPP Garching
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse

def set_fontsize(fig, fontweight='bold', fontsize=12, title_size=None, label_size=None, tick_size=None, legend_size=None, exponent_size=None):
    """
    Set the font size for the title, axis labels, tick labels, and legend of a Matplotlib figure.

    Parameters:
    - fig: Matplotlib figure object (e.g., plt.gcf()).
    - fontweight: Font weight for the text (default is 'bold').
    - fontsize: Default font size to apply if specific sizes are not provided.
    - title_size: Font size for the title (default is `fontsize`).
    - label_size: Font size for the axis labels (default is `fontsize`).
    - tick_size: Font size for the tick labels (default is `fontsize`).
    - legend_size: Font size for the legend (default is `fontsize`).
    - exponent_size: Font size for the exponent in scientific notation (default is `fontsize`).

    Returns:
    - None
    """

    # Use the default fontsize if specific sizes are not provided
    title_size = title_size if title_size is not None else fontsize
    label_size = label_size if label_size is not None else fontsize
    tick_size = tick_size if tick_size is not None else fontsize
    legend_size = legend_size if legend_size is not None else fontsize
    exponent_size = exponent_size if exponent_size is not None else fontsize

    # Iterate over all axes in the figure
    for ax in fig.get_axes():
        # Set title fontsize
        ax.title.set_fontsize(title_size)
        ax.title.set_fontweight(fontweight)

        # Set axis labels fontsize
        ax.xaxis.label.set_fontsize(label_size)
        ax.yaxis.label.set_fontsize(label_size)
        ax.xaxis.label.set_fontweight(fontweight)
        ax.yaxis.label.set_fontweight(fontweight)

        # Set tick labels fontsize
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(tick_size)
            label.set_fontweight(fontweight)

        legend = plt.legend(prop={'weight': fontweight, 'size' : legend_size})
        legend.set_draggable(True)

        # Set exponent fontsize (scientific notation)
        ax.xaxis.get_offset_text().set_fontsize(exponent_size)
        ax.yaxis.get_offset_text().set_fontsize(exponent_size)
        ax.xaxis.get_offset_text().set_fontweight(fontweight)
        ax.yaxis.get_offset_text().set_fontweight(fontweight)


def read_ascii_file(file_path, skiprows, droprows=-1):
    """
    Read data from an ASCII file.

    Args:
        file_path (str): The path to the ASCII file.

    Returns:
        tuple: A tuple containing the column names (as a list) and the data (as a NumPy array).
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract column names from the first line and remove leading/trailing whitespace
    column_names = lines[0].strip().split()
    if skiprows >= 1:
        column_names = lines[skiprows - 1].strip().split()

    # Read data from the remaining lines, converting them to a NumPy array
    data = np.genfromtxt(lines[skiprows:droprows], dtype=float)

    return column_names, data


def plot_ascii_files(file_paths, xscaleType, yscaleType, title, xlabel, ylabel, skiprows, droprows, column_plot_for_X, columns_plot_for_Y, 
                     shift_first_file, color_first_file, pow_for_X, pow_for_Y, multiplier_for_X, multiplier_for_Y, font_size, line_width, xlim, ylim, data_marker, mylegends="", noline = False):
    """
    Plot lines from multiple ASCII files in the same plot.

    Args:
        file_paths (list): List of file paths to the ASCII files.
    """
    # Define line styles (you can customize this list with different line styles)
    linestyles = ['-', '--', '-.', ':']
    # Set the default linewidth for all lines
    mpl.rcParams['lines.linewidth'] = line_width

    plt.figure(figsize=(10, 6))  # Set the figure size

    ilegend = 0
    for i, file_path in enumerate(file_paths):
        column_names, data = read_ascii_file(file_path, skiprows, droprows)
        if column_plot_for_X >= data.shape[1]:
            print("Error: Value of -xc is out of range!!!")
            exit(-1)
        x = data[:, column_plot_for_X]  # Assuming the first column is the x-axis data
        y = data[:, 1:]  # Assuming the remaining columns are y-axis data

        if pow_for_X != 1.0:
            x = pow(x, pow_for_X)
        
        if pow_for_Y != 1.0:
            y = pow(y, pow_for_Y)

        if multiplier_for_X != 1.0:
            x = x * multiplier_for_X
        
        if multiplier_for_Y != 1.0:
            y = y * multiplier_for_Y

        if i == 0:
            x = x + float(shift_first_file)

        linestyle = linestyles[i % len(linestyles)]  # Cycle through linestyles
        if noline:
            linestyle = 'none'

        if columns_plot_for_Y == -1:
            columns_plot_for_Y_mdf = range(y.shape[1])  # Plot all columns
        elif i == 0:
            columns_plot_for_Y_mdf = [j - 1 for j in columns_plot_for_Y]

        for j in columns_plot_for_Y_mdf:
            if j < y.shape[1]:
                if mylegends == "":
                    mylegend = f'{column_names[j + 1]} ({file_path})'
                else:
                    mylegend = mylegends[ilegend]
                    ilegend += 1
                if (i == 0) and (color_first_file != ""):
                    plt.plot(x, y[:, j], data_marker, label=mylegend, linestyle=linestyle, color=color_first_file, markersize=1)
                else:
                    plt.plot(x, y[:, j], data_marker, label=mylegend, linestyle=linestyle, markersize=1)
            else:
                print("Error: Value of -yc is out of range!!!")
                print("file_path: ", file_path)
                exit(-2)

    plt.xlabel(column_names[column_plot_for_X] if xlabel == "" else xlabel, fontsize=font_size)  # Set the x-axis label based on the first file
    plt.ylabel('Y-axis Label' if ylabel == "" else ylabel, fontsize=font_size)  # Set the y-axis label (customize as needed)
    legend = plt.legend()  # Show legends for each line
    legend.set_draggable(True)
    for label in legend.get_texts():
        label.set_fontsize(font_size)
    plt.title(title, fontsize=font_size)  # Set the plot title
    plt.grid(True)  # Show grid lines
    plt.xscale(xscaleType)  # Set linear or log for the x scale
    plt.yscale(yscaleType)  # Set linear or log for the y scale

    # Access the y-axis OffsetText object
    plt.gca().get_yaxis().get_offset_text().set_size(font_size)

    # Axis tick label font size
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)

    if xlim:
        xmin, xmax = xlim
        plt.xlim(xmin, xmax)
    if ylim:
        ymin, ymax = ylim
        plt.ylim(ymin, ymax)

    plt.show()


if __name__ == "__main__":
    # Provide a list of file paths as command-line arguments
    # Example usage: python script.py file1 file2 ... -xs log (or linear) -ys log (or linear) -ti "title name" -xl "xlabel name" -yl "ylabel name"
    parser = argparse.ArgumentParser(description="Plot lines from multiple ASCII files.")
    parser.add_argument("-xs", "--xscale", choices=["linear", "log"], default="linear",
                        help="X-scale (linear or log). Default is linear.")
    parser.add_argument("-ys", "--yscale", choices=["linear", "log"], default="linear",
                        help="Y-scale (linear or log). Default is linear.")
    parser.add_argument("-ti", "--title", default="Multiple ASCII Files Plot", help="Manually define the title.")
    parser.add_argument("-xl", "--xlabel", default="", help="Manually define the xlabel.")
    parser.add_argument("-yl", "--ylabel", default="", help="Manually define the ylabel.")
    parser.add_argument("-dm", "--data_marker", default="", help="Manually define the marker for the data.")
    parser.add_argument("-yc", "--columns_plot_for_Y", type=int, nargs="+", default=-1,
                        help="Columns to be plotted (default is -1 for plotting all) as Y values.")
    parser.add_argument("-xc", "--column_plot_for_X", type=int, default=0,
                        help="Columns to be plotted (default is 0) as X value.")
    parser.add_argument("-sk", "--skiprows", default=1,type=int,
                        help="Define how many rows will be skiped when reading the data.")
    parser.add_argument("-dr", "--droprows", default=-1,type=int,
                        help="Define since which row, data will be droped when reading the data.")
    parser.add_argument("-sf", "--shift_first_file", default=0,
                        help="Shift the data for the first file horizontally (in X direction).")
    parser.add_argument("-co", "--color_first_file", default="",
                        help="Specify the color for the first file.")
    parser.add_argument("-xp", "--pow_for_X", type=float, default=1.0,
                        help="Set power for the X-axis, e.g., r_n = pow(psi_N, 0.5).")
    parser.add_argument("-yp", "--pow_for_Y", type=float, default=1.0,
                        help="Set power for the Y-axis, e.g., f = pow(f, 0.5).")
    parser.add_argument("-xm", "--multiplier_for_X", type=float, default=1.0,
                        help="Multiple the X-data with the specified value.")
    parser.add_argument("-ym", "--multiplier_for_Y", type=float, default=1.0,
                        help="Multiple the Y-data with the specified value.")
    parser.add_argument("file_paths", nargs="+", help="List of file paths to the ASCII files.")
    parser.add_argument("-fz", "--fontsize", default=15,
                        help="Define the fontsize of figure.")
    parser.add_argument("-lw", "--linewidth", default=2,
                        help="Define the linewidth of figure.")
    parser.add_argument('-xlim', nargs=2, type=float, metavar=('xmin', 'xmax'),
                        help='Specify xmin, xmax for the plot.')
    parser.add_argument('-ylim', nargs=2, type=float, metavar=('ymin', 'ymax'),
                        help='Specify ymin, ymax for the plot.')
    parser.add_argument("-legends", type=str, default="", nargs="+", help="Legends for each line.")
    parser.add_argument("-nln", "--noline", type=lambda x: (str(x).lower() in ['true', '1', 'yes']), default=False,
                        help="Plot without lines.")

    args = parser.parse_args()

    plot_ascii_files(args.file_paths, args.xscale, args.yscale, args.title, args.xlabel, args.ylabel,
                     int(args.skiprows), int(args.droprows), args.column_plot_for_X, args.columns_plot_for_Y, args.shift_first_file, args.color_first_file, 
                     args.pow_for_X, args.pow_for_Y, args.multiplier_for_X, args.multiplier_for_Y, args.fontsize, args.linewidth, args.xlim, args.ylim, args.data_marker, args.legends, args.noline)
