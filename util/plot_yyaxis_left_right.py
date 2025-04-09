#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import argparse
import plot_multiple_files as pltmf

# Set global font properties using rcParams
plt.rcParams.update({
    'axes.formatter.limits': (-3, 3),
    'axes.formatter.useoffset': True,  # Use scientific notation
    'font.size': 18,                # Global font size
    'lines.linewidth': 3.0,         # Global line width
    'axes.linewidth': 3.0           # Set line width for both x and y axes
    # 'font.family': 'serif',         # Use serif fonts, which includes Times New Roman
    # 'font.serif': ['Times New Roman'],  # Specifically set Times New Roman
    # 'font.weight': 'bold'           # Font weight for all text elements
    # 'legend.fontsize': 18,          # Font size for the legend
    # 'axes.labelsize': 18,           # Font size for axis labels
    # 'axes.titlesize': 16,           # Font size for title
    # 'xtick.labelsize': 12,          # Font size for x-tick labels
    # 'ytick.labelsize': 12,          # Font size for y-tick labels
})

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot lines from multiple ASCII files.")
    parser.add_argument("file_paths", nargs="+", help="List of file paths to the ASCII files.")
    parser.add_argument('-xlim', nargs=2, type=float, metavar=('xmin', 'xmax'),
                        help='Specify xmin, xmax for the plot.')
    parser.add_argument('-ylim1', nargs=2, type=float, metavar=('ymin', 'ymax'),
                        help='Specify ymin, ymax for the plot of the 1st line.')
    parser.add_argument('-ylim2', nargs=2, type=float, metavar=('ymin', 'ymax'),
                        help='Specify ymin, ymax for the plot of the 2nd line.')
    parser.add_argument("-legends", type=str, default="", nargs="+", help="Legends for each line.")
    parser.add_argument("-xl", "--xlabel", default="", help="Manually define the xlabel.")
    parser.add_argument("-yl1", "--ylabel1", default="", help="Manually define the ylabel of the 1st line.")
    parser.add_argument("-yl2", "--ylabel2", default="", help="Manually define the ylabel of the 2nd line.")
    parser.add_argument("-ti", "--title", default="Multiple ASCII Files Plot", help="Manually define the title.")
    parser.add_argument("-sk", "--skiprows", default=1,type=int,
                        help="Define how many rows will be skiped when reading the data.")
    parser.add_argument("-xs", "--xscale", choices=["linear", "log"], default="linear",
                        help="X-scale (linear or log). Default is linear.")
    parser.add_argument("-ys1", "--yscale1", choices=["linear", "log"], default="linear",
                        help="Y-scale (linear or log). Default is linear of the 1st line.")
    parser.add_argument("-ys2", "--yscale2", choices=["linear", "log"], default="linear",
                        help="Y-scale (linear or log). Default is linear of the 2nd line.")

    parser.add_argument("-xc1", "--column_plot_for_X1", type=int, default=0,
                        help="Columns to be plotted (default is 0) as X value of the 1st line.")
    parser.add_argument("-xc2", "--column_plot_for_X2", type=int, default=0,
                        help="Columns to be plotted (default is 0) as X value of the 2nd line.")
    
    parser.add_argument("-yc1", "--column_plot_for_Y1", type=int, nargs="+", default=1,
                        help="Columns to be plotted (default is 1) as Y values of the 1st line.")                                                
    parser.add_argument("-yc2", "--column_plot_for_Y2", type=int, nargs="+", default=1,
                        help="Columns to be plotted (default is 1) as Y values of the 2nd line.")

    args = parser.parse_args()
        
    column_names, data = pltmf.read_ascii_file(args.file_paths[0], args.skiprows)
    x1 = data[:, args.column_plot_for_X1]
    y1 = data[:, args.column_plot_for_Y2]

    # column_names, data = pltmf.read_ascii_file('/u/hazh/WORKSPACE/20230825_fix_axis_bug_get_all_commits/Base_case2D/' + args.file_paths[0], args.skiprows)
    # x11 = data[:, args.column_plot_for_X1]
    # y11 = data[:, args.column_plot_for_Y2]

    column_names, data = pltmf.read_ascii_file(args.file_paths[1], args.skiprows)
    x2 = data[:, args.column_plot_for_X2]
    y2 = data[:, args.column_plot_for_Y2]

    # column_names, data = pltmf.read_ascii_file('/u/hazh/WORKSPACE/20230825_fix_axis_bug_get_all_commits/Base_case2D/' + args.file_paths[1], args.skiprows)
    # x22 = data[:, args.column_plot_for_X2]
    # y22 = data[:, args.column_plot_for_Y2]

    # Create figure and primary axis
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot first dataset on primary axis (left)
    color = 'tab:blue'
    ax1.set_xlabel(args.xlabel)
    ax1.set_ylabel(args.ylabel1, color=color)
    ax1.plot(x1, y1, color=color, label=args.legends[0])
    # ax1.plot(x11, y11, color=color, linestyle = 'dotted', label='$q_0$ (2D case)' , marker = '')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid()
    # plt.legend()
    plt.yscale(args.yscale1)  # Set linear or log for the y scale
    if args.ylim1:
        ymin, ymax = args.ylim1
        plt.ylim(ymin, ymax)

    # Create secondary axis (right)
    ax2 = ax1.twinx()

    # Plot second dataset on secondary axis
    color = 'tab:red'
    ax2.set_ylabel(args.ylabel2, color=color)
    ax2.plot(x2, y2, linestyle = '--', color=color, label=args.legends[1])
    # ax2.plot(x22, y22, linestyle = ':', color=color, label='$p_0$ (2D case)', marker = '')
    ax2.tick_params(axis='y', labelcolor=color)
    # ax2.grid()
    # plt.legend()
    plt.yscale(args.yscale2)  # Set linear or log for the y scale
    if args.ylim2:
        ymin, ymax = args.ylim2
        plt.ylim(ymin, ymax)

    if args.xlim:
        xmin, xmax = args.xlim
        plt.xlim(xmin, xmax)
    
    plt.xscale(args.xscale)  # Set linear or log for the x scale

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    legend = ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
    legend.set_draggable(True)

    # Title and show
    plt.title(args.title)

    plt.tight_layout()
    plt.show()