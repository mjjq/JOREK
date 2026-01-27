#! /usr/bin/env python3

#
# Purpose: Plot current diffusion rate
#
# Date: 2023-09-18
# Author: Haowei Zhang, IPP Garching
#
# Example:
# ~/SCRATCH/pitagora/20230825_fix_axis_bug_get_all_commits/AUG_36663_3.75_ed22_711_00_ultimate1_eta_zk_fg_q0_gt_1_nmax_4_nflux_80_lft_q0_v4_Dpar0.1_current_source_2_thermalization2_correct_zkperp_etaARAZ_simple
# plot_multiple_files_current_diffusion.py postproc/exprs_averaged_s??[0,2,4,6,8]000.dat -yc 6 -xp 1 -cful true -tm $time2si -xlim 0.028 0.22 -ylim ' -4.2e6' 4.2e6  -xl '$\sqrt{\psi_n}$' -yl 'current diffusion rate [A/m$^2$s]' -ti ''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import plot_q_versus_time
import plot_poincare_all
import math
from scipy.ndimage import gaussian_filter1d


time2si         = 4.607490970E-04  # to [ms]
MU_zero         = 4.e-7*math.pi
central_density = 1.01
central_mass    = 1.0
mass_proton     = 1.67262178e-27
rho_norm        = central_density *1.e20 * central_mass * mass_proton
fact_J_mu_zero  = MU_zero
fact_time       = 1/time2si
fact_eta        = 1.0 / math.sqrt ( MU_zero / rho_norm )
fact_VXB        = math.sqrt( MU_zero * central_mass * central_density * 1.e20 * mass_proton )

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
    data = np.genfromtxt(lines[skiprows:len(lines) + 1 + droprows], dtype=float)

    return column_names, data


def plot_ascii_files(file_paths, xscaleType, yscaleType, title, xlabel, ylabel, skiprows, droprows, column_plot_for_X, columns_plot_for_Y, 
                     shift_first_file, color_first_file, pow_for_X, pow_for_Y, multiplier_for_X, multiplier_for_Y, font_size, line_width, xlim, ylim, data_marker, mylegends="", noline = False, colorful = False, multiplier_for_time = 1.0, abs = False):
    """
    Plot lines from multiple ASCII files in the same plot.

    Args:
        file_paths (list): List of file paths to the ASCII files.
    """
    # Define line styles (you can customize this list with different line styles)
    linestyle_list = ['-', '--', '-.', ':']
    # Set the default linewidth for all lines
    mpl.rcParams['lines.linewidth'] = line_width

    # plt.figure(figsize=(10, 6))  # Set the figure size

    if colorful:
        n_curves = len(file_paths)
        cmap = plt.get_cmap('viridis_r')
        colors = cmap(np.linspace(0, 1, n_curves))
        time_array = []

    ilegend = 0
    lines = []
    for i, file_path in enumerate(file_paths):
        
        column_names, data = read_ascii_file(file_path, skiprows, droprows)

        psiN = data[:,column_names.index("rho_p")]**2
        VxBt_n_eq_zero = data[:,column_names.index("VXBt00")]
        VxBt_n_ge_one = data[:,column_names.index("VXBt_all")]
        gradPsi2 = data[:,column_names.index("gradPsi2")]
        lapPsi = data[:,column_names.index("lapPsi")]
        VxBt_n_all = VxBt_n_ge_one + VxBt_n_eq_zero
        sigma_VxB = 3
        sigma_etaJ = 3
        VxBt_n_all = gaussian_filter1d(VxBt_n_all, sigma_VxB)
        ddrho_VxBt_n_all = np.gradient(VxBt_n_all, psiN)
        d2drho2_VxBt_n_all = np.gradient(ddrho_VxBt_n_all, psiN) * gradPsi2 + ddrho_VxBt_n_all * lapPsi
        d2drho2_VxBt_n_all_si = d2drho2_VxBt_n_all / fact_VXB / MU_zero

        etaJtor = data[:,column_names.index("etaJt")]
        etaSJtor = data[:,column_names.index("etaSJt")]
        etaJtor = gaussian_filter1d(etaJtor, sigma_etaJ)
        etaSJtor = gaussian_filter1d(etaSJtor, sigma_etaJ)
        ddrho_etaJtor = np.gradient(etaJtor, psiN)
        d2drho2_etaJtor = np.gradient(ddrho_etaJtor, psiN) * gradPsi2 + ddrho_etaJtor * lapPsi
        ddrho_etaSJtor = np.gradient(etaSJtor, psiN)
        d2drho2_etaSJtor = np.gradient(ddrho_etaSJtor, psiN) * gradPsi2 + ddrho_etaSJtor * lapPsi
        d2drho2_etaJtor_minus_etaSjtor = d2drho2_etaJtor - d2drho2_etaSJtor
        d2drho2_etaJtor_minus_etaSjtor_si = d2drho2_etaJtor_minus_etaSjtor / fact_eta / fact_J_mu_zero / MU_zero

        if colorful:
            num_time = plot_q_versus_time.extract_number(file_path)
            time = plot_poincare_all.get_time(num_time) * multiplier_for_time
            time_array.append(time)
            print(f'time = {time} [ms]')

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

        if abs:
            y = np.abs(y)


        if columns_plot_for_Y == -1:
            columns_plot_for_Y_mdf = range(y.shape[1])  # Plot all columns
        elif i == 0:
            columns_plot_for_Y_mdf = [j - 1 for j in columns_plot_for_Y]

        if not colorful:
            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        color_plot = colors[i % len(colors)]

        if colorful:
            linestyle_plot = linestyle_list[0]
        else:
            linestyle_plot = linestyle_list[i % len(linestyle_list)]  # Cycle through linestyle_list

        for rankj, j in enumerate(columns_plot_for_Y_mdf):
            if len(file_paths) == 1:
                color_plot = colors[rankj % len(colors)]
                linestyle_plot = linestyle_list[rankj % len(linestyle_list)]  # Cycle through linestyle_list
            if noline:
                linestyle_plot = 'none'
            if j < y.shape[1]:
                if mylegends == "":
                    mylegend = f'{column_names[j + 1]} ({file_path})'
                else:
                    mylegend = mylegends[ilegend]
                    ilegend += 1
                if (i == 0) and (color_first_file != ""):
                    color_plot = color_first_file
                # plt.plot(x, y[:, j], data_marker, label=mylegend, linestyle=linestyle_plot, color=color_plot, markersize=20)
                line1, = plt.plot(x, d2drho2_etaJtor_minus_etaSjtor_si, 'o', label=r'$\langle\nabla^2\left[\eta (J_{\varphi} - S_{j})\right]\rangle/\mu_0$', linestyle='none', color=color_plot, markersize=5)
                line2, = plt.plot(x, -d2drho2_VxBt_n_all_si, data_marker, label=r'$-\langle\nabla^2 \left(\mathbf{v}\times\mathbf{B}\right)_{\varphi}\rangle/\mu_0$', linestyle=linestyle_plot, color=color_plot, markersize=20)
                
                lines.append(line1)
                lines.append(line2)
            else:
                print("Error: Value of -yc is out of range!!!")
                print("file_path: ", file_path)
                exit(-2)

    plt.xlabel(column_names[column_plot_for_X] if xlabel == "" else xlabel, fontsize=font_size)  # Set the x-axis label based on the first file
    plt.ylabel('Y-axis Label' if ylabel == "" else ylabel, fontsize=font_size)  # Set the y-axis label (customize as needed)
    if colorful:
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        from matplotlib.ticker import ScalarFormatter

        time_array = np.array(time_array)
        norm = mcolors.Normalize(vmin=time_array.min(), vmax=time_array.max())
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        # Add the colorbar
        cbar = plt.colorbar(sm, ax=plt.gca())
        cbar.set_label('time [ms]')
        # Set colorbar formatter with custom limits
        formatter = ScalarFormatter()
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1000, 1000))  # Only use scientific notation outside this range
        cbar.ax.yaxis.set_major_formatter(formatter)
        legend = plt.legend(handles=lines[-2:], labels=[r'$\langle\nabla^2\left[\eta (J_{\varphi} - S_{j})\right]\rangle/\mu_0$',
                                               r'$-\langle\nabla^2 \left(\mathbf{v}\times\mathbf{B}\right)_{\varphi}\rangle/\mu_0$'])
        legend.set_draggable(True)
        for label in legend.get_texts():
            label.set_fontsize(font_size)
    else:
        # legend = plt.legend(loc="upper left", bbox_to_anchor=(1, 1))  # Show legends for each line
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

    plt.tight_layout()
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
    parser.add_argument("-fz", "--fontsize", default=18,
                        help="Define the fontsize of figure.")
    parser.add_argument("-lw", "--linewidth", default=3,
                        help="Define the linewidth of figure.")
    parser.add_argument('-xlim', nargs=2, type=float, metavar=('xmin', 'xmax'),
                        help='Specify xmin, xmax for the plot.')
    parser.add_argument('-ylim', nargs=2, type=float, metavar=('ymin', 'ymax'),
                        help='Specify ymin, ymax for the plot.')
    parser.add_argument("-tm", "--multiplier_for_time", type=float, default=1.0,
                        help="Multiple the time-data with the specified value.")
    parser.add_argument("-cful", "--colorful", type=lambda x: (str(x).lower() in ['true', '1', 'yes']), default=False,
                        help="Set different colors for lines.")
    parser.add_argument("-legends", type=str, default="", nargs="+", help="Legends for each line.")
    parser.add_argument("-nln", "--noline", type=lambda x: (str(x).lower() in ['true', '1', 'yes']), default=False,
                        help="Plot without lines.")
    parser.add_argument("-abs", "--abs", type=lambda x: (str(x).lower() in ['true', '1', 'yes']), default=False,
                        help="Plot without abs(Y).")

    args = parser.parse_args()

    plot_ascii_files(args.file_paths, args.xscale, args.yscale, args.title, args.xlabel, args.ylabel,
                     int(args.skiprows), int(args.droprows), args.column_plot_for_X, args.columns_plot_for_Y, args.shift_first_file, args.color_first_file, 
                     args.pow_for_X, args.pow_for_Y, args.multiplier_for_X, args.multiplier_for_Y, args.fontsize, args.linewidth, args.xlim, args.ylim, args.data_marker, args.legends, args.noline, args.colorful, args.multiplier_for_time, args.abs)
