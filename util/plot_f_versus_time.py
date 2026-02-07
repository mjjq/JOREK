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
from tqdm import tqdm
import matplotlib.ticker as ticker
from scipy import integrate

def find_nearest(array, target):
    idx = (np.abs(array - target)).argmin()
    return idx


def plot_radial_integration(times_2d, radial_2d, data_2d_plt, t_start, t_end):
    # Determine which integration function is available
    if hasattr(integrate, 'trapezoid'):
        integrator = integrate.trapezoid
    else:
        integrator = integrate.trapz  # Fallback for older versions

    # Extract 1D vectors
    times_1d = times_2d[:, 0]
    radial_1d = radial_2d[0, :]
    
    # Create mask and slice data
    time_mask = (times_1d >= t_start) & (times_1d <= t_end)
    selected_times = times_1d[time_mask]
    selected_data = data_2d_plt[time_mask, :]
    
    # Perform integration using the detected integrator
    integrated_radial_profile = integrator(selected_data, x=selected_times, axis=0)
    
    # Plotting
    plt.figure(figsize=(8, 5))
    plt.plot(radial_1d, integrated_radial_profile, linewidth=2, color='tab:red')
    plt.title(f'Integrated Profile from $t={t_start}$ to $t={t_end}$')
    plt.xlabel('Radius')
    plt.ylabel('Integrated Signal')
    plt.grid(True)
    plt.show()

    return radial_1d, integrated_radial_profile

def analyze_radial_profile(times_2d, radial_2d, data_2d_plt, t_start, t_end):
    # Compatibility check for SciPy version
    integrator = integrate.trapezoid if hasattr(integrate, 'trapezoid') else integrate.trapz

    # 1. Extract 1D vectors
    times_1d = times_2d[:, 0]
    radial_1d = radial_2d[0, :]
    
    # 2. Mask the time window
    time_mask = (times_1d >= t_start) & (times_1d <= t_end)
    selected_times = times_1d[time_mask]
    selected_data = data_2d_plt[time_mask, :]
    
    if len(selected_times) < 2:
        raise ValueError("The selected time range is too small or contains no data points.")

    # 3. Calculate the Integral
    integrated_profile = integrator(selected_data, x=selected_times, axis=0)
    
    # 4. Calculate the Time-Average
    # We divide by the actual duration of the sliced data
    duration = selected_times[-1] - selected_times[0]
    averaged_profile = integrated_profile / duration
    
    # 5. Plotting the Average Profile
    plt.figure(figsize=(8, 5))
    plt.plot(radial_1d, averaged_profile, linewidth=2, color='darkgreen', label='Time-Average')
    
    plt.title(f'Time-Averaged Radial Profile\n(Interval: {selected_times[0]:.2f} to {selected_times[-1]:.2f})')
    plt.xlabel(r'Radial Coordinate ($r$)')
    plt.ylabel(r'Averaged Value $\langle Data \rangle_t$')
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend()
    plt.show()

    return radial_1d, integrated_profile, averaged_profile

# Example call:
# r, integral, average = analyze_radial_profile(times_2d, radial_2d, data_2d_plt, 0.5, 1.5)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot radial distribution of a variable (y-axis) versus time (x-axis) from postproc files.")
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
    parser.add_argument("-xm", "--multiplier_for_X", type=float, default=1.0,
                        help="Multiple the X-data with the specified value.")
    parser.add_argument("-cm", "--multiplier_for_colored_data", type=float, default=1.0,
                        help="Multiple the colored-data with the specified value.")
    parser.add_argument("-tslc", "--time_slice", type=float, default=-1.e9,
                        help="Specify the time point to plot and save 1D profile vs radial coordinate.")
    parser.add_argument("-rslc", "--radial_slice", type=float, default=-1.e9,
                        help="Specify the raidal coordinate point to plot and save 1D profile vs time.")
    parser.add_argument("-tstart", "--tstart_avg_intg", type=float, default=-1.e9,
                        help="Specify the start time for averaging/integrating data over time.")
    parser.add_argument("-tend", "--tend_avg_intg", type=float, default=-1.e9,
                        help="Specify the end time for averaging/integrating data over time.")

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
    for i, file_name in tqdm(enumerate(sorted_file_names), total = len(sorted_file_names), desc="Loading files"):
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
            # for i, step in enumerate(step_all_1d):
            #     if it == step:
            #         times_1d_array.append(times_all_1d[i])
            #         break
            #     if i == len(step_all_1d) - 1:
            #         print("Did not find step = ", it, " in " + timefile)
    else:
        print("times.dat required, please firstly execute $ plot_live_data.sh -q time ")
        exit(-1)

    # numbers_1d = np.array(numbers_1d_array)
    times_1d = np.array(times_1d_array)
    psi_n_1d = np.array(radial_1d_array)
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

    times_1d = times_1d * args.multiplier_for_X
    data_2d = data_2d * args.multiplier_for_colored_data
    
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

    plt.figure(figsize=(10, 5)) 
    n_level = 40
    if args.clim:
        contour = plt.contourf(times_2d, radial_2d, data_2d_plt, n_level, levels = np.linspace(args.clim[0], args.clim[1], n_level), cmap='PiYG_r')#, extend = 'both')
    else:
        contour = plt.contourf(times_2d, radial_2d, data_2d_plt, n_level, cmap='PiYG_r')

    cbar = plt.colorbar(contour)
    # formatter = ticker.FormatStrFormatter('%.2e')
    # cbar.ax.yaxis.set_major_formatter(formatter)
    # cbar.set_label(r'$\varepsilon_\parallel$ [V/m]')
    # cbar.locator = ticker.MaxNLocator(nbins=5)

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
    
    if args.time_slice != -1.e9:
        id_target = np.abs(times_1d - args.time_slice).argmin()
        data_slice_plt = data_2d[id_target,:]
        file_save = 'postproc/' + variable_name + '_profile_at_t_eq_' + '{:.2f}'.format(times_1d[id_target]) + '.dat'
        data_save=np.column_stack((radial_1d, data_slice_plt))
        np.savetxt(file_save, data_save, delimiter=' ', comments='', header = 'rhop ' + variable_name +'-profile')
        print(variable_name + "-profile is saved to " + file_save)

        plt.figure(figsize=(7, 5))
        plt.plot(radial_1d, data_slice_plt, linewidth = 2)
        plt.xlabel(coord_name if args.ylabel == "" else args.ylabel, fontsize=font_size)
        plt.ylabel(variable_name + ' profile', fontsize=font_size)
        plt.title(variable_name + ' profile @ time = ' + '{:.2f}'.format(times_1d[id_target]), fontsize=font_size)
        plt.xticks(fontsize=font_size)
        plt.yticks(fontsize=font_size)
        plt.xlim(min(radial_1d), max(radial_1d))
        plt.grid()
        plt.tight_layout()
    
    if args.radial_slice != -1.e9:
        id_target = np.abs(radial_1d - args.radial_slice).argmin()
        data_slice_plt = data_2d[:,id_target]
        file_save = 'postproc/' + variable_name + '_vs_time_at_' + coord_name + '_eq_' + '{:.2f}'.format(radial_1d[id_target]) + '.dat'
        data_save=np.column_stack((times_1d, data_slice_plt))
        np.savetxt(file_save, data_save, delimiter=' ', comments='', header = 'time ' + variable_name + '-value')
        print(variable_name + " vs time is saved to " + file_save)
        
        plt.figure(figsize=(10, 5))
        plt.plot(times_1d, data_slice_plt, linewidth = 2)
        plt.xlabel(r'normalized time' if args.xlabel == "" else args.xlabel, fontsize=font_size)
        plt.ylabel(variable_name + ' value', fontsize=font_size)
        plt.title(variable_name + r' vs time @ ' + coord_name + ' = ' + '{:.2f}'.format(radial_1d[id_target]), fontsize=font_size)
        plt.xticks(fontsize=font_size)
        plt.yticks(fontsize=font_size)
        plt.xlim(min(times_1d), max(times_1d))
        plt.grid()
        plt.tight_layout()    

    plt.show()


    if args.tstart_avg_intg >= 0 and args.tend_avg_intg > 0:
        t_start = args.tstart_avg_intg
        t_end   = args.tend_avg_intg
        # r_coord, integrated_profile = plot_radial_integration(times_2d, radial_2d, data_2d, t_start, t_end)
        r_coord, integrated_profile, averaged_profile = analyze_radial_profile(times_2d, radial_2d, data_2d, t_start, t_end)
        file_save = 'postproc/avg_' + variable_name + '_profile_during_t_eq_' + '{:.2f}'.format(t_start) + '_' + '{:.2f}'.format(t_end) + '.dat'
        data_save=np.column_stack((r_coord, averaged_profile))
        np.savetxt(file_save, data_save, delimiter=' ', comments='', header = 'rhop ' + variable_name +'-profile')
        print('avg_' + variable_name + "-profile is saved to " + file_save)


    if False: # plot Jtor over time, $ plot_f_versus_time.py -xm $time2si -fn exprs_averaged -yc 8 -xl "time [ms]" -rp 0.5 -clim 0 3.3 
        rmin = 0.07
        rmax = 0.14
        tmin = 0 # ms
        tid = find_nearest(times_1d, tmin)
        idmin = find_nearest(radial_1d, rmin)
        idmax = find_nearest(radial_1d, rmax)
        print("idmin = ", idmin, ", rmin = ", radial_1d[idmin])
        print("idmax = ", idmax, ", rmax = ", radial_1d[idmax])
        # plot 1D radial averaged quantities over time
        plt.figure(figsize=(10, 5))
        data_1D_intg = np.trapz(data_2d_plt[:,idmin:idmax+1], psi_n_1d[idmin:idmax+1], axis = 1) / (psi_n_1d[idmax] - psi_n_1d[idmin])
        # data_1D_intg = np.mean(data_2d_plt[:,idmin:idmax+1], axis = 1)
        # data_1D_intg = np.max(data_2d_plt[:,idmin:idmax+1], axis = 1)
        plt.plot(times_1d, data_1D_intg, linewidth =2, label='Original data')

        time_extrap = np.linspace(times_1d[tid], 400, 20000)
        degree = 2
        coefficients = np.polyfit(times_1d[tid:], data_1D_intg[tid:], degree)
        fitted_curve = np.poly1d(coefficients)
        # fitted_equation_str = f'{coefficients[0]:.2e}x^2 + {coefficients[1]:.2e}x + {coefficients[2]:.2e}'
        # For higher-degree polynomials, generalize the string representation
        func_coeff_reverse = coefficients[::-1]  # Reverse the order of coefficients
        terms = [f"{coeff:.2e}" if deg == 0 else f"{coeff:.2e}$t^{deg}$" for deg, coeff in enumerate(func_coeff_reverse)]
        fitted_equation_str = "f = " + f"\n   + ".join(terms)
        values_extrap = fitted_curve(time_extrap)
        plt.plot(time_extrap, values_extrap, '--',linewidth =2, label=f'Extrap: {fitted_equation_str}')

        plt.xlabel(r'normalized time' if args.xlabel == "" else args.xlabel, fontsize=font_size)
        plt.ylabel(r'radially averaged ' + variable_name, fontsize=font_size)
        plt.title(r'radially averaged ' + variable_name + ' versus time', fontsize=font_size)
        plt.xticks(fontsize=font_size)
        plt.yticks(fontsize=font_size)
        plt.legend(fontsize=font_size)
        plt.grid()
        
        if args.xlim:
            plt.xlim(xmin, xmax)
        plt.ylim(3.00, 3.35)
        plt.show()
    
    if False: # plot q over time, $ plot_f_versus_time.py -xm $time2si -fn qprofile -yc 1 -xl "time [ms]" -rp 0.5 
        rmin = 0.115
        rmax = 0.185
        tmin = 0 # ms
        tid = find_nearest(times_1d, tmin)
        idmin = find_nearest(radial_1d, rmin)
        idmax = find_nearest(radial_1d, rmax)
        print("idmin = ", idmin, ", rmin = ", radial_1d[idmin])
        print("idmax = ", idmax, ", rmax = ", radial_1d[idmax])
        # plot 1D radial averaged quantities over time
        plt.figure(figsize=(10, 5))
        data_2d_plt = np.abs(data_2d_plt)
        data_1D_intg = np.trapz(data_2d_plt[:,idmin:idmax+1], psi_n_1d[idmin:idmax+1], axis = 1) / (psi_n_1d[idmax] - psi_n_1d[idmin])
        # data_1D_intg = np.mean(data_2d_plt[:,idmin:idmax+1], axis = 1)
        # data_1D_intg = np.max(data_2d_plt[:,idmin:idmax+1], axis = 1)
        # data_1D_intg = np.min(data_2d_plt, axis = 1)
        plt.plot(times_1d, data_1D_intg, linewidth =2, label='Original data')

        time_extrap = np.linspace(times_1d[tid], 400, 20000)
        degree = 2
        coefficients = np.polyfit(times_1d[tid:], data_1D_intg[tid:], degree)
        fitted_curve = np.poly1d(coefficients)
        # fitted_equation_str = f'{coefficients[0]:.2e}x^2 + {coefficients[1]:.2e}x + {coefficients[2]:.2e}'
        # For higher-degree polynomials, generalize the string representation
        func_coeff_reverse = coefficients[::-1]  # Reverse the order of coefficients
        terms = [f"{coeff:.2e}" if deg == 0 else f"{coeff:.2e}$t^{deg}$" for deg, coeff in enumerate(func_coeff_reverse)]
        fitted_equation_str = "f = " + f"\n   + ".join(terms)
        values_extrap = fitted_curve(time_extrap)
        plt.plot(time_extrap, values_extrap, '--',linewidth =2, label=f'Extrap: {fitted_equation_str}')

        plt.xlabel(r'normalized time' if args.xlabel == "" else args.xlabel, fontsize=font_size)
        plt.ylabel(r'radially averaged q', fontsize=font_size)
        plt.title(r'radially averaged q versus time', fontsize=font_size)
        plt.xticks(fontsize=font_size)
        plt.yticks(fontsize=font_size)
        plt.legend(fontsize=font_size)
        plt.grid()
        
        if args.xlim:
            plt.xlim(xmin, xmax)
        plt.ylim(0.94, 1.04)
        plt.show()