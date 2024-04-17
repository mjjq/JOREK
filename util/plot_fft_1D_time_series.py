#! /usr/bin/env python3

#
# Purpose: Perform FFT on 1D time-series data.
#
# Date: 2023-02-27
# Author: Haowei Zhang, IPP Garching
#

import numpy as np
import matplotlib.pyplot as plt
import argparse
import plot_multiple_files as plmf
import testFFT as myFFT
from scipy import interpolate


def proc_ascii_files(file_paths, title, xlabel, ylabel, skiprows, column_plot_for_X, columns_plot_for_Y, 
                     shift_first_file, multiplier_for_X, multiplier_for_Y, font_size, line_width, xlim, ylim, data_marker, noZeroFreq, complex_data, columns_imag_plot_for_Y):
    for i, file_path in enumerate(file_paths):

        column_names, data = plmf.read_ascii_file(file_path, skiprows)
        if column_plot_for_X >= data.shape[1]:
            print("Error: Value of -xc is out of range!!!")
            exit(-1)
        x = data[:, column_plot_for_X]  # Assuming the first column is the x-axis data
        y = data[:, 1:]  # Assuming the remaining columns are y-axis data
        if multiplier_for_X != 1.0:
            x = x * multiplier_for_X
        
        if multiplier_for_Y != 1.0:
            y = y * multiplier_for_Y
        
        if i == 0:
            x = x + float(shift_first_file)

        if columns_plot_for_Y == -1:
            columns_plot_for_Y_mdf = range(y.shape[1])  # Plot all columns
            if(complex_data):
                complex_data = False
                print("WARNING: Disable imaginary data in case of all columns being plotted.")
        elif i == 0:
            columns_plot_for_Y_mdf = [j - 1 for j in columns_plot_for_Y]
            if(complex_data):
                columns_imag_plot_for_Y_mdf = [j - 1 for j in columns_imag_plot_for_Y]

        data_len = len(x) * 10
        fft_window_width = int(data_len / 20)
        time_1d = np.linspace(np.min(x), np.max(x), data_len)
        print("data length              = ", data_len)
        print("fft window width         = ", fft_window_width)
        print("fft window width in time = ", fft_window_width * (np.max(x) - np.min(x))/(data_len - 1.0))

        for icol, j in enumerate(columns_plot_for_Y_mdf):
            if j >= y.shape[1]:
                print("Error: Value of -yc is out of range!!!")
                print("file_path: ", file_path)
                exit(-2)
            if(complex_data):
                f_interp = interpolate.interp1d(x, y[:,j] + 1j * y[:,columns_imag_plot_for_Y_mdf[icol]])
            else:
                f_interp = interpolate.interp1d(x, y[:,j])
            data_1d = f_interp(time_1d)
            time_new, signal_new, times_2d, freq_2d, spec_2d = myFFT.window_fft(time_1d, data_1d, fft_window_width)
            if noZeroFreq:
                spec_2d[:,0] = 0.e0
            
            # spec_2d[np.where(times_2d < 20)] = 0.e0

            # Plot the signal
            plt.figure(figsize=(12, 9))  # Adjust the figure size as needed
            plt.subplot(2, 1, 1)
            if(complex_data):
                plt.plot(x, y[:,j], data_marker, label=f'{column_names[j + 1]} ({file_path})')
                plt.plot(x, y[:,columns_imag_plot_for_Y_mdf[icol]], data_marker, label=f'{column_names[columns_imag_plot_for_Y_mdf[icol] + 1]} ({file_path})')
                plt.plot(time_new, signal_new.real, '', linestyle='-', linewidth = line_width, label='Time Domain Signal (real, interpolated)')
                plt.plot(time_new, signal_new.imag, '', linestyle='--', linewidth = line_width, label='Time Domain Signal (imag, interpolated)')
            else:
                plt.plot(x, y[:,j], data_marker, label=f'{column_names[j + 1]} ({file_path})')
                plt.plot(time_new, signal_new, '', linestyle='--', linewidth = line_width, label='Time Domain Signal (interpolated)')
            if xlim:
                xmin, xmax = xlim
                plt.xlim(xmin, xmax)
            else:
                plt.xlim(min(time_new), max(time_new))
            plt.ylabel('Signal Amplitude', fontsize=font_size)
            plt.grid(True)
            plt.xticks(fontsize=font_size)
            plt.yticks(fontsize=font_size)     
            plt.gca().get_yaxis().get_offset_text().set_size(font_size)
            legend = plt.legend()  # Show legends for each line
            legend.set_draggable(True)
            for label in legend.get_texts():
                label.set_fontsize(10)
            
            plt.subplot(2, 1, 2)
            plt.contourf(times_2d, freq_2d, spec_2d, 20, cmap='jet')
            plt.ylim(0,10)
            plt.xlabel(column_names[column_plot_for_X] if xlabel == "" else xlabel, fontsize=font_size)
            plt.ylabel('Frequency' if ylabel == "" else ylabel, fontsize=font_size)  # Set the y-axis label (customize as needed)
            plt.colorbar(orientation='horizontal', shrink=0.6)
            if xlim:
                xmin, xmax = xlim
                plt.xlim(xmin, xmax)
            else:
                plt.xlim(min(time_new), max(time_new))
            if ylim:
                ymin, ymax = ylim
                plt.ylim(ymin, ymax)     
            plt.xticks(fontsize=font_size)
            plt.yticks(fontsize=font_size)     
            plt.gca().get_yaxis().get_offset_text().set_size(font_size)           
            # plt.title(title, fontsize=font_size)  # Set the plot title
            plt.tight_layout()
    plt.show()
            


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot lines from multiple ASCII files.")
    parser.add_argument("-ti", "--title", default="Multiple ASCII Files Plot", help="Manually define the title.")
    parser.add_argument("-xl", "--xlabel", default="", help="Manually define the xlabel.")
    parser.add_argument("-yl", "--ylabel", default="", help="Manually define the ylabel.")
    parser.add_argument("-dm", "--data_marker", default=".", help="Manually define the marker for the data.")
    parser.add_argument("-yc", "--columns_plot_for_Y", type=int, nargs="+", default=-1,
                        help="Columns (for real data) to be plotted (default is -1 for plotting all) as Y values.")
    parser.add_argument("-xc", "--column_plot_for_X", type=int, default=0,
                        help="Columns to be plotted (default is 0) as X value.")
    parser.add_argument("-sk", "--skiprows", default=1,
                        help="Define how many rows will be skiped when reading the data.")
    parser.add_argument("-sf", "--shift_first_file", default=0,
                        help="Shift the data for the first file horizontally (in X direction).")
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
    parser.add_argument("-cmpx", "--complex_data", default=False,
                        help="Remove the zero frequency component.")
    parser.add_argument("-ycj", "--columns_imag_plot_for_Y", type=int, nargs="+", default=-1,
                        help="Columns (for imaginary data) to be plotted (default is -1 for plotting all) as Y values.")
    parser.add_argument("-noZF", "--noZeroFreq", default=False,
                        help="Remove the zero frequency component.")
    args = parser.parse_args()
    proc_ascii_files(args.file_paths, args.title, args.xlabel, args.ylabel,
                     int(args.skiprows), args.column_plot_for_X, args.columns_plot_for_Y, args.shift_first_file, 
                     args.multiplier_for_X, args.multiplier_for_Y, args.fontsize, args.linewidth, args.xlim, args.ylim, 
                     args.data_marker, args.noZeroFreq, args.complex_data, args.columns_imag_plot_for_Y)