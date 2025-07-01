#! /usr/bin/env python3

import matplotlib.pyplot as plt
import argparse
import random
import numpy as np
import plot_q_versus_time as plot_q_time
import os
import plot_multiple_files as pltmpf

def get_time(num, path = './'):
    timefile = path + 'times.dat'
    time = 0
    if os.path.exists(timefile):
        column_names, data = pltmpf.read_ascii_file(timefile, 1)
        step_all_1d = data[:,0]
        times_all_1d = data[:,1]
        for i, step in enumerate(step_all_1d):
            if num == step:
                time = times_all_1d[i]
    else:
        print("times.dat required, please firstly execute $ plot_live_data.sh -q time ")
        # exit(-1)
        time = 0

    return time

# Define a function to read data from the file and create scatter plots for each group
def read_and_plot_data(filename, xlabel_name, ylabel_name, lines_to_skip=1, font_size = 15, multiplier_for_time = 1.0, colorful = False):
    fig = plt.figure(figsize=(7, 6))
    data = []
    current_group = []
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if i < lines_to_skip:
                continue  # Skip the specified number of lines
            line = line.strip()
            if line:  # If the line is not empty
                parts = line.split()
                if len(parts) == 2:
                    # Assuming each line contains two columns of data
                    x, y = map(float, parts)
                    current_group.append((x, y))
            elif current_group:
                # If an empty line is encountered and we have data in the current group
                data.append(current_group)
                current_group = []

    if current_group:
        data.append(current_group)

    xmin = ymin = 1e5
    xmax = ymax = -1e5
    # Plot each group with a different color
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    for i, group in enumerate(data):
        # group = group[1:]  # Skip the first line
        x, y = zip(*group)
        
        if colorful:
            color = colors[i % len(colors)]
            # color = "#{:02x}{:02x}{:02x}".format(random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
        else:
            color = 'blue'
        plt.scatter(x, y, marker='.', alpha=1.0, s=1, c=color)
        xmin = min(xmin, min(x))
        xmax = max(xmax, max(x))
        ymin = min(ymin, min(y))
        ymax = max(ymax, max(y))

    plt.xlabel(xlabel_name, fontsize = font_size)
    plt.ylabel(ylabel_name, fontsize = font_size)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    # Axis tick label font size
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.grid(True)
    # Adjust the aspect ratio to make axis scales equal
    if xlabel_name == 'R':
        plt.gca().set_aspect('equal', adjustable='box')
    # else:
        # plt.gca().set_aspect('auto')
    title_str = 'Poincaré plot in ' + xlabel_name + '-' + ylabel_name
    if('poinc_R-Z.dat' in filename or 'poinc_rho-theta.dat' in filename):
        plt.title(title_str, fontsize = font_size)
    else:
        num = plot_q_time.extract_number(filename)
        time = get_time(num)
        if multiplier_for_time != 1.0:
            time = time * multiplier_for_time
        
        time_str = "{:.3f}".format(time)
        title_str = title_str + ', '
        title_ext = r'normalized time = ' + time_str if multiplier_for_time == 1.0 else r'time = ' + time_str + '[ms]'
        plt.title(title_str + title_ext, fontsize = font_size)
    plt.tight_layout()
    # plt.legend()
    # plt.show()
    return fig

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot Poincaré plot.")
    parser.add_argument("file_paths", nargs="+",
                        help="filename contains the data.")
    parser.add_argument("-fz", "--fontsize", type=int, default=15,
                        help="Define the fontsize of figure.")
    parser.add_argument("-tm", "--multiplier_for_time", type=float, default=1.0,
                        help="Multiple the time-data with the specified value.")
    parser.add_argument("-cful", "--colorful", type=lambda x: (str(x).lower() in ['true', '1', 'yes']), default=False,
                        help="Set different colors for lines.")
    parser.add_argument("-png", "--png_output", type=lambda x: (str(x).lower() in ['true', '1', 'yes']), default=False,
                        help="Save results into figures.")
    parser.add_argument('-xylim', nargs=4, type=float, metavar=('xmin', 'xmax', 'ymin', 'ymax'),
                        help='Specify xmin, xmax, ymin, ymax for the plot.')
    args = parser.parse_args()

    filenames = args.file_paths

    # Replace 'poinc_R-Z.dat' with the path to your file
    # filename = 'poinc_R-Z.dat'
    lines_to_skip = 1
    for filename in filenames:
        if 'poinc_R-Z' in filename: # 'poinc_R-Z.dat':
            xlabel_name = 'R'
            ylabel_name = 'Z'
        elif 'poinc_rho-theta' in filename: # 'poinc_rho-theta.dat':
            xlabel_name = r'$\rho$'
            ylabel_name = r'$\theta$'
            lines_to_skip = 4
        else:
            xlabel_name = 'X'
            ylabel_name = 'Y'

        fig = read_and_plot_data(filename, xlabel_name, ylabel_name, lines_to_skip, args.fontsize, args.multiplier_for_time, args.colorful)
        
        if(args.png_output):
            figure_name = filename[:-4] + ".png"
            plt.savefig(figure_name)
            print(figure_name + " is saved.")
            plt.close()    
    
    if args.xylim:
        xmin, xmax, ymin, ymax = args.xylim
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)

    if(not args.png_output):
        plt.show()
