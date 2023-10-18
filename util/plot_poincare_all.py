#! /usr/bin/env python3

import matplotlib.pyplot as plt
import argparse
import random


# Define a function to read data from the file and create scatter plots for each group
def read_and_plot_data(filename, xlabel_name, ylabel_name, lines_to_skip=1):
    plt.figure(figsize=(7, 6))
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
    for i, group in enumerate(data):
        # group = group[1:]  # Skip the first line
        x, y = zip(*group)
        # color = "#{:02x}{:02x}{:02x}".format(random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
        color = 'blue'
        plt.scatter(x, y, marker='.', alpha=1.0, s=1, c=color)
        xmin = min(xmin, min(x))
        xmax = max(xmax, max(x))
        ymin = min(ymin, min(y))
        ymax = max(ymax, max(y))

    plt.xlabel(xlabel_name)
    plt.ylabel(ylabel_name)
    plt.title('Poincaré plot in ' + xlabel_name + '-' + ylabel_name)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.grid(True)
    # Adjust the aspect ratio to make axis scales equal
    if xlabel_name == 'R':
        plt.gca().set_aspect('equal', adjustable='box')
    # else:
        # plt.gca().set_aspect('auto')

    plt.tight_layout()
    # plt.legend()
    # plt.show()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot Poincaré plot.")
    parser.add_argument("file_paths", nargs="+",
                        help="filename contains the data.")
    args = parser.parse_args()

    filenames = args.file_paths

    # Replace 'poinc_R-Z.dat' with the path to your file
    # filename = 'poinc_R-Z.dat'
    lines_to_skip = 1
    for filename in filenames:
        if filename == 'poinc_R-Z.dat':
            xlabel_name = 'R'
            ylabel_name = 'Z'
        elif filename == 'poinc_rho-theta.dat':
            xlabel_name = r'$\rho$'
            ylabel_name = r'$\theta$'
            lines_to_skip = 4
        else:
            xlabel_name = 'X'
            ylabel_name = 'Y'

        read_and_plot_data(filename, xlabel_name, ylabel_name, lines_to_skip)

    plt.show()
