#! /usr/bin/env python3

import re
import numpy as np
import h5py
import argparse
import matplotlib.pyplot as plt

def extract_floats(file_path):
    # Define the regex pattern to match float numbers, including negative numbers and those with decimal points
    pattern = r'[-+]?\d*\.\d+'
    # Find all matches in the string
    floats = re.findall(pattern, file_path)
    # Convert the matches to float
    float_numbers = [float(num) for num in floats]
    return float_numbers

def extract_variable_names(file_path):
    variable_names = []

    # Open the file in read mode
    with open(file_path, 'r') as file:
        # Read each line in the file
        for line in file:
            # Check if the line starts with "expressions"
            if line.strip().startswith('expressions'):
                # Split the line by whitespace and extract variable names
                parts = line.strip().split()
                variable_names.extend(parts[1:])  # Exclude 'expressions' from the list

    return variable_names


def plot_3d_dataset(file_path, dataset_name, iphi_plot):

    print("")
    print("Plotting raw data ......")
    
    float_numbers = extract_floats(file_path)
    Rmin = float_numbers[0]
    Rmax = float_numbers[1]
    Zmin = float_numbers[2]
    Zmax = float_numbers[3]
    Phimin  = float_numbers[4]
    Phimax  = float_numbers[5]
    print(f"Rmin: {Rmin}, Rmax: {Rmax}, Zmin: {Zmin}, Zmax: {Zmax}, Phimin: {Phimin}, Phimax: {Phimax}")

    # Example usage:
    variable_names = extract_variable_names("postproc_script_rectangular.txt")
    print(variable_names)

    # Open the .h5 file
    with h5py.File(file_path, 'r') as hdf:
        # Check if the dataset exists in the file
        if dataset_name in hdf:
            data = hdf[dataset_name][:]

            Phi_grid = np.linspace(Phimin, Phimax, hdf["dim"][0])
            R_grid = np.linspace(Rmin, Rmax, hdf["dim"][1])
            Z_grid = np.linspace(Zmin, Zmax, hdf["dim"][2])

            X, Y = np.meshgrid(R_grid, Z_grid)
            
            if iphi_plot >= hdf["dim"][0]:
                print(f'iphi_plot cannot be large than {hdf["dim"][0]}')
                exit(-1)
            # If the data size is 4D (3d raw data), plot it
            if data.ndim == 4:
                fig, axes = plt.subplots(nrows=3, ncols=int((len(variable_names)-1+2)/3), figsize=(15, 10))
                # Flatten the axes array for easier indexing
                axes = axes.flatten()
                for i in range(1,len(variable_names)):
                    ax = axes.flat[i-1]
                    im = ax.contourf(X, Y, data[i,:,:,iphi_plot], 50, cmap='jet') # viridis_r
                    ax.set_title(f'{variable_names[i]}')
                    ax.set_xlabel('R')
                    ax.set_ylabel('Z')
                    ax.margins(x=0)
                    ax.margins(y=0)
                    ax.set_aspect('equal', adjustable='box')
                    fig.colorbar(im, ax=ax, orientation='vertical', shrink=0.5)

                fluxdata = np.loadtxt('fluxsurfaces.dat')
                ax = axes.flat[len(variable_names)-1]
                ax.scatter(fluxdata[:, 0], fluxdata[:, 1], s=1)
                ax.set_xlabel('R')
                ax.set_ylabel('Z')
                ax.set_title(f'Flux Surfaces')
                ax.set_aspect('equal', adjustable='box')

                # Turn off the axis for any remaining subplots
                for ax in axes[len(variable_names):]:
                    ax.axis('off')
                fig.suptitle(f'3D data @ phi = {Phi_grid[iphi_plot]}, time = '+ "{:.2f}".format(hdf["time"][0]*1e3) + " [ms]")
                plt.tight_layout()
                plt.show()
            else:
                print(f"Dataset {dataset_name} is not 4D. It has {data.ndim} dimensions.")
                print(f"Dataset {dataset_name} is {data}")
        else:
            print(f"Dataset {dataset_name} not found in the file.")


def print_datasets_info(file_path):
    # Open the .h5 file
    with h5py.File(file_path, 'r') as hdf:
        print(f"Datasets in the file '{file_path}':")
        # Recursively print all datasets
        def print_info(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"Name: {name}, Shape: {obj.shape}, Data Type: {obj.dtype}")
        
        hdf.visititems(print_info)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot a 3D dataset from an .h5 file.')
    parser.add_argument('file_path', type=str, help='Path to the .h5 file')
    parser.add_argument('-na', '--dataset_name', type=str, default="values" , help='Name of the dataset to plot')
    parser.add_argument('-ip', '--iphi_plot', type=int, default=0, help='The rank of toroidal plane to be plotted.')
    
    args = parser.parse_args()

    # Print datasets information
    print_datasets_info(args.file_path)
    plot_3d_dataset(args.file_path, args.dataset_name, args.iphi_plot)

