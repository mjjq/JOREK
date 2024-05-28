#! /usr/bin/env python3
import argparse
import h5py
import re
import matplotlib.pyplot as plt
import numpy as np

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


def plot_2d_dataset(file_path, dataset_name):
    
    float_numbers = extract_floats(file_path)
    Rmin = float_numbers[0]
    Rmax = float_numbers[1]
    Zmin = float_numbers[2]
    Zmax = float_numbers[3]
    Phi  = float_numbers[4]
    print(f"Rmin: {Rmin}, Rmax: {Rmax}, Zmin: {Zmin}, Zmax: {Zmax}, Phi: {Phi}")

    # Example usage:
    variable_names = extract_variable_names("postproc_script_rectangular.txt")
    print(variable_names)

    # Open the .h5 file
    with h5py.File(file_path, 'r') as hdf:
        # Check if the dataset exists in the file
        if dataset_name in hdf:
            data = hdf[dataset_name][:]

            R_grid = np.linspace(Rmin, Rmax, hdf["dim"][0])
            Z_grid = np.linspace(Zmin, Zmax, hdf["dim"][1])

            X, Y = np.meshgrid(R_grid, Z_grid)
            
            # If the data is 3D, plot it
            if data.ndim == 3:
                for i in range(len(variable_names)):
                    plt.figure()
                    # plt.imshow(data[1], aspect='auto', interpolation='none')
                    plt.contourf(X, Y, data[i], 50, cmap='jet')
                    plt.colorbar()
                    plt.title(f'2D Dataset: {variable_names[i]}')
                    plt.xlabel('R-axis')
                    plt.ylabel('Z-axis')
                    plt.axis('equal')
                plt.show()
            else:
                print(f"Dataset {dataset_name} is not 2D. It has {data.ndim} dimensions.")
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
    parser = argparse.ArgumentParser(description='Plot a 2D dataset from an .h5 file.')
    parser.add_argument('file_path', type=str, help='Path to the .h5 file')
    parser.add_argument('dataset_name', type=str, help='Name of the dataset to plot')
    
    args = parser.parse_args()

    # Print datasets information
    print_datasets_info(args.file_path)
    plot_2d_dataset(args.file_path, args.dataset_name)

