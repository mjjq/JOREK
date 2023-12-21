#! /usr/bin/env python3

import h5py

def print_dataset_shapes(file_path):
    # Open the HDF5 file in read mode
    with h5py.File(file_path, 'r') as file:
        # List all groups and datasets in the file
        print("Keys:", list(file.keys()))

        # Iterate over all items in the file
        def print_shapes(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"Dataset: {name}, Shape: {obj.shape}")

        file.visititems(print_shapes)

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python script.py <h5_file_path>")
        sys.exit(1)

    h5_file_path = sys.argv[1]
    print_dataset_shapes(h5_file_path)
