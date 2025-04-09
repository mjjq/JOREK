#! /usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
import plot_multiple_files as pltm
import math as math
from scipy.interpolate import interp1d

R_axis = 1.713934595E+00
Z_axis = 7.092228919E-02
psin_truncate = 0.800184
loop_voltage = -5.0
fmu0 = 4 * math.pi * 1e-7
Iptot_target = 8.06e5 # A
Iptot_equilibrium = 8.269383806e+05 # A
flux_surface_file = 'q_surface_-3.5366_s000000.txt'
boundary_file = 'jorek_36663_3.75-3.95_RZ_Boundary_ed22.txt'
boundary_truncated_file =  'jorek_36663_3.75-3.95_truncated_RZ_Boundary_ed22.txt' 
averaged_file = 'exprs_averaged_s000000.dat' # for obtaining the boundary psi
intg3D_file = 'int3Dprofile_s000000.dat' # for obtaining the totoal current IP
input_folder = 'input_profiles/'
output_folder = 'input_profiles_truncated/'
postproc_folder = 'postproc/'

def truncate_and_interpolate_profiles(psin_truncate):

    os.makedirs(output_folder, exist_ok=True)  # Create output folder if it doesn't exist
    
    for filename in os.listdir(input_folder):
        if filename.endswith('.txt') and 'RZ_Boundary' not in filename:
            input_path = os.path.join(input_folder, filename)
            output_path = os.path.join(output_folder, filename)
            
            # Read and process data manually
            with open(input_path, 'r') as f:
                lines = f.readlines()
            
            # Extract data
            psin_values = []
            variable_values = []
            for line in lines:
                parts = line.split()
                psin_values.append(float(parts[0]))
                variable_values.append(float(parts[1]))
            
            # Convert to numpy arrays
            psin_values = np.array(psin_values)
            variable_values = np.array(variable_values)
            
            num_points = len(psin_values) + 1
            # Interpolate to have the same number of points
            psin_new = np.linspace(min(psin_values), psin_truncate, num_points)
            variable_new = np.interp(psin_new, psin_values, variable_values)
            
            # renormalize psin from 0 to 1
            psin_renormalize = psin_new / max(psin_new)
            # Save interpolated data
            with open(output_path, 'w') as f:
                for p, v in zip(psin_renormalize, variable_new):
                    f.write(f"{p:.6f}\t{v:.6f}\n")
            
            print(f'Truncated and interpolated file saved: {output_path}')
    

    # Read the file
    with open(input_folder + boundary_file, 'r') as file:
        lines = file.readlines()
    # Extract numerical values from lines
    R_boundary = np.array([float(val) for val in lines[0].split('=')[1].split(',')])
    Z_boundary = np.array([float(val) for val in lines[1].split('=')[1].split(',')])

    column_names, data_flux = pltm.read_ascii_file(postproc_folder+flux_surface_file, 0, droprows = -1)
    R_values = data_flux[:, 0]
    Z_values = data_flux[:, 1]

    angles_rad = np.arctan2(Z_values - Z_axis, R_values - R_axis)
    sorted_indices = np.argsort(angles_rad)
    R_flux_sorted = R_values[sorted_indices]
    Z_flux_sorted = Z_values[sorted_indices]
    angles_rad_uniform = np.linspace(min(angles_rad)+1e-15, max(angles_rad)-1e-15, len(R_boundary))
    Rflux_target = interp1d(angles_rad[sorted_indices], R_flux_sorted, kind='linear', fill_value="extrapolate")
    R_flux_target = Rflux_target(angles_rad_uniform)
    Zflux_target = interp1d(angles_rad[sorted_indices], Z_flux_sorted, kind='linear', fill_value="extrapolate")
    Z_flux_target = Zflux_target(angles_rad_uniform)

    column_names, data_averaged = pltm.read_ascii_file(postproc_folder + averaged_file, 1, droprows = -1)
    psin_averaged = data_averaged[:, 0] 
    psi_averaged = data_averaged[:, 38]
    psi_truncate_func = interp1d(psin_averaged, psi_averaged)
    psi_truncate = psi_truncate_func(psin_truncate)
    psiflux_target = psi_truncate + Z_flux_target * 0.e0

    column_names, data_intg3D = pltm.read_ascii_file(postproc_folder + intg3D_file, 1, droprows = -1)
    psin_data_intg3D = data_intg3D[:, 0] 
    Ip_data_intg3D = data_intg3D[:, 3] / fmu0
    Ip_truncate_func = interp1d(psin_data_intg3D, Ip_data_intg3D)
    Ip_truncate = Ip_truncate_func(psin_truncate)
    Ip_truncate = Ip_truncate * Iptot_target / Iptot_equilibrium

    # Write to file
    
    with open(output_folder + boundary_truncated_file, "w") as file:
        file.write("R_boundary = ")
        file.write(", ".join(map(str, R_flux_target)) + "\n")
        file.write("Z_boundary = ")
        file.write(", ".join(map(str, Z_flux_target)) + "\n")
        file.write("Psi_boundary = ")
        file.write(", ".join(map(str, psiflux_target)) + "\n")

    print(f"Truncated and interpolated file saved: {output_folder + filename}")
    print(f"Truncated and interpolated value of psi at the boundary: {psi_truncate}")
    print(f"Truncated and interpolated value of Ip at the boundary: {Ip_truncate}")
    print(f"Truncated and interpolated value of loop_voltage at the boundary: {Ip_truncate/Iptot_target*loop_voltage}")


    plt.plot(R_boundary, Z_boundary)
    plt.plot(R_flux_target, Z_flux_target)
    plt.xlabel(r'R [m]')
    plt.ylabel(r'Z [m]')
    plt.grid()
    plt.axis('equal')
    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    # psin_truncate = float(input("Enter the truncation value (psin_truncate): "))
    # num_points = int(input("Enter the number of points for interpolation: "))
    truncate_and_interpolate_profiles(psin_truncate)