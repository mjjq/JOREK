#! /usr/bin/env python3

# Try to plot helical magnetic flux surface without field line tracing
# Now it only works for reduced MHD model, as the magnetic field is mainly representded by psi.
# For full MHD, AR and AZ are also influcing the magnetic field, causing some problems.
# Ref: Basic equations for the description of a magnetic island [IslandEquations.pdf], by E. Poli
# Reduced MHD case: /u/hazh/WORKSPACE/tutorial/intear_ntor3_merge_case
# Full MHD case: /u/hazh/SCRATCH/pitagora/2025_0610_scan_visco_eta_Hartmann_Prandtl_correct_heating_rerun_pitagora/00_etax10_viscox10_Hx0.1_Px1_add_particle_source

# example for the reduced MHD case:
# step 0: set_more_outfile_digits_false
# step 1: convert2vtk.sh -si -no0 -only 50 -j 8 ./jorek2vtk intear  # to get the perturbed psi
# step 2: convert2vtk.sh -si -only 50 -j 8 ./jorek2vtk_ntor1 intear # to get the equilibrium psi
# step 3: jorek_post_all 00050                                      # enable si-units, and get the q profile and psi profile.
# step 4: change step_str and q_s in this python script.
# step 5: run with: plot_flux_surface_single_helicity.py
# step 6: compare with: plot_poincare_all.py poincares/poinc_R-Z_s00050.dat

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from scipy.ndimage import gaussian_filter1d
import scipy.interpolate as sci_intp
from scipy.interpolate import interp1d
import vtk
import vtkmodules.util.numpy_support as vtk_nump

def display_values(event, x, y, z):
    if event.inaxes:
        x_data, y_data = event.xdata, event.ydata
        # You can choose a different interpolation method if needed
        z_value = sci_intp.griddata(
            (x, y), z, (x_data, y_data), method='cubic')
        print(
            f"Mouse at (x, y): ({x_data:.2e}, {y_data:.2e}), z: {z_value:.2e}")

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

def plot2D_interpolated(x, y, z, time, clim = [], contour_lines = False):

    # plt.figure(figsize=(9, 7), constrained_layout=True)
    plt.figure(figsize=(5.5, 6), constrained_layout=False)
    # Define a regular grid for the contour plot
    x_grid = np.linspace(min(x), max(x), 400)
    y_grid = np.linspace(min(y), max(y), 400)

    X, Y = np.meshgrid(x_grid, y_grid)

    # Interpolate the z values onto the structured grid (optional)

    # You can choose a different interpolation method if needed
    Z = sci_intp.griddata((x, y), z, (X, Y), method='cubic')

    # Create a contour plot
    # You can choose a different colormap if desired
    if clim:
        cf = plt.contourf(X, Y, Z, 50, levels = np.linspace(clim[0], clim[1], 50), cmap='rainbow')#, extend = 'both')
    else:
        cf = plt.contourf(X, Y, Z, 50, cmap='rainbow')
        
    cs = plt.contour(X, Y, Z, 200, colors='black', linewidths=0.6)
    # Optional: add labels to contour lines
    # plt.clabel(cs, inline=True, fontsize=8)

    # Add a colorbar for reference
    plt.colorbar(cf)

    # Optional: You can add labels to the axes
    plt.xlabel('R')
    plt.ylabel('Z')

    time_str = "{:.2f}".format(time)
    plt.title("Time = " + time_str)

    plt.xlim([min(x), max(x)])
    plt.ylim([min(y), max(y)])
    plt.gca().set_aspect('equal',adjustable='box')

    # plt.axis('equal')
    # plt.gca().set_aspect('equal',adjustable='box')

    # Connect the event handler
    plt.connect('motion_notify_event',
                lambda event: display_values(event, x, y, z))

    # Show the plot
    # plt.show()
    # plt.xticks(fontsize=font_size)
    # plt.yticks(fontsize=font_size)
    # plt.legend(fontsize=font_size)
    plt.tight_layout()

###################################################
#                                                 #
#        calculate psi_t from q and psi_p         #
#                                                 #
###################################################

step_str = '00050'

# integrate q-psi_p to get psi_t
infile_psi = 'postproc/exprs_averaged_s' + step_str + '.dat'
infile_q = 'postproc/qprofile_s' + step_str + '.dat'
q_s = -2/1 # m/n

column_names, data = read_ascii_file(infile_psi, 1)

if column_names[0] == 'rho_p':
    rho_p = data[:, 0]
    psi_n = rho_p ** 2
else:
    psi_n = data[:, 0]
    rho_p = np.sqrt(psi_n)

psi = data[:, column_names.index('Psi')]

Psi_axis = psi[0]
Psi_bnd  = psi[-1]
column_names, data = read_ascii_file(infile_q, 0)
psi_n_tmp = data[:, 0]
q_tmp     = data[:, 1]
qpsi = np.interp(psi_n, psi_n_tmp, q_tmp)
psi_t_mag2bry = cumtrapz(qpsi, psi, initial=0)
q_check = np.gradient(psi_t_mag2bry, psi)

plt.plot(rho_p, qpsi, label='q_new')
plt.plot(np.sqrt(psi_n_tmp), q_tmp, label='q_old')
plt.plot(rho_p, q_check, label='q_check')

plt.legend()
plt.grid()
plt.xlabel(r'$\rho_p$')
plt.ylabel('q')

###################################################
#                                                 #
#        load psi_n=0 and >0 from vtks            #
#                                                 #
###################################################

infile_eq_vtk = 'vtk_iplane1_si/jorek.'+ step_str + '.vtk'
infile_pt_vtk = 'vtk_no0_iplane1_si/jorek.'+ step_str + '.vtk'

reader_eq = vtk.vtkDataSetReader()
reader_eq.SetFileName(infile_eq_vtk)
reader_eq.ReadAllScalarsOn()
reader_eq.Update()
dataset_eq = reader_eq.GetOutput()

reader_pt = vtk.vtkDataSetReader()
reader_pt.SetFileName(infile_pt_vtk)
reader_pt.ReadAllScalarsOn()
reader_pt.Update()
dataset_pt = reader_pt.GetOutput()

scalar_names = [reader_eq.GetScalarsNameInFile(
    i) for i in range(0, reader_eq.GetNumberOfScalarsInFile())]

time = dataset_eq.GetFieldData().GetArray("TIME").GetValue(0)

var = 'Psi'
if var not in scalar_names:
    print("Error input variable name: " + var + ", retry!!!")
    exit(-1)
data_array_eq = dataset_eq.GetPointData().GetArray(var)  # Change the index if needed
data_array_pt = dataset_pt.GetPointData().GetArray(var)  # Change the index if needed

Point_cordinates = reader_eq.GetOutput().GetPoints().GetData()
# Convert the VTK data array to a NumPy array
psi_p_eq = vtk_nump.vtk_to_numpy(data_array_eq)
psi_p_pt = vtk_nump.vtk_to_numpy(data_array_pt)
psi_p_tot = psi_p_eq + psi_p_pt
numpy_coordinates = vtk_nump.vtk_to_numpy(Point_cordinates)

plot2D_interpolated(numpy_coordinates[:, 0], numpy_coordinates[:, 1], psi_p_tot, time, [], False)

###################################################
#                                                 #
#        interpolate psi_t to 2D                  #
#                                                 #
###################################################

# psi_tor = np.interp(psi_p_eq, psi, psi_t_mag2bry)

f_interp = interp1d(psi, psi_t_mag2bry, kind='linear', fill_value="extrapolate")
psi_tor = f_interp(psi_p_eq)

# print(max(psi_p_eq))
# print(min(psi_p_eq))

# print(max(psi))
# print(min(psi))

# print(psi_p_eq.shape)
# print(psi.shape)
# print(psi_t_mag2bry.shape)
plt.figure()
plt.plot(psi, psi_t_mag2bry)
plt.plot(psi_p_eq, psi_tor)

psi_he_tot = psi_p_eq + psi_p_pt - psi_tor / q_s

plot2D_interpolated(numpy_coordinates[:, 0], numpy_coordinates[:, 1], psi_he_tot, time, [], True)

plt.show()
