#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci_intp

ivar_plot = 0 # 0 => q
varNameList = ['q', 'density (m-3)', 'Ti (eV)', 'Te (eV)']
file_path = 'q_profile_3D_ip_1.dat'

def display_values(event, x, y, z):
    if event.inaxes:
        x_data, y_data = event.xdata, event.ydata
        # You can choose a different interpolation method if needed
        z_value = sci_intp.griddata(
            (x, y), z, (x_data, y_data), method='cubic')
        print(
            f"Mouse at (x, y): ({x_data:.2e}, {y_data:.2e}), z: {z_value:.2e}")

if __name__ == "__main__":

    data = np.loadtxt(file_path, skiprows=2)
    # Extract columns: second column for x, third column for y, and fourth+ column for amplitude
    x = data[:, 2]
    y = data[:, 3]
    amplitude = np.abs(data[:, 4 + ivar_plot])

    # Filter out invalid (NaN or infinite) values
    mask = np.isfinite(amplitude)
    x = x[mask]
    y = y[mask]
    amplitude = amplitude[mask]

    mask = x > 0
    x = x[mask]
    y = y[mask]
    amplitude = amplitude[mask]

    if ivar_plot == 0:
        # Filter out data points with amplitude greater than 7
        mask = amplitude <= 10
        x = x[mask]
        y = y[mask]
        amplitude = amplitude[mask]

    # Create the scatter plot
    plt.figure(figsize=(6, 6))
    sc = plt.scatter(x, y, c=amplitude, alpha=1.0, s=0.1, cmap='gist_ncar')#, vmin = 0.8, vmax = 1.05)
    plt.colorbar()
    plt.xlabel('R')
    plt.ylabel('Z')
    plt.title('Scatter plot with ' + varNameList[ivar_plot])
    plt.grid(True)
    plt.axis('equal')

    # Connect the event handler
    plt.connect('motion_notify_event',
                lambda event: display_values(event, x, y, amplitude))
    

    plt.figure(figsize=(6, 6), constrained_layout=True)
    # Define a regular grid for the contour plot
    x_grid = np.linspace(min(x), max(x), 400)
    y_grid = np.linspace(min(y), max(y), 400)

    X, Y = np.meshgrid(x_grid, y_grid)
    Z = sci_intp.griddata((x, y), amplitude, (X, Y), method='linear')

    plt.contourf(X, Y, Z, 20, cmap='viridis_r')
    plt.colorbar()

    if ivar_plot == 0:
        contour_line = plt.contour(X, Y, Z, levels=[0.80, 0.90, 0.99, 1.01, 2.0, 3.0], colors='k', linestyles = 'dotted')
        plt.clabel(contour_line, inline=True, fontsize=8)

    plt.xlabel('R')
    plt.ylabel('Z')
    plt.title('Contourf plot with ' + varNameList[ivar_plot])
    plt.grid(True)
    plt.axis('equal')

    # Connect the event handler
    plt.connect('motion_notify_event',
                lambda event: display_values(event, x, y, amplitude))

    plt.show()
