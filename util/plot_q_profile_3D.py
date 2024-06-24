#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci_intp


def display_values(event, x, y, z):
    if event.inaxes:
        x_data, y_data = event.xdata, event.ydata
        # You can choose a different interpolation method if needed
        z_value = sci_intp.griddata(
            (x, y), z, (x_data, y_data), method='cubic')
        print(
            f"Mouse at (x, y): ({x_data:.2e}, {y_data:.2e}), z: {z_value:.2e}")

if __name__ == "__main__":
    # Load the data
    file_path = 'q_profile_3D.dat'
    data = np.loadtxt(file_path, skiprows=1)

    # Extract columns: second column for x, third column for y, and fourth column for amplitude
    x = data[:, 1]
    y = data[:, 2]
    amplitude = np.abs(data[:, 3])
    # amplitude = data[:, 3]

    # Filter out invalid (NaN or infinite) values
    mask = np.isfinite(amplitude)
    x = x[mask]
    y = y[mask]
    amplitude = amplitude[mask]

    mask = amplitude > 0
    x = x[mask]
    y = y[mask]
    amplitude = amplitude[mask]

    # Filter out data points with amplitude greater than 7
    mask = amplitude <= 10
    x = x[mask]
    y = y[mask]
    amplitude = amplitude[mask]

    # Create the scatter plot
    plt.figure(figsize=(6, 6))
    sc = plt.scatter(x, y, c=amplitude, alpha=1.0, s=0.1, cmap='gist_ncar')
    plt.colorbar(sc, label='q')
    plt.xlabel('R')
    plt.ylabel('Z')
    plt.title('Scatter plot with q')
    plt.grid(True)
    plt.axis('equal')
    # Connect the event handler
    plt.connect('motion_notify_event',
                lambda event: display_values(event, x, y, amplitude))
    # plt.show()



    plt.figure(figsize=(6, 6), constrained_layout=True)
    # Define a regular grid for the contour plot
    x_grid = np.linspace(min(x), max(x), 400)
    y_grid = np.linspace(min(y), max(y), 400)

    X, Y = np.meshgrid(x_grid, y_grid)

    # Interpolate the z values onto the structured grid (optional)
    # You can choose a different interpolation method if needed
    Z = sci_intp.griddata((x, y), amplitude, (X, Y), method='linear')

    # Create a contour plot
    # You can choose a different colormap if desired
    plt.contourf(X, Y, Z, 20, cmap='viridis_r')
    # Add a colorbar for reference
    plt.colorbar()

    contour_line = plt.contour(X, Y, Z, levels=[0.80, 0.90, 0.99, 1.01, 2.0, 3.0], colors='k', linestyles = 'dotted')
    plt.clabel(contour_line, inline=True, fontsize=8)


    # Optional: You can add labels to the axes
    plt.xlabel('R')
    plt.ylabel('Z')
    plt.axis('equal')
    # plt.gca().set_aspect('equal',adjustable='box')

    # Connect the event handler
    plt.connect('motion_notify_event',
                lambda event: display_values(event, x, y, amplitude))

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # plt.legend(fontsize=12)

    plt.show()