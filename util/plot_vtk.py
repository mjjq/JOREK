#! /usr/bin/env python3

#
# Purpose: Plot 2D variables from .vtk files.
#
# Date: 2023-09-20
# Author: Haowei Zhang, IPP Garching
#

import vtk
import vtkmodules.util.numpy_support as vtk_nump
import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.interpolate as sci_intp
import mpl_toolkits.mplot3d as plt3d
# import matplotlib.tri as mtri

def str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in {'false', 'f', '0', 'no', 'n'}:
        return False
    elif value.lower() in {'true', 't', '1', 'yes', 'y'}:
        return True
    raise ValueError(f'{value} is not a valid boolean value')

# Function to display (x, y, z) values at the mouse position
def display_values(event, x, y, z):
    if event.inaxes:
        x_data, y_data = event.xdata, event.ydata
        z_value = sci_intp.griddata((x, y), z, (x_data, y_data), method='cubic')  # You can choose a different interpolation method if needed
        print(f"Mouse at (x, y): ({x_data:.2e}, {y_data:.2e}), z: {z_value:.2e}")

def plot3D_unstructed(x, y, z, title, file_path, time):
    fig = plt.figure(figsize=(9, 7), constrained_layout=True)
    # ax = fig.add_subplot(1,1,1)
    ax = fig.add_subplot(111, projection='3d')
    # Add a colorbar outside the figure
    gcf = ax.plot_trisurf(x, y, z, cmap='jet') # colormap: jet, viridis, etc.
    cbar = plt.colorbar(gcf, pad=0.01, fraction=0.2, location='left', shrink=0.5)  # Adjust pad and fraction as needed

    # Optional: You can add labels to the axes
    ax.set_xlabel('R')
    ax.set_ylabel('Z')
    ax.set_title(title + " @ Time = " + "{:.2f}".format(time) + " of " + file_path)

    # # Remove the Z axis
    # ax.zaxis.line.set_lw(0.)
    # ax.set_zticks([])

    # # Set the background color to be transparent
    # ax.xaxis.pane.fill = False
    # ax.yaxis.pane.fill = False
    # ax.zaxis.pane.fill = False

    # # make the panes transparent
    # ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # # Remove the grid lines for X, Y, and Z axes
    # ax.grid(False)

    # # make the grid lines transparent
    # ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    # ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    # ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

    # # Set the axis limits for X and Y based on their min and max values
    # ax.set_xlim(min(x), max(x))
    # ax.set_ylim(min(y), max(y))

    # # Define the coordinates for the bounding box
    # bbox = [
    #     [min(x), min(y), 0],
    #     [max(x), min(y), 0],
    #     [max(x), max(y), 0],
    #     [min(x), max(y), 0],
    #     [min(x), min(y), 0],  # Close the loop
    # ]
    #
    # # Extract the x, y, and z coordinates
    # bbox_x = [coord[0] for coord in bbox]
    # bbox_y = [coord[1] for coord in bbox]
    # bbox_z = [coord[2] for coord in bbox]
    #
    # # Plot the bounding box
    # ax.plot(bbox_x, bbox_y, bbox_z, color='black')

    # Optional: Customize the viewing angle
    ax.view_init(elev=45, azim=-45)

    # Show the plot
    plt.show()

def plot2D_scatter(x, y, z, title, file_path, time):
    
    plt.figure(figsize=(9, 7), constrained_layout=True)

    # Create a scatter plot
    plt.scatter(x, y, c=z, cmap='jet', marker='o', s=50)  # You can choose a different colormap and marker style

    # Add a colorbar for reference
    plt.colorbar()

    # Optional: You can add labels to the axes
    plt.xlabel('R')
    plt.ylabel('Z')

    plt.title(title + " @ Time = " + "{:.2f}".format(time) + " (" + file_path + ")")

    # plt.xlim([min(x), max(x)])
    # plt.ylim([min(y), max(y)])

    plt.axis('equal')
    # plt.gca().set_aspect('equal',adjustable='box')

    # Connect the event handler
    plt.connect('motion_notify_event', lambda event: display_values(event, x, y, z))

    # Show the plot
    plt.show()    


def plot2D_flow_field(x, y, vx, vy, title, file_path, time):
    
    plt.figure(figsize=(9, 7), constrained_layout=True)
    # Define a regular grid for the contour plot
    x_grid = np.linspace(min(x), max(x), 40)
    y_grid = np.linspace(min(y), max(y), 40)

    X, Y = np.meshgrid(x_grid, y_grid)

    # Interpolate the z values onto the structured grid (optional)
    
    V_X = sci_intp.griddata((x, y), vx, (X, Y), method='cubic')  # You can choose a different interpolation method if needed
    V_Y = sci_intp.griddata((x, y), vy, (X, Y), method='cubic')  # You can choose a different interpolation method if needed
    V = pow(V_X * V_X + V_Y * V_Y, 0.5)

    # Create a quiver plot
    plt.quiver(X, Y, V_X, V_Y, V, cmap='hot')  # You can choose a different colormap if desired

    # Add a colorbar for reference
    plt.colorbar()

    # Optional: You can add labels to the axes
    plt.xlabel('R')
    plt.ylabel('Z')

    plt.title(title + " @ Time = " + "{:.2f}".format(time) + " (" + file_path + ")")

    # plt.xlim([min(x), max(x)])
    # plt.ylim([min(y), max(y)])

    plt.axis('equal')
    # plt.gca().set_aspect('equal',adjustable='box')

    # Show the plot
    plt.show()


def plot2D_interpolated(x, y, z, title, file_path, time):
    
    plt.figure(figsize=(9, 7), constrained_layout=True)
    # Define a regular grid for the contour plot
    x_grid = np.linspace(min(x), max(x), 400)
    y_grid = np.linspace(min(y), max(y), 400)

    X, Y = np.meshgrid(x_grid, y_grid)

    # Interpolate the z values onto the structured grid (optional)
    
    Z = sci_intp.griddata((x, y), z, (X, Y), method='cubic')  # You can choose a different interpolation method if needed

    # Create a contour plot
    plt.contourf(X, Y, Z, 50, cmap='jet')  # You can choose a different colormap if desired

    # Add a colorbar for reference
    plt.colorbar()

    # Optional: You can add labels to the axes
    plt.xlabel('R')
    plt.ylabel('Z')

    plt.title(title + " @ Time = " + "{:.2f}".format(time) + " (" + file_path + ")")

    # plt.xlim([min(x), max(x)])
    # plt.ylim([min(y), max(y)])

    plt.axis('equal')
    # plt.gca().set_aspect('equal',adjustable='box')

    # Connect the event handler
    plt.connect('motion_notify_event', lambda event: display_values(event, x, y, z))

    # Show the plot
    plt.show()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot surfaces from vtk files.")
    parser.add_argument("file_paths", nargs="+", help="List of file paths to the vtk files. To plot all files, please replace the file name like jorek.01000.vtk with jorek.*.vtk")
    parser.add_argument("-va", "--variable", default="rho", help="Plot the specified variable, such as [Psi, AR, AZ, UR, UZ, Uphi, rho, T, psi_norm, B_R, B_Z, B_phi, J_R, J_Z, J_phi, FFprime_equi, F_prof_equi, Grad_P, JxB, V_parallel, etc.].")
    parser.add_argument('-me', "--method", choices=["it", "sc", "sf", "ff"], default="it: interpolate",
                        help="Use the specified method to plot the data, including 1. it: interpolate; 2. sc: scatter; 3. sf: surface; 4. ff: flow field (arrows)")
    # parser.add_argument('-it', "--interpolate", type=str_to_bool, nargs='?', choices=[True, False], default=False,
    #                     help="Interpolate unstructed data in to structed grids? Default is True.")
    
    args = parser.parse_args()

    # print(args.file_paths)

    for ifile, file_path in enumerate(args.file_paths):
        
        # file_path = "vtk_no0_iplane1/jorek.00000.vtk"
        # Create a VTK reader
        reader = vtk.vtkDataSetReader()
        reader.SetFileName(file_path)
        reader.ReadAllScalarsOn()
        reader.Update()

        # Get the VTK dataset
        dataset = reader.GetOutput()

        scalar_names = [reader.GetScalarsNameInFile(i) for i in range(0, reader.GetNumberOfScalarsInFile())]

        print('Variable names:\n', scalar_names)

        # Read the VTK file
        reader.Update()

        # Access the data
        dataset = reader.GetOutput()

        # time = reader.GetOutput().GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())
        # cycle = reader.GetOutput().GetInformation().Get(vtk.vtkDataObject.DATA_TIME_INDEX())
        time  = dataset.GetFieldData().GetArray("TIME").GetValue(0)
        # print('time: ', time)

        # Print dataset information
        # print("Dataset Information:")
        # print("Number of Points:", dataset.GetNumberOfPoints())
        # print("Number of Cells:", dataset.GetNumberOfCells())
        # if ifile == 0:
        #     i = 0
        #     while dataset.GetPointData().GetArrayName(i) != None:
        #         print("Point Data Arrays [", i, "]: ", dataset.GetPointData().GetArrayName(i))
        #         # print("Cell Data Arrays:", dataset.GetCellData().GetArrayName(i))
        #         i += 1

        # Get the 2D data array
        # var = "T"
        var = args.variable

        if var not in scalar_names:
            print("Error input variable name: " + var + ", retry!!!")
            exit(-1)

        data_array = dataset.GetPointData().GetArray(var)  # Change the index if needed

        Point_cordinates = reader.GetOutput().GetPoints().GetData()
        # Convert the VTK data array to a NumPy array
        data = vtk_nump.vtk_to_numpy(data_array)
        numpy_coordinates = vtk_nump.vtk_to_numpy(Point_cordinates)

        print("======== " + str(ifile) + " Plotting " + file_path + " ========")
        if args.method == "sc":
            plot2D_scatter(numpy_coordinates[:,0], numpy_coordinates[:,1], data, var, file_path, time)
        elif args.method == "sf":
            plot3D_unstructed(numpy_coordinates[:,0], numpy_coordinates[:,1], data, var, file_path, time)
        elif args.method == "ff":
            UR_array   = dataset.GetPointData().GetArray("UR")
            UZ_array   = dataset.GetPointData().GetArray("UZ")
            # Uphi_array = dataset.GetPointData().GetArray("Uphi")
            UR   = vtk_nump.vtk_to_numpy(UR_array)
            UZ   = vtk_nump.vtk_to_numpy(UZ_array)
            # Uphi = vtk_nump.vtk_to_numpy(Uphi_array)
            plot2D_flow_field(numpy_coordinates[:,0], numpy_coordinates[:,1], UR, UZ, "2D Flow field", file_path, time)
        else:
            plot2D_interpolated(numpy_coordinates[:,0], numpy_coordinates[:,1], data, var, file_path, time)

        ifile += 1