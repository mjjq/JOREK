# ----------------------------------------------------- #
# Program for reading and plots contributing light      #
# light sources to each point on a lens toghether with  #
# the image plane and their viewing directions          #
# ----------------------------------------------------- #
# Program functions ----------------------------------- #
# Read datasets from HDF5 files
def comput_major_radius(x,y):
  from numpy import power,sqrt
  return sqrt(power(x,2)+power(y,2))

def read_datasets_from_hdf5(filenames,filepath,datasets,separator):
  import h5py
  import numpy as np
  data_list = []
  for filename in filenames:
    data = []
    fhandler = h5py.File("".join([filepath,separator,filename]))
    for setname in datasets:
      data.append(np.array(fhandler[setname]))
    data_list.append(data)  
    fhandler.close()
  return data_list

# Check the consistency of the light data:
# the number of spectra is the same for
# all datasets. Structure of the light data 
# for each light data:
# 0: position of the contributing light sources
# 1: spectral intensities of the contributing light sources
# 2: tokamak limiter (first wall) major radius
# 3: tokamak limiter (first wall) vertical position
def check_equal_light_data(light_data):
  from sys import exit
  num_spectra = light_data[0][1].shape[0]
  test_spectra = [dataset[1].shape[0]==num_spectra for dataset in light_data]
  if(not all(test_spectra)):
    exit('Number of spectra must be the same for all light datasets')

# Plot light and camera data
# structure of the light data for each light data
# 0: position of the contributing light sources
# 1: spectral intensities of the contributing light sources
# 2: tokamak limiter (first wall) major radius
# 3: tokamak limiter (first wall) vertical position
# structure of the camera data for each camera data:
# 0: points on the lens
# 1: vertices of the image planes
# 2: view directions of the image planes
def plot_light_and_camera_data(light_data,camera_data,n_tor=100,markersize=1,\
linewidth=3,fontsize=16):
  import numpy as np
  from mpl_toolkits.mplot3d import Axes3D
  from matplotlib import pyplot as plt
  # extract the number of spectra
  n_spectra = light_data[0][1].shape[0]
  # extract limiter coordinates considered equals for all dataset
  R_limiter = light_data[0][2]
  Z_limiter = light_data[0][3]
  cosphi = np.cos(np.linspace(0,2*np.pi,num=n_tor,dtype=np.float64))
  sinphi = np.sin(np.linspace(0,2*np.pi,num=n_tor,dtype=np.float64))
  R_max_limiter = np.amax(R_limiter)
  R_min_limiter = np.amin(R_limiter)
  x_limiter_min = R_min_limiter*cosphi
  x_limiter_max = R_max_limiter*cosphi
  y_limiter_min = R_min_limiter*sinphi
  y_limiter_max = R_max_limiter*sinphi
  del cosphi,sinphi,R_max_limiter,R_min_limiter
  # find the maximum of the specta
  max_spectrum = 0
  for dataset in light_data:
    max_spectrum = max(max_spectrum,np.amax(dataset[1]))
  # do scatter plots
  for spectra_id in range(n_spectra):
    fig = plt.figure(facecolor='white',edgecolor='white')
    axs = [fig.add_subplot(131),fig.add_subplot(132),fig.add_subplot(133,projection='3d')]
    n_particles = 0
    for dataset in light_data:
      positions = dataset[0]
      spectra = dataset[1]
      n_particles = n_particles + positions.shape[1]
      for time_id,spectrum in enumerate(np.transpose(spectra[spectra_id,:,:],[1,0])):
        mask = np.where(spectrum>0)
        # aggregated top view
        axs[0].scatter(positions[time_id,mask,0],positions[time_id,mask,1],marker='.',\
        s=markersize,c=spectrum[mask],cmap='inferno',\
        vmin=0,vmax=max_spectrum)
        axs[0].set_aspect('equal')
        axs[0].set_facecolor([0,0,0])
        axs[0].set_title('Aggregated top view',fontsize=fontsize,color='red')
        axs[0].set_xlabel('x [m]',fontsize=fontsize,color='red')
        axs[0].set_ylabel('y [m]',fontsize=fontsize,color='red')
        axs[0].tick_params(axis='x',labelsize=fontsize,colors='red')
        axs[0].tick_params(axis='y',labelsize=fontsize,colors='red')
        # aggregated front view
        axs[1].scatter(comput_major_radius(positions[time_id,mask,0],positions[time_id,mask,1]),\
        positions[time_id,mask,2],s=markersize,c=spectrum[mask],marker='.',cmap='inferno',\
        vmin=0,vmax=max_spectrum)
        axs[1].set_aspect('equal')
        axs[1].set_facecolor([0,0,0])
        axs[1].set_title('Aggregated frontal view',fontsize=fontsize,color='red')
        axs[1].set_xlabel('R [m]',fontsize=fontsize,color='red')
        axs[1].set_ylabel('Z [m]',fontsize=fontsize,color='red')
        axs[1].tick_params(axis='x',labelsize=fontsize,colors='red')
        axs[1].tick_params(axis='y',labelsize=fontsize,colors='red')
        # 3D view
        axs[2].scatter3D(positions[time_id,mask,0],positions[time_id,mask,1],positions[time_id,mask,2],\
        s=markersize,c=spectrum[mask],marker='.',cmap='inferno',vmin=0,vmax=max_spectrum)
        axs[2].set_facecolor([0,0,0])
        axs[2].set_aspect('equal')
        axs[2].set_title('3D view',fontsize=fontsize,color='red')
        axs[2].set_xlabel('x [m]',fontsize=fontsize,color='red')
        axs[2].set_ylabel('y [m]',fontsize=fontsize,color='red')
        axs[2].set_zlabel('z [m]',fontsize=fontsize,color='red')
        axs[2].tick_params(axis='x',labelsize=fontsize,colors='red')
        axs[2].tick_params(axis='y',labelsize=fontsize,colors='red')
        axs[2].tick_params(axis='z',labelsize=fontsize,colors='red')
        axs[2].xaxis.set_rotate_label(False)
        axs[2].yaxis.set_rotate_label(False)
        axs[2].zaxis.set_rotate_label(False)
  # Plotting limiter
  axs[0].plot(x_limiter_min,y_limiter_min,color='red',linewidth=linewidth)
  axs[0].plot(x_limiter_max,y_limiter_max,color='red',linewidth=linewidth)
  axs[1].plot(R_limiter,Z_limiter,color='red',linewidth=linewidth)
  del x_limiter_min,y_limiter_min,x_limiter_max,y_limiter_max,R_limiter,Z_limiter
  # Plotting the camera properties
  for dataset in camera_data:
    # Plot points on lens
    lens_point_position_avg = []
    for points in dataset[0]:
      lens_point_position_avg.append(np.array([np.mean(points[:,0]),np.mean(points[:,1]),np.mean(points[:,2])]))
      axs[0].scatter(points[:,0],points[:,1],marker='.',s=markersize,c='green')
      axs[2].scatter3D(points[:,0],points[:,1],points[:,2],marker='.',s=markersize,c='green')
    # Plot image planes
    for plane in dataset[1]:
      point4 = plane[1,:]-plane[0,:]+plane[2,:]
      axs[0].plot([plane[0,0],plane[1,0]],[plane[0,1],plane[1,1]],color='green',linewidth=linewidth)
      axs[0].plot([plane[0,0],plane[2,0]],[plane[0,1],plane[2,1]],color='green',linewidth=linewidth)
      axs[0].plot([plane[1,0],point4[0]],[plane[1,1],point4[1]],color='green',linewidth=linewidth)
      axs[0].plot([plane[2,0],point4[0]],[plane[2,1],point4[1]],color='green',linewidth=linewidth)
      axs[2].plot3D([plane[0,0],plane[1,0]],[plane[0,1],plane[1,1]],[plane[0,2],plane[1,2]],\
      color='green',linewidth=linewidth)
      axs[2].plot3D([plane[0,0],plane[2,0]],[plane[0,1],plane[2,1]],[plane[0,2],plane[2,2]],\
      color='green',linewidth=linewidth)
      axs[2].plot3D([plane[1,0],point4[0]],[plane[1,1],point4[1]],[plane[1,2],point4[2]],\
      color='green',linewidth=linewidth)
      axs[2].plot3D([plane[2,0],point4[0]],[plane[2,1],point4[1]],[plane[2,2],point4[2]],\
      color='green',linewidth=linewidth)
    del point4
    # Plot view directions
    for time_id,direction in enumerate(dataset[2]):
      point = lens_point_position_avg[time_id]
      axs[0].quiver(point[0],point[1],direction[0],direction[1],color='green',\
      linewidth=linewidth)
      axs[2].quiver(point[0],point[1],point[2],direction[0],direction[1],direction[2],\
      color='green',linewidth=linewidth)
    del lens_point_position_avg
  # generating image
  plt.suptitle("".join(['Point light source intensities, spectrum N#:',str(spectra_id)]),\
  fontsize=fontsize)
  plt.show() 

# Main function
# The light source datasets are:
#   contributing_light_positions: position of the contributing
#     lights per lens point
#   contributing_light_intensities: spectral intensity of the
#     contributing lights per lens point
#   limiter_major_radius
#     tokamak limiter (first wall) major radius
#   limiter_vertical_coordinate
#     tokamak limiter (first wall) vertical position
# The camera datasets are:
#   point_on_lens_positions: positions of the points on 
#     the camera lens
#   image_plane_vertices: vertices of the image planes
#   image_plane_directions: viewing directions of each
#     image plane
def load_and_plot_contributing_light_sources(\
light_filenames=[],camera_filenames=[],filepath="",\
light_dataset=["contributing_light_positions","contributing_light_intensities",\
"limiter_major_radius","limiter_vertical_coordinate"],\
camera_dataset=["point_on_lens_positions","image_plane_vertices",\
"image_plane_directions"],separator="/"):
  # Read data from files
  light_data = read_datasets_from_hdf5(light_filenames,filepath,\
  light_dataset,separator)  
  camera_data = read_datasets_from_hdf5(camera_filenames,filepath,\
  camera_dataset,separator)
  # Check consistency of the data
  check_equal_light_data(light_data)
  # Plot data
  plot_light_and_camera_data(light_data,camera_data)

# Run main -------------------------------------------- #
if __name__ == "__main__":
  light_names = ["contributing_light_intesities_0.h5",\
  "contributing_light_intesities_1.h5"]#,"contributing_light_intesities_2.h5",\
  #"contributing_light_intesities_3.h5","contributing_light_intesities_4.h5",\
  #"contributing_light_intesities_5.h5","contributing_light_intesities_6.h5",\
  #"contributing_light_intesities_7.h5","contributing_light_intesities_8.h5",\
  #"contributing_light_intesities_9.h5","contributing_light_intesities_10.h5",\
  #"contributing_light_intesities_11.h5"] 
  camera_names = ["contributing_camera_intesities_0.h5",\
  "contributing_camera_intesities_1.h5","contributing_camera_intesities_2.h5",\
  "contributing_camera_intesities_3.h5","contributing_camera_intesities_4.h5",\
  "contributing_camera_intesities_5.h5","contributing_camera_intesities_6.h5",\
  "contributing_camera_intesities_7.h5","contributing_camera_intesities_8.h5",\
  "contributing_camera_intesities_9.h5","contributing_camera_intesities_10.h5",\
  "contributing_camera_intesities_11.h5"]
  path = "."
  fig=load_and_plot_contributing_light_sources(light_filenames=light_names,\
  camera_filenames=camera_names,filepath=path)
# ----------------------------------------------------- #

