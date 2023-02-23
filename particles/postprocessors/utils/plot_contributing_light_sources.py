# ----------------------------------------------------- #
# Program for reading and plots contributing light      #
# light sources to each point on a lens toghether with  #
# the image plane and their viewing directions          #
# ----------------------------------------------------- #
# Program functions ----------------------------------- #
# Read datasets from HDF5 files
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
# structure of the camera data for each camera data:
# 0: vertices of the image planes
# 1: view directions of the image planes
def plot_light_and_camera_data(light_data,camera_data,markersize=1,fontsize=16):
  import numpy as np
  from mpl_toolkits.mplot3d import Axes3D
  from matplotlib import pyplot as plt
  # extract the number of spectra
  n_spectra = light_data[0][1].shape[0]
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
        axs[0].scatter(positions[time_id,np.where(spectrum>0),0],\
        positions[time_id,np.where(spectrum>0),1],marker='.',s=markersize,\
        c=spectrum[spectrum>0],cmap='inferno',vmin=0,vmax=max_spectrum)
        axs[0].set_aspect('equal')
        axs[0].set_facecolor([0,0,0])
        axs[0].set_title('Aggregated top view',fontsize=fontsize,color='red')
        axs[0].set_xlabel('x [m]',fontsize=fontsize,color='red')
        axs[0].set_ylabel('y [m]',fontsize=fontsize,color='red')
        axs[0].tick_params(axis='x',labelsize=fontsize,colors='red')
        axs[0].tick_params(axis='y',labelsize=fontsize,colors='red')
        major_radius = np.sqrt(np.power(positions[time_id,np.where(spectrum>0),0],2)+\
        np.power(positions[time_id,np.where(spectrum>0),1],2))
        axs[1].scatter(major_radius,positions[time_id,np.where(spectrum>0),2],s=markersize,\
        c=spectrum[spectrum>0],marker='.',cmap='inferno',vmin=0,vmax=max_spectrum)
        axs[1].set_aspect('equal')
        axs[1].set_facecolor([0,0,0])
        axs[1].set_title('Aggregated frontal view',fontsize=fontsize,color='red')
        axs[1].set_xlabel('x [m]',fontsize=fontsize,color='red')
        axs[1].set_ylabel('y [m]',fontsize=fontsize,color='red')
        axs[1].tick_params(axis='x',labelsize=fontsize,colors='red')
        axs[1].tick_params(axis='y',labelsize=fontsize,colors='red')
        axs[2].scatter3D(positions[time_id,np.where(spectrum>0),0],\
        positions[time_id,np.where(spectrum>0),1],positions[time_id,np.where(spectrum>0),2],\
        s=markersize,c=spectrum[spectrum>0],marker='.',cmap='inferno',vmin=0,vmax=max_spectrum)
        axs[2].set_facecolor([0,0,0])
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

  plt.suptitle("".join(['Point light source intensities, spectrum N#:',str(spectra_id)]),\
  fontsize=fontsize)
  plt.show() 

# Main function
# The light source datasets are:
#   contributing_light_positions: position of the contributing
#     lights per lens point
#   contributing_light_intensities: spectral intensity of the
#     contributing lights per lens point
# The camera datasets are:
#   point_on_lens_positions: positions of the points on 
#     the camera lens
#   image_plane_vertices: vertices of the image planes
#   image_plane_directions: viewing directions of each
#     image plane
def load_and_plot_contributing_light_sources(\
light_filenames=[],camera_filenames=[],filepath="",\
light_dataset=["contributing_light_positions","contributing_light_intensities"],\
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

