# --------------------------------------------------------------- #
# Read and plot the camera pixel and filter intensities HDF5 file
# --------------------------------------------------------------- #
# Method used for plotting 4d arrays as 2d images
# inputs:
#   frames_spectra: (nx,ny,n_spectra,n_times) array to plot
#   x_positions:    (nx,n_times)(optional) x positions for scaling
#                   default: no scaling is used
#   y_positions:    (ny,n_times)(optional) y positions for scaling
#                   default: no scaling is used
#   title:          (string) figure title, default: empty
def imshow_4d(frames_spectra,x_positions=[],y_positions=[],title=""):
  from matplotlib import pyplot as plt
  for frame_id,frame in enumerate(frames_spectra):
    if(not ((len(x_positions)==0) and (len(y_positions)==0))):
      dx = 0.5*(x_positions[1,frame_id]-x_positions[0,frame_id])
      dy = 0.5*(y_positions[1,frame_id]-y_positions[0,frame_id])
      ext = [x_positions[0,frame_id]-dx, x_positions[-1,frame_id]+dx,\
      y_positions[0,frame_id]-dy,y_positions[-1,frame_id]+dy]
    else:
      ext = None
    for spectrum_id,spectrum in enumerate(frame):
      plt.figure()
      plt.imshow(spectrum,extent=ext)
      plt.colorbar()
      plt.title("".join([title,' spectrum N# ',str(spectrum_id+1),\
      ' time N#: ',str(frame_id+1)]))   
  return  

# main function
def load_and_plot_images(filename="",filepath=".",image_datasetname="",\
pixel_x_postion_datasetname="",pixel_y_postion_datasetname="",\
separator="/"):
  import h5py
  import numpy as np
  from matplotlib import pyplot as plt
  # Initialisation
  # Read image
  fhandler = h5py.File("".join([filepath,separator,filename]),'r')
  pixel_filter_intensities = np.array(fhandler[image_datasetname])
  x_positions = np.array(fhandler[pixel_x_postion_datasetname])
  y_positions = np.array(fhandler[pixel_y_postion_datasetname])
  fhandler.close()
  # Reorder pixel and filter intensities for plotting
  x_positions = np.transpose(x_positions,axes=[1,0])
  y_positions = np.transpose(y_positions,axes=[1,0])
  pixel_filter_intensities = np.transpose(pixel_filter_intensities,axes=[3,4,0,1,2])
  pixel_intensities = pixel_filter_intensities[0]
  filter_intensities = pixel_filter_intensities[1]
  # Plot pixel and filter intensities
  imshow_4d(pixel_intensities,x_positions=x_positions,y_positions=y_positions,title="Image for" )
  imshow_4d(filter_intensities,x_positions=x_positions,y_positions=y_positions,title="Filter for" )
  plt.show()
  return

# Run main ------------------------------------------------------ #
if __name__ == "__main__":
  load_and_plot_images(filename="pixel_filter_intensities.h5",\
  filepath=".",image_datasetname="pixel_filter_intensities",\
  pixel_x_postion_datasetname="x_pixel_coordinates",\
  pixel_y_postion_datasetname="y_pixel_coordinates")

# --------------------------------------------------------------- #
