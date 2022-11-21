# --------------------------------------------------------------- #
# Read and plot the camera pixel and filter intensities HDF5 file
# --------------------------------------------------------------- #

import h5py
import numpy as np
from matplotlib import pyplot as plt

# Data ---------------------------------------------------------- #
filepath = "."
filename = "pixel_filter_intensities.h5"
datasetname = "pixel_filter_intensities"
separator = "/"

# Initialisation ------------------------------------------------ #
# Read image
fhandler = h5py.File("".join([filepath,separator,filename]),'r')
pixel_filter_intensities = np.array(fhandler[datasetname])
fhandler.close()

# Reorder pixel and filter intensities for plotting
pixel_filter_intensities = np.transpose(pixel_filter_intensities,axes=[3,4,0,1,2])
pixel_intensities = pixel_filter_intensities[0]
filter_intensities = pixel_filter_intensities[1]

# Plot pixel and filter intensities ----------------------------- #
for frame_id,frame in enumerate(pixel_intensities):
  for spectrum_id,spectrum in enumerate(frame):
    plt.figure()
    plt.imshow(spectrum)
    plt.colorbar()
    plt.title("".join(['Image for spectrum N# ',str(spectrum_id+1),\
    ' time N#: ',str(frame_id+1)]))

for img_filter_id,img_filter in enumerate(filter_intensities):
  for spectrum_id,spectrum in enumerate(img_filter):
    plt.figure()
    plt.imshow(spectrum)
    plt.colorbar()
    plt.title("".join(['Filter for spectrum N# ',str(spectrum_id+1),\
    ' time N#: ',str(img_filter_id+1)]))

plt.show()

# --------------------------------------------------------------- #
