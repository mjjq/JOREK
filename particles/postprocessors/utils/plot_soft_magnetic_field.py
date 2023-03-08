# Plot the SOFT-compatible JOREK MHD field
# main function
def load_and_plot_images(filename="",filepath=".",separator="/",\
fontsize=32,linewidth=3,markersize=3):
  import h5py
  import numpy as np
  from matplotlib import pyplot as plt
  # Initialisation
  separatrix = np.array([])
  wall = np.array([])
  # Read image
  fhandler = h5py.File("".join([filepath,separator,filename]),'r')
  major_radius = np.array(fhandler['r'])
  vertical_position = np.array(fhandler['z'])
  BR = np.transpose(np.array(fhandler['Br']))
  BZ = np.transpose(np.array(fhandler['Bz']))
  Bphi = np.transpose(np.array(fhandler['Bphi']))
  poloidal_flux = np.transpose(np.array(fhandler['Psi']))
  if('wall' in fhandler):
    wall = np.transpose(np.array(fhandler['wall']))
  if('separatrix' in fhandler):
    lcfs = np.transpose(np.array(fhandler['separatrix']))
  axis = np.array(fhandler['maxis'])
  fhandler.close()

  # Plot magnetic field
  xlabels = [str(np.around(value,decimals=2)) for value in major_radius]
  ylabels = [str(np.around(value,decimals=2)) for value in vertical_position]
  extension = [major_radius[0],major_radius[-1],vertical_position[0],vertical_position[-1]]
  fig = plt.figure(facecolor='white',edgecolor='white')
  ax = [fig.add_subplot(221),fig.add_subplot(222),fig.add_subplot(223),fig.add_subplot(224)]
  im0=ax[0].imshow(BR,extent=extension)
  ax[0].set_aspect('equal')
  ax[0].set_xlabel('R [m]',fontsize=fontsize)
  ax[0].set_ylabel('Z [m]',fontsize=fontsize)
  if(len(lcfs)!=0):
    ax[0].plot(lcfs[0,:],lcfs[1,:],color='r',linewidth=linewidth)
  if(len(wall)!=0):
    ax[0].plot(wall[0,:],wall[1,:],color='r',linewidth=linewidth)
  ax[0].plot(axis[0],axis[1],marker='o',markersize=markersize,markerfacecolor='r',markeredgecolor='r')
  ax[0].set_title('Radial magnetic field',fontsize=fontsize)
  ax[0].tick_params(axis='x',labelsize=fontsize)
  ax[0].tick_params(axis='y',labelsize=fontsize)
  fig.colorbar(im0,ax=ax[0])
  im1=ax[1].imshow(BZ,extent=extension)
  ax[1].set_aspect('equal')
  ax[1].set_xlabel('R [m]',fontsize=fontsize)
  ax[1].set_ylabel('Z [m]',fontsize=fontsize)
  if(len(wall)!=0):
    ax[1].plot(wall[0,:],wall[1,:],color='r',linewidth=linewidth)
  if(len(lcfs)!=0):
    ax[1].plot(lcfs[0,:],lcfs[1,:],color='r',linewidth=linewidth)
  ax[1].plot(axis[0],axis[1],marker='o',markersize=markersize,markerfacecolor='r',markeredgecolor='r')
  ax[1].set_title('Vertical magnetic field',fontsize=fontsize)
  ax[1].tick_params(axis='x',labelsize=fontsize)
  ax[1].tick_params(axis='y',labelsize=fontsize)
  fig.colorbar(im1,ax=ax[1])
  im2=ax[2].imshow(Bphi,extent=extension)
  ax[2].set_aspect('equal')
  ax[2].set_xlabel('R [m]',fontsize=fontsize)
  ax[2].set_ylabel('Z [m]',fontsize=fontsize)
  if(len(lcfs)!=0):
    ax[2].plot(lcfs[0,:],lcfs[1,:],color='r',linewidth=linewidth)
  if(len(wall)!=0):
    ax[2].plot(wall[0,:],wall[1,:],color='r',linewidth=linewidth)
  ax[2].plot(axis[0],axis[1],marker='o',markersize=markersize,markerfacecolor='r',markeredgecolor='r')
  ax[2].set_title('Toroidal magnetic field',fontsize=fontsize)
  ax[2].tick_params(axis='x',labelsize=fontsize)
  ax[2].tick_params(axis='y',labelsize=fontsize)
  fig.colorbar(im2,ax=ax[2])
  im3=ax[3].imshow(poloidal_flux,extent=extension)
  ax[3].set_aspect('equal')
  ax[3].set_xlabel('R [m]',fontsize=fontsize)
  ax[3].set_ylabel('Z [m]',fontsize=fontsize)
  if(len(lcfs)!=0):
    ax[3].plot(lcfs[0,:],lcfs[1,:],color='r',linewidth=linewidth)
  if(len(wall)!=0):
    ax[3].plot(wall[0,:],wall[1,:],color='r',linewidth=linewidth)
  ax[3].plot(axis[0],axis[1],marker='o',markersize=markersize,markerfacecolor='r',markeredgecolor='r')
  ax[3].set_title('Poloidal flux',fontsize=fontsize)
  ax[3].tick_params(axis='x',labelsize=fontsize)
  ax[3].tick_params(axis='y',labelsize=fontsize)
  fig.colorbar(im3,ax=ax[3])
  plt.show()
  return

# Run main ------------------------------------------------------ #
if __name__ == "__main__":
  load_and_plot_images(filename="magnetic_field_jorek_to_soft.h5",filepath=".")

# --------------------------------------------------------------- #

