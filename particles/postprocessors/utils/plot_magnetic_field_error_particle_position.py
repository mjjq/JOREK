# read hdf5
def read_datasets_from_hdf5(filename,filepath,datasets,separator):
  import h5py
  import numpy as np
  data = {}
  fhandler = h5py.File("".join([filepath,separator,filename]))
  for setname in datasets:
    data[setname] = np.array(fhandler[setname])
  fhandler.close()
  return data

# extract a subdictionary of all keys containing a string
def extract_subdict_from_string(string,dictionary):
  subdict = dict()
  for key,data in dictionary.items():
    if(string in key):
      subdict[key] = data
  return data

# plot scatter 2d
def plot_scatter_2d(ax,xy,values,title='',xlab='x',ylab='y',fontsize=12,\
markertype='.',markersize='1',colormap='inferno',aspect='equal'):
  from numpy import amin,amax
  mask = ((xy[:,0]==0.)&(xy[:,1]==0.))==False
  ax.scatter(xy[mask,0],xy[mask,1],s=markersize,\
  marker=markertype,c=values[mask],cmap=colormap,\
  vmin=amin(values[mask],vmax=amax(values[mask])))
  if(len(aspect)!=0):
    ax.set_aspect(aspect)
  ax.set_title(title,fontsize=fontsize)
  ax.set_xlabel(xlab,fontsize=fontsize)
  ax.set_ylabel(ylab,fontsize=fontsize)
  ax.tick_params(axis='x',labelsize=fontsize)
  ax.tick_params(axis='y',labelsize=fontsize)
  ax.grid()

# plot scatter 3d
def plot_scatter_3d(ax,xyz,values,title='',xlab='x',ylab='y',zlab='z',\
fontsize=12,markertype='.',markersize='1',markercolor='b',\
aspect='equal',colormap='inferno',rotate=False):
  from numpy import amin,amax
  mask = ((xyz[:,0]==0.)&(xyz[:,1]==0.)&(xyz[:,2]==0.))==False
  ax.scatter3D(xyz[mask,0],xyz[mask,1],xyz[mask,2],\
  s=markersize,marker=markertype,c=values[mask],cmap=colormap,\
  vmin=amin(values[mask]),vmax=amax(values[mask]))
  if(len(aspect)!=0):
    ax.set_aspect(aspect)
  ax.set_title(title,fontsize=fontsize)
  ax.set_xlabel(xlab,fontsize=fontsize)
  ax.set_ylabel(ylab,fontsize=fontsize)
  ax.set_zlabel(zlab,fontsize=fontsize)
  ax.tick_params(axis='x',labelsize=fontsize)
  ax.tick_params(axis='y',labelsize=fontsize)
  ax.tick_params(axis='z',labelsize=fontsize)
  ax.xaxis.set_rotate_label(rotate)
  ax.yaxis.set_rotate_label(rotate)
  ax.zaxis.set_rotate_label(rotate)

# compute the number of rows of a subplot given the number
# of plots and subplot columns
def compute_number_of_subplot_rows(nplots,ncols):
  from numpy import ceil 
  nrows = int(ceil(nplots/ncols)) 
  return nrows
 
# plot the magnetic field errors in the RZ plane
def plot_RZ_error_fields(RZPhi,error_dict,fontsize=12,markertype='.',markersize=1,\
colormap='inferno',aspect='equal',ncols=3):
  from matplotlib.pyplot import figure,show
  nerrors = len(error_dict)
  nrows = compute_number_of_subplot_rows(nerrors,ncols)
  fig     = figure(facecolor='white',edgecolor='white')
  print([int("".join([str(nrows),str(ncols),str(ids+1)])) for ids in range(nerrors)])
  axs = [fig.add_subplot(int("".join([str(nrows),str(ncols),str(ids+1)]))) for ids in range(nerrors)]
  for errorid,errorkey in enumerate(error_dict.keys()):
    plot_scatter_2d(axs[errorid],RZPhi[:,0:2],error_dict[errorkey],\
    title=().join(['Error field: ',errorkey]),xlab='R [m]',ylab='Z [m]',\
    fontsize=fontsize, markertype=markertype,markersize=markersize,\
    colormap=colormap,aspect=aspect)
  show()
    
# plot the magnetic field errors in 3D
def plot_xyz_error_fields(pos,error_dict,fontsize=12,markertype='.',markersize=1,\
colormap='inferno',aspect='equal',rotate=False,ncols=3):
  from matplotlib.pyplot import figure,show
  nerrors = len(error_dict)
  nrows = compute_number_of_subplot_rows(nerrors,ncols)
  fig = figure(facecolor='white',edgecolor='white')
  axs = [fig.add_subplot(int("".join([str(nrows),str(ncols),str(ids+1)])),\
  projection='3d') for ids in range(nerrors)]
  for errorid,errorkey in enumerate(error_dict.keys()):
    plot_scatter_3d(axs[errorid],pos,error_dict[errorkey],\
    title=().join(['Error field: ',errorkey]),xlab='x [m]',ylab='y [m]',zlab='z [m]',\
    fontsize=fontsize,markertype=markertype,markersize=markersize,markercolor=markercolor,\
    aspect=aspect,colormap=colormap,rotate=rotate)
  show()

# plot the magnetic field errors
def read_plot_magnetic_field_error(datasets,filename,filepath='.',separator='/',\
  error_str='error',rzphi_str='RZPhi',pos_str='x',fontsize=12,markertype='.',\
  markersize=1,colormap='inferno',aspect=True,rotate=False,ncols=3):
  data_dict = read_datasets_from_hdf5(filename,filepath,datasets,separator)
  # extract error dictionary
  error_dict = extract_subdict_from_string(error_str,data_dict)
  # plot magnetic field errors in the RZ plane
  plot_RZ_error_fields(data_dict[rzphi_str],error_dict,fontsize=fontsize,\
  markertype=markertype,markersize=markersize,colormap=colormap,aspect=aspect,ncols=ncols)
  # plot 3d magnetic field error
  plot_xyz_error_fields(data_dict[pos_str],error_dict,fontsize=fontsize,\
  markertype=markertype,markersize=markersize,colormap=colormap,aspect=aspect,\
  rotate=rotate,ncols=ncols)

# input parser
def generate_argument_parser():
  import argparse
  parser = argparse.ArgumentParser(description='plot magnetic field error at particle positions')
  parser.add_argument('--filename','-f',type=str,required=False,action='store',\
  dest='filename',default='soft_jorek_magnetic_field_error.h5',\
  help='name of the hdf5 to be read, default: soft_jorek_magnetic_field_error.h5')
  parser.add_argument('--filepath','-fp',type=str,required=False,\
  dest='filepath',action='store',default='.',help='path to the hdf5 to be loaded, default: .')
  parser.add_argument('--separator','-sep',type=str,required=False,\
  dest='separator',action='store',default='/',help='file separator, default: /')
  parser.add_argument('--fontsize','-fsize',type=int,required=False,\
  dest='fontsize',action='store',default=16,help='fontsize for videos, default: 16')
  parser.add_argument('--markersize','-msize',type=int,required=False,\
  dest='markersize',action='store',default=1,help='marker size, default: 1')
  parser.add_argument('--markertype','-mtype',type=str,required=False,\
  dest='markertype',action='store',default='.',help='marker type, default: . (dot)')
  parser.add_argument('--datasets','-dset',type=str,nargs='*',\
  required=False,dest='datasets',action='store',default=['RZPhi','x','error_BR','error_BZ',\
  'error_Bphi','error_Babs'],\
  help='list of the name of the datasets to read, default: [RZPhi,x,error_BR,error_BZ],error_Bphi,error_Babs')
  parser.add_argument('--colormap','-cmap',type=str,required=False,\
  dest='colormap',action='store',default='inferno',help='plot colormap, default: inferno')
  parser.add_argument('--n_plot_columns','-ncols',type=int,required=False,\
  dest='ncols',action='store',default=3,help='number of plot columns, default: 3')
  parser.add_argument('--aspect_equal','-asp',type=bool,required=False,\
  dest='aspect',action='store',default=True,help='use aspect equal for plotting, default: True')
  parser.add_argument('--rotate_3d_plot','-rot3d',type=bool,required=False,\
  dest='rotate',action='store',default=False,help='rotate a 3d plot, default: False')
  parser.add_argument('--error_string','-estr',type=str,required=False,\
  dest='error_str',action='store',default='error',help='string identifying error dataset, default: error')
  parser.add_argument('--cylindrical_string','-cylstr',type=str,required=False,\
  dest='rzphi_str',action='store',default='RZPhi',help='string identifying cylindrical positions, default: RZPhi')
  parser.add_argument('--cartesian_string','-cartstr',type=str,required=False,\
  dest='pos_str',action='store',default='x',help='string identifying cartesian positions, default: x')
  return parser.parse_args()

if __name__ == '__main__':
  args = generate_argument_parser()
  read_plot_magnetic_field_error(args.datasets,args.filename,filepath=args.filepath,\
  separator=args.separator,error_str=args.error_str,rzphi_str=args.rzphi_str,pos_str=args.pos_str,\
  fontsize=args.fontsize,markertype=args.markertype,markersize=args.markersize,colormap=args.colormap,\
  aspect=args.aspect,ncols=args.ncols,rotate=args.rotate)
