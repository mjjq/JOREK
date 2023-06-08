# small python script used for plotting the variance and 
# standard deviation convergence from synchrotron images
def read_datasets_from_hdf5(filenames,filepath,separator):
  from h5py import File
  from numpy import array,transpose
  stddev   = []
  variance = []
  average  = []
  for filename in filenames:
    path = "".join([filepath,separator,filename])
    print('Reading file: ',path)
    fhandler = File(path)
    stddev.append(transpose(array(fhandler['standard_deviation']),[2,1,0]))
    variance.append(transpose(array(fhandler['variance']),[2,1,0]))
    average.append(transpose(array(fhandler['average']),[2,1,0]))
    fhandler.close()
  return stddev,variance,average

# compute the number of rows of a subplot given the number
# of plots and subplot columns
def compute_number_of_subplot_rows(nplots,ncols):
  from numpy import ceil
  nrows = int(ceil(nplots/ncols))
  return nrows

# Method used for plotting 4d arrays as 2d image
def imshow_3d(fig,ax,frame,x_positions=[],y_positions=[],title="",fontsize=12):
  from numpy import newaxis
  from matplotlib.pyplot import figure,imshow,colorbar
  from matplotlib.pyplot import title as tit
  # plot with extension
  if(not ((len(x_positions)==0) and (len(y_positions)==0))):
    dx = 0.5*(x_positions[1]-x_positions[0])
    dy = 0.5*(y_positions[1]-y_positions[0])
    ext = [x_positions[0]-dx, x_positions[-1]+dx,\
    y_positions[0]-dy,y_positions[-1]+dy]
  else:
    ext = None
  for spectrum_id,spectrum in enumerate(frame):
    im=ax.imshow(spectrum,extent=ext)
    fig.colorbar(im,ax=ax)
    ax.set_title("".join([title,' spectrum N# ',str(spectrum_id+1)]),fontsize=fontsize)
  return

# Plot convergence for 3d arrays of the form [nx,ny,nerror]
def plot_convergence_3d(ax,x,values,accept_rate=2e-1,title='',fontsize=12,xlabel='x',\
ylabel='error',loglog=True,linewidth=3,markersize=10,markertype='s'):
  from numpy.random import rand
  if(loglog):
    for valuesyid,valuesy in enumerate(values):
      print('Doing pixel row: ',valuesyid,' accept rate: ',accept_rate)
      for valuesx in valuesy:
        if(all(valuesx>0e0) and rand()<accept_rate):
          ax.loglog(x,valuesx,linewidth=linewidth,markersize=markersize,marker=markertype)
  else:
    for valuesyid,valuesy in enumerate(values):
      print('Doing pixel row: ',valuesyid,'accept rate: ',accept_rate)
      for valuesx in valuesy:
        if(all(valuesx>0e0) and rand()<accept_rate):
          ax.plot(x,valuesx,linewidth=linewidth,markersize=markersize,maker=markertype)
  ax.set_title(title,fontsize=fontsize)
  ax.set_xlabel(xlabel,fontsize=fontsize)
  ax.set_ylabel(ylabel,fontsize=fontsize)
  ax.tick_params(axis='x',labelsize=fontsize)
  ax.tick_params(axis='y',labelsize=fontsize)
  ax.grid()

# plot multiple subplots
def plot_frame_lists(data,titleroot,fontsize=12,ncols=3):
  from numpy import array
  from matplotlib.pyplot import figure,show
  nplots = len(data)
  nrows = compute_number_of_subplot_rows(nplots,ncols)
  if(nplots<=ncols):
    ncols = nplots
    nrows = 1
  # plot variance
  fig = figure(facecolor='white',edgecolor='white')
  axs = [fig.add_subplot(int("".join([str(nrows),str(ncols),str(ids+1)]))) for ids in range(nplots)]
  for datid,dat in enumerate(data):
    imshow_3d(fig,axs[datid],dat,title=titleroot,fontsize=fontsize) 
  show()

# plot the variance and standard deviation values per pixel
def plot_per_pixel_quantities(x,data,title='',accept_rate=0.2,fontsize=12,xlabel='x',\
ylabel='error',loglog=True,linewidth=3,markersize=10,markertype='s',ncols=3):
  from numpy import array,transpose
  from matplotlib.pyplot import figure,show
  fig = figure(facecolor='white',edgecolor='white')
  data = transpose(array(data),[1,2,3,0])
  nrows = compute_number_of_subplot_rows(data.shape[0],ncols)
  if(data.shape[0]<=ncols):
    ncols = data.shape[0]
    nrows = 1
  axs = [fig.add_subplot(int("".join([str(nrows),str(ncols),str(ids+1)]))) for ids in range(data.shape[0])]
  for dataid,dat in enumerate(data):
    plot_convergence_3d(axs[dataid],x,dat,accept_rate=accept_rate,\
    title="".join([title,' for spectrum id: ',str(dataid)]),\
    fontsize=fontsize,xlabel=xlabel,ylabel=ylabel,loglog=loglog,linewidth=linewidth,\
    markersize=markersize,markertype=markertype)
  show()

# plot average,variance and standard deviation
def plot_average_variance_stddev(average,variance,stddev,fontsize=12,ncols=3):
  plot_frame_lists(average,'Image intensity mean: ',fontsize=fontsize,ncols=ncols)
  plot_frame_lists(variance,'Image intensity variance: ',fontsize=fontsize,ncols=ncols)
  plot_frame_lists(stddev,'Image intensity std deviation: ',fontsize=fontsize,ncols=ncols)

# plot convergence w.r.t. x of the variance and the standard deviation
def plot_variance_stdded_covergence_per_pixel(x,variance,stddev,accept_rate=0.2,\
nplots=10,fontsize=12,xlabel='x',loglog=True,linewidth=3,markersize=10,markertype='s',ncols=3,):
  plot_per_pixel_quantities(x,stddev,accept_rate=accept_rate,title='Standard deviation convergence',\
  fontsize=fontsize,xlabel=xlabel,ylabel='$\sigma$',loglog=loglog,\
  linewidth=linewidth,markersize=markersize,markertype=markertype,ncols=ncols)
  plot_per_pixel_quantities(x,variance,accept_rate=accept_rate,title='Variance convergence',\
  fontsize=fontsize,xlabel=xlabel,ylabel='$\sigma^2$',loglog=loglog,\
  linewidth=linewidth,markersize=markersize,markertype=markertype,ncols=ncols)

# main functions
def read_and_plot_data(x,filenames,filepath,separator,accept_rate=0.2,fontsize=12,\
xlabel='x',loglog=True,linewidth=3,markersize=10,markertype='s',ncols=3):
  stddev,variance,average = read_datasets_from_hdf5(filenames,filepath,separator)
  plot_average_variance_stddev(average,variance,stddev,ncols=ncols,fontsize=fontsize)
  plot_variance_stdded_covergence_per_pixel(x,variance,stddev,accept_rate=accept_rate,\
  fontsize=fontsize,xlabel=xlabel,loglog=True,linewidth=linewidth,markersize=markersize,\
  markertype=markertype,ncols=ncols)

# argument parser
def generate_argumet_parser():
  from numpy import float64
  import argparse
  parser = argparse.ArgumentParser(description='plot per pixel light intensity avrage, variance and standard deviation')
  parser.add_argument('--conv_values','-cval',type=float64,nargs='*',required=True,action='store',\
  dest='x',help='list of values used for convergence study')
  parser.add_argument('--filename','-f',type=str,nargs='*',required=True,action='store',\
  dest='filenames',help='name of the files containing the light intensity average, variance and std deviation')
  parser.add_argument('--filepath','-fp',type=str,required=False,\
  dest='filepath',action='store',default='.',help='path to the hdf5 to be loaded, default: .')
  parser.add_argument('--separator','-sep',type=str,required=False,\
  dest='separator',action='store',default='/',help='file separator, default: /')
  parser.add_argument('--fontsize','-fsize',type=int,required=False,\
  dest='fontsize',action='store',default=12,help='fontsize for videos, default: 12')
  parser.add_argument('--markersize','-msize',type=int,required=False,\
  dest='markersize',action='store',default=10,help='marker size, default: 10')
  parser.add_argument('--markertype','-mtype',type=str,required=False,\
  dest='markertype',action='store',default='s',help='marker type, default: s (square)')
  parser.add_argument('--linewidth','-lw',type=int,required=False,\
  dest='linewidth',action='store',default=3,help='line width, default: 3')
  parser.add_argument('--n_plot_columns','-ncols',type=int,required=False,\
  dest='ncols',action='store',default=2,help='number of plot columns, default: 2')
  parser.add_argument('--xlabel','-xlab',type=str,required=False,dest='xlabel',action='store',\
  default='N# particles',help='name of the convergence parameter, default: N# particles')
  parser.add_argument('--loglog_plot','-logp',type=bool,required=False,\
  dest='loglog',action='store',default=True,help='plot convergence using logarithm plot, default: True')
  parser.add_argument('--pixel_accept_rate','-parp',type=float64,required=False,\
  dest='accept_rate',action='store',default=0.2,help='pixel accept rate for plotting, default: 0.2')
  return parser.parse_args()

if __name__ == '__main__':
  args = generate_argumet_parser() 
  read_and_plot_data(args.x,args.filenames,args.filepath,args.separator,accept_rate=args.accept_rate,\
  xlabel=args.xlabel,loglog=args.loglog,linewidth=args.linewidth,markersize=args.markersize,\
  markertype=args.markertype,fontsize=args.fontsize,ncols=args.ncols)
