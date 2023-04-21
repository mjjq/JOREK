# ----------------------------------------------------- #
# The function plot the particle properties read by
# a JOREK restart file
# ----------------------------------------------------- #
# Read datasets from HDF5 files
def cylindrical_to_cartesian(RZphi):
  from numpy import cos,sin
  return np.array([RZphi[0]*cos(-RZphi[2]),RZphip[0]*sin(-RZphi[2]),RZphi[1]])

# compute the radial volumes for cylindrical coordinates
def compute_cylindrical_radial_volumes(r):
  from numpy import power
  return 5e-1*(power(r[2:],2)-power(r[1:-2],2))
# compute cartesian volumes
def computes_cartesian_volumes(x):
  return x[2,:]-x[1:-2]
# compute the radial volumes for spherical coordinates
def compute_spherical_radial_volumes(r):
  from numpy import power
  return (power(r[2:],3)-power(r[1:-2],3))/3e0
# compute spherical azimuthal volume
def compute_spherical_azimuthal_volumes(theta):
  from numpy import cos
  return cos(theta[1:-2])-cos(theta[2:])

# Return the correct y labels and titles in velocity space
def define_histogram_labels_titles(key,ptype):
  if(key=='x'):
    titles = ['Major radius','Vertical coordinate','Toroidal angle']
    labels = ['R [m]','Z [m]','\phi [r]']
  elif(key=='v'):
    if(ptype=='particle_kinetic_relativistic'):
      titles = ['x-cartesian momentum','y-cartesian momentum','z-cartesian momentum']
      labels = ['px [AMUm/s]','py [AMUm/s]','pz [AMUm/s]']
    elif(ptype=='particle_gc_relativistic'):
      titles = ['Parallel momentum','Magnetic moment']
      labels = ['p_par [AMUm/s]','mu [AMUm^2/Cs]']
    else:
      titles = ['velocity 1','velocity 2','velocity 3']
      labels = ['v1','v2','v3']
  return titles,labels

# Read a jorek particle restart file. Zero ended bytes structre 
# ('S'-type) are identified and transformed in strings
def read_jorek_particle_restart_file(filename,filepath,separator):
  import h5py
  from numpy import array,float64,transpose,floor_divide,pi
  groups = []
  fhandler = h5py.File("".join([filepath,separator,filename]))
  for group in fhandler['groups'].values():
    particles = {}
    for part_data_name,part_data in group.items():
      particles[part_data_name] = array(part_data)
      if((part_data_name=='x') or (part_data_name=='v')):
        particles[part_data_name] = transpose(particles[part_data_name])
      if(part_data_name=='x'):
        particles[part_data_name][2,:] = particles[part_data_name][2,:]-\
        2e0*pi*floor_divide(particles[part_data_name][2,:],2e0*pi)
      if('S' in str(particles[part_data_name].dtype)):
        particles[part_data_name] = str(particles[part_data_name])
        particles[part_data_name] = particles[part_data_name][2:-1]
    groups.append(particles)
  sim_time = float64(fhandler['time'])
  fhandler.close()
  return groups,sim_time

# generate 1d histograms for a set of positions (physical or velocity space)
# results are stored in a list having elements of the form: 
# [histogram,histogram_edges]
def create_phase_space_1d_histograms(data_array,p_weights,bins=[]):
  from numpy import histogram
  hists = []
  for ids,data in enumerate(data_array):
    histo,edges = histogram(data,bins=bins[ids],weights=p_weights)
    hists.append([histo,edges])
  return hists

# generate 2d histograms for a set of positions (physical or velocity space)
# results are stored in a list of the form 
def create_phase_space_2d_histograms():

# generate 2d histograms for a set of positions (physical or velocity space)
# given a values along a third axis results are stored in a list of the form 
def create_phase_space_2d_histograms_slice():

# plot 1d histograms
def plot_1d_histograms(hists,titles,xlabels,ylabels,fontsize=18,ncols=3):
  from numpy import ceil as npceil
  from matplotlib.pyplot import subplots,stairs
  nrows = int(npceil(len(hists)/ncols))
  if(ncols>len(hists)):
    ncols = len(hists)
  fig,axs = subplots(nrows=nrows,ncols=ncols,facecolor='white',edgecolor='white')
  for histo_id,histo in enumerate(hists):
    axs[histo_id].stairs(histo[0],edges=histo[1],fill=True)
    axs[histo_id].set_title(titles[histo_id],fontsize=fontsize)
    axs[histo_id].set_xlabel(xlabels[histo_id],fontsize=fontsize)
    axs[histo_id].set_ylabel(ylabels[histo_id],fontsize=fontsize)
    axs[histo_id].tick_params(axis='x',labelsize=fontsize)
    axs[histo_id].tick_params(axis='y',labelsize=fontsize)
    axs[histo_id].grid()
  return fig,axs

# generate and plot 1d histograms from the jorek particle distribution
def computes_1d_histograms_jorek_particles(groups,key,bins=[100,100,100],\
ylabels=['','',''],n_cols=3,fontsize=18):
  from matplotlib.pyplot import show
  for group in groups: 
    # extract titles and x lables
    titles,xlabels = define_histogram_labels_titles(key,group['type'])
    # compute histograms physical space
    hists = create_phase_space_1d_histograms(group[key],group['weight'],bins=bins)
    fig,axs = plot_1d_histograms(hists,titles,xlabels,ylabels,fontsize=fontsize,ncols=n_cols)
    show()

# main function
def read_analyse_plot_jorek_restart(filename,filepath,separator,\
pos_bins=[100,100,100],n_cols_pos=3,fontsize=18):
  # read the jorek particle restart data
  p_groups,sim_time = read_jorek_particle_restart_file(filename,filepath,separator)
  # plot 1d histograms
  computes_1d_histograms_jorek_particles(p_groups,'x',bins=pos_bins,\
  ylabels=['Nphys','Nphys','Nphys'],n_cols=n_cols_pos,fontsize=fontsize)
  computes_1d_histograms_jorek_particles(p_groups,'v',bins=pos_bins,\
  ylabels=['Nphys','Nphys','Nphys'],n_cols=n_cols_pos,fontsize=fontsize)
 
# argument parser 
def generate_argument_parser():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='read and plot a JOREK particle restart')
  parser.add_argument('--filename','-f',type=str,action='store',required=False,\
  dest='filename',default='part_restart.h5',\
  help='name of the particle restart file, default: part_restart.h5')
  parser.add_argument('--filepath','-fpath',type=str,action='store',required=False,\
  dest='filepath',default='.',help='path of the file to be read, default: .')
  parser.add_argument('--separator','-sep',type=str,action='store',required=False,\
  dest='separator',default='/',help='file separator, default: /')
  parser.add_argument('--pos_bins','-pbins',nargs='*',action='store',required=False,\
  default=[100,100,100],dest='pos_bins',\
  help='number of or method for computing the bins in physical space, default: auto')
  parser.add_argument('--n_cols_subplots_pos','-ncsp',type=int,action='store',required=False,\
  dest='n_cols_pos',default=3,help='number of columns for physical space subplots, default: 3')
  parser.add_argument('--fontsize','-fsize',type=int,action='store',required=False,\
  dest='fontsize',default=18,help='plot font size, default: 18')
  return parser.parse_args()

# Run main -------------------------------------------- #
if __name__ == "__main__":
  args = generate_argument_parser()
  read_analyse_plot_jorek_restart(args.filename,args.filepath,args.separator,\
  pos_bins=args.pos_bins,n_cols_pos=args.n_cols_pos,fontsize=args.fontsize)

# ----------------------------------------------------- #
