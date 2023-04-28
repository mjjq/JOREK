# ----------------------------------------------------- #
# The function plots the particle properties read by
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

# Return the correct x labels and titles in phase space
def define_histogram_labels_titles(key,ptype):
  if(key=='x'):
    titles = ['Major radius','Vertical coordinate','Toroidal angle']
    labels = ['R [m]','Z [m]','phi [r]']
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

# Return the correct x,y labels, titles and aspect ratio in phase space
def define_histogram_labels_titles_aspectratio_2d(key,ptype):
  if(key=='x'):
    titles  = ['Major radius - vertical coordinate','Major radius - toroidal angle',\
    'Vertical coordinate - toroidal angle']
    xlabels = ['R [m]','R [m]','Z [m]']
    ylabels = ['Z [m]','phi [r]','phi [r]']
    aspectratio = [True,False,False]
  elif(key=='v'):
    if(ptype=='particle_kinetic_relativistic'):
      titles  = ['Momenta: px-py','Momenta: px-pz','Momenta: py-pz']
      xlabels = ['px [AMUm/s]','px [AMUm/s]','py [AMUm/s]']
      ylabels = ['py [AMUm/s]','pz [AMUm/s]','pz [AMUm/s]']
      aspectratio = [True,True,True]
    elif(ptype=='particle_gc_relativistic'):
      titles  = ['Parallel momentum - magnetic moment']
      xlabels = ['p_par [AMUm/s]']
      ylabels = ['mu [AMUm^2/Cs]']
      aspectratio = [False]
    else:
      titles  = ['velocities v1-v2','velocities v1-v3','velocities v2-v3']
      xlabels = ['v1','v1','v2']
      ylabels = ['v2','v3','v3']
      aspectratio = [False,False,False]
  return titles,xlabels,ylabels,aspectratio

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

# compute the amount of uninitialised particles w.r.t. the total number of particles
def compute_uninitialised_particles(groups,deathflags,deathvalues):
  for group_id,group in enumerate(groups):
    for deathflag_id,deathflag in enumerate(deathflags):
      values = group[deathflag]
      n_dead_particles = len(values[values<=deathvalues[deathflag_id]]) 
      print("group id: ",group_id," number of particles with ",\
      deathflag,"<=",deathvalues[deathflag_id],1e2*float(n_dead_particles)/len(values),"%, ",n_dead_particles)

# generate 1d histograms for a set of positions (physical or velocity space)
# results are stored in a list having elements of the form: 
# [histogram,histogram_edges]
def create_phase_space_1d_histograms(data_array,p_weights,flags,deathvalue,bins=[]):
  from numpy import histogram
  hists = []
  for ids,data in enumerate(data_array):
    histo,edges = histogram(data[flags>deathvalue],bins=bins[ids],weights=p_weights[flags>deathvalue])
    hists.append([histo,edges])
  return hists

# generate 2d histograms for a set of positions (physical or velocity space)
# results are stored in a list of the form 
def create_phase_space_2d_histograms(data_array,p_weights,flags,deathvalue,bins2d=[]):
  from numpy import histogram2d,amin,amax
  n_histograms = 0
  histos = []
  for id1,data1 in enumerate(data_array):
    local_histos = []
    for id2,data2 in enumerate(data_array[id1+1:]):
      histo,xedges,yedges = histogram2d(data1[flags>deathvalue],data2[flags>deathvalue],\
      bins=[bins2d[id1],bins2d[id1+id2+1]],weights=p_weights[flags>deathvalue])
      local_histos.append([histo,xedges,yedges])
    n_histograms = n_histograms + len(local_histos)
    histos.append(local_histos)
  return histos,n_histograms

# plot 1d histograms
def plot_1d_histograms(hists,titles,xlabels,ylabels,fontsize=18,ncols=3):
  from numpy import ceil as npceil
  from matplotlib.pyplot import subplots,stairs
  n_histos = len(hists)
  nrows = max(int(npceil(n_histos/ncols)),1)
  if(ncols>n_histos):
    ncols = n_histos
  fig,axs = subplots(nrows=nrows,ncols=ncols,facecolor='white',edgecolor='white')
  if(n_histos==1):
    axs = [axs]
  for histo_id,histo in enumerate(hists):
    axs[histo_id].stairs(histo[0],edges=histo[1],fill=True)
    axs[histo_id].set_title(titles[histo_id],fontsize=fontsize)
    axs[histo_id].set_xlabel(xlabels[histo_id],fontsize=fontsize)
    axs[histo_id].set_ylabel(ylabels[histo_id],fontsize=fontsize)
    axs[histo_id].tick_params(axis='x',labelsize=fontsize)
    axs[histo_id].tick_params(axis='y',labelsize=fontsize)
    axs[histo_id].grid()
  return fig,axs

# plot 2d histograms
def plot_2d_histograms(hists2d,n_histos,titles,xlabels,ylabels,aspectequal,\
fontsize=18,ncols=3,colormap='inferno'):
  from numpy import ceil as npceil
  from numpy import amax
  from matplotlib.pyplot import subplots,pcolormesh
  count = 0
  nrows = max(int(npceil(n_histos/ncols)),1)
  if(ncols>n_histos):
    ncols = n_histos
  fig,axs = subplots(nrows=nrows,ncols=ncols,facecolor='white',edgecolor='white')
  if(n_histos==1): 
    axs = [axs]
  for histos_id,histos in enumerate(hists2d):
    for histo_id,histo in enumerate(histos):
      im = axs[count].pcolormesh(histo[1],histo[2],histo[0],\
      cmap=colormap,vmin=0.,vmax=amax(histo[0]),edgecolors='none',shading='flat')
      axs[count].set_title(titles[histos_id+histo_id],fontsize=fontsize)
      axs[count].set_xlabel(xlabels[histos_id+histo_id],fontsize=fontsize)
      axs[count].set_ylabel(ylabels[histos_id+histo_id],fontsize=fontsize)
      axs[count].tick_params(axis='x',labelsize=fontsize)
      axs[count].tick_params(axis='y',labelsize=fontsize)
      if(aspectequal[count]):
        axs[count].set_aspect('equal',adjustable='datalim')
      fig.colorbar(im,ax=axs[count])
      count = count + 1
  return fig,axs
      
# generate and plot 1d histograms from the jorek particle distribution
def computes_1d_histograms_jorek_particles(groups,key,deathvalue,bins=[100,100,100],\
ylabels=['','',''],n_cols=3,fontsize=18):
  from matplotlib.pyplot import show
  for group in groups: 
    # extract titles and x lables
    titles,xlabels = define_histogram_labels_titles(key,group['type'])
    # compute histograms physical space and plot it
    hists = create_phase_space_1d_histograms(group[key],group['weight'],group['i_elm'],\
    deathvalue,bins=bins)
    fig,axs = plot_1d_histograms(hists,titles,xlabels,ylabels,fontsize=fontsize,ncols=n_cols)
    show()

# generate and plot 1d histograms from the jorek particle distribution
def computes_2d_histograms_jorek_particles(groups,key,deathvalue=0,bins=[100,100,100],\
n_cols=3,fontsize=18,colormap='inferno'):
  from matplotlib.pyplot import show
  for group in groups:
    # extract titles, xlabels, ylabels and plot aspect ratio
    titles,xlabels,ylabels,aspectratio = define_histogram_labels_titles_aspectratio_2d(key,group['type'])
    # compute histograms physical space and plot it
    hists2d,n_histos = create_phase_space_2d_histograms(group[key],group['weight'],group['i_elm'],\
    deathvalue,bins2d=bins)
    fig,axs = plot_2d_histograms(hists2d,n_histos,titles,xlabels,ylabels,aspectratio,\
    fontsize=fontsize,ncols=n_cols,colormap=colormap)
    show()
  
# main function
def read_analyse_plot_jorek_restart(filename,filepath,separator,deathvalue=0,bins1d=[100,100,100],
bins2d=[100,100,100],n_cols=3,fontsize=18,colormap='inferno'):
  # read the jorek particle restart data
  p_groups,sim_time = read_jorek_particle_restart_file(filename,filepath,separator)
  # 0d analysis
  compute_uninitialised_particles(p_groups,['i_elm','weight'],[0,0.])
  # plot 1d histograms
  computes_1d_histograms_jorek_particles(p_groups,'x',deathvalue,bins=bins1d,\
  ylabels=['Nphys','Nphys','Nphys'],n_cols=n_cols,fontsize=fontsize)
  computes_1d_histograms_jorek_particles(p_groups,'v',deathvalue,bins=bins1d,\
  ylabels=['Nphys','Nphys','Nphys'],n_cols=n_cols,fontsize=fontsize)
  # plot 2d histograms
  computes_2d_histograms_jorek_particles(p_groups,'x',deathvalue,bins=bins2d,\
  n_cols=n_cols,fontsize=fontsize,colormap=colormap)
  computes_2d_histograms_jorek_particles(p_groups,'v',deathvalue,bins=bins2d,\
  n_cols=n_cols,fontsize=fontsize,colormap=colormap)
  
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
  parser.add_argument('--bins1d','-pbins1d',nargs='*',action='store',required=False,\
  default=[1000,1000,1000],dest='bins1d',\
  help='number of or method for computing the bins for 1D histograms, default: 100')
  parser.add_argument('--bins2d','-pbins2d',nargs='*',action='store',required=False,\
  default=[1000,1000,1000],dest='bins2d',\
  help='number of or method for computing the bins for 2D histograms, default: 100')
  parser.add_argument('--n_cols_subplots','-ncsp',type=int,action='store',required=False,\
  dest='n_cols',default=3,help='number of columns for histogram subplots, default: 3')
  parser.add_argument('--fontsize','-fsize',type=int,action='store',required=False,\
  dest='fontsize',default=18,help='plot font size, default: 18')
  parser.add_argument('--colormap','-cmap',type=str,action='store',required=False,\
  dest='colormap',default='inferno',help='colormap of 2D histograms, default: inferno')
  parser.add_argument('--deathvalue','-dval',type=int,action='store',required=False,\
  dest='deathvalue',default=0,\
  help='flag value below which a particle is set to inactive, default: 0') 
  return parser.parse_args()

# Run main -------------------------------------------- #
if __name__ == "__main__":
  args = generate_argument_parser()
  read_analyse_plot_jorek_restart(args.filename,args.filepath,args.separator,\
  deathvalue=args.deathvalue,bins1d=args.bins1d,bins2d=args.bins2d,\
  n_cols=args.n_cols,fontsize=args.fontsize,colormap=args.colormap)

# ----------------------------------------------------- #
