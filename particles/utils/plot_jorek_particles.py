# ----------------------------------------------------- #
# The function plot the particle properties read by
# a JOREK restart file
# ----------------------------------------------------- #
# Read datasets from HDF5 files
def cylindrical_to_cartesian(RZphi):
  from numpy import cos,sin
  return np.array([RZphi[0]*cos(-RZphi[2]),RZphip[0]*sin(-RZphi[2]),RZphi[1]])

# Read a jorek particle restart file
def read_jorek_particle_restart_file(filename,filepath,separator):
  import h5py
  from numpy import array,float64
  groups = []
  fhandler = h5py.File("".join([filepath,separator,filename]))
  for group in fhandler['groups'].values():
    particles = {}
    for part_data_name,part_data in group.items():
      particles[part_data_name] = array(part_data)
    groups.append(particles)
  sim_time = float64(fhandler['time'])
  fhandler.close()
  return groups,sim_time

# main function
def read_analyse_plot_jorek_restart(filename,filepath,separator):
  # read the jorek particle restart data
  p_groups,sim_time = read_jorek_particle_restart_file(filename,filepath,separator)
  
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
  return parser.parse_args()

# Run main -------------------------------------------- #
if __name__ == "__main__":
  args = generate_argument_parser()
  read_analyse_plot_jorek_restart(args.filename,args.filepath,args.separator)

# ----------------------------------------------------- #
