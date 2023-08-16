# ------------------------------------------------ #
# Python script to be used for comparing the       #
# histograms of two JOREK particle populations.    #
# The JOREK particles should be stored in two      #
# different HDF5 files.                            #
# ------------------------------------------------ #

# Convert JOREK particle phase space data into table
def convert_jorek_particle_velocities_charges_to_table(group):
  from numpy import array,transpose,zeros,float64
  velocity = array([]);  charge = array([]);
  if('q' in group):
    charge = transpose(array(group['q']))
  if('v' in group):
    velocity = transpose(array(group['v']))
  elif('mu' in group):
    if('Vpar' in group):
      dummy_array = array(group['Vpar'])
    elif('E' in group):
      dummy_array = array(group['E'])
    mu       = array(group['mu'])
    velocity = zeros((2,mu.size),dtype=float64)
    velocity[0,:] = dummy_array; velocity[1,:] = mu;   
  return velocity,charge

# Read a JOREK particle restart file. Zero ended
# bytes structure ('S'-type) are identified and
# transformed in strings. The JOREK particle data 
# are homogenized in a N-dimensional table.
def read_jorek_particle_restart_file(filename,filepath,separator):
  from h5py import File
  from numpy import array,float64,transpose,zeros,floor_divide,pi
  particle_tables = {}; mass = {}; weights = {};
  fhandler = File("".join([filepath,separator,filename]))
  for group_number,group in fhandler['groups'].items():
    positions          = transpose(array(group['x']))
    velocities,charges = convert_jorek_particle_velocities_charges_to_table(group)
    ndim = positions.shape[0]
    if(velocities.shape[0]==positions.shape[1]):
      ndim = ndim+1
    else:
      ndim = ndim + velocities.shape[0]
    if(charges.size>0):
      ndim = ndim + 1
    particles = zeros((ndim,positions.shape[1]),dtype=float64)
    particles[0:positions.shape[0],:] = positions
    if(velocities.size>0):
      particles[positions.shape[0]:positions.shape[0]+velocities.shape[0],:] = velocities
    if(charges.size>0):
      particles[-1,:] = charges
    particle_tables[group_number] = particles
    mass[group_number]    = float64(group['mass'])
    weights[group_number] = array(group['weight'])
  sim_time = float64(fhandler['time'])
  fhandler.close()
  return particle_tables,weights,mass,sim_time

# argument parser
def generate_argument_parser():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='read and compute metrics for statistical coherency of particle histograms')
  parser.add_argument('--filename','-f',type=str,action='store',required=False,\
  dest='filename',default='part_restart.h5',\
  help='name of the particle restart file, default: part_restart.h5')
  parser.add_argument('--filepath','-fpath',type=str,action='store',required=False,\
  dest='filepath',default='.',help='path of the file to be read, default: .')
  parser.add_argument('--separator','-sep',type=str,action='store',required=False,\
  dest='separator',default='/',help='file separator, default: /')
  return parser.parse_args()

# Run main --------------------------------------- #
if __name__ == "__main__":
  args = generate_argument_parser()
  particle_tables,weights,mass,sim_time = read_jorek_particle_restart_file(\
  args.filename,args.filepath,args.separator) 
  print(particle_tables['001'][6,:])
  print(weights['001'].shape)
  print(mass)
  print(sim_time)

# ------------------------------------------------ #
