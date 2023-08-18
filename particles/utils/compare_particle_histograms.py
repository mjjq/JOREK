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
  particle_tables = []; mass = []; weights = [];
  fhandler = File("".join([filepath,separator,filename]))
  for group in fhandler['groups'].values():
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
    particle_tables.append(particles)
    mass.append(float64(group['mass']))
    weights.append(array(group['weight']))
  sim_time = (float64(fhandler['time']))
  fhandler.close()
  return particle_tables,weights,array(mass),sim_time

# compute the range of coordinates for each dimensions of each particle table
def compute_particle_coordinates_range(particle_tables):
  from numpy import amin,amax,array
  coordinate_ranges = []
  for particles in particle_tables:
    min_coord = amin(particles,axis=1)
    max_coord = amax(particles,axis=1)
    coordinate_ranges.append(array([array(min_coord),array(max_coord)]))
  return coordinate_ranges

# compute the multivariate histograms for all particle tables    
def compute_particle_histograms(particle_tables,weights,coordinate_ranges,nbins):
  from numpy import histogramdd,transpose
  histograms = []
  for particle_id,particles in enumerate(particle_tables):
    histogram = histogramdd(particles,bins=nbins,range=coordinate_ranges[particle_id],\
    density=None,weights=weights[particle_id])
    histograms.append(histogram)
  return histograms

# main function
def compute_error_between_histograms(test_filename,ref_filename,filepath,separator,nbins):
  test_tables,test_weights,test_mass,test_sim_time = read_jorek_particle_restart_file(\
  test_filename,filepath,separator)
  ref_tables,ref_weights,ref_mass,ref_sim_time = read_jorek_particle_restart_file(\
  ref_filename,filepath,separator)

# argument parser
def generate_argument_parser():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='read and compute metrics for statistical coherency of particle histograms')
  parser.add_argument('--filename','-f',type=str,action='store',required=False,\
  dest='filename',default='part_restart.h5',\
  help='name of the particle restart file, default: part_restart.h5')
  parser.add_argument('--ref-filename','-rf',type=str,action='store',required=False,\
  dest='ref_filename',default='part_restart_original.h5',\
  help='name of the reference particle restart file, default: part_restart_original.h5') 
  parser.add_argument('--filepath','-fpath',type=str,action='store',required=False,\
  dest='filepath',default='.',help='path of the file to be read, default: .')
  parser.add_argument('--separator','-sep',type=str,action='store',required=False,\
  dest='separator',default='/',help='file separator, default: /')
  parser.add_argument('--bins','-b',type=int,nargs='*',action='store',required=False,\
  default=[10,10,10,10,10,10,10],dest='nbins',help='number of bins for each particle coordinate, default: 10')
  return parser.parse_args()

# Run main --------------------------------------- #
if __name__ == "__main__":
  args = generate_argument_parser()
  compute_error_between_histograms(args.filename,args.ref_filename,args.filepath,\
  args.separator,args.nbins)

# ------------------------------------------------ #
