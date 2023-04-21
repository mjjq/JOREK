# ----------------------------------------------------------- #
# The function plots the particle diagnostics produce by a    #
# JOREK particle diagnostics                                  #
# ----------------------------------------------------------- #

# read a jorek particle diagnistics file. Zero ended bytes objects
# (S-type) are identified and transformed in strings
def read_jorek_particle_diagnostics_file(filename,filepath,separator):
  from h5py import File as h5pyfile
  from numpy import array
  groups = []
  fhandler = h5pyfile("".join([filepath,separator,filename]))
  for group in fhandler['groups'].values():
    diagnostics = {}
    for diag_name,diag in group.items():
      diagnostics[diag_name] = array(diag)
      if('S' in str(diagnostics[diag_name].dtype)):
        diagnostics[diag_name] = str(diagnostics[diag_name])
        diagnostics[diag_name] = diagnostics[diag_name][2:-1]
    groups.append(diagnostics)
  psi_axis = array(fhandler['psi_axis'])
  psi_bnd  = array(fhandler['psi_sep'])
  fhandler.close()
  return groups,psi_axis,psi_bnd

# main function
def read_analyse_plot_jorek_diagnostics(filename,filepath,separator):
  # read jorek particle diagnostics file
  p_groups,psi_axis,psi_bnd = read_jorek_particle_diagnostics_file(filename,filepath,separator)

# argument parser
def generate_argument_parser():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='read and plot JOREK particle diagnostics')
  parser.add_argument('--filename','-f',type=str,action='store',required=False,\
  dest='filename',default='part_diag.h5',\
  help='name of the particle diagnostic file, default: part_diag.h5')
  parser.add_argument('--filepath','-fpath',type=str,action='store',required=False,\
  dest='filepath',default='.',help='path to the particle diagnostic file, default: .')
  parser.add_argument('--separator','-sep',type=str,action='store',required=False,\
  dest='separator',default='/',help='file separator, default: /')
  return parser.parse_args()

# Run main -------------------------------------------------- #
if __name__ == "__main__":
  args = generate_argument_parser()
  read_analyse_plot_jorek_diagnostics(args.filename,args.filepath,args.separator)
# ----------------------------------------------------------- #


