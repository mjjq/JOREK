# --------------------------------------------- #
# The python script generates a new fast camera #
# example file modifying the lines requested by #
# the used. The main application is to          #
# automatically generate example files with the #
# right input parameters for running non        #
# regression tests.                             #
# --------------------------------------------- #

# Argument parser -------------------------------
# outputs:
#   parser: (namespace) namespace containing the 
#           inputs as attributes
def generate_argument_parser():
  from argparse import ArgumentParser
  parser = ArgumentParser(description=\
  "Automatically generate Fast camera example files with user's inputs")
  parser.add_argument('--example_dir','-ed',type=str,required=False,\
  action='store',dest='example_dir',default='./particles/postprocessors/examples',\
  help='Directory path to the examples, default: ./particles/postprocessors/examples')
  parser.add_argument('--example_name','-efn',type=str,required=False,\
  action='store',dest='example_name',default='camera_RE_gyroaverage_synchrotron_example',\
  help='Filename of the example to be modified, default: camera_RE_gyroaverage_synchrotron_example')
  return parser.parse_args()

# Execute script --------------------------------
if __name__ == "__main__":
  args = generate_argument_parser()
# -----------------------------------------------
