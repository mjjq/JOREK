# --------------------------------------------------------------- #
# Script aimed to the automatic generation, compilation and       #
# execution of unit tests. Unit test results are synthesized in   #
# a XML file compatible with the Atlassian Bamboo CI/CD installed #
# on the ITER HPC infrastructure.                                 #
# --------------------------------------------------------------- #

# Argument parsers ---------------------------------------------- #
def generate_argument_parser():
  import argparse
  parser = argparse.ArgumentParser(\
  description='generate, compile and execute JOREK unit tests')
  parser.add_argument('--directories','-d',type=str,nargs='*',\
  required=False,action='store',dest='test_dirs',default=[\
  './particles/tests','./non_regression_tests/unit_tests'],\
  help='relative paths of the directories containing unit tests')
  return parser.parse_args()

# Execute script ------------------------------------------------ #
if __name__ == '__main__':
  args = generate_argument_parser()
  print('Test directories: ',args.test_dirs)
# End-of-the-scripts -------------------------------------------- #
