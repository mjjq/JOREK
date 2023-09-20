# --------------------------------------------------------------- #
# Script aimed to the automatic generation, compilation and       #
# execution of unit tests. Unit test results are synthesized in   #
# a XML file compatible with the Atlassian Bamboo CI/CD installed #
# on the ITER HPC infrastructure.                                 #
# --------------------------------------------------------------- #

# create a dictionary of empty list from keys
# inputs:
#   keys_list:  (list,strings) list of keys
# outputs:
#   list_dict:  (dict,list) dictionary of empty lists
def create_list_dictionary_from_keys(keys_list):
  list_dict = {}
  for key in list(set(keys_list)):
    list_dict[key] = []
  return list_dict

# Find unit test modules in a given folder. Unit test modules are
# identified via the prefix mod_, the suffix _test and the
# extension .f90
# inputs:
#   test_dir:      (string) path to the search directory
#   test_prefix:   (string) prefix identifying a test module
#   test_suffix:   (string) suffix identifying a test module
#   test_ext:      (string) file extension identifying a test module
#   test_parallel: (list,string) list of the parallelization types 
# outputs:
#   test_modules: (dict) dictionary containing the filenames
#                 of the test modules
def find_unit_test_modules(test_dir,test_prefix,test_suffix,\
test_ext,test_parallel):
  from pathlib import Path
  # initialize test module dictionary
  test_parallel.append('serial')
  test_modules = create_list_dictionary_from_keys(test_parallel)
  test_parallel.remove('serial')
  # loop posix path for all test modules in test_dir and 
  # store them in the test_modules dictionary as a 
  # function of their parallelism
  test_dir_path = Path(test_dir)
  for posix_path in  test_dir_path.glob(\
  "".join([test_prefix,'*',test_suffix,'.',test_ext])):
   for key in test_parallel:
     if(key in posix_path.name):
       test_modules[key].append(posix_path)
     else:
       test_modules['serial'].append(posix_path)
  # return the module dictionary   
  return test_modules 

# Argument parsers ---------------------------------------------- #
# outpus:
#   parser: (namespace) namespace having the inputs as attributes
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
  import sys
  args = generate_argument_parser()
  test_modules = find_unit_test_modules(args.test_dirs[0],'mod_','_test',\
  'f90',['mpi'])
  print(test_modules)
  # exit program with success
  sys.exit(0)
# End-of-the-scripts -------------------------------------------- #
