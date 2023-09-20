# --------------------------------------------------------------- #
# Script aimed to the automatic generation, compilation and       #
# execution of unit tests. Unit test results are synthesized in   #
# a XML file compatible with the Atlassian Bamboo CI/CD installed #
# on the ITER HPC infrastructure.                                 #
# --------------------------------------------------------------- #

# remove file if requested
# inputs:
#   filename: (string) filename of file to be removed
#   remove:   (boolean) if true file is removed
def remove_file(filename,remove):
  from os import system
  if(remove):
    system(''.join(['rm -f ',filename]))

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

# write the pre-fruit basket calls in the unit test file
# inputs:
#   driver_path:        (path)   path to the unit test driver file
#   test_name:          (string) basename of the unit test
#   test_basket_prefix: (string) prefix of the basket subroutine
#   driver_suffix:      (string) suffix of the unit test driver
#   test_prefix:        (string) prefix of the unit test module
#   test_suffix:        (string) suffix of the unit test module
# outputs: 
def write_test_driver_serial(driver_path,test_name,test_basket_prefix,\
driver_suffix,test_prefix,test_suffix):
  with driver_path.open(mode='w') as driver:
    driver.write(''.join(['program ',test_name,driver_suffix,'\n']))
    driver.write('use fruit\n')
    driver.write(''.join(['use ',''.join([test_prefix,test_name,\
    test_suffix,', only: ',test_basket_prefix,test_name,'\n'])]))
    driver.write('  implicit none\n')
    driver.write('\n')
    driver.write('  ! init fruit suite\n')
    driver.write(' call init_fruit_xml\n')
    driver.write(' call init_fruit\n')
    driver.write('\n')
    driver.write(''.join(['  ! run ',test_name,' test basket\n']))
    driver.write(''.join(['  call ',test_basket_prefix,test_name,'\n']))
    driver.write('\n')
    driver.write('  ! write test summary and finalize test suit\n')
    driver.write('  call fruit_summary_xml\n')
    driver.write('  call fruit_summary\n')
    driver.write('  call fruit_finalize\n')
    driver.write(''.join(['end program ',test_name,driver_suffix,'\n']))

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
#   test_modules: (dict) dictionary containing the path associated
#                 to the test modules
def find_unit_test_modules(test_dir,test_prefix,test_suffix,\
test_ext,test_parallel):
  from pathlib import Path
  # initialize test module dictionary
  test_parallel.append('serial')
  test_modules = create_list_dictionary_from_keys(test_parallel)
  test_parallel.remove('serial')
  # loop on the paths for all test modules in test_dir 
  # and store them in the test_modules dictionary as a 
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

# generate a unit test driver from a test module path
# inputs:
#   test_path:          (path) path posix of the test module
#   test_dir:           (string) unit test directory
#   test_parallel:      (string) type of unit test parallelism
#   test_basket_prefix: (string) prefix of the unit test basket
#   test_ext:           (string) test extension
#   driver_suffix:      (string) suffix of the unit test driver
# outputs:
#   driver_path:        (path) path posix of the test driver
def generate_unit_test_driver(test_path,test_parallel,test_basket_prefix,\
test_prefix,test_suffix,test_ext,driver_suffix):
  from pathlib import Path
  # create the unit test driver path
  test_name = test_path.name.replace(test_prefix,'').replace(\
  ''.join([test_suffix,'.',test_ext]),'')
  driver_path = Path(''.join([test_name,driver_suffix,'.',test_ext]))
  # force remove old file and create new one
  driver_path.unlink(missing_ok=True)
  driver_path.touch()
  # write unit test driver as a function of the parallelism
  if(test_parallel=='serial'):
    write_test_driver_serial(driver_path,test_name,test_basket_prefix,\
    driver_suffix,test_prefix,test_suffix)
  return driver_path

# compile the unit test driver using 8 threads
# inputs:
#   driver_path: (path) path posix of the test driver
def compile_unit_test_driver(driver_path):
  import os
  exec_name = driver_path.name.replace(driver_path.suffix,'')
  os.system('make cleanall')
  os.system(''.join(['rm -f ',exec_name]))
  os.system(''.join(['make -j8 ',exec_name]))

# run the unit test driver using 2 mpi tasks maximum 
# and 2 omp threads maximum
# inputs:
#   driver_path:   (path) path posix of the test driver
#   test_parallel: (string) type of unit test parallelism
def run_unit_test_driver(driver_path,test_parallel):
  import os
  # set the number of OMP threads 
  os.system('export OMP_NUM_THREADS=2')
  exec_name = driver_path.name.replace(driver_path.suffix,'')
  if(test_parallel=='serial'):
    os.system(''.join(['./',exec_name]))
  # remove executable
  os.system(''.join(['rm -f ',exec_name]))

# execute unit test: generate, compile, run unit test driver
#   test_path:          (path) path posix of the test module
#   test_dir:           (string) unit test directory
#   test_parallel:      (string) type of unit test parallelism
#   test_basket_prefix: (string) prefix of the unit test basket
#   test_ext:           (string) test extension
#   driver_suffix:      (string) suffix of the unit test driver
#   remove_driver:      (boolean) if true the driver file
def execute_unit_test(test_path,test_parallel,test_basket_prefix,\
test_prefix,test_suffix,test_ext,driver_suffix,remove_driver):
  # generate the unit test driver
  driver_path = generate_unit_test_driver(test_path,test_parallel,\
  test_basket_prefix,test_prefix,test_suffix,test_ext,driver_suffix) 
  # compile the unit test driver
  compile_unit_test_driver(driver_path)
  # execute the unit test driver
  run_unit_test_driver(driver_path,test_parallel)
  # remove the driver
  remove_file(driver_path.name,remove_driver)

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
  execute_unit_test(test_modules['serial'][0],'serial',\
  'run_fruit_','mod_','_test','f90','_test_driver',True)
  # exit program with success
  sys.exit(0)
# End-of-the-scripts -------------------------------------------- #
