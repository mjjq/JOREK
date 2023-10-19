# --------------------------------------------------------------- #
# Script aimed to the automatic generation, compilation and       #
# execution of unit tests. Unit test results are synthesized in   #
# a XML file compatible with the Atlassian Bamboo CI/CD installed #
# on the ITER HPC infrastructure.                                 #
# --------------------------------------------------------------- #

# general purpose routines -------------------------------------- #
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

# routines for handling unit test drivers ----------------------- #

# write the driver file for fruit serial unit tests
# inputs:
#   driver_path:        (path)   path to the unit test driver file
#   test_name:          (string) basename of the unit test
#   test_basket_prefix: (string) prefix of the basket subroutine
#   driver_suffix:      (string) suffix of the unit test driver
#   test_prefix:        (string) prefix of the unit test module
#   test_suffix:        (string) suffix of the unit test module
#   log_fruit_summary:  (bool) if true, the fruit summary is logged 
def write_test_driver_serial(driver_path,test_name,test_basket_prefix,\
driver_suffix,test_prefix,test_suffix,log_fruit_summary):
  with driver_path.open(mode='w') as driver:
    driver.write(''.join(['program ',test_name,driver_suffix,'\n']))
    driver.write('use fruit\n')
    driver.write(''.join(['use ',''.join([test_prefix,test_name,\
    test_suffix,', only: ',test_basket_prefix,test_name,'\n'])]))
    driver.write('  implicit none\n')
    driver.write('\n')
    driver.write('  ! init fruit suite\n')
    driver.write('  call init_fruit_xml\n')
    driver.write('  call init_fruit\n')
    driver.write('  call fruit_hide_dots')
    driver.write('\n')
    driver.write(''.join(['  ! run ',test_name,' test basket\n']))
    driver.write(''.join(['  call ',test_basket_prefix,test_name,'\n']))
    driver.write('\n')
    driver.write('  ! write test summary and finalize test suit\n')
    driver.write('  call fruit_summary_xml\n')
    if(log_fruit_summary):
      driver.write('  call fruit_summary\n')
    driver.write('  call fruit_finalize\n')
    driver.write(''.join(['end program ',test_name,driver_suffix,'\n']))

# write the driver file for fruit mpi-enabled unit tests
# inputs:
#   driver_path:        (path)   path to the unit test driver file
#   test_name:          (string) basename of the unit test
#   test_basket_prefix: (string) prefix of the basket subroutine
#   driver_suffix:      (string) suffix of the unit test driver
#   test_prefix:        (string) prefix of the unit test module
#   test_suffix:        (string) suffix of the unit test module
#   log_fruit_summary:  (bool) if true, the fruit summary is logged
def write_test_driver_parallel(driver_path,test_name,test_basket_prefix,\
driver_suffix,test_prefix,test_suffix,log_fruit_summary):
  with driver_path.open(mode='w') as driver:
     driver.write(''.join(['program ',test_name,driver_suffix,'\n'])) 
     driver.write('use fruit\n')
     driver.write('use fruit_mpi\n')   
     driver.write('use mod_mpi_tools, only: init_mpi_threads\n')
     driver.write('use mod_mpi_tools, only: finalize_mpi_threads\n')
     driver.write(''.join(['use ',''.join([test_prefix,test_name,\
     test_suffix,', only: ',test_basket_prefix,test_name,'\n'])]))
     driver.write('  implicit none\n')
     driver.write('  integer :: rank,n_tasks,ifail\n')
     driver.write('\n')
     driver.write('  ! init the mpi communicator\n')
     driver.write('  call init_mpi_threads(rank,n_tasks,ifail)\n')
     driver.write('\n')
     driver.write('  ! init the fruit suit\n')
     driver.write('  call init_fruit\n')
     driver.write('  call fruit_init_mpi_xml(rank)\n')
     driver.write('  call fruit_hide_dots')
     driver.write('\n')
     driver.write(''.join(['  ! run ',test_name,' test basket\n']))
     driver.write(''.join(['  call ',test_basket_prefix,test_name,\
     '(rank,n_tasks,ifail)','\n']))
     driver.write('\n')  
     driver.write('  ! write test summary and finalize test suit\n')
     if(log_fruit_summary):
       driver.write('  call fruit_summary_mpi(n_tasks,rank)\n')
     driver.write('  call fruit_summary_mpi_xml(n_tasks,rank)\n')
     driver.write('  call fruit_finalize_mpi(n_tasks,rank)\n')
     driver.write('\n')
     driver.write('  ! finalize MPI communicator\n')
     driver.write('  call finalize_mpi_threads(ifail)\n')
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
#   test_basket_prefix: (string) prefix of the unit test basket
#   test_ext:           (string) test extension
#   driver_suffix:      (string) suffix of the unit test driver
#   log_fruit_summary:  (bool) if true, the fruit summary is logged 
# outputs:
#   driver_path:        (path) path posix of the test driver
def generate_unit_test_driver(test_path,test_basket_prefix,\
test_prefix,test_suffix,test_ext,driver_suffix,log_fruit_summary):
  from pathlib import Path
  # create the unit test driver path
  test_name = test_path.name.replace(test_prefix,'').replace(\
  ''.join([test_suffix,'.',test_ext]),'')
  driver_path = Path(''.join([test_name,driver_suffix,'.',test_ext]))
  # force remove old file and create new one
  driver_path.unlink(missing_ok=True)
  driver_path.touch()
  # write unit test driver as a function of the parallelism
  if('mpi' in test_name):
    write_test_driver_parallel(driver_path,test_name,test_basket_prefix,\
    driver_suffix,test_prefix,test_suffix,log_fruit_summary)
  else:
    write_test_driver_serial(driver_path,test_name,test_basket_prefix,\
    driver_suffix,test_prefix,test_suffix,log_fruit_summary)
  return driver_path

# compile the unit test driver using 8 threads
# inputs:
#   driver_path: (path) path posix of the test driver
def compile_unit_test_driver(driver_path):
  from os import system
  exec_name = driver_path.name.replace(driver_path.suffix,'')
  system('make cleanall')
  system(''.join(['rm -f ',exec_name]))
  system(''.join(['make -j8 ',exec_name]))

# find the right launcher for the application
# inputs:
#   exec_name: (string) name of the application
def find_launcher_type(exec_name):
  if('mpi' in exec_name):
    return 'mpi'
  else:
    return 'serial'

# run the unit test driver using 2 mpi tasks maximum 
# and 2 omp threads maximum
# inputs:
#   driver_path:   (path) path posix of the test driver
#   launchers:     (dict(string)) launchers to be invoked for 
#                  executing a unit test application
def run_unit_test_driver(driver_path,launchers):
  from os import system,environ
  # set the number of OMP threads
  omp_num_threads_old = str(1)
  if('OMP_NUM_THREADS' in environ):
    omp_num_threads_old = environ['OMP_NUM_THREADS']
  system('export OMP_NUM_THREADS=2')
  exec_name = driver_path.name.replace(driver_path.suffix,'')
  system(''.join([launchers[find_launcher_type(exec_name)],exec_name]))
  # remove executable
  system(''.join(['rm -f ',exec_name]))
  # restore original number of omp threads
  system(''.join(['export OMP_NUM_THREADS=',omp_num_threads_old]))

# routines for handling FRUIT XML files ------------------------- #

# Read multiple fruit files and store their path in list. Discard
# temporary result file.
# inputs:
#   result_dir:    (string) directory containing the FRUIT results
#   result_prefix: (string) root of the FRUIT result filenames
#   result_ext:    (string) extension of the FRUIT result filenames
#   remove_result: (bool) if True temporary result files are removed  
# outputs:
#   result_paths:  (list,path) path list of all FRUIT result files
def find_fruit_result_file(result_dir,result_prefix,\
result_ext,remove_result):
  from pathlib import Path
  result_dir_path = Path(result_dir)
  result_paths = []
  for result in result_dir_path.glob("".join(\
  [result_prefix,'*','.',result_ext])):
    if('_tmp' not in result.name):
      result_paths.append(result)
    else:
      remove_file(result.name,remove_result)
  return result_paths

# Read the file containing the test results as run by FRUIT.
# The total number of successes and failures is returned.
# inputs:
#   result_path:   (path) path of the FRUIT result file
#   result_map:    (dict) dictionary containing the association
#                         between the fruit errors,tests,failures,id
#                         fields and the index of the integer in list
#   remove_result: (bool) if True result files are removed
# outputs:
#   n_successes: (int) number of test ending in success
#   n_failures:  (int) number of test ending in failures
#   n_errors:    (int) number of test ending in errors
#   test_id:     (int) test id
def read_fruit_result_file(result_path,result_map,remove_result):
  from re import findall
  with result_path.open(mode='r') as result:
    lines = result.readlines()
  # extract the line containing the results
  for line in lines:
    if('failures' in line):
      results = [int(result) for result in findall(r'\d+',line)]
      break
  # remove result file if required
  remove_file(result_path.name,remove_result)
  # return results
  return results[result_map['tests']]-results[result_map['failures']],\
  results[result_map['failures']],results[result_map['errors']],\
  results[result_map['id']]

# Read and reduce all fruit results (more than one result file can
# exist for the same FRUIT test if run using MPI)
# inputs:
#   result_dir:     (string) directory containing the FRUIT results
#   result_prefix:  (string) root of the FRUIT result filenames
#   result_ext:     (string) extension of the FRUIT result filenames
#   result_map:     (dict) dictionary containing the association
#                          between the fruit errors,tests,failures,id
#                          fields and the index of the integer in list
#   remove_results: (boolean) if True result files are removed
# outputs:
#   n_successes_tot: (int) total number of test ending in success
#   n_failures_tot:  (int) total number of test ending in failures
#   n_errors_tot:    (int) total number of test ending in errors
def read_reduce_fruit_results(result_dir,result_prefix,result_ext,\
result_map,remove_results):
  # initialisations
  n_successes_tot=0; n_failures_tot=0; n_errors_tot=0;
  # find all FRUIT result paths
  result_path_list = find_fruit_result_file(result_dir,\
  result_prefix,result_ext,remove_results)
  # read and reduce fruit results
  for result_path in result_path_list:
    n_successes,n_failures,n_errors,test_id = read_fruit_result_file(\
   result_path,result_map,remove_results)
    n_successes_tot = n_successes_tot + n_successes
    n_failures_tot = n_failures_tot + n_failures
    n_errors_tot = n_errors_tot + n_errors    
  return n_successes_tot,n_failures_tot,n_errors_tot

# global unit test routines ------------------------------------- #

# check results, if there is a non-zero number of failures
# or a non zero number of erros, the script return the 
# error exit code of 1. The success exit code of 0 is 
# returned otherwise
# inputs: 
#   n_failures: (integer) total number of failed tests
#   n_errors:   (integer) total number of errors
def check_results(n_failures,n_errors):
  if(n_failures!=0):
    print(''.join(['Unit test failed: N# failures: ',\
    str(n_failures),' N# errors: ',str(n_errors)]))
    return 1
  if(n_errors!=0):
    print(''.join(['Unit test failed: N# failures: ',\
    str(n_failures),' N# errors: ',str(n_errors)]))
    return 1
  else: 
    print('Unit test successfully completed!')
    return 0

# execute unit test: generate, compile, run unit test driver,
# read and reduce fruit result for the unit test
#   test_path:          (path) path posix of the test module
#   test_basket_prefix: (string) prefix of the unit test basket
#   test_ext:           (string) test extension
#   driver_suffix:      (string) suffix of the unit test driver
#   result_dir:         (string) directory containing the FRUIT results
#   result_prefix:      (string) root of the FRUIT result filenames
#   result_ext:         (string) extension of the FRUIT result filenames
#   result_map:         (dict) dictionary containing the association
#                              between the fruit errors,tests,failures,id
#                              fields and the index of the integer in list
#   remove_driver:      (boolean) if true the driver file
#   remove_results:     (boolean) if True result files are removed
#   log_fruit_summary:  (bool) if true, the fruit summary is logged 
def execute_unit_test(test_path,test_basket_prefix,test_prefix,\
test_suffix,test_ext,driver_suffix,launchers,result_dir,\
result_prefix,result_ext,result_map,remove_driver,remove_results,\
log_fruit_summary):
  # generate the unit test driver
  driver_path = generate_unit_test_driver(test_path,test_basket_prefix,\
  test_prefix,test_suffix,test_ext,driver_suffix,log_fruit_summary) 
  # compile the unit test driver
  compile_unit_test_driver(driver_path)
  # execute the unit test driver
  run_unit_test_driver(driver_path,launchers)
  # remove the driver
  remove_file(driver_path.name,remove_driver)
  # find, read and reduce all FRUIT results and remove result files
  n_successes,n_failures,n_errors = read_reduce_fruit_results(\
  result_dir,result_prefix,result_ext,result_map,remove_results)
  return n_successes,n_failures,n_errors

# log unit tests results
# inputs: 
#   unit_test_log: (list(characters)) path of the unit tests to be logged
#   log_header:    (character) header of the unit test log to be printed
def log_unit_test_results(unit_test_log,log_header):
  if(len(unit_test_log)>0):
    print(log_header)
    for unit_test in unit_test_log:
      print("".join(['  ',unit_test]))

# execute the overall suite of unit tests
# inputs:
#   test_dirs:          (list)(string) list of paths to the search directories
#   test_parallel:      (list,string) list of the parallelization types 
#   test_basket_prefix: (string) prefix of the unit test basket
#   test_ext:           (string) test extension
#   driver_suffix:      (string) suffix of the unit test driver
#   result_dir:         (string) directory containing the FRUIT results
#   result_prefix:      (string) root of the FRUIT result filenames
#   result_ext:         (string) extension of the FRUIT result filenames
#   result_map:         (dict) dictionary containing the association
#                              between the fruit errors,tests,failures,id
#                              fields and the index of the integer in list
#   remove_driver:      (boolean) if true the driver file
#   remove_results:     (boolean) if True result files are removed
#   log_fruit_summary:  (bool) if true, the fruit summary is logged 
# outputs:
#   exit_code:          (integer) 0 if all tests terminated successfully
#                                 1 otherwise
def execute_all_unit_tests(test_dirs,test_parallel,test_basket_prefix,\
test_prefix,test_suffix,test_ext,driver_suffix,launchers,result_dir,\
result_prefix,result_ext,result_map,remove_driver,remove_results,\
log_fruit_summary):
  # initialise the failure and error counters
  n_failures=0; n_errors=0; failed_tests=[]; error_tests=[];
  # loop over all test directories
  for test_dir in test_dirs:
    # find all unit test modules in test_dir
    test_modules =  find_unit_test_modules(test_dir,test_prefix,\
    test_suffix,test_ext,test_parallel)
    # execute all unit tests in the test modules
    for tests in test_modules.values():
      for test in tests:
        # execute a unit test
        n_successes_loc,n_failures_loc,n_errors_loc = \
        execute_unit_test(test,test_basket_prefix,test_prefix,\
        test_suffix,test_ext,driver_suffix,launchers,result_dir,\
        result_prefix,result_ext,result_map,remove_driver,\
        remove_results,log_fruit_summary)
        # store the name of failed tests
        if(n_failures_loc!=0): 
          failed_tests.append("".join([test_dir,'/',test.name]))
        # store the name of tests with errors
        if(n_errors_loc!=0):
          error_tests.append("".join([test_dir,'/',test.name]))
        # reduce the total number of successes, failures and errors
        n_failures  = n_failures + n_failures_loc
        n_errors    = n_errors + n_errors_loc
  # log tests with failures and/or errors
  log_unit_test_results(failed_tests,'Unit test module with failed tests:')
  log_unit_test_results(error_tests,'Unit test module with errors:')
  # check the validity of the overall result
  return check_results(n_failures,n_errors)
     
# Argument parsers ---------------------------------------------- #
# outpus:
#   parser: (namespace) namespace having the inputs as attributes
def generate_argument_parser():
  from argparse import ArgumentParser 
  parser = ArgumentParser(\
  description='generate, compile and execute JOREK unit tests')
  parser.add_argument('--directories','-d',type=str,nargs='*',\
  required=False,action='store',dest='test_dirs',default=['./particles/tests'],\
  help='relative paths of the directories containing unit tests')
  parser.add_argument('--parallelisms','-p',type=str,nargs='*',\
  required=False,action='store',dest='test_parallel',\
  default=['mpi'],help='type of parallelism of the unit tests,default: mpi')
  parser.add_argument('--fruit-basket-prefix','-fbp',type=str,\
  required=False,action='store',dest='test_basket_prefix',default='run_fruit_',\
  help='prefix of the fruit basket to be run,default: run_fruit_')
  parser.add_argument('--test-prefix','-tp',type=str,required=False,\
  action='store',dest='test_prefix',default='mod_',\
  help='prefix of the unit test module file, default: mod_')
  parser.add_argument('--test-suffix','-ts',type=str,required=False,\
  action='store',dest='test_suffix',default='_test',\
  help='suffix of the unit test module file, default: _test') 
  parser.add_argument('--test-extension','-te',type=str,required=False,\
  action='store',dest='test_ext',default='f90',\
  help='extension of the unit test module file, default: f90')
  parser.add_argument('--test-driver-suffix','-ds',type=str,required=False,\
  action='store',dest='driver_suffix',default='_test_driver',\
  help='suffix of unit test driver file, default: _test_driver')
  parser.add_argument('--result-dir','-rd',type=str,required=False,\
  action='store',dest='result_dir',default='.',\
  help='folder of the unit test results, default: .')
  parser.add_argument('--result-prefix','-rp',type=str,required=False,\
  action='store',dest='result_prefix',default='result',\
  help='prefix of the test result file, default: result')
  parser.add_argument('--result-extension','-re',type=str,required=False,\
  action='store',dest='result_ext',default='xml',\
  help='extension of the test result file, default: xml')
  parser.add_argument('--remove-driver','-rmd',type=bool,required=False,\
  action='store',dest='remove_drivers',default=True,\
  help='if true the test drivers are removed after execution, default: true')
  parser.add_argument('--remove-results','-rmr',type=bool,required=False,\
  action='store',dest='remove_results',default=True,\
  help='if true the test results are removed after execution, default: true')
  parser.add_argument('--log-fruit-summary','-lfs',type=bool,required=False,\
  action='store',dest='log_fruit_summary',default=False,\
  help='if true the fruit summary is logged, default: false')
  parser.add_argument('--fruit-result-map','-frm',type=int,nargs='*',\
  required=False,action='store',dest='list_fruit_result_map',default=[0,1,2,3],\
  help='relative position of the fruit error,tests,failures,id as read by result file, default: [0,1,2,3]')
  parser.add_argument('--launchers','-l',type=str,nargs='*',\
  required=False,action='store',dest='list_launchers',default=['./','mpirun -np 2 '],\
  help='launcher to be used for executing a test in the order: serial, mpi, default: [./,mpirun -np 2 ]')  
  return parser.parse_args()

# Execute script ------------------------------------------------ #
if __name__ == '__main__':
  from sys import exit as sysexit
  # parse the inputs
  args = generate_argument_parser()
  fruit_result_map = {'errors':args.list_fruit_result_map[0],\
  'tests':args.list_fruit_result_map[1],'failures':args.list_fruit_result_map[2],\
  'id':args.list_fruit_result_map[3]}
  launchers = {'serial':args.list_launchers[0],'mpi':args.list_launchers[1]}
  # run all unit test suites
  print(args.remove_drivers)
  exit_code = execute_all_unit_tests(args.test_dirs,args.test_parallel,\
  args.test_basket_prefix,args.test_prefix,args.test_suffix,args.test_ext,\
  args.driver_suffix,launchers,args.result_dir,args.result_prefix,\
  args.result_ext,fruit_result_map,args.remove_drivers,args.remove_results,\
  args.log_fruit_summary) 
  # exit with the appropriate exit code
  sysexit(exit_code)

# End-of-the-scripts -------------------------------------------- #
