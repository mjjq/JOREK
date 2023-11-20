# --------------------------------------------------- #
# The procedures contained in the module are used     #
# for modifying a fortran file from a input           #
# dictionary line by line. Example of its usage is    #
# the modification of particle example files for      #
# non regression testing. The method is write safe    #
# in the sens that modifications are applied only     #
# to a copy of the source file but not on the source. #
# --------------------------------------------------- #

# Class containing procedures used for modifying
# FORTRAN files
class ModifyFortranFile():
  # Class constructor
  # inputs:
  #   filename: (string) the name of the file
  #   filedir:  (string) directory containing the file 
  def __init__(self,source_file,source_dir='.', \
    dest_dir='.',modify_name='_modified',separator='/'): 
    from pathlib import Path
    self.separator = separator
    self.source_path = Path("".join([source_dir,separator,source_file]))
    self.dest_path   = Path("".join([dest_dir,separator,\
    self.source_path.stem,modify_name,self.source_path.suffix]))

  # Class destructor, reset attributes to zero for safety
  def __del__(self):
    self.separator    = None
    self.source_path  = None
    self.dest_path    = None

  # Override class printer
  def __str__(self):
    return "".join(['Class for modiying FORTRAN ocdes source: ',\
    self.source_path.as_posix(),', destination: ',self.dest_path.as_posix()])

  # Copy file from source to destination
  def copy_file(self):
    # define substring to search for
    from shutil import copyfile
    copyfile(str(self.source_path),str(self.dest_path))

  # convert a variable in a fortran string
  def convert_variable_fortran_string(self,variable):
    if(variable is dict):
      print('Warning dictionaries are not converted in fortran string!')
      return ''
    fortran_string = str(variable)
    fortran_string = fortran_string.replace('False','.false.')
    fortran_string = fortran_string.replace('True','.true.')
    fortran_string = fortran_string.replace('(','[')
    fortran_string = fortran_string.replace(')',']')
    fortran_string = fortran_string.replace('{','[')
    fortran_string = fortran_string.replace('}',']')
    return fortran_string

  # Modify variables in destination files from dictionary.
  # variables are identified looking for 'key =' substrings
  # while reading the file. Warning: the substitution 
  # operation is performed for each matching line!
  # inputs: 
  #   dict_variables: (dict) dictionary containing the variables
  #                   to be substituted
  def modify_variables(self,dict_variables):
    from string import whitespace
    import re
    # create dictionary of compiled regular expression
    dict_compile_reg    = dict((key,re.compile(re.escape(key+\
    "(\s+|)=.*?(;|$)")) for key in dict_variables.keys())
    dict_substitute_str = dict((key,"".join([key,'=',self.convert_variable_fortran_string(\
    variable),';'])) for key,variable in dict_variables.items()) 
    with self.dest_path.open(mode='r+') as file:
      for line in file:
        replaced_line = line
        if(any(key in line for key in dict_variables)):
          for key,substitute in dict_substitute_str.items():
            replaced_line = dict_compile_reg[key].sub(substitute,replaced_line)
        print(replaced_line)

# Test main
if __name__ == '__main__':
  source_file    = 'ex1.f90'
  source_dir     = './particles/examples'
  dest_dir       = './particles'
  modify_name    = 'mode'
  separator      = '+'
  dict_variables = {'p%x':1,'p%v':[1,2,3]}

  # Constructor 
  fortran_modifier = ModifyFortranFile(source_file,source_dir=source_dir,\
  dest_dir=dest_dir,modify_name=modify_name,separator=separator)  
  print(fortran_modifier)
  del fortran_modifier
  fortran_modifier = ModifyFortranFile(source_file,source_dir=source_dir)  
  print(fortran_modifier)
  fortran_modifier.copy_file()
  fortran_modifier.modify_variables(dict_variables)
  del fortran_modifier 
