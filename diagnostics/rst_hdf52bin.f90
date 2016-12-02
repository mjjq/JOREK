!> Program to convert a JOREK2 HDF5 restart file into BINARY restart file
program RST_convert_hdf52bin

  use data_structure
  use phys_module
  ! Argument parsing
  use cla
  use mod_import_restart
  use mod_export_restart

  implicit none

  type (type_node_list)   , pointer :: node_list
  type (type_element_list), pointer :: element_list

  integer :: ierr, i

  character(len=80) :: filein, fileout
  logical :: verbose, file_exists

  allocate(node_list)
  allocate(element_list)

#ifndef USE_HDF5
#error " Should be compiled with -DUSE_HDF5"
#endif

  ! Parse command line arguments
  call cla_init
  call pcla_register('filename', 'name of the restart hdf5 file to convert',  cla_char, 'jorek_restart.h5')    
  call cla_register('-v','--verbose','enable verbose output', cla_flag,'v')
  call cla_validate("rst_hdf52bin")
  call cla_get('filename',filein)
  verbose = cla_key_present('--verbose')

  ! Create output filename
  fileout = filein(1:index(filein,'.h5',.true.)) // 'rst' ! .true. searches backwards

  ! --- Check for presence of the restart file
  inquire(file=filein, exist=file_exists)
  if (.not. file_exists) then
    write(*,*) "File " // trim(filein) // " not found", ""
    call cla_help('rst_hdf52bin')
    call exit(1)
  endif

  ! --- Initialize mode and mode_type arrays
  call det_modes()

  ! --- Read the restart HDF5 file
  if (verbose) write (6,*) " =============> rst_hdf52bin for filename = ",filein
  call import_hdf5_restart(node_list, element_list, filein, ierr)

  index_now = index_start
  t_now     = t_start
  
  visco     = visco_rst
  visco_par = visco_par_rst

  eta       = eta_rst

  ! -- Write the BINARY restart file
  if (verbose) write (6,*) " =============> rst_hdf52bin, write BIN file = ",fileout
  call export_binary_restart(node_list, element_list, fileout)

end program RST_convert_hdf52bin
