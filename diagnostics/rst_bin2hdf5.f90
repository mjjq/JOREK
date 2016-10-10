!> Program to convert a JOREK2 BINARY restart file into HDF5 restart file
program RST_convert_bin2hdf5

  use data_structure
  use phys_module
  ! Argument parsing
  use cla
  use import_restart
  use export_restart

  implicit none

  type (type_node_list)   , pointer :: node_list
  type (type_element_list), pointer :: element_list

  integer :: ierr, i

  character(len=80) :: filein, fileout
  logical :: verbose, file_exists

#ifndef USE_HDF5
#error " Should be compiled with -DUSE_HDF5"
#endif

  ! Parse command line arguments
  call cla_init
  call pcla_register('filename', 'name of the restart file to convert',  cla_char, 'jorek_restart.rst')    
  call cla_register('-v','--verbose','enable verbose output', cla_flag,'v')
  call cla_validate("rst_bin2hdf5")
  call cla_get('filename',filein)
  verbose = cla_key_present('--verbose')

  ! Create output filename
  fileout = filein(1:index(filein,'.rst',.true.)) // 'h5' ! .true. searches backwards

  allocate(node_list)
  allocate(element_list)

  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  rst_format = 0

  ! --- Check for presence of the restart file
  inquire(file=filein, exist=file_exists)
  if (.not. file_exists) then
    write(*,*) "File " // trim(filein) // " not found", ""
    call cla_help('rst_bin2hdf5')
    call exit(1)
  endif

  ! --- Read the restart binary file
  if (verbose) write (6,*) " =============> rst_bin2hdf5 for filename = ",filein
  call import_binary_restart(node_list, element_list, filein, rst_format, ierr)

  index_now = index_start
  t_now     = t_start

  visco     = visco_rst
  visco_par = visco_par_rst

  eta       = eta_rst

  ! -- Write the HDF5 restart file
  if (verbose) write (6,*) " =============> rst_bin2hdf5, write HDF5 file = ",fileout
  call export_hdf5_restart(node_list, element_list, fileout)

end program RST_convert_bin2hdf5
