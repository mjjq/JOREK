!> Program to convert a JOREK2 HDF5 restart file into BINARY restart file
program RST_convert_hdf52bin

  use data_structure
  use phys_module

  implicit none

  type (type_node_list)   , pointer :: node_list
  type (type_element_list), pointer :: element_list

  integer :: ierr, i

  character*5 :: index
  character*16 :: fileout_h5
  character*17 :: fileout_bin

  allocate(node_list)
  allocate(element_list)

#ifndef USE_HDF5
#error " Should be compiled with -DUSE_HDF5"
#endif

  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  rst_format = 0

  ! --- Read the restart HDF5 file
  fileout_h5 = "jorek_restart.h5"
  write (6,*) " =============> rst_hdf52bin for filename = ",fileout_h5
  call import_hdf5_restart(node_list, element_list, fileout_h5, rst_format, ierr)

  index_now = index_start
  t_now     = t_start
  
  visco     = visco_rst
  visco_par = visco_par_rst

  eta       = eta_rst

  ! -- Write the BINARY restart file
  fileout_bin = "jorek_restart.rst"
  write (6,*) " =============> rst_hdf52bin, write BIN file = ",fileout_bin
  call export_binary_restart(node_list, element_list, fileout_bin)

end program RST_convert_hdf52bin
