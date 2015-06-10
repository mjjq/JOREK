!> Program to convert a JOREK2 HDF5 restart file into BINARY restart file
program RST_convert_hdf52bin

  use data_structure
  use phys_module

  implicit none

  type (type_node_list)   , pointer :: node_list
  type (type_element_list), pointer :: element_list

  integer :: ierr, i

  character*5 :: index
  character*13 :: fileout_h5
  character*14 :: fileout_bin

  allocate(node_list)
  allocate(element_list)

#ifndef USE_HDF5
#error " Should be compiled with -DUSE_HDF5"
#endif

  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  rst_format = 0

  ! --- Define the file name
  open (10, file = "file.out")
  read (10,*) index

  write (fileout_h5,'(A5,A5,A3)') 'jorek',index,".h5"
  write (fileout_bin,'(A5,A5,A4)') 'jorek',index,'.rst'

  ! --- Read the restart HDF5 file
  write (6,*) " =============> rst_hdf52bin for filename = ",fileout_h5
  call import_hdf5_restart(node_list, element_list, fileout_h5, rst_format, ierr)

  index_now = index_start
  t_now     = t_start
  
  visco     = visco_rst
  visco_par = visco_par_rst

  eta       = eta_rst

  ! -- Write the BINARY restart file
  write (6,*) " =============> rst_hdf52bin, write BIN file = ",fileout_bin
  call export_binary_restart(node_list, element_list, fileout_bin)

end program RST_convert_hdf52bin
