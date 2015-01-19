!> Program to convert a JOREK2 BINARY restart file into HDF5 restart file
program RST_convert_bin2hdf5

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

  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  rst_format = 0

  ! --- Define the file name
  open (10, file = "file.out")
  read (10,*) index

  write (fileout_bin,'(A5,A5,A4)') 'jorek',index,'.rst'
  write (fileout_h5,'(A5,A5,A3)') 'jorek',index,".h5"

  ! --- Read the restart binary file
  write (6,*) " =============> rst_bin2hdf5 for filename = ",fileout_bin
  call import_binary_restart(node_list, element_list, fileout_bin, rst_format, ierr)

  index_now = index_start
  t_now     = t_start

  visco     = visco_rst
  visco_par = visco_par_rst

  eta       = eta_rst

  ! -- Write the HDF5 restart file
  write (6,*) " =============> rst_bin2hdf5, write HDF5 file = ",fileout_h5
  call export_hdf5_restart(node_list, element_list, fileout_h5)

end program RST_convert_bin2hdf5
