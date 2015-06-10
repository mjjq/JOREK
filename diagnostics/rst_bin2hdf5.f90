!> Program to convert a JOREK2 BINARY restart file into HDF5 restart file
program RST_convert_bin2hdf5

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

  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  rst_format = 0

!CP  ! --- Define the file name
!CP  open (10, file = "file.out")
!CP  read (10,*) index
!CP
!CP  write (fileout_bin,'(A5,A5,A4)') 'jorek',index,'.rst'
!CP  write (fileout_h5,'(A5,A5,A3)') 'jorek',index,".h5"

  ! --- Read the restart binary file
  fileout_bin = "jorek_restart.rst"
  write (6,*) " =============> rst_bin2hdf5 for filename = ",fileout_bin
  call import_binary_restart(node_list, element_list, fileout_bin, rst_format, ierr)

  index_now = index_start
  t_now     = t_start

  visco     = visco_rst
  visco_par = visco_par_rst

  eta       = eta_rst

  ! -- Write the HDF5 restart file
  fileout_h5 = "jorek_restart.h5"
  write (6,*) " =============> rst_bin2hdf5, write HDF5 file = ",fileout_h5
  call export_hdf5_restart(node_list, element_list, fileout_h5)

end program RST_convert_bin2hdf5
