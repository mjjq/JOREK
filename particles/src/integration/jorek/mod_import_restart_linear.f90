module mod_import_restart_linear
contains
!> At step i_step, look for next file and use the appropriate method to import
!> The next file can be up to 10 istep away, and if it is more than 2 away
!> deltas are recalculated to fix the import
subroutine import_next_restart(node_list,element_list, istep, istep_out, rst_format)
  use data_structure
  use mpi
  implicit none

  ! --- Routine parameters
  type(type_node_list),    intent(inout) :: node_list !< previously imported node list at istep
  type(type_element_list), intent(inout) :: element_list !< Previously imported element list at istep
  integer,                 intent(in)    :: istep !< Current timestep
  integer,                 intent(out)   :: istep_out !< timestep of next restart file found
  integer,                 intent(in)    :: rst_format !< Format of the restart file

  ! --- Internal variables
  character*17 :: restart_file
  logical :: file_exists
  integer :: i, ierr

  ! Check if file for next iteration exists
  write(restart_file,'(A,i5.5,A)') 'jorek', istep+1, '.rst'
  inquire(file=restart_file, exist=file_exists)
  if (file_exists) then
    ! If so, import it and we're done
    call import_binary_restart(node_list,element_list,restart_file,rst_format,ierr)
    istep_out = istep+1
  else
    ! If not, keep looping (up to 10) to find one, and use the merge import
    ! This assumes that the current node_list contains the values
    ! at time istep (but does not need to contain the deltas, these are calculated)
    do i=istep+2,istep+10
      write(restart_file,'(A,i5.5,A)') 'jorek', i, '.rst'
      inquire(file=restart_file, exist=file_exists)
      if (file_exists) then
        call import_merge_restart(node_list,element_list,restart_file,rst_format,ierr)
        istep_out=i
        exit ! the loop
      endif
    enddo
  endif

  if (ierr .ne. 0) then
    write(*,*) "Error reading restart file", restart_file
    call MPI_ABORT(MPI_COMM_WORLD, ierr)
  endif
end subroutine import_next_restart



!> Import a binary restart file and merges it with the values currently known
!> This can then be used to interpolate linearly between any two restart files
subroutine import_merge_restart(node_list,element_list, restart_file, format_rst, ierr)
  use data_structure
  use phys_module
  implicit none

  ! --- Routine parameters
  type(type_node_list),    intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  character(len=*),        intent(in)    :: restart_file !< Filename of new restart file to import
  integer,                 intent(out)   :: ierr
  integer,                 intent(in)    :: format_rst !< Restart file format

  ! --- Internal variables
  real*8, allocatable, dimension(:,:,:,:) :: values
  integer :: inode
  real*8 :: tstart_old

  ! Save the old values to calculate the new deltas
  allocate(values(n_tor,n_order+1,n_var,node_list%n_nodes))
  do inode=1,node_list%n_nodes
    values(:,:,:,inode) = node_list%node(inode)%values(:,:,:)
  enddo
  tstart_old = t_start
  
  ! Import new values
  call import_binary_restart(node_list,element_list, restart_file, format_rst, ierr)

  ! Calculate deltas as values_new - values_old
  do inode=1,node_list%n_nodes
    node_list%node(inode)%deltas = node_list%node(inode)%values - values(:,:,:,inode)
  enddo

  ! Set timestep to time between restart files
  tstep = t_start - tstart_old

  deallocate(values)
end subroutine import_merge_restart
end module mod_import_restart_linear
