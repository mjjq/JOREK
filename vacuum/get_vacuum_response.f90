subroutine get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list, coils)
!-------------------------------------------------------------------
! Determines the vacuum response for an ideal or resistive wall
!-------------------------------------------------------------------

  use data_structure
  use phys_module
  use vacuum_response_module
  
  implicit none
    
  include 'mpif.h'
  
  integer,                     intent(in) :: my_id            ! MPI thread number of current thread
  type(type_node_list),        intent(in) :: node_list        ! List of boundary nodes
  type(type_bnd_element_list), intent(in) :: bnd_elm_list     ! List of boundary elements
  type(type_bnd_node_list),    intent(in) :: bnd_node_list    ! List of boundary nodes
  logical                                 :: coils            ! import coil contributions or not
  integer :: ierr
  
  
  n_dof_bnd = 2*bnd_node_list%n_bnd_nodes                     ! degrees of freedom on the boundary
  
  write(*,*) '************************************'
  write(*,*) '*       get_vacuum_response        *'
  write(*,*) '************************************'
  
  ! --- Output some information about the boundary.
  230 format(A,' = ',I8)
  231 format(A,' = ',L8)
  write(*,230) 'n_bnd_elements',bnd_elm_list%n_bnd_elements
  write(*,230) 'n_bnd_nodes   ',bnd_node_list%n_bnd_nodes
  write(*,230) 'n_dof_bnd     ', n_dof_bnd
  write(*,231) 'resistive_wall', resistive_wall
  if ( .not. resistive_wall ) write(*,231) 'use_starwall  ', use_starwall
  
  ! --- Write out the boundary information for STARWALL.
  if (my_id .eq. 0) call export_boundary(node_list, bnd_elm_list, bnd_node_list)
    
  ! --- Resistive wall
  if ( resistive_wall ) then
    
    ! --- Get the STARWALL response matrices
    call resistive_wall_starwall(my_id,node_list,bnd_elm_list,bnd_node_list)
  
  ! --- Ideal wall 
  else
    
    ! --- Allocate the ideal wall vacuum response matrix
    if ( allocated(vacuum_response) ) deallocate(vacuum_response)
    allocate(vacuum_response(n_dof_bnd,n_dof_bnd,n_tor))
    vacuum_response = 0.
    
    ! --- Get the STARWALL response matrix
    if ( use_starwall ) then
      call ideal_wall_starwall(my_id,node_list,bnd_elm_list,bnd_node_list)
      
    ! --- Let JOREK determine the response (works only in special cases; for testing)
    else
!      call ideal_wall(my_id,node_list,bnd_elm_list,bnd_node_list)      
    end if
    
    ! --- Send the vacuum response matrix.
    call MPI_bcast(vacuum_response, n_dof_bnd*n_dof_bnd*n_tor, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
  end if
  
  if (coils) call import_external_fields()
    
end subroutine get_vacuum_response
