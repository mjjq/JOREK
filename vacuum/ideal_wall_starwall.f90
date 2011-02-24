subroutine ideal_wall_starwall(my_id,node_list,bnd_elm_list,bnd_node_list)
!-------------------------------------------------------------------
! Reads the ideal wall vacuum response matrices written out by STARWALL
!-------------------------------------------------------------------

  use data_structure
  use vacuum_response
  use phys_module
  
  implicit none

  integer,                     intent(in) :: my_id            ! MPI thread number of current thread
  type(type_node_list),        intent(in) :: node_list        ! List of boundary nodes
  type(type_bnd_element_list), intent(in) :: bnd_elm_list     ! List of boundary elements
  type(type_bnd_node_list),    intent(in) :: bnd_node_list    ! List of boundary nodes

  real*8  :: TWOPI
  integer :: dim(2)

  TWOPI = 8.d0 * atan(1.d0)
  
  write(*,*) '@@> ideal_wall_starwall'
  
  dim(:) = response_index(bnd_node_list%n_bnd_nodes, n_tor, 2)
  call read_response_matrix( vac_response, dim, 'starwall_m_id' )
  vac_response(:,:) = vac_response(:,:) * TWOPI

  write(*,*) '@@< ideal_wall_starwall'

end subroutine ideal_wall_starwall
