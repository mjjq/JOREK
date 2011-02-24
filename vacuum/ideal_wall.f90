subroutine ideal_wall(my_id, node_list, boundary_list, bnd_node_list)
!-------------------------------------------------------------------
! routine calculates the vacuum reponse matrix of a vacuum region
! up to an ideally conducting wall
!-------------------------------------------------------------------
  
  use data_structure
  use vacuum_response
  use mumps_module
  
  implicit none
  
  include 'mpif.h'
  
  integer,                     intent(in) :: my_id            ! MPI thread number of current thread
  type(type_node_list),        intent(in) :: node_list        ! List of boundary nodes
  type(type_bnd_element_list), intent(in) :: boundary_list    ! List of boundary elements
  type(type_bnd_node_list),    intent(in) :: bnd_node_list    ! List of boundary nodes
  
  !### CURRENTLY NOT IMPLEMENTED ###

end subroutine ideal_wall
