subroutine resistive_wall_starwall(my_id, node_list, boundary_list)

  use data_structure
  use vacuum_response_module
  use phys_module
  
  
  implicit none


  integer,                     intent(in) :: my_id            ! MPI thread number of current thread
  type(type_node_list),        intent(in) :: node_list        ! List of boundary nodes
  type(type_bnd_element_list), intent(in) :: boundary_list    ! List of boundary elements


  character(len=128) :: file_response_starwall = 'vacuum_response_starwall'

  integer :: dim
  integer :: i_n, j_n
  integer :: itor, jtor
  integer :: inode, jnode
  integer :: ibas, jbas
  integer :: icossin, jcossin ! 0: cos, 1: sin
  integer :: iindex, jindex
  integer :: rn_response, cn_response
  integer :: iindex2, jindex2
  real*8  :: response
  real*8  :: TWOPI
  
  TWOPI=8.*atan(1.)

  open(42, FILE='###')
  
  ! READ MATRICES

  close(42)
  
  !####

end subroutine resistive_wall_starwall

