module mod_fix_axis_nodes
contains

!*******************************************************************************
!* Subroutine: fix_axis_nodes                                                  *
!*******************************************************************************
!*                                                                             *
!* Add condition for the axis directly in the matrix.                          *
!* This is aimed at stabilising numerical noise on the axis.                   *
!*                                                                             *
!* Parameters:                                                                 *
!*   node_list    - List of nodes                                              *
!*   element_list - List of all elements                                       *
!*   local_elms   - List of local elements                                     *
!*   n_local_elms - Number of local elements                                   *
!*   index_min    - Minimal index of local elements                            *
!*   index_max    - Maximal index of local elements (                          *
!*                                                                             *
!*******************************************************************************
subroutine fix_nodes_on_axis(node_list, element_list, local_elms, n_local_elms, index_min, index_max)

  use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
  use data_structure
  use global_distributed_matrix
  use mod_locate_irn_jcn

  implicit none

  ! Subroutine parameters
  integer                  :: local_elms(*)
  integer                  :: n_local_elms
  integer                  :: index_min, index_max
  type (type_node_list)    :: node_list
  type (type_element_list) :: element_list

  ! Internal parameters
  real*8  :: zbig
  integer :: i, in, iv, inode, k
  integer :: index_large_i, index_node, index_node2, ielm
  integer :: ijA_position,ijA_position2, ilarge2


  zbig = 1.d12
  do i=1, n_local_elms

    ielm = local_elms(i)

    do iv=1, n_vertex_max

      inode = element_list%element(ielm)%vertex(iv)

      if (node_list%node(inode)%axis_node) then

        do in=1, n_tor
          do k=1, n_var

            ! --- For t-derivative
            index_node = node_list%node(inode)%index(3)
            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
              call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
              index_large_i = n_tor * n_var * (index_node - 1)
              ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in
              irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
              jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
              A_glob(ilarge2)   = zbig
            end if

            ! --- For cross st-derivative
            index_node = node_list%node(inode)%index(4)
            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
              call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
              index_large_i = n_tor * n_var * (index_node - 1)
              ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in
              irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
              jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
              A_glob(ilarge2)   = zbig
            end if

          enddo
        enddo

      endif

    enddo  ! n_vertex
  enddo ! n_elements

  return

end subroutine fix_nodes_on_axis

end module mod_fix_axis_nodes
