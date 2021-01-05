!> Add condition for the axis directly in the matrix.
!> This is aimed at stabilising numerical noise on the axis.
module mod_fix_axis_nodes
contains

subroutine fix_nodes_on_axis(node_list, element_list, local_elms, n_local_elms, index_min, index_max, & 
  ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat, i_tor_min, i_tor_max )

  use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
  use data_structure
  use mod_locate_irn_jcn
  use mod_integer_types

  implicit none

  ! Subroutine parameters
  integer,                            intent(in)    :: local_elms(*)         !< List of local elements
  integer,                            intent(in)    :: n_local_elms          !< Number of local elements
  integer,                            intent(in)    :: index_min, index_max  !< Min/max index of local elements
  type (type_node_list),              intent(in)    :: node_list             !< List of nodes
  type (type_element_list),           intent(in)    :: element_list          !< List of all elements
  integer(kind=int_all), allocatable, intent(in)    :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:)
  integer,                            intent(in)    :: i_tor_min, i_tor_max
  integer(kind=int_all), allocatable, intent(inout) :: irn(:), jcn(:)
  real*8,                allocatable, intent(inout) :: A_mat(:)
  ! Internal parameters
  real*8                :: zbig
  integer               :: i, in, iv, inode, k, ielm, ilarge2, n_tor_local
  integer               :: index_node, index_node2
  integer(kind=int_all) :: index_large_i
  integer(kind=int_all) :: ijA_position,ijA_position2

  n_tor_local = (i_tor_max - i_tor_min + 1)

  zbig = 1.d12
  do i=1, n_local_elms

    ielm = local_elms(i)

    do iv=1, n_vertex_max

      inode = element_list%element(ielm)%vertex(iv)

      if (node_list%node(inode)%axis_node) then

        do in=i_tor_min, i_tor_max
          do k=1, n_var

            ! --- For t-derivative
            index_node = node_list%node(inode)%index(3)
            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
              call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
              index_large_i = n_tor_local * n_var * (index_node - 1)
              ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                + (k-1)*n_tor_local + in - i_tor_min + 1
              irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
              jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
              A_mat(ilarge2)   = zbig
            end if

            ! --- For cross st-derivative
            index_node = node_list%node(inode)%index(4)
            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
              call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
              index_large_i = n_tor_local * n_var * (index_node - 1)
              ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                + (k-1)*n_tor_local + in - i_tor_min + 1
              irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
              jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
              A_mat(ilarge2)   = zbig
            end if

          enddo
        enddo

      endif

    enddo  ! n_vertex
  enddo ! n_elements

  return

end subroutine fix_nodes_on_axis

!> Add condition for the axis directly in the matrix.
!> This is aimed at applying C0 continuity on the grid axis.
subroutine penalize_third_dof_on_axis(node_list, element_list, local_elms, n_local_elms, index_min, index_max, &
  ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat, i_tor_min, i_tor_max )

  use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
  use data_structure
  use mod_locate_irn_jcn

  implicit none

  ! Subroutine parameters
  integer,                   intent(in)    :: local_elms(*)         !< List of local elements
  integer,                   intent(in)    :: n_local_elms          !< Number of local elements
  integer,                   intent(in)    :: index_min, index_max  !< Min/max index of local elements
  type (type_node_list),     intent(in)    :: node_list             !< List of nodes
  type (type_element_list),  intent(in)    :: element_list          !< List of all elements
  integer, allocatable,      intent(in)    :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:)
  integer,                   intent(in)    :: i_tor_min, i_tor_max
  integer, allocatable,      intent(inout) :: irn(:), jcn(:)
  real*8,  allocatable,      intent(inout) :: A_mat(:)
  ! Internal parameters
  real*8  :: zbig
  integer :: i, in, iv, inode, k
  integer :: index_large_i, index_node, index_node2, ielm
 integer :: ijA_position,ijA_position2, ilarge2, n_tor_local

  n_tor_local = (i_tor_max - i_tor_min + 1)

  zbig = 1.d12
  do i=1, n_local_elms

    ielm = local_elms(i)

    do iv=1, n_vertex_max

      inode = element_list%element(ielm)%vertex(iv)

      if (node_list%node(inode)%axis_node) then

        do in=i_tor_min, i_tor_max
          do k=1, n_var

            index_node = node_list%node(inode)%index(3)
            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
              call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
              index_large_i = n_tor_local * n_var * (index_node - 1)
              ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                + (k-1)*n_tor_local + in - i_tor_min + 1
              irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
              jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
              A_mat(ilarge2)   = zbig
            end if

          enddo
        enddo

      endif

    enddo  ! n_vertex
  enddo ! n_elements

  return

end subroutine penalize_third_dof_on_axis

end module mod_fix_axis_nodes
