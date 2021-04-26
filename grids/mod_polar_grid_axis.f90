subroutine identify_axis_elements(node_list,element_list)
use data_structure
use phys_module, only: treat_axis, treat_axis2

implicit none
integer :: ie, iv, j1, j2 ,j3 ,j4
type(type_node_list),    intent(inout) :: node_list
type(type_element_list), intent(inout) :: element_list

do ie = 1, element_list%n_elements
  j1 = element_list%element(ie)%vertex(1)
  j2 = element_list%element(ie)%vertex(2)
  j3 = element_list%element(ie)%vertex(3)
  j4 = element_list%element(ie)%vertex(4)

  element_list%element(ie)%axis_element = .false.

  ! The first and fourth vertex is on the grid-axis.
  if ( (treat_axis .or. treat_axis2) .and. ( node_list%node(j1)%axis_node .and. node_list%node(j4)%axis_node) ) then
     element_list%element(ie)%axis_element = .true.
  endif
enddo

end subroutine identify_axis_elements
