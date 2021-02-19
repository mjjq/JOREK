subroutine identify_axis_elements(node_list,element_list)
use data_structure
use phys_module, only: treat_axis, n_flux

implicit none
integer :: ie, iv, j
type(type_node_list),    intent(inout) :: node_list
type(type_element_list), intent(inout) :: element_list

do ie = 1, element_list%n_elements
  do iv = 1, 4 
    j = element_list%element(ie)%vertex(iv)
    element_list%element(ie)%axis_element = .false.
    if (treat_axis .and. node_list%node(j)%axis_node) then
       element_list%element(ie)%axis_element = .true.
    endif
  enddo
enddo

end subroutine identify_axis_elements
