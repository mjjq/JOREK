subroutine export_grid(node_list,element_list)
!**************************************************************************
! export the grid and values of the nodes/elements                        *
!**************************************************************************
use mod_parameters
use data_structure

implicit none

type(type_node_list)    :: node_list
type(type_element_list) :: element_list
integer :: i, j

write(*,*) '****************************************'
write(*,*) '* export grid/values                   *'
write(*,*) '****************************************'

open(21,file='bezier.grid',form='unformatted')

write(21) n_dim, n_var
write(21) node_list%n_nodes

do i=1, node_list%n_nodes
 write(21) node_list%node(i)%x
 write(21) node_list%node(i)%boundary
 do j=1, n_var
   write(21) node_list%node(i)%values(:,:,j)
 enddo
enddo
write(21) element_list%n_elements
write(21) element_list%element(1:element_list%n_elements)

close(21)

write(*,*) ' n_nodes    : ',node_list%n_nodes
write(*,*) ' n_elements : ',element_list%n_elements
write(*,*) ' n_var      : ',n_var
write(*,*) '*********** end export *****************'

return
end
