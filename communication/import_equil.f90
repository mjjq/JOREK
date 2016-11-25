subroutine import_equil(node_list,element_list)
!**************************************************************************
! import the grid and values of the nodes/elements                        *
!**************************************************************************
use mod_parameters
use data_structure

implicit none

type(type_node_list)    :: node_list
type(type_element_list) :: element_list
integer :: n_dim_tmp, n_var_tmp, i, j

write(*,*) '****************************************'
write(*,*) '* import grid/values                   *'
write(*,*) '****************************************'

open(22,file='bezier.grid',form='unformatted')

read(22) n_dim_tmp, n_var_tmp
read(22) node_list%n_nodes

do i=1, node_list%n_nodes
 read(22) node_list%node(i)%x
 read(22) node_list%node(i)%boundary
 node_list%node(i)%values = 0.d0
 do j=1, n_var_tmp
   read(22) node_list%node(i)%values(:,:,j)
 enddo
enddo
read(22) element_list%n_elements
read(22) element_list%element(1:element_list%n_elements)

close(22)

write(*,*) ' n_nodes    : ',node_list%n_nodes
write(*,*) ' n_elements : ',element_list%n_elements
write(*,*) ' n_var      : ',n_var_tmp
write(*,*) '*********** end import *****************'

return
end