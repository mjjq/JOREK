subroutine import_boundary(node_List,boundary_list) 

use data_structure 
implicit none

type(type_bnd_element_list) :: boundary_list
type(type_node_list)     :: node_list

integer :: i,iv1,iv2,idir1,idir2

write(*,*) '******************************'
write(*,*) '* import boundary            *'
write(*,*) '******************************'
 
open(22,file='boundary.txt')
 
read(*,*) boundary_list%n_bnd_elements, bnd_node_list%n_bnd_nodes

do i=1,boundary_list%n_bnd_elements
    iv1   = boundary_list%bnd_element(i)%vertex(1)
    iv2   = boundary_list%bnd_element(i)%vertex(2)
    idir1 = boundary_list%bnd_element(i)%direction(1,2)
    idir2 = boundary_list%bnd_element(i)%direction(2,2)
    read(22,'(3i6,12e16.8)') i,iv1,iv2, node_list%node(iv1)%x(1,1:2),          &
                                        node_list%node(iv1)%x(idir1,1:2),      &
                                        boundary_list%bnd_element(i)%size(1:2,2), &	
                                        node_list%node(iv2)%x(1,1:2),          &
                                        node_list%node(iv2)%x(idir2,1:2),      &
                                        boundary_list%bnd_element(i)%size(1:2,2) 
enddo

do i=1, bnd_node_list%n_bnd_nodes
  read(22,'(4i7)') i,bnd_node_list%bnd_node(i)%index_starwall,bnd_node_list%bnd_node(i)%n_dof
enddo

close(22)
return
end
