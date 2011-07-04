subroutine export_restart(node_list,element_list,filename)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
character*(*)            :: filename
integer                  :: i

open(21,file=filename,form='unformatted')

write(21) n_tor
write(21) node_list%n_nodes,element_list%n_elements
write(21) node_list%n_dof

do i=1,node_list%n_nodes
  write(21) node_list%node(i)%x
  write(21) node_list%node(i)%values
  write(21) node_list%node(i)%deltas
  write(21) node_list%node(i)%index
  write(21) node_list%node(i)%boundary
  write(21) node_list%node(i)%parents		     
  write(21) node_list%node(i)%parent_elem			      
  write(21) node_list%node(i)%ref_lambda
  write(21) node_list%node(i)%ref_mu		     
  write(21) node_list%node(i)%constrained
enddo

write(21) element_list%element(1:element_list%n_elements)
write(21) tstep,eta,visco,visco_par
write(21) index_now
write(21) t_now
if (index_now .gt. 0) then
  write(21) xtime(1:index_now)
  write(21) energies(:,:,1:index_now)
endif
if (use_pellet) then
  write(21) pellet_particles, pellet_R, pellet_Z
endif

close(21)

return
end

