!> Export the current simulation state as a restart file that can be read back into JOREK or into
!! a diagnostic program by the routine import_restart.
subroutine export_restart(node_list,element_list,filename)

use data_structure
use phys_module
use vacuum, only: export_restart_vacuum

implicit none

! --- Routine parameters
type(type_node_list),    intent(in)    :: node_list
type(type_element_list), intent(in)    :: element_list
character*(*),           intent(in)    :: filename

! --- Local variables
integer :: i

open(21, file=filename, form='unformatted', status='replace', action='write')

write(21) n_tor
write(21) node_list%n_nodes,element_list%n_elements
write(21) node_list%n_dof

do i=1,node_list%n_nodes
  write(21) node_list%node(i)%x
  write(21) node_list%node(i)%values
  write(21) node_list%node(i)%deltas
#ifdef fullmhd
  write(21) node_list%node(i)%psi_eq               !< equilibrium flux at the nodes
  write(21) node_list%node(i)%Fprof_eq             !< equilibrium profile R*B_phi at the nodes
#endif
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

#ifdef USE_HDF5
  write(21) h5_nbsave_all
#endif

if (index_now .gt. 0) then
  write(21) xtime(1:index_now)
  write(21) energies(:,:,1:index_now)
endif

if (use_pellet) then
  write(21) pellet_particles, pellet_R, pellet_Z
endif

call export_restart_vacuum(21, freeboundary, resistive_wall, index_now)

close(21)

return
end subroutine export_restart

