subroutine import_restart(node_list,element_list)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use tr_module 
use data_structure
use phys_module
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
integer :: i, j, n_tor_tmp
real*8  :: growth_mag, growth_kin, eta_tmp, visco_tmp, visco_par_tmp, amplitude

open(21,file='jorek_restart.rst', form='unformatted', status='old', action='read')

read(21) n_tor_tmp

if (n_tor_tmp .gt. n_tor) write(*,*) ' IMPORT WARNING : reducing number of harmonics!', n_tor, n_tor_tmp
if (n_tor_tmp .lt. n_tor) write(*,*) ' IMPORT WARNING : increasing number of harmonics!',n_tor,n_tor_tmp

n_tor_tmp = min(n_tor,n_tor_tmp)

write(*,'(A,i5,A)') ' importing ',n_tor_tmp,' harmonics'

read(21) node_list%n_nodes,element_list%n_elements
read(21) node_list%n_dof

do i=1,node_list%n_nodes
  read(21) node_list%node(i)%x
  read(21) node_list%node(i)%values(1:n_tor_tmp,:,:)
  read(21) node_list%node(i)%deltas(1:n_tor_tmp,:,:)
  read(21) node_list%node(i)%index
  read(21) node_list%node(i)%boundary
  read(21) node_list%node(i)%parents			     
  read(21) node_list%node(i)%parent_elem			      
  read(21) node_list%node(i)%ref_lambda
  read(21) node_list%node(i)%ref_mu		     
  read(21) node_list%node(i)%constrained
enddo

read(21) element_list%element(1:element_list%n_elements)
read(21) tstep,eta_tmp,visco_tmp,visco_par_tmp
read(21) index_start
read(21) t_start

write(*,*) ' IMPORT CHECK nstep = ',nstep

if (index_start .ge. 1) then

  if (allocated(xtime)) call tr_deallocate(xtime,"xtime")
  call tr_allocate(xtime,1,index_start+nstep,"xtime")

  if (allocated(energies)) call tr_deallocate(energies,"energies")
  call tr_allocate(energies,1,n_tor,1,2,1,index_start+nstep,"energies")

  energies = 0.d0

  read(21) xtime(1:index_start)
  read(21) energies(1:n_tor_tmp,:,1:index_start)

endif

if (use_pellet) then
  read(21)  pellet_particles, pellet_R, pellet_Z
endif

close(21)

write(*,*) '************* restart ******************'
write(*,'(A19,i6,f14.6,A)') ' *  restart time : ',index_start,t_start,' *'
write(*,*) '****************************************'


do i=2,index_start

 Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,i)/energies(n_tor,1,i-1))) &
             / (xtime(i)-xtime(i-1))
 Growth_kin  = 0.5d0*log(abs(energies(n_tor,2,i)/energies(n_tor,2,i-1))) &
             / (xtime(i)-xtime(i-1))

! write(*,'(i7,f10.3,200e14.6)') i,xtime(i),energies(1:n_tor,:,i),growth_mag,growth_kin
! write(*,'(i7,f10.3,200e14.6)') i,xtime(i),energies(1:n_tor,:,i)

enddo

amplitude = 1.d-10

if (n_tor_tmp .lt. n_tor) then ! initialise new harmonics (only density and temperature, to be improved)
  do i=1,node_list%n_nodes
    node_list%node(i)%values(n_tor_tmp+1:n_tor,:,1:4)= 0.d0
    do j=n_tor_tmp+1, n_tor
      node_list%node(i)%values(j,:,5)= amplitude * node_list%node(i)%values(1,:,5)
      node_list%node(i)%values(j,:,6)= amplitude * node_list%node(i)%values(1,:,6)
    enddo
  enddo
endif

write(*,*) '********* end restart ******************'

!call add_pellet(node_list,element_list,25.d0,0.06d0,0.03d0,3.78d0,0.14d0)

return
end
