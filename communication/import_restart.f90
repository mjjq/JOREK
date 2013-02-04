!> Imports a restart file written out by the routine export_restart.
subroutine import_restart(node_list, element_list, filename, format_rst, error)

use tr_module 
use data_structure
use phys_module
use vacuum, only: import_restart_vacuum

implicit none

! --- Routine parameters
type(type_node_list),    intent(inout) :: node_list
type(type_element_list), intent(inout) :: element_list
character(len=*),        intent(in)    :: filename
integer,                 intent(out)   :: error
integer,                 intent(in)    :: format_rst  ! format of restart file

! --- Local variables
integer              :: i, j, m, k, n_tor_tmp
real*8               :: growth_mag, growth_kin, eta_tmp, visco_tmp, visco_par_tmp, amplitude
integer, allocatable :: mode_tmp(:)
real*8,  allocatable :: values_tmp(:,:,:), deltas_tmp(:,:,:)

error = 0

write(*,*) 'Importing restart file "', trim(filename), '".'
write(*,*) '  Using format : ',rst_format

open(21,file=trim(filename), form='unformatted', status='old', action='read', iostat=error)
if ( error /= 0 ) then
  write(*,*) '...failed!'
  return
end if

read(21) n_tor_tmp

allocate(mode_tmp(n_tor_tmp), values_tmp(n_tor_tmp,n_order+1,n_var), deltas_tmp(n_tor_tmp,n_order+1,n_var))

if (format_rst == 1) then
  read(21) mode_tmp
  write(*,*) ' NEW format (1) : ',mode_tmp
elseif (format_rst == 0) then
  write(*,*) ' mode : ',mode
  mode_tmp(1) = mode(1) ! temporary, remove
  if (n_tor .eq. n_tor_tmp) mode_tmp = mode
  write(*,*) ' OLD format (0) : '
  write(*,'(A,32i4)') ' previous modenumbers : ',mode_tmp
  write(*,'(A,32i4)') ' new mode numbers     : ',mode
else
  write(*,'(A,i3)') ' restart file format not supported : ',format_rst
endif

if (n_tor_tmp .gt. n_tor) write(*,'(3(a,i4))') ' IMPORT WARNING : Reducing number of harmonics from', n_tor_tmp, ' to', n_tor, '!'
if (n_tor_tmp .lt. n_tor) write(*,'(3(a,i4))') ' IMPORT WARNING : Increasing number of harmonics from', n_tor_tmp, ' to', n_tor, '!'

write(*,'(A,i5,A)') ' Importing ',n_tor_tmp,' harmonics'

read(21) node_list%n_nodes,element_list%n_elements
read(21) node_list%n_dof

do i=1,node_list%n_nodes

  read(21) node_list%node(i)%x

  read(21) values_tmp
  read(21) deltas_tmp

#ifdef fullmhd
  read(21) node_list%node(i)%psi_eq               !< equilibrium flux at the nodes
  read(21) node_list%node(i)%Fprof_eq             !< equilibrium profile R*B_phi at the nodes
#endif
  read(21) node_list%node(i)%index
  read(21) node_list%node(i)%boundary
  read(21) node_list%node(i)%parents
  read(21) node_list%node(i)%parent_elem
  read(21) node_list%node(i)%ref_lambda
  read(21) node_list%node(i)%ref_mu
  read(21) node_list%node(i)%constrained

  node_list%node(i)%values = 0.d0
  node_list%node(i)%deltas = 0.d0

  do m=1,n_tor_tmp,2
    do k=1, n_tor,2
      if (mode_tmp(m) .eq. mode(k)) then
        if ((m .eq. 1) .and. (k.eq.1)) then
          node_list%node(i)%values(k,:,:) = values_tmp(m,:,:)
          node_list%node(i)%deltas(k,:,:) = deltas_tmp(m,:,:)
        else
          node_list%node(i)%values(k-1,:,:) = values_tmp(m-1,:,:)
          node_list%node(i)%deltas(k-1,:,:) = deltas_tmp(m-1,:,:)
          node_list%node(i)%values(k,:,:)   = values_tmp(m,:,:)
          node_list%node(i)%deltas(k,:,:)   = deltas_tmp(m,:,:)
        endif
      endif
    enddo
  enddo

enddo

read(21) element_list%element(1:element_list%n_elements)
read(21) tstep,eta_tmp,visco_tmp,visco_par_tmp
read(21) index_start
read(21) t_start
 
#ifdef USE_HDF5
  read(21) h5_nbsave_all
#endif

if (index_start .ge. 1) then

  if (allocated(xtime)) call tr_deallocate(xtime,"xtime",CAT_UNKNOWN)
  call tr_allocate(xtime,1,index_start+nstep,"xtime",CAT_UNKNOWN)

  if (allocated(energies)) call tr_deallocate(energies,"energies",CAT_UNKNOWN)
  call tr_allocate(energies,1,n_tor,1,2,1,index_start+nstep,"energies",CAT_UNKNOWN)

  energies = 0.d0

  read(21) xtime(1:index_start)
  read(21) energies(1:n_tor_tmp,:,1:index_start)

endif

if (use_pellet) then
  read(21)  pellet_particles, pellet_R, pellet_Z
endif

call import_restart_vacuum(21, freeboundary, resistive_wall, index_start)

close(21)

write(*,*) '************* restart ******************'
write(*,'(A19,i6,f14.6,A)') ' *  restart time : ',index_start,t_start,' *'
#ifdef USE_HDF5
  write(*,'(A19,f14.6,A)') ' *  HDF5 files written : ',h5_nbsave_all,' *'
#endif
write(*,*) '****************************************'


do i=2,index_start
if ( (energies(n_tor,1,i).ne.0.) .and. (energies(n_tor,1,i-1).ne.0.)) then
   Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,i)/energies(n_tor,1,i-1))) &
        / (xtime(i)-xtime(i-1))
else
    Growth_mag  = 0.
endif
if ( (energies(n_tor,2,i).ne.0.) .and. (energies(n_tor,2,i-1).ne.0.)) then
 Growth_kin  = 0.5d0*log(abs(energies(n_tor,2,i)/energies(n_tor,2,i-1))) &
             / (xtime(i)-xtime(i-1))
else
    Growth_kin  = 0.
endif

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
end subroutine import_restart

