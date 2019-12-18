!> Writes out the q-profile
subroutine q_profile(node_list,element_list,surface_list,psi_axis,psi_bnd,psi_xpoint,Z_xpoint)

use tr_module 
use data_structure
use phys_module

implicit none

! --- Gaussian points between (-1.,1.) for Gauss-integration
real*8, parameter :: xgs(4) = (/-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053 /)
real*8, parameter :: wgs(4) = (/ 0.347854845137454,  0.652145154862546, 0.652145154862546,  0.347854845137454 /)

! --- Input parameters.
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
type (type_surface_list), intent(in)    :: surface_list
real*8,                   intent(in)    :: psi_axis
real*8,                   intent(in)    :: psi_bnd        !< Psi value at plasma boundary (x-point or limiter)
real*8,                   intent(in)    :: psi_xpoint(2)
real*8,                   intent(in)    :: Z_xpoint(2)

! --- Local variables
real*8, allocatable :: q(:),rad(:)
integer :: i

write(*,*) '**********************************'
write(*,*) '*        q-profile               *'
write(*,*) '**********************************'

call tr_allocate(q,1,surface_list%n_psi,"q",CAT_GRID)
call tr_allocate(rad,1,surface_list%n_psi,"rad",CAT_GRID)

call determine_q_profile(node_list,element_list,surface_list,psi_axis,psi_xpoint,Z_xpoint,q,rad)

if ( write_ps ) call lplot(2,1,1,surface_list%psi_values(2),q(2),surface_list%n_psi-1,1,'q-profile',9,'flux',4,'q',1)

! --- Write out the q-profile versus the poloidal flux to "qprofile.dat".
! data written is psi_n, q, and the flux surface average of minor radius=sqrt((R-R_geo)^2+(Z-Z_geo)^2)
open(42, file='qprofile.dat', action='write', status='replace')
do i=2, surface_list%n_psi
   write(42,*) (surface_list%psi_values(i)-psi_axis)/(psi_bnd-psi_axis), q(i),rad(i)
end do
close(42)

call tr_deallocate(q,"q",CAT_GRID)
call tr_deallocate(rad,"rad",CAT_GRID)

return
end subroutine q_profile
