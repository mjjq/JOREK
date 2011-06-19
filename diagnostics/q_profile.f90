subroutine q_profile(node_list,element_list,surface_list,psi_axis,psi_xpoint,Z_xpoint)
!----------------------------------------------------------------------
! subroutine calculats the q-profile from the flux surface representation
! (adapted from helena20)
!---------------------------------------------------------------------
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
real*8,                   intent(in)    :: psi_xpoint
real*8,                   intent(in)    :: Z_xpoint

! --- Internal variables
integer :: n_int
integer :: i_elm, j, k, n1, n2, n3
real*8  :: t,rr1, rr2, drr1, drr2, ss1, ss2, dss1, dss2, ri, si, dri, dsi, dl
real*8  :: RRgi, dRRgi_dr, dRRgi_ds, ZZgi, dZZgi_dr, dZZgi_ds, dRRgi_dt, dZZgi_dt
real*8  :: PSgi, dPSgi_dr, dPSgi_ds, PSI_R, PSI_Z, RZJAC, grad_psi, psi_n
real*8  :: sum_dl, B_tot2, PI
real*8  :: dRRgi_drs,dRRgi_drr,dRRgi_dss, dZZgi_drs,dZZgi_drr,dZZgi_dss, dPSgi_drs,dPSgi_drr,dPSgi_dss
real*8, allocatable :: q(:)
integer :: i,m, ig, ip

write(*,*) '**********************************'
write(*,*) '*        q-profile               *'
write(*,*) '**********************************'

call tr_allocate(q,1,surface_list%n_psi,"q")

PI = 2.d0 *asin(1.d0)

write(*,*) '   i     psi           q          sum_dl'

do i=2, surface_list%n_psi

  q(i)   = 0.d0
  sum_dl = 0.d0

  do k=1, surface_list%flux_surfaces(i)%n_pieces

    do ig = 1, 4

      t = xgs(ig)

      rr1  = surface_list%flux_surfaces(i)%s(1,k)
      drr1 = surface_list%flux_surfaces(i)%s(2,k)
      rr2  = surface_list%flux_surfaces(i)%s(3,k)
      drr2 = surface_list%flux_surfaces(i)%s(4,k)

      ss1  = surface_list%flux_surfaces(i)%t(1,k)
      dss1 = surface_list%flux_surfaces(i)%t(2,k)
      ss2  = surface_list%flux_surfaces(i)%t(3,k)
      dss2 = surface_list%flux_surfaces(i)%t(4,k)

      call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
      call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

      i_elm = surface_list%flux_surfaces(i)%elm(k)

      call interp(node_list,element_list,i_elm,1,1,ri,si,PSgi,dPSgi_dr,dPSgi_ds,dPSgi_drs,dPSgi_drr,dPSgi_dss)

      call interp_RZ(node_list,element_list,i_elm,ri,si,RRgi,dRRgi_dr,dRRgi_ds,dRRgi_drs,dRRgi_drr,dRRgi_dss, &
                                                        ZZgi,dZZgi_dr,dZZgi_ds,dZZgi_drs,dZZgi_drr,dZZgi_dss)
                                                        
      ! --- Make sure that for flux surfaces at Psi_N < 1, the surface integral is carried out only
      !     over the flux surface segments of the plasma region.
      !     I.e., ignore flux surface segments in the private flux region below the x-point.
      if ( xpoint .and. ((PSgi < psi_xpoint) .and. (ZZgi < z_xpoint)) ) cycle

      dRRgi_dt = dRRgi_dr * dri + dRRgi_ds * dsi
      dZZgi_dt = dZZgi_dr * dri + dZZgi_ds * dsi

      dl = sqrt(dRRgi_dt**2 + dZZgi_dt**2)

      RZjac  = DRRgi_dr * dZZgi_ds - dRRgi_ds * dZZgi_dr

      PSI_R = (   dPSgi_dr * dZZgi_ds - dPSgi_ds * dZZgi_dr ) / RZjac
      PSI_Z = ( - dPSgi_dr * dRRgi_ds + dPSgi_ds * dRRgi_dr ) / RZjac

      grad_psi = sqrt(PSI_R * PSI_R + PSI_Z * PSI_Z)

      B_tot2 =  (F0 / RRgi)**2 + (grad_psi/ RRgi)**2

      sum_dl = sum_dl +  wgs(ig) * dl

      q(i) = q(i) +  wgs(ig) / (RRgi * grad_psi) * dl

    enddo

  enddo

  q(i) = F0 * q(i) / (2.d0 * PI)

  write(*,'(i5,3es13.5)') i, surface_list%psi_values(i), q(i), sum_dl

enddo

! --- Write out the q-profile to "qprofile.dat".
open(42, file='qprofile.dat', action='write', status='replace')
do i=2, surface_list%n_psi
  write(42,'(2ES13.5)') (surface_list%psi_values(i)-psi_axis)/(psi_xpoint-psi_axis), q(i)
end do
close(42)


!----------------------------------- values on axis
!q(1) =  PI / sqrt(CRR_axis*CZZ_Axis)

call lplot(2,1,1,surface_list%psi_values(2),q(2),surface_list%n_psi-1,1,'q-profile',9,'flux',4,'q',1)

call tr_deallocate(q,"q")

return
end subroutine q_profile
