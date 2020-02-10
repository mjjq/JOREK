subroutine find_theta_surface(node_list,element_list,surface_list,j_surf,theta,R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)
!---------------------------------------------------------------------------
! subroutine finds a theta value on a specific surface
!---------------------------------------------------------------------------

use constants
use data_structure
use mod_interp, only: interp_RZ

implicit none

! --- Routine parameters
type (type_node_list),    intent(in)    :: node_list       !< node list
type (type_element_list), intent(in)    :: element_list    !< element list
type (type_surface_list), intent(in)    :: surface_list    !< surface list
integer,                  intent(in)    :: j_surf          !< index of surface we are intersecting in surface_list
real*8,                   intent(in)    :: theta           !< theta-angle of the line we are intersecting surface with
real*8,                   intent(in)    :: R_axis          !< R_axis
real*8,                   intent(in)    :: Z_axis          !< Z_axis
integer,                  intent(out)   :: i_elm_find(*)   !< the elm indices where we found intersections
real*8,                   intent(out)   :: s_find(*)       !< the local s-value of each element where we found intersections
real*8,                   intent(out)   :: t_find(*)       !< the local t-value of each element where we found intersections
integer,                  intent(out)   :: i_find          !< number of intersections we found

! --- local variables
integer :: i, k, i_elm, ifail, iterMax
real*8  :: ri, dri, si ,dsi, tht_in
real*8  :: rr1, drr1, rr2, drr2, ss1, dss1, ss2, dss2, t, t2, t3, tht1, tht2, tht_out
real*8  :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8  :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8  :: RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss
real*8  :: ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss
real*8  :: dRRg1_dt, dZZg1_dt, dRRg2_dt, dZZg2_dt, RZ1, RZ2, dRZ1, dRZ2, RZ0, A0, A1, A2, A3

!write(*,'(A,2e16.8)') ' finding THETA on surface : ',surface_list%psi_values(j_surf),theta
!write(*,'(A,2e16.8)') ' R_axis, Z_axis : ',R_axis, Z_axis

tht_in = theta
if (tht_in .lt. 0.d0)    tht_in = tht_in + 2.d0 * PI
if (tht_in .gt. 2.d0*PI) tht_in = tht_in - 2.d0 * PI

i_find = 0
do k=1,surface_list%flux_surfaces(j_surf)%n_pieces

  rr1  = surface_list%flux_surfaces(j_surf)%s(1,k);   ss1  = surface_list%flux_surfaces(j_surf)%t(1,k)
  drr1 = surface_list%flux_surfaces(j_surf)%s(2,k);   dss1 = surface_list%flux_surfaces(j_surf)%t(2,k)
  rr2  = surface_list%flux_surfaces(j_surf)%s(3,k);   ss2  = surface_list%flux_surfaces(j_surf)%t(3,k)
  drr2 = surface_list%flux_surfaces(j_surf)%s(4,k);   dss2 = surface_list%flux_surfaces(j_surf)%t(4,k)

  i_elm = surface_list%flux_surfaces(j_surf)%elm(k)

  call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                      ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
  call interp_RZ(node_list,element_list,i_elm,rr2,ss2,RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss, &
                                                      ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss)

  tht1 = atan2(ZZg1-Z_axis,RRg1-R_axis)
  if (tht1 .lt. 0.d0) tht1 = tht1 + 2.d0*PI
  tht2 = atan2(ZZg2-Z_axis,RRg2-R_axis)
  if (tht2 .lt. 0.d0) tht2 = tht2 + 2.d0*PI

!  write(*,*) k,tht1,tht2

  dRRg1_dt = dRRg1_dr * drr1 + dRRg1_ds * dss1
  dZZg1_dt = dZZg1_dr * drr1 + dZZg1_ds * dss1
  dRRg2_dt = dRRg2_dr * drr2 + dRRg2_ds * dss2
  dZZg2_dt = dZZg2_dr * drr2 + dZZg2_ds * dss2

!  RZ1  = RRg1     * tan(theta) - ZZg1
!  RZ2  = RRg2     * tan(theta) - ZZg2
!  dRZ1 = dRRg1_dt * tan(theta) - dZZg1_dt
!  dRZ2 = dRRg2_dt * tan(theta) - dZZg2_dt

  RZ1  = RRg1     * sin(theta) - ZZg1     * cos(theta)
  RZ2  = RRg2     * sin(theta) - ZZg2     * cos(theta)
  dRZ1 = dRRg1_dt * sin(theta) - dZZg1_dt * cos(theta)
  dRZ2 = dRRg2_dt * sin(theta) - dZZg2_dt * cos(theta)

  RZ0  = R_axis  * sin(theta) - Z_axis * cos(theta)

  a3 = (      RZ1 + dRZ1 -      RZ2 + dRZ2 )/4.d0
  a2 = (          - dRZ1            + dRZ2 )/4.d0
  a1 = (-3.d0*RZ1 - dRZ1 + 3.d0*RZ2 - dRZ2 )/4.d0
  a0 = ( 2.d0*RZ1 + dRZ1 + 2.d0*RZ2 - dRZ2 )/4.d0 - RZ0

  call SOLVP3(a0,a1,a2,a3,t,t2,t3,ifail)

!  if(   ( (tht1 .lt. theta) .and. (tht2 .gt. theta) ) &
!    .or.( (tht1 .gt. theta) .and. (tht2 .lt. theta) ) ) then
  if (abs(t) .lt. 1.01d0) then
    
    call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
    call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

    call interp_RZ(node_list,element_list,i_elm,ri,si,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                      ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

    tht_out = atan2(ZZg1-Z_axis,RRg1-R_axis)
    if (tht_out .lt. 0.d0) tht_out = tht_out + 2.d0 * PI
!------------------------------------------------------ now step to improve quality of solution
    iterMax =5
    do i=1,iterMax
      dRRg1_dt = dRRg1_dr * dri + dRRg1_ds * dsi
      dZZg1_dt = dZZg1_dr * dri + dZZg1_ds * dsi

      t = t - ((RRg1-R_axis)*sin(theta)-(ZZg1-Z_axis)*cos(theta)) / (dRRg1_dt * sin(theta) - dZZg1_dt * cos(theta))

      if (abs(t) .le. 1.d0) then

        call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
        call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

        call interp_RZ(node_list,element_list,i_elm,ri,si,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                      ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

        tht_out = atan2(ZZg1-Z_axis,RRg1-R_axis)
        if (tht_out .lt. 0.d0) tht_out = tht_out + 2.d0 * PI

!        write(*,'(A,3e16.8)') ' find_theta check : ',tht_out,theta,(RRg1-R_axis)*sin(theta)-(ZZg1-Z_axis)*cos(theta)

        if (abs(tht_in - tht_out) .lt. 1.D-4) then

          s_find(i_find+1)     = ri
          t_find(i_find+1)     = si
          i_elm_find(i_find+1) = i_elm
          i_find               = i_find + 1
          exit

        endif

      endif
    enddo

  endif

enddo

!if (i_find .eq. 0) write(*,*) ' WARNING : no theta surface found'

return
end subroutine find_theta_surface
