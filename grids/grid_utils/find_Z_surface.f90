subroutine find_Z_surface(node_list,element_list,surface_list,j_surf,Z_find,i_elm_find,s_find,t_find,st_find,i_find)
!---------------------------------------------------------------------------
! subroutine finds a Z value on a specific surface
!---------------------------------------------------------------------------

use data_structure
use mod_interp, only: interp_rz

implicit none

! --- Routine parameters
type (type_node_list),    intent(in)    :: node_list       !< node list
type (type_element_list), intent(in)    :: element_list    !< element list
type (type_surface_list), intent(in)    :: surface_list    !< surface list
integer,                  intent(in)    :: j_surf          !< index of surface we are intersecting in surface_list
real*8,                   intent(in)    :: Z_find          !< Z-value of horizontal line we are intersecting surface with
integer,                  intent(out)   :: i_elm_find(*)   !< the elm indices where we found intersections
real*8,                   intent(out)   :: s_find(*)       !< the local s-value of each element where we found intersections
real*8,                   intent(out)   :: t_find(*)       !< the local t-value of each element where we found intersections
real*8,                   intent(out)   :: st_find(*)      !< the surface local parameter at which we found intersection on surface piece (-1<=st<=+1)
integer,                  intent(out)   :: i_find          !< number of intersections we found

! --- local variables
integer :: i, k, i_elm, ifail, iterMax
real*8  :: ri, dri, si ,dsi
real*8  :: rr1, drr1, rr2, drr2, ss1, dss1, ss2, dss2, t, t2, t3, delta_t
real*8  :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8  :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8  :: RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss
real*8  :: ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss
real*8  :: dRRg1_dt, dZZg1_dt, dRRg2_dt, dZZg2_dt, RZ1, RZ2, dRZ1, dRZ2, RZ0, A0, A1, A2, A3


!write(*,*) ' finding Z on surface : ',surface_list%psi_values(j_surf),Z_find

i_find = 0
i_elm_find(1:8) = 0

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
!  write(*,'(A,i8,4f15.4)') ' find_Z check 0 : ',k,RRg1,ZZg1,RRg2,ZZg2
  dRRg1_dt = dRRg1_dr * drr1 + dRRg1_ds * dss1
  dZZg1_dt = dZZg1_dr * drr1 + dZZg1_ds * dss1
  dRRg2_dt = dRRg2_dr * drr2 + dRRg2_ds * dss2
  dZZg2_dt = dZZg2_dr * drr2 + dZZg2_ds * dss2

  RZ1  = + ZZg1
  RZ2  = + ZZg2
  dRZ1 = + dZZg1_dt
  dRZ2 = + dZZg2_dt

  RZ0  = Z_find

  a3 = (   RZ1 + dRZ1 -   RZ2 + dRZ2 )/4.d0
  a2 = (       - dRZ1         + dRZ2 )/4.
  a1 = (-3.d0*RZ1 - dRZ1 + 3.d0*RZ2 - dRZ2 )/4.d0
  a0 = ( 2.d0*RZ1 + dRZ1 + 2.d0*RZ2 - dRZ2 )/4.d0 - RZ0

  call SOLVP3(a0,a1,a2,a3,t,t2,t3,ifail)

  if (abs(t) .le. 1.05d0) then

    call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
    call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

    call interp_RZ(node_list,element_list,i_elm,ri,si,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                      ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)


!    write(*,'(A,2i8,4e16.8)') ' find_Z check 1 : ',i_find+1,i_elm,RRg1,t,ZZg1,ZZg1 - Z_find

!------------------------------------------------------ one step to improve quality of solution
    iterMax =2
    do i=1,iterMax
      dZZg1_dt = dZZg1_dr * dri + dZZg1_ds * dsi

      if (dZZg1_dt .ne. 0.d0)  delta_t = - (ZZg1-Z_find) / (dZZg1_dt)

      if (delta_t .le. 0.25) t = t + delta_t

      if (abs(t) .le. 1.000001d0) then

        call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
        call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

        call interp_RZ(node_list,element_list,i_elm,ri,si,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                      ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)


!        write(*,'(A,2i8,4e16.8)') ' find_Z check 2 : ',i_find+1,i_elm,RRg1,t,ZZg1,ZZg1 - Z_find

        if (abs(Z_find - ZZg1) .lt. 1.D-6) then
          s_find(i_find+1) = ri
          t_find(i_find+1) = si
          if (s_find(i_find+1) .lt. 0.d0) s_find(i_find+1) = 0.d0
          if (s_find(i_find+1) .gt. 1.d0) s_find(i_find+1) = 1.d0
          if (t_find(i_find+1) .lt. 0.d0) t_find(i_find+1) = 0.d0
          if (t_find(i_find+1) .gt. 1.d0) t_find(i_find+1) = 1.d0
          i_elm_find(i_find+1) = i_elm
          st_find(i_find+1) = t
          i_find = i_find + 1
          exit
        endif

      endif
    enddo

  endif

enddo

!if (i_find .eq. 0) write(*,*) ' WARNING : no Z surface found'

return
end subroutine find_Z_surface
