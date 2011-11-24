!> Subroutine defines the new grid_points from crossing of polar and radial coordinate lines
subroutine define_new_grid_points(node_list, element_list, flux_list, &
                                   xcase, R_xpoint, Z_xpoint, psi_xpoint, n_grids, stpts, sigmas, nwpts)

use constants
use tr_module 
use data_structure
use grid_xpoint_data

implicit none

! --- Routine parameters
type (type_surface_list)    , intent(inout) :: flux_list
type (type_node_list)       , intent(inout) :: node_list
type (type_element_list)    , intent(inout) :: element_list
type (type_strategic_points), intent(inout) :: stpts
type (type_new_points)      , intent(inout) :: nwpts
integer,                      intent(inout) :: n_grids(10)
integer,                      intent(in)    :: xcase
real*8,                       intent(in)    :: sigmas(16)
real*8,                       intent(in)    :: psi_xpoint(2), R_xpoint(2), Z_xpoint(2)

! --- local variables
real*8, allocatable :: s_tmp(:), theta_sep(:)
real*8, allocatable :: xp(:),yp(:)
integer             :: i, j, k, m, i_elm_find(8), i_find, i_elm_axis
integer             :: i_sep1, i_sep2, i_max, n_loop, n_start, i_surf
integer             :: n_psi, n_tht_2, n_tht_mid, n_tht_mid2
integer	            :: n_flux, n_tht,   n_open,   n_outer,   n_inner    
integer	            :: n_private,   n_up_priv,   n_leg,   n_up_leg
integer             :: npl, ifail, my_id
real*8              :: delta, ss, tmp1, tmp2
real*8              :: R_cub1d(4), Z_cub1d(4)
real*8              :: psi_axis, R_axis, Z_axis, s_axis, t_axis
real*8              :: tht_x1, tht_x2
real*8              :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8              :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8              :: PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss
real*8              :: s_find(8), t_find(8)
real*8              :: R_beg, R_end
real*8              :: Z_beg, Z_end
real*8              :: SIG_theta
real*8              :: SIG_leg_0, SIG_leg_1
real*8              :: SIG_up_leg_0, SIG_up_leg_1
real*8              :: SIG_0, SIG_1

write(*,*) '*****************************************'
write(*,*) '* X-point grid : Define new grid points *'
write(*,*) '*****************************************'
write(*,*) '                 Define extrapolation points'


my_id = 1 ! Just don't want the printout...
call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

SIG_theta    = sigmas(2) 
SIG_leg_0    = sigmas(8) ; SIG_leg_1    = sigmas(9) 
SIG_up_leg_0 = sigmas(10); SIG_up_leg_1 = sigmas(11)

n_flux    = n_grids(1); n_tht     = n_grids(2)
n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
n_private = n_grids(6); n_up_priv = n_grids(7)
n_leg	  = n_grids(8); n_up_leg  = n_grids(9)

nwpts%k_cross = 0



!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!*********************************** First part: find extrapolation points  *********************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!



!-------------------------------------------------------------------------------------------!
!------- Extrapolation points for first part of the grid (everything except legs) ----------!
!-------------------------------------------------------------------------------------------!

!-------------------------------- define the polar coordinate
n_psi   = n_flux + n_open + n_outer + n_inner + n_private + n_up_priv + 1   ! this includes the magnetic axis
n_tht_2 = n_tht + 2*n_leg + 2*n_up_leg

call tr_allocate(theta_sep,1,n_tht_2,"theta_sep",CAT_GRID)
  
if (xcase .ne. 3) then 
  if (xcase .eq. 1) tht_x1 = atan2(Z_xpoint(1)-Z_axis,R_xpoint(1)-R_axis)
  if (xcase .eq. 2) tht_x1 = atan2(Z_xpoint(2)-Z_axis,R_xpoint(2)-R_axis)
  n_tht_mid = n_tht/2
  
  call tr_allocate(s_tmp,1,n_tht,"s_tmp",CAT_GRID)
  s_tmp = 0
  call meshac2(n_tht,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta,0.8d0,1.0d0)

  do j=1,n_tht
    theta_sep(j) = tht_x1 + 2.d0 * PI * s_tmp(j)
    if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
    if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
  enddo
  call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  
else
  if (psi_xpoint(1) .le. psi_xpoint(2)) then
    tht_x1 = atan2(Z_xpoint(1)-Z_axis,R_xpoint(1)-R_axis)
    tht_x2 = atan2(Z_xpoint(2)-Z_axis,R_xpoint(2)-R_axis)
  else
    tht_x2 = atan2(Z_xpoint(1)-Z_axis,R_xpoint(1)-R_axis)
    tht_x1 = atan2(Z_xpoint(2)-Z_axis,R_xpoint(2)-R_axis)
  endif
  if (tht_x1 .lt. 0.d0) tht_x1 = tht_x1 + 2.d0 * PI
  if (tht_x2 .lt. 0.d0) tht_x2 = tht_x2 + 2.d0 * PI
  !write(*,'(A,2f)') ' angles : ',tht_x1,tht_x2
  
  ! Spread out points evenly (outer angle between tht_x1 and tht_x2 is usually bigger than inner angle)
  if (psi_xpoint(1) .le. psi_xpoint(2)) then
    n_tht_mid = int(n_tht * (2.d0*PI - (tht_x1 - tht_x2)) / (2.d0*PI))
    ! Make sure n_tht_mid is odd and save it to n_grids for later use
    if(mod(n_tht_mid,2) .eq. 0) n_tht_mid = n_tht_mid + 1
    n_grids(10) = n_tht_mid
    
    call tr_allocate(s_tmp,1,n_tht_mid,"s_tmp",CAT_GRID)
    s_tmp = 0
    call meshac2(n_tht_mid,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta,0.8d0,1.0d0)

    do j=1,n_tht_mid
      theta_sep(j) = (tht_x1-2.d0*PI) + (tht_x2-(tht_x1-2.d0*PI)) * s_tmp(j)
      if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
      if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  
    n_tht_mid2 = n_tht-n_tht_mid
    call tr_allocate(s_tmp,1,n_tht_mid2,"s_tmp",CAT_GRID)
    s_tmp = 0
    call meshac2(n_tht_mid2,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta,0.8d0,1.0d0)

    do j=n_tht_mid+1,n_tht
      theta_sep(j) = tht_x2 + (tht_x1-tht_x2) * s_tmp(j-n_tht_mid)
      if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
      if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  else
    n_tht_mid = int(n_tht * (tht_x2 - tht_x1) / (2.d0*PI))
    ! Make sure n_tht_mid is odd and save it to n_grids for later use
    if(mod(n_tht_mid,2) .eq. 0) n_tht_mid = n_tht_mid + 1
    n_grids(10) = n_tht_mid
    
    call tr_allocate(s_tmp,1,n_tht_mid,"s_tmp",CAT_GRID)
    s_tmp = 0
    call meshac2(n_tht_mid,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta,0.8d0,1.0d0)
    do j=1,n_tht_mid
      theta_sep(j) = tht_x1 + (tht_x2-tht_x1) * s_tmp(j)
      if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
      if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)

    n_tht_mid2 = n_tht-n_tht_mid
    call tr_allocate(s_tmp,1,n_tht_mid2,"s_tmp",CAT_GRID)
    s_tmp = 0
    call meshac2(n_tht_mid2,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta,0.8d0,1.0d0)
    do j=n_tht_mid+1,n_tht
      theta_sep(j) = (tht_x2-2.d0*PI) + (tht_x1-(tht_x2-2.d0*PI)) * s_tmp(j-n_tht_mid)
      if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
      if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  
  endif      
endif

!------------------------------------- find crossing with separatrix
i_sep1  = n_flux
i_sep2  = n_flux + n_open
if (xcase .ne. 3) then
  do j=2,n_tht-1
    call find_theta_surface(node_list,element_list,flux_list,i_sep1,theta_sep(j),R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)
    if(i_find .eq. 0) return
    
    call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                   RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                   ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    if( ((ZZg1 .lt. Z_xpoint(1)) .and. (xcase .eq. 1)) .or. ((ZZg1 .gt. Z_xpoint(2)) .and. (xcase .eq. 2)) ) then
      call interp_RZ(node_list,element_list,i_elm_find(2),s_find(2),t_find(2),&
    		     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,  &
    		     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    endif

    nwpts%R_sep(j) = RRg1
    nwpts%Z_sep(j) = ZZg1
  enddo
else
  do j=2,n_tht-1
    if((j .ne. n_tht_mid) .and. (j .ne. n_tht_mid+1)) then
      if (theta_sep(j) .ge. pi) then
        if (psi_xpoint(1) .le. psi_xpoint(2)) then
          call find_theta_surface(node_list,element_list,flux_list,i_sep1,theta_sep(j),R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)
        else
          call find_theta_surface(node_list,element_list,flux_list,i_sep2,theta_sep(j),R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)
        endif
        if(i_find .eq. 0) return

        call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
        	       RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,	&
        	       ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
        if(ZZg1 .lt. Z_xpoint(1)) then
          call interp_RZ(node_list,element_list,i_elm_find(2),s_find(2),t_find(2),&
        	         RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,	&
        	         ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
	endif
        
	nwpts%R_sep(j) = RRg1
        nwpts%Z_sep(j) = ZZg1
      else
        if (psi_xpoint(1) .le. psi_xpoint(2)) then
          call find_theta_surface(node_list,element_list,flux_list,i_sep2,theta_sep(j),R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)
        else
          call find_theta_surface(node_list,element_list,flux_list,i_sep1,theta_sep(j),R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)
        endif
        if(i_find .eq. 0) return

        call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
        	       RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,	&
        	       ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
        if(ZZg1 .gt. Z_xpoint(2)) then
          call interp_RZ(node_list,element_list,i_elm_find(2),s_find(2),t_find(2),&
        	         RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,	&
        	         ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
	endif

	nwpts%R_sep(j) = RRg1
        nwpts%Z_sep(j) = ZZg1    
      endif
    endif
  enddo
endif

if (xcase .eq. 1) then
  nwpts%R_sep(1)	     = R_xpoint(1) ! this one is known - safer...
  nwpts%Z_sep(1)	     = Z_xpoint(1) ! this one is known - safer...
  nwpts%R_sep(n_tht)         = R_xpoint(1) ! this one is known - safer...
  nwpts%Z_sep(n_tht)         = Z_xpoint(1) ! this one is known - safer...
endif
if (xcase .eq. 2) then
  nwpts%R_sep(1)	     = R_xpoint(2) ! this one is known - safer...
  nwpts%Z_sep(1)	     = Z_xpoint(2) ! this one is known - safer...
  nwpts%R_sep(n_tht)         = R_xpoint(2) ! this one is known - safer...
  nwpts%Z_sep(n_tht)         = Z_xpoint(2) ! this one is known - safer...
endif
if (xcase .eq. 3) then
  if (psi_xpoint(1) .le. psi_xpoint(2)) then
    nwpts%R_sep(1)	     = R_xpoint(1) ! this one is known - safer...
    nwpts%Z_sep(1)	     = Z_xpoint(1) ! this one is known - safer...
    nwpts%R_sep(n_tht)       = R_xpoint(1) ! this one is known - safer...
    nwpts%Z_sep(n_tht)       = Z_xpoint(1) ! this one is known - safer...
    nwpts%R_sep(n_tht_mid)   = R_xpoint(2) ! this one is known - safer...
    nwpts%Z_sep(n_tht_mid)   = Z_xpoint(2) ! this one is known - safer...
    nwpts%R_sep(n_tht_mid+1) = R_xpoint(2) ! this one is known - safer...
    nwpts%Z_sep(n_tht_mid+1) = Z_xpoint(2) ! this one is known - safer...
  else
    nwpts%R_sep(1)	     = R_xpoint(2) ! this one is known - safer...
    nwpts%Z_sep(1)	     = Z_xpoint(2) ! this one is known - safer...
    nwpts%R_sep(n_tht)       = R_xpoint(2) ! this one is known - safer...
    nwpts%Z_sep(n_tht)       = Z_xpoint(2) ! this one is known - safer...
    nwpts%R_sep(n_tht_mid)   = R_xpoint(1) ! this one is known - safer...
    nwpts%Z_sep(n_tht_mid)   = Z_xpoint(1) ! this one is known - safer...
    nwpts%R_sep(n_tht_mid+1) = R_xpoint(1) ! this one is known - safer...
    nwpts%Z_sep(n_tht_mid+1) = Z_xpoint(1) ! this one is known - safer...
  endif
endif

!------------------------------------ find crossing with last fluxsurface 
do j=1,n_tht

  if (nwpts%Z_sep(j) .le. Z_axis) then
  
    if ( (xcase .eq. 1) .or. ((xcase .eq. 3) .and. (psi_xpoint(1) .le. psi_xpoint(2))) ) then        
      if (j .gt. n_tht_mid) then
        nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_LowerInnerLeg - Z_xpoint(1)) * ((nwpts%Z_sep(j) - Z_axis)/(Z_xpoint(1) - Z_axis))**2
	i_max = n_flux + n_open + n_outer + n_inner
      else
        nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_LowerOuterLeg - Z_xpoint(1)) * ((nwpts%Z_sep(j) - Z_axis)/(Z_xpoint(1) - Z_axis))**2
	i_max = n_flux + n_open + n_outer
      endif      
      call find_Z_surface(node_list,element_list,flux_list,i_max,nwpts%Z_max(j),i_elm_find,s_find,t_find,i_find)    
    elseif ( (xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) ) then      
      if (j .gt. n_tht_mid) then
        nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_LowerOuterLeg - Z_xpoint(1)) * ((nwpts%Z_sep(j) - Z_axis)/(Z_xpoint(1) - Z_axis))**2
	i_max = n_flux + n_open + n_outer
      else
        nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_LowerInnerLeg - Z_xpoint(1)) * ((nwpts%Z_sep(j) - Z_axis)/(Z_xpoint(1) - Z_axis))**2
	i_max = n_flux + n_open + n_outer + n_inner
      endif      
      call find_Z_surface(node_list,element_list,flux_list,i_max,nwpts%Z_max(j),i_elm_find,s_find,t_find,i_find)    
    else    
      i_max = n_flux + n_open
      call find_theta_surface(node_list,element_list,flux_list,i_max,theta_sep(j),R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)    
    endif

    call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1), &
        	   RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
        	   ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

    if(    (xcase .eq. 2)                                                                             &
      .or. (     ( (xcase .eq. 1) .or. ((xcase .eq. 3) .and. (psi_xpoint(1) .le. psi_xpoint(2)) ) ) &
           .and. (    ( (RRg1 .gt. R_xpoint(1)) .and. (j.lt.n_tht_mid) ) &
                 .or. ( (RRg1 .lt. R_xpoint(1)) .and. (j.gt.n_tht_mid) ) )                            ) &
      .or. (     ( (xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) )                        &
           .and. (    ( (RRg1 .lt. R_xpoint(1)) .and. (j.le.n_tht_mid) ) &
                 .or. ( (RRg1 .gt. R_xpoint(1)) .and. (j.gt.n_tht_mid) ) )                            ) ) then

      nwpts%R_max(j)  	         = RRg1
      nwpts%Z_max(j)  	         = ZZg1
      nwpts%RR_new(i_max+1,j)    = RRg1
      nwpts%ZZ_new(i_max+1,j)    = ZZg1
      nwpts%ielm_flux(i_max+1,j) = i_elm_find(1)
      nwpts%s_flux(i_max+1,j)    = s_find(1)
      nwpts%t_flux(i_max+1,j)    = t_find(1)
      nwpts%t_tht(i_max+1,j)     = 1.d0
      nwpts%k_cross(i_max+1,j)   = 3
      
    else

      call interp_RZ(node_list,element_list,i_elm_find(2),s_find(2),t_find(2),&
    		     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
    		     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      nwpts%R_max(j)  	         = RRg1
      nwpts%Z_max(j)  	         = ZZg1
      nwpts%RR_new(i_max+1,j)    = RRg1
      nwpts%ZZ_new(i_max+1,j)    = ZZg1
      nwpts%ielm_flux(i_max+1,j) = i_elm_find(2)
      nwpts%s_flux(i_max+1,j)    = s_find(2)
      nwpts%t_flux(i_max+1,j)    = t_find(2)
      nwpts%t_tht(i_max+1,j)     = 1.d0
      nwpts%k_cross(i_max+1,j)   = 3

    endif
    
    if (     ( (xcase .eq. 1) .or. ((xcase .eq. 3) .and. (psi_xpoint(1) .le. psi_xpoint(2)) ) ) &
       .and. ((j .eq. 1) .or. (j .eq. n_tht))                                                   ) then
      nwpts%R_max(1)  	          = stpts%RLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%Z_max(1)  	          = stpts%ZLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,1)     = stpts%RLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,1)     = stpts%ZLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%R_max(n_tht)	  = stpts%RLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%Z_max(n_tht)	  = stpts%ZLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht) = stpts%RLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht) = stpts%ZLimit_LowerInnerLeg ! this one is known - safer...
    endif
    if (     ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) )  &
       .and. ((j .eq. 1) .or. (j .eq. n_tht))                            ) then
      nwpts%R_max(n_tht_mid)	        = stpts%RLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%Z_max(n_tht_mid)	        = stpts%ZLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht_mid)   = stpts%RLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht_mid)   = stpts%ZLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%R_max(n_tht_mid+1)  	= stpts%RLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%Z_max(n_tht_mid+1)  	= stpts%ZLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht_mid+1) = stpts%RLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht_mid+1) = stpts%ZLimit_LowerOuterLeg ! this one is known - safer...
    endif

  else
        
    if (    ( (j .gt. n_tht_mid) .and. (xcase .eq. 3) .and. (psi_xpoint(1) .le. psi_xpoint(2)) ) & 
       .or. ( (j .lt. n_tht_mid) .and. (xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) ) & 
       .or. ( (j .lt. n_tht_mid) .and. (xcase .eq. 2) )                                          ) then
      nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_UpperInnerLeg - Z_xpoint(2)) * ((nwpts%Z_sep(j) - Z_axis)/(Z_xpoint(2) - Z_axis))**2
      i_max = n_flux + n_open + n_outer + n_inner
      call find_Z_surface(node_list,element_list,flux_list,i_max,nwpts%Z_max(j),i_elm_find,s_find,t_find,i_find)
    endif 
    if (    ( (j .le. n_tht_mid) .and. (xcase .eq. 3) .and. (psi_xpoint(1) .le. psi_xpoint(2)) ) &
       .or. ( (j .gt. n_tht_mid) .and. (xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) ) & 
       .or. ( (j .gt. n_tht_mid) .and. (xcase .eq. 2) )                                          ) then
      nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_UpperOuterLeg - Z_xpoint(2)) * ((nwpts%Z_sep(j) - Z_axis)/(Z_xpoint(2) - Z_axis))**2
      i_max = n_flux + n_open + n_outer
      call find_Z_surface(node_list,element_list,flux_list,i_max,nwpts%Z_max(j),i_elm_find,s_find,t_find,i_find)
    endif
    if (xcase .eq. 1) then
      call find_theta_surface(node_list,element_list,flux_list,i_max,theta_sep(j),R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)
    endif

    call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
    		   RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
    		   ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)


    if(    (xcase .eq. 1)                                                                              &
      .or. (     ( (xcase .eq. 3) .and. (psi_xpoint(1) .le. psi_xpoint(2)) )                           &
           .and. (    ( (RRg1 .ge. R_xpoint(2)) .and. (j .le. n_tht_mid) ) &
                 .or. ( (RRg1 .lt. R_xpoint(2)) .and. (j .gt. n_tht_mid) )   )                     ) &
      .or. (     ( (xcase .eq. 2) .or. ( (xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) ) )   &
           .and. (    ( (RRg1 .ge. R_xpoint(2)) .and. (j .gt. n_tht_mid) ) &
                 .or. ( (RRg1 .lt. R_xpoint(2)) .and. (j .lt. n_tht_mid) ) ) 			       ) ) then

      nwpts%R_max(j)  	         = RRg1
      nwpts%Z_max(j)  	         = ZZg1
      nwpts%RR_new(i_max+1,j)    = RRg1
      nwpts%ZZ_new(i_max+1,j)    = ZZg1
      nwpts%ielm_flux(i_max+1,j) = i_elm_find(1)
      nwpts%s_flux(i_max+1,j)    = s_find(1)
      nwpts%t_flux(i_max+1,j)    = t_find(1)
      nwpts%t_tht(i_max+1,j)     = 1.d0
      nwpts%k_cross(i_max+1,j)   = 3
      
    else

      call interp_RZ(node_list,element_list,i_elm_find(2),s_find(2),t_find(2),&
    		     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
    		     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      nwpts%R_max(j)  	         = RRg1
      nwpts%Z_max(j)  	         = ZZg1
      nwpts%RR_new(i_max+1,j)    = RRg1
      nwpts%ZZ_new(i_max+1,j)    = ZZg1
      nwpts%ielm_flux(i_max+1,j) = i_elm_find(2)
      nwpts%s_flux(i_max+1,j)    = s_find(2)
      nwpts%t_flux(i_max+1,j)    = t_find(2)
      nwpts%t_tht(i_max+1,j)     = 1.d0
      nwpts%k_cross(i_max+1,j)   = 3

    endif
    
    if (     ( (xcase .eq. 2) .or. ( (xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) ) ) & 
       .and. ( (j .eq. 1) .or. (j .eq. n_tht_2) )                                                ) then
      nwpts%R_max(n_tht)	    = stpts%RLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%Z_max(n_tht)	    = stpts%ZLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht)   = stpts%RLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht)   = stpts%ZLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%R_max(1)  	            = stpts%RLimit_UpperInnerLeg ! this one is known - safer...
      nwpts%Z_max(1)  	            = stpts%ZLimit_UpperInnerLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,1)       = stpts%RLimit_UpperInnerLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,1)       = stpts%ZLimit_UpperInnerLeg ! this one is known - safer...
    endif
    
    if (     ( (xcase .eq. 3) .and. (psi_xpoint(1) .le. psi_xpoint(2)) ) & 
       .and. ( (j .eq. n_tht_mid) .or. (j .eq. n_tht_mid+1) )            ) then
      nwpts%R_max(n_tht_mid)	        = stpts%RLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%Z_max(n_tht_mid)	        = stpts%ZLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht_mid)   = stpts%RLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht_mid)   = stpts%ZLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%R_max(n_tht_mid+1)	        = stpts%RLimit_UpperInnerLeg ! this one is known - safer...
      nwpts%Z_max(n_tht_mid+1)	        = stpts%ZLimit_UpperInnerLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht_mid+1) = stpts%RLimit_UpperInnerLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht_mid+1) = stpts%ZLimit_UpperInnerLeg ! this one is known - safer...
    endif

  endif

enddo





!--------------------------------------------------------------------------!
!------- Extrapolation points for second part of the grid (Legs) ----------!
!--------------------------------------------------------------------------!

!------------------------------ Intersections with private surfaces
do i=1,4
  
  if (xcase .ne. 2) then 
    if (i .eq. 1) then
      n_loop  = n_leg
      n_start = n_tht
      i_surf  = n_flux+n_open+n_outer+n_inner+n_private
      R_beg   = stpts%RRightCorn_LowerInnerLeg
      R_end   = stpts%RMiddle_LowerPrivate
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_leg
      n_start = n_tht + n_leg
      i_surf  = n_flux+n_open+n_outer+n_inner+n_private
      R_beg   = stpts%RLeftCorn_LowerOuterLeg
      R_end   = stpts%RMiddle_LowerPrivate
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 3) then
      if(xcase .eq. 3) then
        n_loop  = n_up_leg
        n_start = n_tht + n_leg + n_leg
        i_surf  = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv
        R_beg   = stpts%RRightCorn_UpperInnerLeg
        R_end   = stpts%RMiddle_UpperPrivate
        SIG_0   = SIG_up_leg_0
        SIG_1   = SIG_up_leg_1
      else
        exit
      endif
    endif
    if (i .eq. 4) then
      if(xcase .eq. 3) then
        n_loop  = n_up_leg
        n_start = n_tht + n_leg + n_leg + n_up_leg
        i_surf  = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv
        R_beg   = stpts%RLeftCorn_UpperOuterLeg
        R_end   = stpts%RMiddle_UpperPrivate
        SIG_0   = SIG_up_leg_0
        SIG_1   = SIG_up_leg_1
      else
        exit
      endif
    endif
  else
    if (i .eq. 1) then
      n_loop  = n_up_leg
      n_start = n_tht
      i_surf  = n_flux+n_open+n_up_priv
      R_beg   = stpts%RRightCorn_UpperInnerLeg
      R_end   = stpts%RMiddle_UpperPrivate
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_up_leg
      n_start = n_tht + n_up_leg
      i_surf  = n_flux+n_open+n_up_priv
      R_beg   = stpts%RLeftCorn_UpperOuterLeg
      R_end   = stpts%RMiddle_UpperPrivate
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    if (i .eq. 3) exit
    if (i .eq. 4) exit
  endif
  
  call tr_allocate(s_tmp,1,n_loop,"s_tmp",CAT_GRID)
  s_tmp = 0
  call meshac2(n_loop,s_tmp,0.d0,1.d0,SIG_0,SIG_1,0.6d0,1.0d0)
  do j=1,n_loop

    nwpts%R_min(n_start + j) = R_beg + (R_end-R_beg) * s_tmp(j)

    call find_R_surface(node_list,element_list,flux_list,i_surf,nwpts%R_min(n_start+j),i_elm_find,s_find,t_find,i_find)

    do k=1,i_find

      call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),&
    		     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
    		     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ( (xcase .eq. 1) .and. (ZZg1 .le. Z_xpoint(1)) ) exit
      if ( (xcase .eq. 2) .and. (ZZg1 .ge. Z_xpoint(2)) ) exit
      if ( (xcase .eq. 3) .and. (ZZg1 .le. Z_xpoint(1)) .and. (i .le. 2) ) exit
      if ( (xcase .eq. 3) .and. (ZZg1 .ge. Z_xpoint(2)) .and. (i .ge. 3) ) exit

    enddo

    nwpts%Z_min(n_start + j) = ZZg1

  enddo
  call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  if ( (xcase .ne. 2) .and. (i .eq. 1) ) nwpts%Z_min(n_start+1) = stpts%ZRightCorn_LowerInnerLeg ! this one is known - safer...
  if ( (xcase .ne. 2) .and. (i .eq. 2) ) nwpts%Z_min(n_start+1) = stpts%ZLeftCorn_LowerOuterLeg  ! this one is known - safer...
  if ( (xcase .eq. 3) .and. (i .eq. 3) ) nwpts%Z_min(n_start+1) = stpts%ZRightCorn_UpperInnerLeg ! this one is known - safer...
  if ( (xcase .eq. 3) .and. (i .eq. 4) ) nwpts%Z_min(n_start+1) = stpts%ZLeftCorn_UpperOuterLeg  ! this one is known - safer...
  if ( (xcase .eq. 2) .and. (i .eq. 1) ) nwpts%Z_min(n_start+1) = stpts%ZRightCorn_UpperInnerLeg ! this one is known - safer...
  if ( (xcase .eq. 2) .and. (i .eq. 2) ) nwpts%Z_min(n_start+1) = stpts%ZLeftCorn_UpperOuterLeg  ! this one is known - safer...
enddo

!------------------------------ Intersections with open surfaces
do i=1,4
  
  if (xcase .ne. 2) then 
    if (i .eq. 1) then
      n_loop  = n_leg
      n_start = n_tht
      i_surf  = n_flux+n_open+n_outer+n_inner
      Z_beg   = stpts%ZLeftCorn_LowerInnerLeg
      Z_end   = stpts%ZLimit_LowerInnerLeg
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_leg
      n_start = n_tht + n_leg
      i_surf  = n_flux+n_open+n_outer
      Z_beg   = stpts%ZRightCorn_LowerOuterLeg
      Z_end   = stpts%ZLimit_LowerOuterLeg
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 3) then
      if(xcase .eq. 3) then
        n_loop  = n_up_leg
        n_start = n_tht + n_leg + n_leg
        i_surf  = n_flux+n_open+n_outer+n_inner
        Z_beg   = stpts%ZLeftCorn_UpperInnerLeg
        Z_end   = stpts%ZLimit_UpperInnerLeg
        SIG_0   = SIG_up_leg_0
        SIG_1   = SIG_up_leg_1
      else
        exit
      endif
    endif
    if (i .eq. 4) then
      if(xcase .eq. 3) then
        n_loop  = n_up_leg
        n_start = n_tht + n_leg + n_leg + n_up_leg
        i_surf  = n_flux+n_open+n_outer
        Z_beg   = stpts%ZRightCorn_UpperOuterLeg
        Z_end   = stpts%ZLimit_UpperOuterLeg
        SIG_0   = SIG_up_leg_0
        SIG_1   = SIG_up_leg_1
      else
        exit
      endif
    endif
  else
    if (i .eq. 1) then
      n_loop  = n_up_leg
      n_start = n_tht
      i_surf  = n_flux+n_open
      Z_beg   = stpts%ZLeftCorn_UpperInnerLeg
      Z_end   = stpts%ZLimit_UpperInnerLeg
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_up_leg
      n_start = n_tht + n_up_leg
      i_surf  = n_flux+n_open
      Z_beg   = stpts%ZRightCorn_UpperOuterLeg
      Z_end   = stpts%ZLimit_UpperOuterLeg
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    if (i .eq. 3) exit
    if (i .eq. 4) exit
  endif
  
  call tr_allocate(s_tmp,1,n_loop,"s_tmp",CAT_GRID)
  s_tmp = 0
  call meshac2(n_loop,s_tmp,0.d0,1.d0,SIG_0,SIG_1,0.6d0,1.0d0)
  do j=1,n_loop

    nwpts%Z_max(n_start + j) = Z_beg + (Z_end-Z_beg) * s_tmp(j)

    call find_Z_surface(node_list,element_list,flux_list,i_surf,nwpts%Z_max(n_start+j),i_elm_find,s_find,t_find,i_find)

    do k=1,i_find

      call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),&
    		     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
    		     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ( (xcase .ne. 2) .and. (i .eq. 1) .and. (RRg1 .le. R_xpoint(1)) ) exit
      if ( (xcase .ne. 2) .and. (i .eq. 2) .and. (RRg1 .ge. R_xpoint(1)) ) exit
      if ( (xcase .ne. 2) .and. (i .eq. 3) .and. (RRg1 .le. R_xpoint(2)) ) exit
      if ( (xcase .ne. 2) .and. (i .eq. 4) .and. (RRg1 .ge. R_xpoint(2)) ) exit
      if ( (xcase .eq. 2) .and. (i .eq. 1) .and. (RRg1 .le. R_xpoint(2)) ) exit
      if ( (xcase .eq. 2) .and. (i .eq. 2) .and. (RRg1 .ge. R_xpoint(2)) ) exit

    enddo

    nwpts%R_max(n_start + j) = RRg1

  enddo
  call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  if ( (xcase .eq. 2) .and. (i .eq. 1) ) nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht)       ! this one is known - safer...
  if ( (xcase .eq. 2) .and. (i .eq. 2) ) nwpts%R_max(n_start+n_loop) = nwpts%R_max(1)           ! this one is known - safer...
  if ( (xcase .eq. 3) .and. (i .eq. 3) ) nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht_mid+1) ! this one is known - safer...
  if ( (xcase .eq. 3) .and. (i .eq. 4) ) nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht_mid)   ! this one is known - safer...
  if ( (xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) .and. (i .eq. 1) ) &
                                         nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht_mid)   ! this one is known - safer...
  if ( (xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) .and. (i .eq. 2) ) &
                                         nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht_mid+1) ! this one is known - safer...
  if ( (xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) .and. (i .eq. 3) ) &
                                         nwpts%R_max(n_start+n_loop) = nwpts%R_max(1)           ! this one is known - safer...
  if ( (xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) .and. (i .eq. 4) ) &
                                         nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht)       ! this one is known - safer...
  if ( (xcase .eq. 2) .and. (i .eq. 1) ) nwpts%R_max(n_start+n_loop) = nwpts%R_max(1)           ! this one is known - safer...
  if ( (xcase .eq. 2) .and. (i .eq. 2) ) nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht)       ! this one is known - safer...
enddo

!------------------------------ Intersections with separatrices
do i=1,4
  
  if (xcase .ne. 2) then 
    if (i .eq. 1) then
      n_loop  = n_leg
      n_start = n_tht
      i_surf  = n_flux
      if((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) i_surf  = n_flux+n_open
      R_beg   = stpts%RStrike_LowerInnerLeg
      R_end   = R_xpoint(1)
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_leg
      n_start = n_tht + n_leg
      i_surf  = n_flux
      if((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) i_surf  = n_flux+n_open
      R_beg   = stpts%RStrike_LowerOuterLeg
      R_end   = R_xpoint(1)
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 3) then
      if(xcase .eq. 3) then
        n_loop  = n_up_leg
        n_start = n_tht + n_leg + n_leg
        i_surf  = n_flux+n_open
        if((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) i_surf  = n_flux
        R_beg   = stpts%RStrike_UpperInnerLeg
        R_end   = R_xpoint(2)
        SIG_0   = SIG_up_leg_0
        SIG_1   = SIG_up_leg_1
      else
        exit
      endif
    endif
    if (i .eq. 4) then
      if(xcase .eq. 3) then
        n_loop  = n_up_leg
        n_start = n_tht + n_leg + n_leg + n_up_leg
        i_surf  = n_flux+n_open
        if((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) i_surf  = n_flux
        R_beg   = stpts%RStrike_UpperOuterLeg
        R_end   = R_xpoint(2)
        SIG_0   = SIG_up_leg_0
        SIG_1   = SIG_up_leg_1
      else
        exit
      endif
    endif
  else
    if (i .eq. 1) then
      n_loop  = n_up_leg
      n_start = n_tht
      i_surf  = n_flux
      R_beg   = stpts%RStrike_UpperInnerLeg
      R_end   = R_xpoint(2)
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_up_leg
      n_start = n_tht + n_up_leg
      i_surf  = n_flux
      R_beg   = stpts%RStrike_UpperOuterLeg
      R_end   = R_xpoint(2)
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    if (i .eq. 3) exit
    if (i .eq. 4) exit
  endif
  
  call tr_allocate(s_tmp,1,n_loop,"s_tmp",CAT_GRID)
  s_tmp = 0
  call meshac2(n_loop,s_tmp,0.d0,1.d0,SIG_0,SIG_1,0.6d0,1.0d0)
  do j=1,n_loop

    nwpts%R_sep(n_start + j) = R_beg + (R_end-R_beg) * s_tmp(j)

    call find_R_surface(node_list,element_list,flux_list,i_surf,nwpts%R_sep(n_start+j),i_elm_find,s_find,t_find,i_find)

    do k=1,i_find

      call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),&
    		     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
    		     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      
      if ( (xcase .eq. 1) .and. (ZZg1 .le. Z_xpoint(1)) ) exit
      if ( (xcase .eq. 2) .and. (ZZg1 .ge. Z_xpoint(2)) ) exit
      if ( (xcase .eq. 3) .and. (ZZg1 .le. Z_xpoint(1)) .and. (i .le. 2) ) exit
      if ( (xcase .eq. 3) .and. (ZZg1 .ge. Z_xpoint(2)) .and. (i .ge. 3) ) exit

    enddo

    nwpts%Z_sep(n_start + j) = ZZg1

  enddo
  call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  if ( (xcase .ne. 2) .and. (i .le. 2) ) nwpts%Z_sep(n_start+n_loop) = Z_xpoint(1) ! this one is known - safer...
  if ( (xcase .eq. 3) .and. (i .ge. 3) ) nwpts%Z_sep(n_start+n_loop) = Z_xpoint(2) ! this one is known - safer...
  if   (xcase .eq. 2)                    nwpts%Z_sep(n_start+n_loop) = Z_xpoint(2) ! this one is known - safer...
enddo

!call lincol(2)
!call lplot6(1,1,nwpts%R_max,nwpts%Z_max,-(n_tht_3),' ')
!call lincol(0)
!call lplot6(1,1,nwpts%R_sep,nwpts%Z_sep,-(n_tht_3),' ')







!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!*********************************** Second part: find crossings of lines ***********************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!
write(*,*) '                 Find crossings between coordinate lines'





!------------------------------ Construct polar coordinate lines
do j=1,n_tht

  delta = 0.1

  if ( (j .eq. 1) .or. (j .eq. n_tht) )       delta = 0.d0
  if ( (j .eq. 2) .or. (j .eq. n_tht - 1) )   delta = 0.05d0
  if ( ((j .eq. n_tht_mid)   .or. (j .eq. n_tht_mid+1)) .and. (xcase .eq. 3) ) delta = 0.d0
  if ( ((j .eq. n_tht_mid-1) .or. (j .eq. n_tht_mid+2)) .and. (xcase .eq. 3) ) delta = 0.05d0

  nwpts%R_polar(1,1,j) = R_axis
  nwpts%R_polar(1,4,j) = delta * R_axis + (1.d0 - delta) * nwpts%R_sep(j)
  nwpts%R_polar(1,2,j) = ( 2.d0 * nwpts%R_polar(1,1,j)  +         nwpts%R_polar(1,4,j) ) / 3.d0
  nwpts%R_polar(1,3,j) = (        nwpts%R_polar(1,1,j)  +  2.d0 * nwpts%R_polar(1,4,j) ) / 3.d0

  nwpts%Z_polar(1,1,j) = Z_axis
  nwpts%Z_polar(1,4,j) = delta * Z_axis + (1.d0 - delta) * nwpts%Z_sep(j)
  nwpts%Z_polar(1,2,j) = ( 2.d0 * nwpts%Z_polar(1,1,j)  +         nwpts%Z_polar(1,4,j) ) / 3.d0
  nwpts%Z_polar(1,3,j) = (        nwpts%Z_polar(1,1,j)  +  2.d0 * nwpts%Z_polar(1,4,j) ) / 3.d0

  nwpts%R_polar(3,1,j) = nwpts%R_max(j)
  nwpts%R_polar(3,4,j) = delta * nwpts%R_max(j) + (1.d0 - delta) * nwpts%R_sep(j)
  nwpts%R_polar(3,2,j) = ( 2.d0 * nwpts%R_polar(3,1,j)  +         nwpts%R_polar(3,4,j) ) / 3.d0
  nwpts%R_polar(3,3,j) = (        nwpts%R_polar(3,1,j)  +  2.d0 * nwpts%R_polar(3,4,j) ) / 3.d0

  nwpts%Z_polar(3,1,j) = nwpts%Z_max(j)
  nwpts%Z_polar(3,4,j) = delta * nwpts%Z_max(j) + (1.d0 - delta) * nwpts%Z_sep(j)
  nwpts%Z_polar(3,2,j) = ( 2.d0 * nwpts%Z_polar(3,1,j)  +         nwpts%Z_polar(3,4,j) ) / 3.d0
  nwpts%Z_polar(3,3,j) = (        nwpts%Z_polar(3,1,j)  +  2.d0 * nwpts%Z_polar(3,4,j) ) / 3.d0

  nwpts%R_polar(2,1,j) = nwpts%R_polar(1,4,j)
  nwpts%R_polar(2,4,j) = nwpts%R_polar(3,4,j)
  nwpts%R_polar(2,2,j) = ( nwpts%R_polar(2,1,j) +  2.d0 * nwpts%R_sep(j) ) / 3.d0
  nwpts%R_polar(2,3,j) = ( nwpts%R_polar(2,4,j) +  2.d0 * nwpts%R_sep(j) ) / 3.d0

  nwpts%Z_polar(2,1,j) = nwpts%Z_polar(1,4,j)
  nwpts%Z_polar(2,4,j) = nwpts%Z_polar(3,4,j)
  nwpts%Z_polar(2,2,j) = ( nwpts%Z_polar(2,1,j) +  2.d0 * nwpts%Z_sep(j) ) / 3.d0
  nwpts%Z_polar(2,3,j) = ( nwpts%Z_polar(2,4,j) +  2.d0 * nwpts%Z_sep(j) ) / 3.d0

enddo

if(xcase .ne. 2) then
  do j=1,2*n_leg

    n_start = n_tht
    delta = 0.2
    if ( (j .eq. n_leg)   .or. (j .eq. 2*n_leg) )    delta = 0.d0
    if ( (j .eq. n_leg-1) .or. (j .eq. 2*n_leg-1) )  delta = 0.05d0

    nwpts%R_polar(1,1,n_start+j) = nwpts%R_min(n_start+j)
    nwpts%R_polar(1,4,n_start+j) = delta * nwpts%R_min(n_start+j) + (1.d0 - delta) * nwpts%R_sep(n_start+j)
    nwpts%R_polar(1,2,n_start+j) = ( 2.d0 * nwpts%R_polar(1,1,n_start+j)  +	    nwpts%R_polar(1,4,n_start+j) ) / 3.d0
    nwpts%R_polar(1,3,n_start+j) = (        nwpts%R_polar(1,1,n_start+j)  +  2.d0 * nwpts%R_polar(1,4,n_start+j) ) / 3.d0

    nwpts%Z_polar(1,1,n_start+j) = nwpts%Z_min(n_start+j)
    nwpts%Z_polar(1,4,n_start+j) = delta * nwpts%Z_min(n_start+j) + (1.d0 - delta) * nwpts%Z_sep(n_start+j)
    nwpts%Z_polar(1,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(1,1,n_start+j)  +	    nwpts%Z_polar(1,4,n_start+j) ) / 3.d0
    nwpts%Z_polar(1,3,n_start+j) = (        nwpts%Z_polar(1,1,n_start+j)  +  2.d0 * nwpts%Z_polar(1,4,n_start+j) ) / 3.d0

    nwpts%R_polar(3,1,n_start+j) = nwpts%R_max(n_start+j)
    nwpts%R_polar(3,4,n_start+j) = delta * nwpts%R_max(n_start+j) + (1.d0 - delta) * nwpts%R_sep(n_start+j)
    nwpts%R_polar(3,2,n_start+j) = ( 2.d0 * nwpts%R_polar(3,1,n_start+j)  +	    nwpts%R_polar(3,4,n_start+j) ) / 3.d0
    nwpts%R_polar(3,3,n_start+j) = (        nwpts%R_polar(3,1,n_start+j)  +  2.d0 * nwpts%R_polar(3,4,n_start+j) ) / 3.d0

    nwpts%Z_polar(3,1,n_start+j) = nwpts%Z_max(n_start+j)
    nwpts%Z_polar(3,4,n_start+j) = delta * nwpts%Z_max(n_start+j) + (1.d0 - delta) * nwpts%Z_sep(n_start+j)
    nwpts%Z_polar(3,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(3,1,n_start+j)  +         nwpts%Z_polar(3,4,n_start+j) ) / 3.d0
    nwpts%Z_polar(3,3,n_start+j) = (        nwpts%Z_polar(3,1,n_start+j)  +  2.d0 * nwpts%Z_polar(3,4,n_start+j) ) / 3.d0

    nwpts%R_polar(2,1,n_start+j) = nwpts%R_polar(1,4,n_start+j)
    nwpts%R_polar(2,4,n_start+j) = nwpts%R_polar(3,4,n_start+j)
    nwpts%R_polar(2,2,n_start+j) = ( nwpts%R_polar(2,1,n_start+j) +  2.d0 * nwpts%R_sep(n_start+j) ) / 3.d0
    nwpts%R_polar(2,3,n_start+j) = ( nwpts%R_polar(2,4,n_start+j) +  2.d0 * nwpts%R_sep(n_start+j) ) / 3.d0

    nwpts%Z_polar(2,1,n_start+j) = nwpts%Z_polar(1,4,n_start+j)
    nwpts%Z_polar(2,4,n_start+j) = nwpts%Z_polar(3,4,n_start+j)
    nwpts%Z_polar(2,2,n_start+j) = ( nwpts%Z_polar(2,1,n_start+j) +  2.d0 * nwpts%Z_sep(n_start+j) ) / 3.d0
    nwpts%Z_polar(2,3,n_start+j) = ( nwpts%Z_polar(2,4,n_start+j) +  2.d0 * nwpts%Z_sep(n_start+j) ) / 3.d0

  enddo
endif

if(xcase .ne. 1) then
  do j=1,2*n_up_leg

    delta = 0.2
    if ( (j .eq. n_up_leg)   .or. (j .eq. 2*n_up_leg) )    delta = 0.d0
    if ( (j .eq. n_up_leg-1) .or. (j .eq. 2*n_up_leg-1) )  delta = 0.05d0
    if(xcase .eq. 2) n_start = n_tht
    if(xcase .eq. 3) n_start = n_tht+2*n_leg

    nwpts%R_polar(1,1,n_start+j) = nwpts%R_min(n_start+j)
    nwpts%R_polar(1,4,n_start+j) = delta * nwpts%R_min(n_start+j) + (1.d0 - delta) * nwpts%R_sep(n_start+j)
    nwpts%R_polar(1,2,n_start+j) = ( 2.d0 * nwpts%R_polar(1,1,n_start+j)  +	    nwpts%R_polar(1,4,n_start+j) ) / 3.d0
    nwpts%R_polar(1,3,n_start+j) = (        nwpts%R_polar(1,1,n_start+j)  +  2.d0 * nwpts%R_polar(1,4,n_start+j) ) / 3.d0

    nwpts%Z_polar(1,1,n_start+j) = nwpts%Z_min(n_start+j)
    nwpts%Z_polar(1,4,n_start+j) = delta * nwpts%Z_min(n_start+j) + (1.d0 - delta) * nwpts%Z_sep(n_start+j)
    nwpts%Z_polar(1,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(1,1,n_start+j)  +	    nwpts%Z_polar(1,4,n_start+j) ) / 3.d0
    nwpts%Z_polar(1,3,n_start+j) = (        nwpts%Z_polar(1,1,n_start+j)  +  2.d0 * nwpts%Z_polar(1,4,n_start+j) ) / 3.d0

    nwpts%R_polar(3,1,n_start+j) = nwpts%R_max(n_start+j)
    nwpts%R_polar(3,4,n_start+j) = delta * nwpts%R_max(n_start+j) + (1.d0 - delta) * nwpts%R_sep(n_start+j)
    nwpts%R_polar(3,2,n_start+j) = ( 2.d0 * nwpts%R_polar(3,1,n_start+j)  +	    nwpts%R_polar(3,4,n_start+j) ) / 3.d0
    nwpts%R_polar(3,3,n_start+j) = (        nwpts%R_polar(3,1,n_start+j)  +  2.d0 * nwpts%R_polar(3,4,n_start+j) ) / 3.d0

    nwpts%Z_polar(3,1,n_start+j) = nwpts%Z_max(n_start+j)
    nwpts%Z_polar(3,4,n_start+j) = delta * nwpts%Z_max(n_start+j) + (1.d0 - delta) * nwpts%Z_sep(n_start+j)
    nwpts%Z_polar(3,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(3,1,n_start+j)  +         nwpts%Z_polar(3,4,n_start+j) ) / 3.d0
    nwpts%Z_polar(3,3,n_start+j) = (        nwpts%Z_polar(3,1,n_start+j)  +  2.d0 * nwpts%Z_polar(3,4,n_start+j) ) / 3.d0

    nwpts%R_polar(2,1,n_start+j) = nwpts%R_polar(1,4,n_start+j)
    nwpts%R_polar(2,4,n_start+j) = nwpts%R_polar(3,4,n_start+j)
    nwpts%R_polar(2,2,n_start+j) = ( nwpts%R_polar(2,1,n_start+j) +  2.d0 * nwpts%R_sep(n_start+j) ) / 3.d0
    nwpts%R_polar(2,3,n_start+j) = ( nwpts%R_polar(2,4,n_start+j) +  2.d0 * nwpts%R_sep(n_start+j) ) / 3.d0

    nwpts%Z_polar(2,1,n_start+j) = nwpts%Z_polar(1,4,n_start+j)
    nwpts%Z_polar(2,4,n_start+j) = nwpts%Z_polar(3,4,n_start+j)
    nwpts%Z_polar(2,2,n_start+j) = ( nwpts%Z_polar(2,1,n_start+j) +  2.d0 * nwpts%Z_sep(n_start+j) ) / 3.d0
    nwpts%Z_polar(2,3,n_start+j) = ( nwpts%Z_polar(2,4,n_start+j) +  2.d0 * nwpts%Z_sep(n_start+j) ) / 3.d0

  enddo
endif


call lincol(3)

npl = 11
call tr_allocate(xp,1,npl,"xp",CAT_GRID)
call tr_allocate(yp,1,npl,"yp",CAT_GRID)
do j=1,n_tht_2

  do m=1,n_pieces

    do k=1,npl
      ss = -1. + 2.*float(k-1)/float(npl-1)

      R_cub1d = (/ nwpts%R_polar(m,1,j), 3.d0/2.d0 *(nwpts%R_polar(m,2,j)-nwpts%R_polar(m,1,j)), &
                   nwpts%R_polar(m,4,j), 3.d0/2.d0 *(nwpts%R_polar(m,4,j)-nwpts%R_polar(m,3,j))  /)
      Z_cub1d = (/ nwpts%Z_polar(m,1,j), 3.d0/2.d0 *(nwpts%Z_polar(m,2,j)-nwpts%Z_polar(m,1,j)), &
                   nwpts%Z_polar(m,4,j), 3.d0/2.d0 *(nwpts%Z_polar(m,4,j)-nwpts%Z_polar(m,3,j)) /)

      call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),ss,xp(k), tmp1)
      call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),ss,yp(k), tmp2)
    enddo
    call lincol(3)
    write(51,*) ' .1 setlinewidth'
    call lplot6(1,1,xp,yp,-npl,' ')
    write(51,*) ' stroke'

  enddo

  call lincol(0)

enddo
call tr_deallocate(xp,"xp",CAT_GRID)
call tr_deallocate(yp,"yp",CAT_GRID)


!----------------------------------- find grid_points from crossing of coordinate lines

do j=1, n_tht          ! the magnetic axis

  nwpts%RR_new(1,j)    = R_axis
  nwpts%ZZ_new(1,j)    = Z_axis
  nwpts%ielm_flux(1,j) = i_elm_axis
  nwpts%s_flux(1,j)    = s_axis
  nwpts%t_flux(1,j)    = t_axis
  nwpts%t_tht(1,j)     = -1.d0          ! expressed in cubic Hermite (-1<t<+1)

enddo

nwpts%k_cross(1,:) = 1

do i=1,n_flux+n_open+n_outer+n_inner-1        ! The main part (without privates) ! We avoid last surface, points are already known...
  do j=1, n_tht
    do k=1,n_pieces       ! 3 line pieces per coordinate line

      R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
                   nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
      Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
                   nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

      call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
                       nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail)


!      if((ifail .eq. 0) .and. (i .eq. n_flux_2+n_open_2+n_outer_2+n_inner_2)) write(*,*)'diff = ',nwpts%RR_new(i+1,j)-nwpts%R_max(j),nwpts%ZZ_new(i+1,j)-nwpts%Z_max(j)
      if (ifail .eq. 0) then
        nwpts%k_cross(i+1,j) = k
        exit
      endif

    enddo

    if (ifail .ne. 0) then
      if (i .eq. n_flux+n_open+n_outer+n_inner) then
        nwpts%k_cross(i+1,j) = 3
	nwpts%RR_new(i+1,j)  = nwpts%R_max(j)
	nwpts%ZZ_new(i+1,j)  = nwpts%Z_max(j)
        write(*,*) ' WARNING node not found for last openflux surface -> using RZ_max '
      else
        write(*,*) ' WARNING node not found for central grid (without legs) : ',ifail,i,j,theta_sep(j)
      endif
    endif

  enddo
enddo
!----------------------------------- Print a python file that plots the bound points
open(100,file='plot_bound_points.py')
  write(100,'(A)')	   '#!/usr/bin/env python'
  write(100,'(A)')	   'import numpy as N'
  write(100,'(A)')	   'import pylab'
  write(100,'(A)')	   'def main():'
  write(100,'(A,i6,A)')     ' r = N.zeros(',(n_flux+1)*n_tht,')'
  write(100,'(A,i6,A)')     ' z = N.zeros(',(n_flux+1)*n_tht,')'
  do i=1,n_flux+1
  do j=1,n_tht
    write(100,'(A,i6,A,f)') ' r[',(i-1)*n_tht+j-1,'] = ',nwpts%RR_new(i+1,j)
    write(100,'(A,i6,A,f)') ' z[',(i-1)*n_tht+j-1,'] = ',nwpts%ZZ_new(i+1,j)
  enddo
  enddo
  write(100,'(A,i6,A)')     ' for i in range (0,',(n_flux+1)*n_tht,'):'
  write(100,'(A)')	   '  pylab.plot(r[i:i+1],z[i:i+1], "r.")'
  write(100,'(A)')	   ' pylab.axis("equal")'
  write(100,'(A)')	   ' pylab.show()'
  write(100,'(A)')	   ' '
  write(100,'(A)')	   'main()'
close(100)

if(xcase .ne. 3) then
  do i=n_flux,n_psi-1          ! With the private parts
    do j=n_tht+1, n_tht_2
      do k=1,n_pieces	    ! 3 line pieces per coordinate line

    	R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
    		     nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
    	Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
    		     nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

    	call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
    			   nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail)

    	if (ifail .eq. 0) then
    	  nwpts%k_cross(i+1,j) = k
    	  exit
    	endif

      enddo

      if ( (ifail .ne. 0) .and. (xcase .eq. 1) ) write(*,*) ' WARNING node not found for lower part of grid : ',ifail,i,j
      if ( (ifail .ne. 0) .and. (xcase .eq. 2) ) write(*,*) ' WARNING node not found for upper part of grid : ',ifail,i,j

    enddo
  enddo
else
  do i=n_flux,n_flux+n_open        ! The sandwich parts
    if(psi_xpoint(2) .lt. psi_xpoint(1))then
      n_start = n_tht+2*n_leg+1
      n_loop  = n_tht_2
    else
      n_start = n_tht+1
      n_loop  = n_tht_2-2*n_up_leg
    endif
    do j=n_start,n_loop
      do k=1,n_pieces	    ! 3 line pieces per coordinate line

    	R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
    		     nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
    	Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
    		     nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

    	call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
    			   nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail)

    	if (ifail .eq. 0) then
    	  nwpts%k_cross(i+1,j) = k
    	  exit
    	endif

      enddo

      if (ifail .ne. 0) write(*,*) ' WARNING node not found for the sandwich part of grid : ',ifail,i,j

    enddo
  enddo
  do i=n_flux+n_open,n_psi-n_private-n_up_priv-1          ! Outer/Inner parts
    do j=n_tht+1, n_tht_2
      do k=1,n_pieces	    ! 3 line pieces per coordinate line

    	R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
    		     nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
    	Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
    		     nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

    	call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
    			   nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail)

    	if (ifail .eq. 0) then
    	  nwpts%k_cross(i+1,j) = k
    	  exit
    	endif

      enddo

      if (ifail .ne. 0) write(*,*) ' WARNING node not found for the outer/inner part of the legs : ',ifail,i,j

    enddo
  enddo
  do i=n_psi-n_private-n_up_priv,n_psi-n_up_priv-1          ! Lower private
    do j=n_tht+1, n_tht_2-2*n_up_leg
      do k=1,n_pieces	    ! 3 line pieces per coordinate line

    	R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
    		     nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
    	Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
    		     nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

    	call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
    			   nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail)

    	if (ifail .eq. 0) then
    	  nwpts%k_cross(i+1,j) = k
    	  exit
    	endif

      enddo

      if (ifail .ne. 0) write(*,*) ' WARNING node not found lower private part of the grid : ',ifail,i,j

    enddo
  enddo
  do i=n_psi-n_up_priv,n_psi-1          !Upper private 
    do j=n_tht+2*n_leg+1, n_tht_2
      do k=1,n_pieces	    ! 3 line pieces per coordinate line

    	R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
    		     nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
    	Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
    		     nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

    	call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
    			   nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail)

    	if (ifail .eq. 0) then
    	  nwpts%k_cross(i+1,j) = k
    	  exit
    	endif

      enddo

      if (ifail .ne. 0) write(*,*) ' WARNING node not found upper private part of the grid : ',ifail,i,j

    enddo
  enddo
endif




!----------------------------------- Print a python file that plots the bound points
!open(100,file='plot_bound_points.py')
!  write(100,'(A)')	    '#!/usr/bin/env python'
!  write(100,'(A)')	    'import numpy as N'
!  write(100,'(A)')	    'import pylab'
!  write(100,'(A)')	    'def main():'
!  write(100,'(A,i6,A)')     ' r = N.zeros(',2*n_tht_2+2*n_leg+2*n_up_leg,')'
!  write(100,'(A,i6,A)')     ' z = N.zeros(',2*n_tht_2+2*n_leg+2*n_up_leg,')'
!  do j=1,n_tht_2
!    write(100,'(A,i6,A,f)') ' r[',j-1,'] = ',nwpts%R_sep(j)
!    write(100,'(A,i6,A,f)') ' z[',j-1,'] = ',nwpts%Z_sep(j)
!  enddo
!  do j=1,n_tht_2
!    write(100,'(A,i6,A,f)') ' r[',n_tht_2+j-1,'] = ',nwpts%R_max(j)
!    write(100,'(A,i6,A,f)') ' z[',n_tht_2+j-1,'] = ',nwpts%Z_max(j)
!  enddo
!  do j=1,2*n_leg
!    write(100,'(A,i6,A,f)') ' r[',2*n_tht_2+j-1,'] = ',nwpts%R_min(j)
!    write(100,'(A,i6,A,f)') ' z[',2*n_tht_2+j-1,'] = ',nwpts%Z_min(j)
!  enddo
!  do j=1,2*n_up_leg
!    write(100,'(A,i6,A,f)') ' r[',2*n_tht_2+2*n_leg+j-1,'] = ',nwpts%R_min(2*n_leg+j)
!    write(100,'(A,i6,A,f)') ' z[',2*n_tht_2+2*n_leg+j-1,'] = ',nwpts%Z_min(2*n_leg+j)
!  enddo
!  write(100,'(A,i6,A)')     ' for i in range (0,',2*n_tht_2+2*n_leg+2*n_up_leg,'):'
!  write(100,'(A)')	    '  pylab.plot(r[i:i+1],z[i:i+1], "r.")'
!  write(100,'(A)')	    ' pylab.axis("equal")'
!  write(100,'(A)')	    ' pylab.show()'
!  write(100,'(A)')	    ' '
!  write(100,'(A)')	    'main()'
!close(100)

!----------------------------------- Print a python file that plots the extrapolation points
!open(101,file='plot_extra_points.py')
!  write(101,'(A)')	    '#!/usr/bin/env python'
!  write(101,'(A)')	    'import numpy as N'
!  write(101,'(A)')	    'import pylab'
!  write(101,'(A)')	    'def main():'
!  write(101,'(A,i6,A)')     ' r = N.zeros(',3*n_tht_2,')'
!  write(101,'(A,i6,A)')     ' z = N.zeros(',3*n_tht_2,')'
!  do j=1,n_tht
!    write(101,'(A,i6,A,f)') ' r[',3*(j-1),'] = ',R_axis
!    write(101,'(A,i6,A,f)') ' z[',3*(j-1),'] = ',Z_axis
!    write(101,'(A,i6,A,f)') ' r[',3*(j-1)+1,'] = ',nwpts%R_sep(j)
!    write(101,'(A,i6,A,f)') ' z[',3*(j-1)+1,'] = ',nwpts%Z_sep(j)
!    write(101,'(A,i6,A,f)') ' r[',3*(j-1)+2,'] = ',nwpts%RR_new(n_psi,j)
!    write(101,'(A,i6,A,f)') ' z[',3*(j-1)+2,'] = ',nwpts%ZZ_new(n_psi,j)
!  enddo
!  do j=n_tht+1,n_tht_2
!    write(101,'(A,i6,A,f)') ' r[',3*(j-1),'] = ',nwpts%R_min(j)
!    write(101,'(A,i6,A,f)') ' z[',3*(j-1),'] = ',nwpts%Z_min(j)
!    write(101,'(A,i6,A,f)') ' r[',3*(j-1)+1,'] = ',nwpts%R_sep(j)
!    write(101,'(A,i6,A,f)') ' z[',3*(j-1)+1,'] = ',nwpts%Z_sep(j)
!    write(101,'(A,i6,A,f)') ' r[',3*(j-1)+2,'] = ',nwpts%R_max(j)
!    write(101,'(A,i6,A,f)') ' z[',3*(j-1)+2,'] = ',nwpts%Z_max(j)
!  enddo
!  write(101,'(A,i6,A)')     ' for i in range (0,',n_tht_2,'):'
!  write(101,'(A)')	    '  pylab.plot(r[3*i:3*i+3],z[3*i:3*i+3], "r")'
!  write(101,'(A)')	    '  pylab.plot(r[3*i:3*i+3],z[3*i:3*i+3], "r+")'
!  write(101,'(A)')	    ' pylab.axis("equal")'
!  write(101,'(A)')	    ' pylab.show()'
!  write(101,'(A)')	    ' '
!  write(101,'(A)')	    'main()'
!close(101)

!----------------------------------- Print a python file that plots the new grid points
!open(102,file='plot_new_points.py')
!  write(102,'(A)')	   '#!/usr/bin/env python'
!  write(102,'(A)')	   'import numpy as N'
!  write(102,'(A)')	   'import pylab'
!  write(102,'(A)')	   'def main():'
!  write(102,'(A,i6,A)')     ' r = N.zeros(',(n_psi)*n_tht_2,')'
!  write(102,'(A,i6,A)')     ' z = N.zeros(',(n_psi)*n_tht_2,')'
!  do i=1,n_psi
!    do j=1,n_tht_2
!    write(102,'(A,i6,A,f)') ' r[',(i-1)*(n_tht_2)+j-1,'] = ',nwpts%RR_new(i,j)
!    write(102,'(A,i6,A,f)') ' z[',(i-1)*(n_tht_2)+j-1,'] = ',nwpts%ZZ_new(i,j)
!    enddo
!  enddo
!  write(102,'(A,i6,A,i6,A)')' pylab.plot(r[0:',(n_psi)*n_tht_2,'],z[0:',(n_psi)*n_tht_2,'], "r.")'
!  write(102,'(A)')	   ' pylab.axis("equal")'
!  write(102,'(A)')	   ' pylab.show()'
!  write(102,'(A)')	   ' '
!  write(102,'(A)')	   'main()'
!close(102)

return
end subroutine define_new_grid_points
