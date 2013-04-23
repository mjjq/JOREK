subroutine find_strategic_points(node_list, element_list, flux_list, &
                                 xcase, R_xpoint, Z_xpoint, psi_xpoint, R_axis, Z_axis, n_grids, stpts)
!----------------------------------------------------------------------------------------
! subroutine finds all the strategic points on the legs (Leg corners, strike points etc.)
!----------------------------------------------------------------------------------------

use constants
use tr_module 
use data_structure
use grid_xpoint_data

implicit none

! --- Routine parameters
type (type_surface_list),     intent(inout) :: flux_list
type (type_node_list),        intent(inout) :: node_list
type (type_element_list),     intent(inout) :: element_list
type (type_strategic_points), intent(inout) :: stpts
integer,                      intent(in)    :: n_grids(10) 
integer,                      intent(in)    :: xcase  
real*8,                       intent(in)    :: psi_xpoint(2), R_xpoint(2), Z_xpoint(2), R_axis, Z_axis


! --- local variables
integer  :: i, k, l, i_elm, i_surf, i_find, i_elm_find(8) , i_max
integer  :: n_flux,   n_open,   n_outer,   n_inner,   n_private,   n_up_priv  
real*8   :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8   :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8   :: rr1, ss1, s_find(8), t_find(8)
real*8   :: tht_x, tht_x1, tht_x2
real*8   :: angle_LowerCorner 
real*8   :: angle_UpperCorner

n_flux    = n_grids(1)
n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
n_private = n_grids(6); n_up_priv = n_grids(7)

write(*,*) '****************************************'
write(*,*) '* X-point grid : Find strategic points *'
write(*,*) '****************************************'


!--------------------------------------------------------------------------------------------
!-------------------------------- Define all strategic points -------------------------------
!--------------------------------------------------------------------------------------------
if(xcase .eq. 2) then
  stpts%RLeftCorn_LowerInnerLeg  = 0.d0;    stpts%ZLeftCorn_LowerInnerLeg  = 0.d0	   
  stpts%RRightCorn_LowerInnerLeg = 0.d0;    stpts%ZRightCorn_LowerInnerLeg = 0.d0	   
  stpts%RLeftCorn_LowerOuterLeg  = 0.d0;    stpts%ZLeftCorn_LowerOuterLeg  = 0.d0	   
  stpts%RRightCorn_LowerOuterLeg = 0.d0;    stpts%ZRightCorn_LowerOuterLeg = 0.d0	   
  stpts%RStrike_LowerInnerLeg	 = 0.d0;     stpts%ZStrike_LowerInnerLeg   = 0.d0 	     
  stpts%RStrike_LowerOuterLeg	 = 0.d0;     stpts%ZStrike_LowerOuterLeg   = 0.d0 	     
else
  stpts%RLeftCorn_LowerInnerLeg  = 999.d0;  stpts%ZLeftCorn_LowerInnerLeg  = 1.d10	
  stpts%RRightCorn_LowerInnerLeg = 1.d10;   stpts%ZRightCorn_LowerInnerLeg = 999.d0	
  stpts%RLeftCorn_LowerOuterLeg  = -1.d10;  stpts%ZLeftCorn_LowerOuterLeg  = 999.d0	
  stpts%RRightCorn_LowerOuterLeg = 999.d0;  stpts%ZRightCorn_LowerOuterLeg = 1.d10	
  stpts%RStrike_LowerInnerLeg	 = 999.d0;   stpts%ZStrike_LowerInnerLeg   = 1.d10	
  stpts%RStrike_LowerOuterLeg	 = 999.d0;   stpts%ZStrike_LowerOuterLeg   = 1.d10	
endif

if(xcase .eq. 1) then
  stpts%RLeftCorn_UpperInnerLeg  = 0.d0;    stpts%ZLeftCorn_UpperInnerLeg  = 0.d0	   
  stpts%RRightCorn_UpperInnerLeg = 0.d0;    stpts%ZRightCorn_UpperInnerLeg = 0.d0	   
  stpts%RLeftCorn_UpperOuterLeg  = 0.d0;    stpts%ZLeftCorn_UpperOuterLeg  = 0.d0	   
  stpts%RRightCorn_UpperOuterLeg = 0.d0;    stpts%ZRightCorn_UpperOuterLeg = 0.d0	   
  stpts%RStrike_UpperInnerLeg	 = 0.d0;     stpts%ZStrike_UpperInnerLeg   = 0.d0 	
  stpts%RStrike_UpperOuterLeg	 = 0.d0;     stpts%ZStrike_UpperOuterLeg   = 0.d0 	
else
  stpts%RLeftCorn_UpperInnerLeg  = 999.d0;  stpts%ZLeftCorn_UpperInnerLeg  = -1.d10
  stpts%RRightCorn_UpperInnerLeg = 1.d10;   stpts%ZRightCorn_UpperInnerLeg = 999.d0
  stpts%RLeftCorn_UpperOuterLeg  = -1.d10;  stpts%ZLeftCorn_UpperOuterLeg  = 999.d0
  stpts%RRightCorn_UpperOuterLeg = 999.d0;  stpts%ZRightCorn_UpperOuterLeg = -1.d10
  stpts%RStrike_UpperInnerLeg	 = 999.d0;   stpts%ZStrike_UpperInnerLeg   = -1.d10
  stpts%RStrike_UpperOuterLeg	 = 999.d0;   stpts%ZStrike_UpperOuterLeg   = -1.d10
endif

if(xcase .ne. 3) then
  stpts%RSecondStrike_InnerLeg   = 0.d0;    stpts%ZSecondStrike_InnerLeg   = 0.d0       
  stpts%RSecondStrike_OuterLeg   = 0.d0;    stpts%ZSecondStrike_OuterLeg   = 0.d0       
else
  if (psi_xpoint(1) .lt. psi_xpoint(2)) then
    stpts%RSecondStrike_InnerLeg = 999.d0;  stpts%ZSecondStrike_InnerLeg   = 1.d10   
    stpts%RSecondStrike_OuterLeg = 999.d0;  stpts%ZSecondStrike_OuterLeg   = 1.d10   
  else
    stpts%RSecondStrike_InnerLeg = 999.d0;  stpts%ZSecondStrike_InnerLeg   = -1.d10
    stpts%RSecondStrike_OuterLeg = 999.d0;  stpts%ZSecondStrike_OuterLeg   = -1.d10  
  endif
endif

stpts%RMiddle_LowerPrivate	 = 999.d0;  stpts%ZMiddle_LowerPrivate	   = 1.d10
stpts%RMiddle_UpperPrivate	 = 999.d0;  stpts%ZMiddle_UpperPrivate	   = -1.d1



!--------------------------------------------------------------------------------------------
!-------------------------------- Now find all the points -----------------------------------
!--------------------------------------------------------------------------------------------

if (xcase .ne. 3) then
  ! ---------------------------------- The last open flux surface (SOL boundary)
  i_surf = n_flux+n_open 
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces    
    do l=1,3,2
      
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
        						  ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if (xcase .eq. 1) then
        if ((ZZg1 .lt. stpts%ZLeftCorn_LowerInnerLeg) .and. (RRg1 .lt. R_xpoint(1))) then
          stpts%RLeftCorn_LowerInnerLeg = RRg1
          stpts%ZLeftCorn_LowerInnerLeg = ZZg1
        endif
        if ((ZZg1 .lt. stpts%ZRightCorn_LowerOuterLeg) .and. (RRg1 .gt. R_xpoint(1))) then
          stpts%RRightCorn_LowerOuterLeg = RRg1
          stpts%ZRightCorn_LowerOuterLeg = ZZg1
        endif
      else
        if ((ZZg1 .gt. stpts%ZLeftCorn_UpperInnerLeg) .and. (RRg1 .lt. R_xpoint(2))) then
          stpts%RLeftCorn_UpperInnerLeg = RRg1
          stpts%ZLeftCorn_UpperInnerLeg = ZZg1
        endif
        if ((ZZg1 .gt. stpts%ZRightCorn_UpperOuterLeg) .and. (RRg1 .gt. R_xpoint(2))) then
          stpts%RRightCorn_UpperOuterLeg = RRg1
          stpts%ZRightCorn_UpperOuterLeg = ZZg1
        endif
      endif
      
    enddo
  enddo

else
  ! ---------------------------------- The last open flux surface (SOL boundary) on outer board (LFS)
  i_surf = n_flux+n_open+n_outer  
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces    
    do l=1,3,2
      
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
        						  ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ((ZZg1 .lt. stpts%ZRightCorn_LowerOuterLeg) .and. (RRg1 .gt. R_xpoint(1))) then
        stpts%RRightCorn_LowerOuterLeg = RRg1
        stpts%ZRightCorn_LowerOuterLeg = ZZg1
      endif
      if ((ZZg1 .gt. stpts%ZRightCorn_UpperOuterLeg) .and. (RRg1 .gt. R_xpoint(2))) then
        stpts%RRightCorn_UpperOuterLeg = RRg1
        stpts%ZRightCorn_UpperOuterLeg = ZZg1
      endif
      
    enddo
  enddo

  ! ---------------------------------- The last open flux surface (SOL boundary) on inner board (HFS)
  i_surf = n_flux+n_open+n_outer+n_inner 
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces    
    do l=1,3,2
      
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
        						  ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ((ZZg1 .lt. stpts%ZLeftCorn_LowerInnerLeg) .and. (RRg1 .lt. R_xpoint(1))) then
        stpts%RLeftCorn_LowerInnerLeg = RRg1
        stpts%ZLeftCorn_LowerInnerLeg = ZZg1
      endif
      if ((ZZg1 .gt. stpts%ZLeftCorn_UpperInnerLeg) .and. (RRg1 .lt. R_xpoint(2))) then
        stpts%RLeftCorn_UpperInnerLeg = RRg1
        stpts%ZLeftCorn_UpperInnerLeg = ZZg1
      endif
      
    enddo
  enddo
endif

! ---------------------------------- The last open flux surface (Private boundary) under lower X-point 
if (xcase .ne. 2) then
  i_surf = n_flux+n_open+n_outer+n_inner+n_private  
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces
    do l=1,3,2
    
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
  	  						  ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ((RRg1 .lt. stpts%RRightCorn_LowerInnerLeg) .and. (ZZg1 .lt. Z_xpoint(1))) then
    	stpts%RRightCorn_LowerInnerLeg = RRg1
    	stpts%ZRightCorn_LowerInnerLeg = ZZg1
      endif
      if ((RRg1 .gt. stpts%RLeftCorn_LowerOuterLeg) .and. (ZZg1 .lt. Z_xpoint(1))) then
    	stpts%RLeftCorn_LowerOuterLeg = RRg1
    	stpts%ZLeftCorn_LowerOuterLeg = ZZg1
      endif
    
    enddo
  enddo
endif

! ---------------------------------- The last open flux surface (Private boundary) above upper X-point 
if (xcase .ne. 1) then
  i_surf = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv 
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces
    do l=1,3,2
    
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
        						  ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ((RRg1 .lt. stpts%RRightCorn_UpperInnerLeg) .and. (ZZg1 .gt. Z_xpoint(2))) then
        stpts%RRightCorn_UpperInnerLeg = RRg1
        stpts%ZRightCorn_UpperInnerLeg = ZZg1
      endif
      if ((RRg1 .gt. stpts%RLeftCorn_UpperOuterLeg) .and. (ZZg1 .gt. Z_xpoint(2))) then
        stpts%RLeftCorn_UpperOuterLeg = RRg1
        stpts%ZLeftCorn_UpperOuterLeg = ZZg1
      endif
    
    enddo
  enddo
endif  

! ---------------------------------- Find line from axis to lower X-point and get intersection ZMiddle_LowerPrivate with last private surface
if (xcase .ne. 2) then
  tht_x = atan2(Z_xpoint(1)-Z_axis,R_xpoint(1)-R_axis)
  if (tht_x .lt. 0.d0) tht_x = tht_x + 2.d0 * PI
  i_surf = n_flux+n_open+n_outer+n_inner+n_private
  call find_theta_surface(node_list,element_list,flux_list,i_surf,tht_x,R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)

  do i=1,i_find
    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                                            ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    stpts%RMiddle_LowerPrivate = RRg1
    stpts%ZMiddle_LowerPrivate = ZZg1
    if (stpts%ZMiddle_LowerPrivate .le. Z_xpoint(1)) exit
  enddo
endif

! ---------------------------------- Find line from axis to upper X-point and get intersection ZMiddle_UpperPrivate with last private surface
if (xcase .ne. 1) then
  tht_x = atan2(Z_xpoint(2)-Z_axis,R_xpoint(2)-R_axis)
  if (tht_x .lt. 0.d0) tht_x = tht_x + 2.d0 * PI
  i_surf = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv
  call find_theta_surface(node_list,element_list,flux_list,i_surf,tht_x,R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)

  do i=1,i_find
    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                                            ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    stpts%RMiddle_UpperPrivate = RRg1
    stpts%ZMiddle_UpperPrivate = ZZg1
    if (stpts%ZMiddle_UpperPrivate .ge. Z_xpoint(2)) exit
  enddo
endif

! ---------------------------------- Find lower strike points
if (xcase .ne. 2) then
  if ( (xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) ) then
    i_surf = n_flux + n_open
  else
    i_surf = n_flux
  endif
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces
    do l=1,3,2
    
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
    							  ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ((ZZg1 .lt. stpts%ZStrike_LowerInnerLeg) .and. (RRg1 .lt. R_xpoint(1))) then
    	stpts%RStrike_LowerInnerLeg = RRg1
    	stpts%ZStrike_LowerInnerLeg = ZZg1
      endif
      if ((ZZg1 .lt. stpts%ZStrike_LowerOuterLeg) .and. (RRg1 .gt. R_xpoint(1))) then
    	stpts%RStrike_LowerOuterLeg = RRg1
    	stpts%ZStrike_LowerOuterLeg = ZZg1
      endif

    enddo
  enddo
endif

! ---------------------------------- Find upper strike points
if (xcase .ne. 1) then
  if ( (xcase .eq. 3) .and. (psi_xpoint(1) .lt. psi_xpoint(2)) ) then
    i_surf = n_flux + n_open
  else
    i_surf = n_flux
  endif
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces
    do l=1,3,2
    
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
  							  ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ((ZZg1 .gt. stpts%ZStrike_UpperInnerLeg) .and. (RRg1 .lt. R_xpoint(2))) then
  	stpts%RStrike_UpperInnerLeg = RRg1
  	stpts%ZStrike_UpperInnerLeg = ZZg1
      endif
      if ((ZZg1 .gt. stpts%ZStrike_UpperOuterLeg) .and. (RRg1 .gt. R_xpoint(2))) then
  	stpts%RStrike_UpperOuterLeg = RRg1
  	stpts%ZStrike_UpperOuterLeg = ZZg1
      endif

    enddo
  enddo
endif

! ---------------------------------- Find strike points of second separatrix
if (xcase .eq. 3) then
  i_surf = n_flux+n_open
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces
    do l=1,3,2
    
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
  							  ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if (psi_xpoint(1) .lt. psi_xpoint(2)) then
        if ((ZZg1 .lt. stpts%ZSecondStrike_InnerLeg) .and. (RRg1 .lt. R_xpoint(1))) then
     	  stpts%RSecondStrike_InnerLeg = RRg1
     	  stpts%ZSecondStrike_InnerLeg = ZZg1
        endif
        if ((ZZg1 .lt. stpts%ZSecondStrike_OuterLeg) .and. (RRg1 .gt. R_xpoint(1))) then
     	  stpts%RSecondStrike_OuterLeg = RRg1
     	  stpts%ZSecondStrike_OuterLeg = ZZg1
        endif
      else
        if ((ZZg1 .gt. stpts%ZSecondStrike_InnerLeg) .and. (RRg1 .lt. R_xpoint(2))) then
     	  stpts%RSecondStrike_InnerLeg = RRg1
     	  stpts%ZSecondStrike_InnerLeg = ZZg1
        endif
        if ((ZZg1 .gt. stpts%ZSecondStrike_OuterLeg) .and. (RRg1 .gt. R_xpoint(2))) then
     	  stpts%RSecondStrike_OuterLeg = RRg1
     	  stpts%ZSecondStrike_OuterLeg = ZZg1
        endif
      endif
      
    enddo
  enddo
endif

!----------------------------------- Define the lines separating the central and upper/lower parts of the grid
if (xcase .ne. 2) tht_x1 = atan2(Z_xpoint(1)-Z_axis,R_xpoint(1)-R_axis)
if (xcase .eq. 2) tht_x1 = atan2(Z_xpoint(2)-Z_axis,R_xpoint(2)-R_axis)
if (xcase .eq. 3) tht_x2 = atan2(Z_xpoint(2)-Z_axis,R_xpoint(2)-R_axis)
if (tht_x1 .lt. 0.d0) tht_x1 = tht_x1 + 2.d0*PI
if (tht_x2 .lt. 0.d0) tht_x2 = tht_x2 + 2.d0*PI

if(xcase .ne. 2) then
  
  angle_LowerCorner       = atan2(stpts%ZLeftCorn_LowerInnerLeg-Z_xpoint(1),stpts%RLeftCorn_LowerInnerLeg-R_xpoint(1))
  stpts%angle_LowerLeft   = tht_x1 + 1.5d0*PI
  stpts%angle_LowerRight  = tht_x1 + 0.5d0*PI
  !!! --- Depending on the equilibrium, it may be better to have a horizontal line... (eg. near double-null at JET)
  !!!stpts%angle_LowerLeft   = PI
  !!!stpts%angle_LowerRight  = 0.d0
  !!! --- Depending on the equilibrium, it may be better to have a horizontal line... (eg. near double-null at JET)
  if (angle_LowerCorner .lt. 0.d0)          angle_LowerCorner       = angle_LowerCorner       + 2.d0*PI
  if (stpts%angle_LowerLeft   .gt. 2.d0*PI) stpts%angle_LowerLeft   = stpts%angle_LowerLeft   - 2.d0*PI
  if (stpts%angle_LowerRight  .gt. 2.d0*PI) stpts%angle_LowerRight  = stpts%angle_LowerRight  - 2.d0*PI
  
  if (stpts%angle_LowerLeft .gt. angle_LowerCorner) then   ! check if Limit_LowerInnerLeg is above ZLeftCorn_LowerInnerLeg : if not adjust angle
    stpts%angle_LowerLeft  = angle_LowerCorner - 0.1
    stpts%angle_LowerRight = stpts%angle_LowerLeft - PI
    if (stpts%angle_LowerRight  .lt. 0.d0) stpts%angle_LowerRight  = stpts%angle_LowerRight + 2.d0*PI
  endif
  
endif

if(xcase .ne. 1) then

  angle_UpperCorner       = atan2(stpts%ZLeftCorn_UpperInnerLeg-Z_xpoint(2),stpts%RLeftCorn_UpperInnerLeg-R_xpoint(2))
  stpts%angle_UpperLeft   = tht_x2 + 0.5d0*PI; if(xcase .eq. 2) stpts%angle_UpperLeft   = tht_x1 + 0.5d0*PI
  stpts%angle_UpperRight  = tht_x2 + 1.5d0*PI; if(xcase .eq. 2) stpts%angle_UpperRight  = tht_x1 + 1.5d0*PI
  !!! --- Depending on the equilibrium, it may be better to have a horizontal line... (eg. near double-null at JET)
  !!!stpts%angle_UpperLeft   = PI
  !!!stpts%angle_UpperRight  = 0.d0
  !!! --- Depending on the equilibrium, it may be better to have a horizontal line... (eg. near double-null at JET)
  if (angle_UpperCorner .lt. 0.d0)          angle_UpperCorner       = angle_UpperCorner + 2.d0*PI
  if (stpts%angle_UpperLeft   .gt. 2.d0*PI) stpts%angle_UpperLeft   = stpts%angle_UpperLeft   - 2.d0*PI
  if (stpts%angle_UpperRight  .gt. 2.d0*PI) stpts%angle_UpperRight  = stpts%angle_UpperRight  - 2.d0*PI
  
  if (stpts%angle_UpperLeft .lt. angle_UpperCorner) then   ! check if Limit_UpperInnerLeg is below ZLeftCorn_UpperInnerLeg : if not adjust angle
    stpts%angle_UpperLeft  = angle_UpperCorner + 0.1
    stpts%angle_UpperRight = stpts%angle_UpperLeft + PI
    if (stpts%angle_UpperRight  .gt. 2.d0*PI) stpts%angle_UpperRight  = stpts%angle_UpperRight  - 2.d0*PI
  endif
  
endif

if(xcase .ne. 2) then
  
  i_max = n_flux + n_open + n_outer
  call find_theta_surface(node_list,element_list,flux_list,i_max,stpts%angle_LowerRight,R_xpoint(1),Z_xpoint(1),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_LowerOuterLeg = RRg1
  stpts%ZLimit_LowerOuterLeg = ZZg1

  i_max = n_flux + n_open + n_outer + n_inner
  call find_theta_surface(node_list,element_list,flux_list,i_max,stpts%angle_LowerLeft,R_xpoint(1),Z_xpoint(1),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_LowerInnerLeg = RRg1
  stpts%ZLimit_LowerInnerLeg = ZZg1

endif

if(xcase .ne. 1) then
  
  i_max = n_flux + n_open + n_outer
  call find_theta_surface(node_list,element_list,flux_list,i_max,stpts%angle_UpperRight,R_xpoint(2),Z_xpoint(2),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_UpperOuterLeg = RRg1
  stpts%ZLimit_UpperOuterLeg = ZZg1

  i_max = n_flux + n_open + n_outer + n_inner
  call find_theta_surface(node_list,element_list,flux_list,i_max,stpts%angle_UpperLeft,R_xpoint(2),Z_xpoint(2),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_UpperInnerLeg = RRg1
  stpts%ZLimit_UpperInnerLeg = ZZg1

endif


! ---------------------------------- And print the output
write(*,'(A)')                  ' _________________________________________________________'
write(*,'(A)')                  '|                                                         |'
write(*,'(A)')  		'| LEG POINTS (R,Z)                                        |'
write(*,'(A)')  		'|_________________________________________________________|'
write(*,'(A)')                  '|                                                         |'

if (xcase .ne. 2) then
  write(*,'(A)')                '| Lower Legs : -------------------------------------------|'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Upper corner of the Lower Inner Leg  : (',stpts%RLimit_LowerInnerLeg,    ', ', stpts%ZLimit_LowerInnerLeg,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left corner  of the Lower Inner Leg  : (',stpts%RLeftCorn_LowerInnerLeg, ', ', stpts%ZLeftCorn_LowerInnerLeg, ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Strike point of the Lower Inner Leg  : (',stpts%RStrike_LowerInnerLeg,   ', ', stpts%ZStrike_LowerInnerLeg,   ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right corner of the Lower Inner Leg  : (',stpts%RRightCorn_LowerInnerLeg,', ', stpts%ZRightCorn_LowerInnerLeg,') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Middle of the Lower Private surface  : (',stpts%RMiddle_LowerPrivate,    ', ', stpts%ZMiddle_LowerPrivate ,   ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left corner  of the Lower Outer Leg  : (',stpts%RLeftCorn_LowerOuterLeg, ', ', stpts%ZLeftCorn_LowerOuterLeg, ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Strike point of the Lower Outer Leg  : (',stpts%RStrike_LowerOuterLeg,   ', ', stpts%ZStrike_LowerOuterLeg,   ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right corner of the Lower Outer Leg  : (',stpts%RRightCorn_LowerOuterLeg,', ', stpts%ZRightCorn_LowerOuterLeg,') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Upper corner of the Lower Outer Leg  : (',stpts%RLimit_LowerOuterLeg,    ', ', stpts%ZLimit_LowerOuterLeg,    ') |'
endif

if (xcase .ne. 1) then
  write(*,'(A)')                '| Upper Legs : -------------------------------------------|'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Lower corner of the Upper Inner Leg  : (',stpts%RLimit_UpperInnerLeg,    ', ', stpts%ZLimit_UpperInnerLeg,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left corner  of the Upper Inner Leg  : (',stpts%RLeftCorn_UpperInnerLeg, ', ', stpts%ZLeftCorn_UpperInnerLeg, ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Strike point of the Upper Inner Leg  : (',stpts%RStrike_UpperInnerLeg,   ', ', stpts%ZStrike_UpperInnerLeg,   ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right corner of the Upper Inner Leg  : (',stpts%RRightCorn_UpperInnerLeg,', ', stpts%ZRightCorn_UpperInnerLeg,') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Middle of the Upper Private surface  : (',stpts%RMiddle_UpperPrivate,    ', ', stpts%ZMiddle_UpperPrivate ,   ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left corner  of the Upper Outer Leg  : (',stpts%RLeftCorn_UpperOuterLeg, ', ', stpts%ZLeftCorn_UpperOuterLeg, ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Strike point of the Upper Outer Leg  : (',stpts%RStrike_UpperOuterLeg,   ', ', stpts%ZStrike_UpperOuterLeg,   ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right corner of the Upper Outer Leg  : (',stpts%RRightCorn_UpperOuterLeg,', ', stpts%ZRightCorn_UpperOuterLeg,') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Lower corner of the Upper Outer Leg  : (',stpts%RLimit_UpperOuterLeg,    ', ', stpts%ZLimit_UpperOuterLeg,    ') |'
endif

if (xcase .eq. 3) then
  write(*,'(A)')                '| Secondary Strike Points : ------------------------------|'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left  Strike point of 2nd separatrix : (',stpts%RSecondStrike_InnerLeg,  ', ', stpts%ZSecondStrike_InnerLeg,  ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right Strike point of 2nd separatrix : (',stpts%RSecondStrike_OuterLeg,  ', ', stpts%ZSecondStrike_OuterLeg,  ') |'
endif
write(*,'(A)')                  '|_________________________________________________________|'


return
end subroutine find_strategic_points
