subroutine define_flux_values(node_list, element_list, flux_list, sep_list, &
                              xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_axis, n_grids, sigmas)
  !-----------------------------------------------------------------------
  ! subroutine defines the flux values of the flux surfaces on which the
  ! finite element grid will be aligned
  !-----------------------------------------------------------------------
  
  use tr_module 
  use data_structure
  use grid_xpoint_data
  
  implicit none
  
  ! --- Routine parameters
  type (type_surface_list), intent(inout) :: flux_list, sep_list
  type (type_node_list),    intent(inout) :: node_list
  type (type_element_list), intent(inout) :: element_list
  integer,		    intent(in)    :: n_grids(10), xcase
  real*8,		    intent(in)    :: sigmas(16)
  real*8,		    intent(in)    :: psi_axis, R_xpoint(2), Z_xpoint(2)
  real*8				  :: psi_xpoint(2)
  
  ! --- local variables
  real*8, allocatable :: s_tmp(:)
  integer	      :: i, j, nPieces, i_elm
  integer	      :: n_flux,      n_open
  integer	      :: n_outer,     n_inner
  integer	      :: n_private,   n_up_priv 
  integer	      :: n_leg,       n_up_leg  
  real*8	      :: SIG_closed, SIG_open, SIG_outer, SIG_inner, SIG_private, SIG_up_priv
  real*8	      :: SIG_leg_0, SIG_leg_1, SIG_up_leg_0, SIG_up_leg_1
  real*8	      :: dPSI_open, dPSI_outer, dPSI_inner, dPSI_private, dPSI_up_priv
  real*8	      :: bgf_open, bgf_closed
  real*8	      :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
  real*8	      :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
  real*8	      :: rr, ss, drr, dss, tt
  real*8	      :: rr1, ss1, drr1, dss1
  real*8	      :: rr2, ss2, drr2, dss2
  real*8	      :: psi_bnd, psi_bnd2
  logical	      :: xpoint
  
  SIG_closed   = sigmas(1) 
  SIG_open     = sigmas(3) ; SIG_outer    = sigmas(4) ; SIG_inner = sigmas(5)  
  SIG_private  = sigmas(6) ; SIG_up_priv  = sigmas(7) 
  SIG_leg_0    = sigmas(8) ; SIG_leg_1    = sigmas(9) 
  SIG_up_leg_0 = sigmas(10); SIG_up_leg_1 = sigmas(11)
  dPSI_open    = sigmas(12); dPSI_outer   = sigmas(13); dPSI_inner = sigmas(14)
  dPSI_private = sigmas(15); dPSI_up_priv = sigmas(16)

  n_flux    = n_grids(1)
  n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
  n_private = n_grids(6); n_up_priv = n_grids(7)
  n_leg     = n_grids(8); n_up_leg  = n_grids(9)
  
  bgf_open   = 0.6d0
  bgf_closed = 0.2d0
  
  write(*,*) '*************************************'
  write(*,*) '* X-point grid : Define Flux Values *'
  write(*,*) '*************************************'
  
  !-------------------------------- Define psi_bnd and psi_bn2
  xpoint = .true.
  if(xcase .eq. 1) psi_bnd = psi_xpoint(1)
  if(xcase .eq. 2) psi_bnd = psi_xpoint(2)
  if(xcase .eq. 3) then
    if(psi_xpoint(2) .lt. psi_xpoint(1)) then
      psi_bnd  = psi_xpoint(2)
      psi_bnd2 = psi_xpoint(1)
    else
      psi_bnd  = psi_xpoint(1)
      psi_bnd2 = psi_xpoint(2)  
    endif
    ! If we have a symmetric double-null, force the single separatrix
    if (abs(psi_xpoint(1)-psi_xpoint(2)) .lt. symmetric_threshold) then
      psi_xpoint(1) = (psi_xpoint(1)+psi_xpoint(2))/2.d0
      psi_xpoint(2) = psi_xpoint(1)
      psi_bnd  = psi_xpoint(1)
      psi_bnd2 = psi_bnd  
    endif
  endif
  
  !-------------------------------- Closed flux surfaces
  call tr_allocate(s_tmp,1,n_flux+1,"s_tmp",CAT_GRID)
  s_tmp = 0
  j	= 0
  call meshac2(n_flux+1,s_tmp,1.d0,9999.d0,SIG_closed,9999.d0,bgf_closed,1.0d0)
  do i=1,n_flux
    flux_list%psi_values(i+j) = psi_axis + (psi_bnd - psi_axis) * s_tmp(i+1)**2
  enddo
  flux_list%psi_values(n_flux) = psi_bnd
  call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  
  !-------------------------------- Open flux surfaces (in case of single-null)
  !-------------------------------- OR Sandwich flux surfaces (in case of double-null) - in between the two separatrices
  if (psi_xpoint(1) .ne. psi_xpoint(2)) then ! Ignore in case of symmetric double-null
    call tr_allocate(s_tmp,1,n_open+1,"s_tmp",CAT_GRID)
    s_tmp = 0
    j	  = n_flux
    if(xcase .ne. 3) then
      call meshac2(n_open+1,s_tmp,0.d0,9999.d0,SIG_open,9999.d0,bgf_open,1.0d0)
      do i=1,n_open
  	flux_list%psi_values(i+j) = psi_axis + (psi_bnd - psi_axis) * (1.d0 + dPSI_open*s_tmp(i+1))**2
      enddo
    else
      call meshac2(n_open+1,s_tmp,0.d0,1.d0,SIG_open,SIG_open,0.8d0,1.0d0)
      do i=1,n_open
  	flux_list%psi_values(i+j) = psi_bnd + (psi_bnd2 - psi_bnd) * s_tmp(i+1)
      enddo
      flux_list%psi_values(n_flux+n_open) = psi_bnd2
    endif
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  endif
  
  !-------------------------------- Outer open flux surfaces (in case of double-null)
  if(xcase .eq. 3) then
    call tr_allocate(s_tmp,1,n_outer+1,"s_tmp",CAT_GRID)
    s_tmp = 0
    j	  = n_flux+n_open
    call meshac2(n_outer+1,s_tmp,0.d0,9999.d0,SIG_outer,9999.d0,0.6d0,1.0d0)
    do i=1,n_outer
      flux_list%psi_values(i+j) = psi_axis + (psi_bnd2 - psi_axis) * (1.d0 + dPSI_outer*s_tmp(i+1))**2
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  endif
  
  !-------------------------------- Inner open flux surfaces (in case of double-null)
  if(xcase .eq. 3) then
    call tr_allocate(s_tmp,1,n_inner+1,"s_tmp",CAT_GRID)
    s_tmp = 0
    j	  = n_flux+n_open+n_outer
    call meshac2(n_inner+1,s_tmp,0.d0,9999.d0,SIG_inner,9999.d0,0.6d0,1.0d0)
    do i=1,n_inner
      flux_list%psi_values(i+j) = psi_axis + (psi_bnd2 - psi_axis) * (1.d0 + dPSI_inner*s_tmp(i+1))**2
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  endif
  
  !-------------------------------- Lower private flux surfaces
  if(xcase .ne. 2) then
    call tr_allocate(s_tmp,1,n_private+1,"s_tmp",CAT_GRID)
    s_tmp = 0
    j	  = n_flux+n_open+n_outer+n_inner
    call meshac2(n_private+1,s_tmp,0.d0,9999.d0,SIG_private,9999.d0,0.6d0,1.0d0)
    do i=1,n_private
      flux_list%psi_values(i+j) = psi_axis + (psi_xpoint(1) - psi_axis) * (1.d0 - dPSI_private*s_tmp(i+1))**2
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  endif
  
  !-------------------------------- Upper private flux surfaces
  if(xcase .ne. 1) then
    call tr_allocate(s_tmp,1,n_up_priv+1,"s_tmp",CAT_GRID)
    s_tmp = 0
    j	  = n_flux+n_open+n_outer+n_inner+n_private
    call meshac2(n_up_priv+1,s_tmp,0.d0,9999.d0,SIG_up_priv,9999.d0,0.6d0,1.0d0)
    do i=1,n_up_priv
      flux_list%psi_values(i+j) = psi_axis + (psi_xpoint(2) - psi_axis) * (1.d0 - dPSI_up_priv*s_tmp(i+1))**2
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  endif
  
  !-------------------------------- Plot separatrices with different colours
  sep_list%psi_values(1) = flux_list%psi_values(n_flux)
  sep_list%psi_values(2) = flux_list%psi_values(n_flux+n_open)
  sep_list%psi_values(3) = flux_list%psi_values(n_flux+n_open+n_private+n_up_priv)
  if(xcase .eq. 3) then
    sep_list%psi_values(1) = psi_xpoint(1)
    sep_list%psi_values(2) = psi_xpoint(2)
    sep_list%psi_values(3) = flux_list%psi_values(n_flux+n_open+n_outer)
    sep_list%psi_values(4) = flux_list%psi_values(n_flux+n_open+n_outer+n_inner)
    sep_list%psi_values(5) = flux_list%psi_values(n_flux+n_open+n_outer+n_inner+n_private)
    sep_list%psi_values(6) = flux_list%psi_values(n_flux+n_open+n_outer+n_inner+n_private+n_up_priv)
  endif
  
  call find_flux_surfaces(xpoint,xcase,node_list,element_list,flux_list)
  call find_flux_surfaces(xpoint,xcase,node_list,element_list,sep_list)  
  if(xcase .eq. 3) then
    do i=1,6
      nPieces = sep_list%flux_surfaces(i)%n_pieces
      sep_list%flux_surfaces(i)%n_pieces = 0
      do j=1,nPieces
  	i_elm = sep_list%flux_surfaces(i)%elm(j)
  	rr1   = sep_list%flux_surfaces(i)%s(1,j)
  	drr1  = sep_list%flux_surfaces(i)%s(2,j)
  	rr2   = sep_list%flux_surfaces(i)%s(3,j)
  	drr2  = sep_list%flux_surfaces(i)%s(4,j)
  
  	ss1   = sep_list%flux_surfaces(i)%t(1,j)
  	dss1  = sep_list%flux_surfaces(i)%t(2,j)
  	ss2   = sep_list%flux_surfaces(i)%t(3,j)
  	dss2  = sep_list%flux_surfaces(i)%t(4,j)
  	tt    = 0.d0
  	call CUB1D(rr1, drr1, rr2, drr2, tt, rr, drr)
  	call CUB1D(ss1, dss1, ss2, dss2, tt, ss, dss)
  	call interp_RZ(node_list,element_list,i_elm,rr,ss,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
  							  ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
  	if(i .eq. 1) then
  	  if( (ZZg1 .lt. Z_xpoint(2)) .or. (psi_xpoint(1) .gt. psi_xpoint(2)) ) then
  	    sep_list%flux_surfaces(i)%n_pieces = sep_list%flux_surfaces(i)%n_pieces + 1
  	    sep_list%flux_surfaces(i)%elm(sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%elm(j)
  	    sep_list%flux_surfaces(i)%s(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%s(:,j)
  	    sep_list%flux_surfaces(i)%t(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%t(:,j)
  	  endif
  	endif
  	if(i .eq. 2) then
  	  if( (ZZg1 .gt. Z_xpoint(1)) .or. (psi_xpoint(2) .gt. psi_xpoint(1)) ) then
  	    sep_list%flux_surfaces(i)%n_pieces = sep_list%flux_surfaces(i)%n_pieces + 1
  	    sep_list%flux_surfaces(i)%elm(sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%elm(j)
  	    sep_list%flux_surfaces(i)%s(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%s(:,j)
  	    sep_list%flux_surfaces(i)%t(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%t(:,j)
  	  endif
  	endif
  	if(i .eq. 3) then
  	  if( (RRg1 .gt. R_xpoint(1)) .and. (RRg1 .gt. R_xpoint(2)) ) then
  	    sep_list%flux_surfaces(i)%n_pieces = sep_list%flux_surfaces(i)%n_pieces + 1
  	    sep_list%flux_surfaces(i)%elm(sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%elm(j)
  	    sep_list%flux_surfaces(i)%s(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%s(:,j)
  	    sep_list%flux_surfaces(i)%t(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%t(:,j)
  	  endif
  	endif	   
  	if(i .eq. 4) then
  	  if( (RRg1 .lt. R_xpoint(1)) .and. (RRg1 .lt. R_xpoint(2)) ) then
  	    sep_list%flux_surfaces(i)%n_pieces = sep_list%flux_surfaces(i)%n_pieces + 1
  	    sep_list%flux_surfaces(i)%elm(sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%elm(j)
  	    sep_list%flux_surfaces(i)%s(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%s(:,j)
  	    sep_list%flux_surfaces(i)%t(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%t(:,j)
  	  endif
  	endif	   
  	if(i .eq. 5) then
  	  if(ZZg1 .lt. Z_xpoint(1)) then
  	    sep_list%flux_surfaces(i)%n_pieces = sep_list%flux_surfaces(i)%n_pieces + 1
  	    sep_list%flux_surfaces(i)%elm(sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%elm(j)
  	    sep_list%flux_surfaces(i)%s(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%s(:,j)
  	    sep_list%flux_surfaces(i)%t(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%t(:,j)
  	  endif
  	endif	   
  	if(i .eq. 6) then
  	  if(ZZg1 .gt. Z_xpoint(2)) then
  	    sep_list%flux_surfaces(i)%n_pieces = sep_list%flux_surfaces(i)%n_pieces + 1
  	    sep_list%flux_surfaces(i)%elm(sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%elm(j)
  	    sep_list%flux_surfaces(i)%s(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%s(:,j)
  	    sep_list%flux_surfaces(i)%t(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%t(:,j)
  	  endif
  	endif	   
      enddo
    enddo
  else
    nPieces = sep_list%flux_surfaces(3)%n_pieces
    sep_list%flux_surfaces(3)%n_pieces = 0
    do j=1,nPieces
      i_elm = sep_list%flux_surfaces(3)%elm(j)
      rr1   = sep_list%flux_surfaces(3)%s(1,j)
      drr1  = sep_list%flux_surfaces(3)%s(2,j)
      rr2   = sep_list%flux_surfaces(3)%s(3,j)
      drr2  = sep_list%flux_surfaces(3)%s(4,j)
  
      ss1   = sep_list%flux_surfaces(3)%t(1,j)
      dss1  = sep_list%flux_surfaces(3)%t(2,j)
      ss2   = sep_list%flux_surfaces(3)%t(3,j)
      dss2  = sep_list%flux_surfaces(3)%t(4,j)
      tt    = 0.d0
      call CUB1D(rr1, drr1, rr2, drr2, tt, rr, drr)
      call CUB1D(ss1, dss1, ss2, dss2, tt, ss, dss)
      call interp_RZ(node_list,element_list,i_elm,rr,ss,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
  							ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      if( (ZZg1 .lt. Z_xpoint(1)) .and. (xcase .eq. 1) ) then
  	sep_list%flux_surfaces(3)%n_pieces = sep_list%flux_surfaces(3)%n_pieces + 1
  	sep_list%flux_surfaces(3)%elm(sep_list%flux_surfaces(3)%n_pieces) = sep_list%flux_surfaces(3)%elm(j)
  	sep_list%flux_surfaces(3)%s(:,sep_list%flux_surfaces(3)%n_pieces) = sep_list%flux_surfaces(3)%s(:,j)
  	sep_list%flux_surfaces(3)%t(:,sep_list%flux_surfaces(3)%n_pieces) = sep_list%flux_surfaces(3)%t(:,j)
      endif
      if( (ZZg1 .gt. Z_xpoint(2)) .and. (xcase .eq. 2) ) then
  	sep_list%flux_surfaces(3)%n_pieces = sep_list%flux_surfaces(3)%n_pieces + 1
  	sep_list%flux_surfaces(3)%elm(sep_list%flux_surfaces(3)%n_pieces) = sep_list%flux_surfaces(3)%elm(j)
  	sep_list%flux_surfaces(3)%s(:,sep_list%flux_surfaces(3)%n_pieces) = sep_list%flux_surfaces(3)%s(:,j)
  	sep_list%flux_surfaces(3)%t(:,sep_list%flux_surfaces(3)%n_pieces) = sep_list%flux_surfaces(3)%t(:,j)
      endif
    enddo
  endif
  
  return
end subroutine define_flux_values
  
  
  
  

  
  
  
  
subroutine redefine_flux_values(node_list, element_list, surface_list, xcase, n_grids, &
                                R_xpoint, Z_xpoint, psi_xpoint, Z_axis, psi_axis,      &
                                n_int_max, n_target, R_target, Z_target, index_target, &
				n_surf_max, n_int_surf, index_int_surf)
  !------------------------------------------------------------------------
  ! subroutine redefines the flux values of the flux surfaces to add values
  ! for the wall corners inside the grid
  !------------------------------------------------------------------------
  
  use tr_module 
  use data_structure
  use phys_module, only : n_limiter, R_limiter, Z_limiter
  
  implicit none
  
  ! --- Routine parameters
  type (type_surface_list), intent(inout)	:: surface_list
  type (type_node_list),    intent(in)		:: node_list
  type (type_element_list), intent(in)		:: element_list
  real*8,                   intent(in)		:: R_xpoint(2), Z_xpoint(2), psi_xpoint(2), Z_axis, psi_axis
  integer,		    intent(inout)	:: n_grids(10), xcase
  integer,		    intent(inout)	:: n_int_max, n_target, index_target(n_int_max,4)
  real*8,		    intent(inout)	:: R_target(n_int_max), Z_target(n_int_max)
  integer,                  intent(inout)	:: n_surf_max
  integer,                  intent(inout)	:: n_int_surf(n_surf_max) ! number of intersections for each surface
  integer,                  intent(inout)	:: index_int_surf(n_surf_max,n_int_max) ! index of intersections on surface
  
  ! --- local variables
  integer		:: i, j, k, i_surf, istart, iend
  integer		:: n_flux,      n_open
  integer		:: n_outer,     n_inner
  integer		:: n_private,   n_up_priv 
  integer		:: n_leg,       n_up_leg  
  integer, parameter	:: LowerLeft =1
  integer, parameter	:: LowerRight=2
  integer, parameter	:: UpperLeft =3
  integer, parameter	:: UpperRight=4
  integer		:: n_pieces
  integer		:: i_beg(4), i_end(4)
  integer		:: i_elm, i_int
  real*8		:: ss
  real*8		:: tt
  real*8		:: R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt
  real*8		:: Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt
  real*8		:: Rmin_lower,Rmax_lower
  real*8		:: Rmin_upper,Rmax_upper
  character*256		:: filename
  logical		:: debug
  
  write(*,*) '*******************************************'
  write(*,*) '* X-point grid : Add Flux Values For Wall *'
  write(*,*) '*******************************************'
  
  ! --- Debug?
  debug = .true.
  
  n_flux    = n_grids(1)
  n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
  n_private = n_grids(6); n_up_priv = n_grids(7)
  n_leg     = n_grids(8); n_up_leg  = n_grids(9)
  
  ! --- First find the beginings of targets (in private region)
  Rmax_lower = -1.d10
  Rmin_lower = 1.d10
  Rmax_upper = -1.d10
  Rmin_upper = 1.d10
  do i_surf = 1,surface_list%n_psi
    do i = 1,n_int_surf(i_surf)
      i_int = index_int_surf(i_surf,i)
      ! --- LowerLeft
      if (index_target(i_int,4) .eq. LowerLeft) then
        if (R_target(i_int) .gt. Rmax_lower) then
	  Rmax_lower = R_target(i_int)
	  i_beg(LowerLeft) = index_target(i_int,3)
	endif
      endif
      ! --- LowerRight
      if (index_target(i_int,4) .eq. LowerRight) then
        if (R_target(i_int) .lt. Rmin_lower) then
	  Rmin_lower = R_target(i_int)
	  i_beg(LowerRight) = index_target(i_int,3)
	endif
      endif
      ! --- UpperLeft
      if (index_target(i_int,4) .eq. UpperLeft) then
        if (R_target(i_int) .gt. Rmax_upper) then
	  Rmax_upper = R_target(i_int)
	  i_beg(UpperLeft) = index_target(i_int,3)
	endif
      endif
      ! --- UpperRight
      if (index_target(i_int,4) .eq. UpperRight) then
        if (R_target(i_int) .lt. Rmin_upper) then
	  Rmin_upper = R_target(i_int)
	  i_beg(UpperRight) = index_target(i_int,3)
	endif
      endif
    enddo
  enddo
  
  ! --- Then find end of targets in SOL regions
  if (xcase .ne. 3) then
    do i_surf = n_flux+1,n_flux+n_open
      do i = 1,n_int_surf(i_surf)
        i_int = index_int_surf(i_surf,i)
        if (xcase .eq. 1) then
	  if (index_target(i_int,4) .eq. LowerLeft)  i_end(LowerLeft)  = index_target(i_int,3)
          if (index_target(i_int,4) .eq. LowerRight) i_end(LowerRight) = index_target(i_int,3)
	endif
        if (xcase .eq. 2) then
	  if (index_target(i_int,4) .eq. UpperLeft)  i_end(UpperLeft)  = index_target(i_int,3)
          if (index_target(i_int,4) .eq. UpperRight) i_end(UpperRight) = index_target(i_int,3)
	endif
      enddo
    enddo
  else
    do i_surf = n_flux+n_open+1,n_flux+n_open+n_outer
      do i = 1,n_int_surf(i_surf)
        i_int = index_int_surf(i_surf,i)
        if (index_target(i_int,4) .eq. LowerRight) i_end(LowerRight) = index_target(i_int,3)
        if (index_target(i_int,4) .eq. UpperRight) i_end(UpperRight) = index_target(i_int,3)
      enddo
    enddo
    do i_surf = n_flux+n_open+n_outer+1,n_flux+n_open+n_outer+n_inner
      do i = 1,n_int_surf(i_surf)
        i_int = index_int_surf(i_surf,i)
        if (index_target(i_int,4) .eq. LowerLeft) i_end(LowerLeft) = index_target(i_int,3)
        if (index_target(i_int,4) .eq. UpperLeft) i_end(UpperLeft) = index_target(i_int,3)
      enddo
    enddo
  endif

  ! --- Now we need to add flux surfaces at target corners
  if (xcase .ne. 2) then
    istart = i_beg(LowerLeft)
    iend   = i_end(LowerLeft)
    if (istart .gt. iend) then
      istart = i_end(LowerLeft)
      iend   = i_beg(LowerLeft)
    endif
    do i=istart+1,iend
      call add_flux_surface(node_list, element_list, surface_list, xcase, n_grids, &
                            R_xpoint,Z_xpoint, psi_xpoint, Z_axis, psi_axis, i,    &
                            n_int_max, n_target, R_target, Z_target, index_target, &
			    n_surf_max, n_int_surf, index_int_surf)
      ! --- The grid values need to be updated after we added flux surfaces
      n_flux	= n_grids(1)
      n_open	= n_grids(3); n_outer	= n_grids(4); n_inner = n_grids(5)
      n_private = n_grids(6); n_up_priv = n_grids(7)
      n_leg	= n_grids(8); n_up_leg  = n_grids(9)
    enddo
    istart = i_beg(LowerRight)
    iend   = i_end(LowerRight)
    if (istart .gt. iend) then
      istart = i_end(LowerRight)
      iend   = i_beg(LowerRight)
    endif
    do i=istart+1,iend
      call add_flux_surface(node_list, element_list, surface_list, xcase, n_grids, &
                            R_xpoint,Z_xpoint, psi_xpoint, Z_axis, psi_axis, i,    &
                            n_int_max, n_target, R_target, Z_target, index_target, &
			    n_surf_max, n_int_surf, index_int_surf)
      ! --- The grid values need to be updated after we added flux surfaces
      n_flux	= n_grids(1)
      n_open	= n_grids(3); n_outer	= n_grids(4); n_inner = n_grids(5)
      n_private = n_grids(6); n_up_priv = n_grids(7)
      n_leg	= n_grids(8); n_up_leg  = n_grids(9)
    enddo
  endif
  if (xcase .ne. 1) then
    istart = i_beg(UpperLeft)
    iend   = i_end(UpperLeft)
    if (istart .gt. iend) then
      istart = i_end(UpperLeft)
      iend   = i_beg(UpperLeft)
    endif
    do i=istart+1,iend
      call add_flux_surface(node_list, element_list, surface_list, xcase, n_grids, &
                            R_xpoint,Z_xpoint, psi_xpoint, Z_axis, psi_axis, i,    &
                            n_int_max, n_target, R_target, Z_target, index_target, &
			    n_surf_max, n_int_surf, index_int_surf)
      ! --- The grid values need to be updated after we added flux surfaces
      n_flux	= n_grids(1)
      n_open	= n_grids(3); n_outer	= n_grids(4); n_inner = n_grids(5)
      n_private = n_grids(6); n_up_priv = n_grids(7)
      n_leg	= n_grids(8); n_up_leg  = n_grids(9)
    enddo
    istart = i_beg(UpperRight)
    iend   = i_end(UpperRight)
    if (istart .gt. iend) then
      istart = i_end(UpperRight)
      iend   = i_beg(UpperRight)
    endif
    do i=istart+1,iend
      call add_flux_surface(node_list, element_list, surface_list, xcase, n_grids, &
                            R_xpoint,Z_xpoint, psi_xpoint, Z_axis, psi_axis, i,    &
                            n_int_max, n_target, R_target, Z_target, index_target, &
			    n_surf_max, n_int_surf, index_int_surf)
      ! --- The grid values need to be updated after we added flux surfaces
      n_flux	= n_grids(1)
      n_open	= n_grids(3); n_outer	= n_grids(4); n_inner = n_grids(5)
      n_private = n_grids(6); n_up_priv = n_grids(7)
      n_leg	= n_grids(8); n_up_leg  = n_grids(9)
    enddo
  endif
  
  ! --- Re-find end of targets in SOL regions, but this time save it to the intersections
  if (xcase .ne. 3) then
    do i_surf = n_flux+1,n_flux+n_open
      do i = 1,n_int_surf(i_surf)
        i_int = index_int_surf(i_surf,i)
        if (xcase .eq. 1) then
	  if (index_target(i_int,4) .eq. LowerLeft)  i_end(LowerLeft)  = i_int
          if (index_target(i_int,4) .eq. LowerRight) i_end(LowerRight) = i_int
	endif
        if (xcase .eq. 2) then
	  if (index_target(i_int,4) .eq. UpperLeft)  i_end(UpperLeft)  = i_int
          if (index_target(i_int,4) .eq. UpperRight) i_end(UpperRight) = i_int
	endif
      enddo
    enddo
  else
    do i_surf = n_flux+n_open+1,n_flux+n_open+n_outer
      do i = 1,n_int_surf(i_surf)
        i_int = index_int_surf(i_surf,i)
        if (index_target(i_int,4) .eq. LowerRight) i_end(LowerRight) = i_int
        if (index_target(i_int,4) .eq. UpperRight) i_end(UpperRight) = i_int
      enddo
    enddo
    do i_surf = n_flux+n_open+n_outer+1,n_flux+n_open+n_outer+n_inner
      do i = 1,n_int_surf(i_surf)
        i_int = index_int_surf(i_surf,i)
        if (index_target(i_int,4) .eq. LowerLeft) i_end(LowerLeft) = i_int
        if (index_target(i_int,4) .eq. UpperLeft) i_end(UpperLeft) = i_int
      enddo
    enddo
  endif

  ! --- If the end of the target is not on the last flux surface, we need to add an artificial intersection
  if (xcase .ne. 2) then
    ! --- LowerLeft
    i_surf   = n_flux+n_open
    if (xcase .eq. 3) i_surf = n_flux+n_open+n_outer+n_inner
    if (index_target(i_end(LowerLeft),1) .ne. i_surf) then
      n_pieces = 1
      i_elm    = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
      ss       = surface_list%flux_surfaces(i_surf)%s(1,n_pieces)
      tt       = surface_list%flux_surfaces(i_surf)%t(1,n_pieces)
      call interp_RZ(node_list,element_list,i_elm,ss,tt,      &
    		     R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
    		     Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
      if (R .gt. R_xpoint(1)) then
    	n_pieces = surface_list%flux_surfaces(i_surf)%n_pieces
    	i_elm	 = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
    	ss	 = surface_list%flux_surfaces(i_surf)%s(3,n_pieces)
    	tt	 = surface_list%flux_surfaces(i_surf)%t(3,n_pieces)
    	call interp_RZ(node_list,element_list,i_elm,ss,tt,	&
    		       R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
    		       Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
      endif
      n_target = n_target + 1
      R_target(n_target) = R
      Z_target(n_target) = Z
      index_target(n_target,1) = i_surf
      index_target(n_target,2) = n_pieces
      index_target(n_target,3) = -1 ! --- This is dangerous, but this is an artificial intersection, so it has no index on the wall
      index_target(n_target,4) = LowerLeft
    endif
    ! --- LowerRight
    i_surf   = n_flux+n_open
    if (xcase .eq. 3) i_surf = n_flux+n_open+n_outer
    if (index_target(i_end(LowerRight),1) .ne. i_surf) then
      n_pieces = 1
      i_elm = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
      ss    = surface_list%flux_surfaces(i_surf)%s(1,n_pieces)
      tt    = surface_list%flux_surfaces(i_surf)%t(1,n_pieces)
      call interp_RZ(node_list,element_list,i_elm,ss,tt,      &
    		     R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
    		     Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
      if (R .lt. R_xpoint(1)) then
    	n_pieces = surface_list%flux_surfaces(i_surf)%n_pieces
    	i_elm	 = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
    	ss	 = surface_list%flux_surfaces(i_surf)%s(3,n_pieces)
    	tt	 = surface_list%flux_surfaces(i_surf)%t(3,n_pieces)
    	call interp_RZ(node_list,element_list,i_elm,ss,tt,	&
    		       R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
    		       Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
      endif
      n_target = n_target + 1
      R_target(n_target) = R
      Z_target(n_target) = Z
      index_target(n_target,1) = i_surf
      index_target(n_target,2) = n_pieces
      index_target(n_target,3) = -1 ! --- This is dangerous, but this is an artificial intersection, so it has no index on the wall
      index_target(n_target,4) = LowerRight
    endif
  endif
  if (xcase .ne. 1) then
    ! --- UpperLeft
    i_surf   = n_flux+n_open
    if (xcase .eq. 3) i_surf = n_flux+n_open+n_outer+n_inner
    if (index_target(i_end(UpperLeft),1) .ne. i_surf) then
      n_pieces = 1
      i_elm = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
      ss    = surface_list%flux_surfaces(i_surf)%s(1,n_pieces)
      tt    = surface_list%flux_surfaces(i_surf)%t(1,n_pieces)
      call interp_RZ(node_list,element_list,i_elm,ss,tt,      &
    		     R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
    		     Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
      if (R .gt. R_xpoint(2)) then
    	n_pieces = surface_list%flux_surfaces(i_surf)%n_pieces
    	i_elm	 = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
    	ss	 = surface_list%flux_surfaces(i_surf)%s(3,n_pieces)
    	tt	 = surface_list%flux_surfaces(i_surf)%t(3,n_pieces)
    	call interp_RZ(node_list,element_list,i_elm,ss,tt,	&
    		       R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
    		       Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
      endif
      n_target = n_target + 1
      R_target(n_target) = R
      Z_target(n_target) = Z
      index_target(n_target,1) = i_surf
      index_target(n_target,2) = n_pieces
      index_target(n_target,3) = -1 ! --- This is dangerous, but this is an artificial intersection, so it has no index on the wall
      index_target(n_target,4) = UpperLeft
    endif
    ! --- UpperRight
    i_surf   = n_flux+n_open
    if (xcase .eq. 3) i_surf = n_flux+n_open+n_outer
    if (index_target(i_end(UpperRight),1) .ne. i_surf) then
      n_pieces = 1
      i_elm = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
      ss    = surface_list%flux_surfaces(i_surf)%s(1,n_pieces)
      tt    = surface_list%flux_surfaces(i_surf)%t(1,n_pieces)
      call interp_RZ(node_list,element_list,i_elm,ss,tt,      &
    		     R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
    		     Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
      if (R .lt. R_xpoint(2)) then
    	n_pieces = surface_list%flux_surfaces(i_surf)%n_pieces
    	i_elm	 = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
    	ss	 = surface_list%flux_surfaces(i_surf)%s(3,n_pieces)
    	tt	 = surface_list%flux_surfaces(i_surf)%t(3,n_pieces)
    	call interp_RZ(node_list,element_list,i_elm,ss,tt,	&
    		       R,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt, &
    		       Z,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
      endif
      n_target = n_target + 1
      R_target(n_target) = R
      Z_target(n_target) = Z
      index_target(n_target,1) = i_surf
      index_target(n_target,2) = n_pieces
      index_target(n_target,3) = -1 ! --- This is dangerous, but this is an artificial intersection, so it has no index on the wall
      index_target(n_target,4) = UpperRight
    endif
  endif
  
      
  if (debug) then
    filename = 'plot_added_flux_surfaces.py'
    call print_py_plot_prepare_plot(filename)
    call print_py_plot_ordered_flux_surfaces(filename, node_list, element_list, surface_list)
    call print_py_plot_wall(filename)
    call print_py_plot_points(filename,n_target,R_target,Z_target)
    call print_py_plot_finish_plot(filename)
  endif
  
  
  return
end subroutine redefine_flux_values
  
  
  
  
  
  
  
  
  

subroutine add_flux_surface(node_list, element_list, surface_list, xcase, n_grids, &
                            R_xpoint,Z_xpoint,psi_xpoint, Z_axis,psi_axis, i_wall, &
                            n_int_max, n_target, R_target, Z_target, index_target, &
			    n_surf_max, n_int_surf, index_int_surf)
  !------------------------------------------------------------------------
  ! subroutine redefines the flux values of the flux surfaces to add values
  ! for the wall corners inside the grid
  !------------------------------------------------------------------------
  
  use tr_module 
  use data_structure
  use grid_xpoint_data
  use phys_module, only : n_limiter, R_limiter, Z_limiter
  
  implicit none
  
  ! --- Routine parameters
  type (type_surface_list), intent(inout)	:: surface_list
  type (type_node_list),    intent(in)		:: node_list
  type (type_element_list), intent(in)		:: element_list
  real*8,                   intent(in)		:: R_xpoint(2), Z_xpoint(2), psi_xpoint(2), Z_axis, psi_axis
  integer,                  intent(in)		:: i_wall
  integer,		    intent(inout)	:: n_grids(10), xcase
  integer,		    intent(inout)	:: n_int_max, n_target, index_target(n_int_max,4)
  real*8,		    intent(inout)	:: R_target(n_int_max), Z_target(n_int_max)
  integer,                  intent(inout)	:: n_surf_max
  integer,                  intent(inout)	:: n_int_surf(n_surf_max) ! number of intersections for each surface
  integer,                  intent(inout)	:: index_int_surf(n_surf_max,n_int_max) ! index of intersections on surface
  
  ! --- local variables
  type (type_surface_list)	:: surface_list_tmp
  type (type_surface_list)	:: surface_list_single
  real*8			:: psi_values_tmp(n_surf_max)
  integer			:: n_int_surf_tmp(n_surf_max) ! number of intersections for each surface
  integer			:: index_int_surf_tmp(n_surf_max,n_int_max) ! index of intersections on surface
  integer			:: i, j, i_surf, i_int, i_add, istart, iend, indent, ifail
  integer			:: n_flux,      n_open
  integer			:: n_outer,     n_inner
  integer			:: n_private,   n_up_priv 
  integer			:: n_leg,       n_up_leg  
  integer, parameter		:: LowerLeft =1
  integer, parameter		:: LowerRight=2
  integer, parameter		:: UpperLeft =3
  integer, parameter		:: UpperRight=4
  real*8			:: R,R_out,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt
  real*8			:: Z,Z_out,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt
  real*8			:: s_out, t_out
  integer			:: i_elm_out
  real*8			:: psi,dpsi_ds,dpsi_dt,dpsi_dst,dpsi_dss,dpsi_dtt
  integer			:: location
  logical, parameter		:: xpoint = .true.
  
  
  write(*,*) '*************************************************************'
  write(*,*) '* X-point grid : Adding Flux Value For Wall Point',i_wall
  write(*,*) '*************************************************************'
  
  ! --- Grid values
  n_flux    = n_grids(1)
  n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
  n_private = n_grids(6); n_up_priv = n_grids(7)
  n_leg     = n_grids(8); n_up_leg  = n_grids(9)

  ! --- First find the psi-value of the wall corner
  R = R_limiter(i_wall)
  Z = Z_limiter(i_wall)
  call find_RZ(node_list,element_list,R,Z,R_out,Z_out,i_elm_out,s_out,t_out,ifail)
  if (ifail .ne. 0) then
    write(*,*)'Warning! Failed to find RZ position of wall corner ',i_wall,R,Z
    write(*,*)'         on grid when adding flux surface. Aborting...'
    return
  endif
  call interp_RZ(node_list,element_list,i_elm_out,s_out,t_out, &
    		 R_out,dRR_ds,dRR_dt,dRR_dst,dRR_dss,dRR_dtt,  &
        	 Z_out,dZZ_ds,dZZ_dt,dZZ_dst,dZZ_dss,dZZ_dtt)
  call interp(node_list,element_list,i_elm_out,1,1,s_out,t_out,psi,dpsi_ds,dpsi_dt,dpsi_dst,dpsi_dss,dpsi_dtt)

  ! --- Determine which region of the grid the wall corner belongs to
  if ( (xcase .ne. 2) .and. (psi .lt. psi_xpoint(1)) &
                      .and. (Z_out .lt. Z_axis)       			) location = private
  if ( (xcase .ne. 1) .and. (psi .lt. psi_xpoint(2)) &
                      .and. (Z_out .gt. Z_axis)       			) location = upper_private
  if ( (xcase .eq. 1) .and. (psi .gt. psi_xpoint(1)) &
                      .and. (Z_out .lt. Z_axis)       			) location = SOL
  if ( (xcase .eq. 2) .and. (psi .gt. psi_xpoint(2)) &
                      .and. (Z_out .gt. Z_axis)       			) location = SOL
  if ( (xcase .eq. 3) .and. (psi .gt. psi_xpoint(1)) &
                      .and. (psi .lt. psi_xpoint(2)) 			) location = sandwich
  if ( (xcase .eq. 3) .and. (psi .lt. psi_xpoint(1)) &
                      .and. (psi .gt. psi_xpoint(2)) 			) location = sandwich
  if ( (xcase .eq. 3) .and. (psi .gt. max(psi_xpoint(1),psi_xpoint(2))) &
                      .and. (Z_out .lt. Z_axis) &
		      .and. (R_out .gt. R_xpoint(1)) 			) location = outer
  if ( (xcase .eq. 3) .and. (psi .gt. max(psi_xpoint(1),psi_xpoint(2))) &
                      .and. (Z_out .lt. Z_axis) &
		      .and. (R_out .lt. R_xpoint(1)) 			) location = inner
  if ( (xcase .eq. 3) .and. (psi .gt. max(psi_xpoint(1),psi_xpoint(2))) &
                      .and. (Z_out .gt. Z_axis) &
		      .and. (R_out .gt. R_xpoint(2)) 			) location = outer
  if ( (xcase .eq. 3) .and. (psi .gt. max(psi_xpoint(1),psi_xpoint(2))) &
                      .and. (Z_out .gt. Z_axis) &
		      .and. (R_out .lt. R_xpoint(2)) 			) location = inner
  
  ! --- Determine loop parameters depending on location
  if (location .eq. private) then
    istart = n_flux+n_open+n_outer+n_inner+n_private
    iend   = n_flux+n_open+n_outer+n_inner+1
    indent = -1
  endif
  if (location .eq. upper_private) then
    istart = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv
    iend   = n_flux+n_open+n_outer+n_inner+n_private+1
    indent = -1
  endif
  if ((location .eq. SOL) .or. (location .eq. sandwich)) then
    istart = n_flux+1
    iend   = n_flux+n_open
    indent = +1
  endif
  if (location .eq. outer) then
    istart = n_flux+n_open+1
    iend   = n_flux+n_open+n_outer
    indent = +1
  endif
  if (location .eq. inner) then
    istart = n_flux+n_open+n_outer+1
    iend   = n_flux+n_open+n_outer+n_inner
    indent = +1
  endif
  
  ! --- Loop over flux surfaces to determine where we need to add surface
  do i_surf=istart,iend,indent
    if (psi .lt. surface_list%psi_values(i_surf)) then
      i_add = i_surf
      exit
    endif
  enddo
  
  ! --- Find surface (we create an artificial list of surfaces)
  surface_list_single%n_psi = 1
  allocate(surface_list_single%psi_values(1))
  surface_list_single%psi_values(1) = psi
  call find_flux_surfaces(xpoint,xcase,node_list,element_list,surface_list_single)
  
  ! --- Record psi-values of temporary flux-surfaces individually
  surface_list_single%flux_surfaces(1)%psi = psi
  
  ! --- Clean up temporary surfaces
  call clean_single_surface(node_list,element_list,surface_list_single%flux_surfaces(1),location,psi_xpoint,R_xpoint,Z_xpoint)
  
  ! --- Order temporary flux surfaces
  call reorder_flux_surfaces(node_list, element_list, surface_list_single, ifail)
  if (ifail .ne. 0) write(*,*)'Warning! reorder_flux_surfaces failed inside add_flux_surface:',ifail
  
  ! --- Allocate temporary surface list for copy
  surface_list_tmp%n_psi = surface_list%n_psi+1
  allocate(surface_list_tmp%psi_values   (surface_list_tmp%n_psi))
  allocate(surface_list_tmp%flux_surfaces(surface_list_tmp%n_psi))
  
  ! --- Copy all surfaces up to added surface
  do i_surf=1,i_add-1
    ! --- Copy surface
    surface_list_tmp%psi_values(i_surf)				= surface_list%psi_values(i_surf)
    surface_list_tmp%flux_surfaces(i_surf)%psi			= surface_list%flux_surfaces(i_surf)%psi
    surface_list_tmp%flux_surfaces(i_surf)%n_pieces		= surface_list%flux_surfaces(i_surf)%n_pieces	
    surface_list_tmp%flux_surfaces(i_surf)%n_parts		= surface_list%flux_surfaces(i_surf)%n_parts	
    surface_list_tmp%flux_surfaces(i_surf)%parts_index(:)	= surface_list%flux_surfaces(i_surf)%parts_index(:)
    do j=1,surface_list%flux_surfaces(i_surf)%n_pieces
      surface_list_tmp%flux_surfaces(i_surf)%elm(j) 		= surface_list%flux_surfaces(i_surf)%elm(j)
      surface_list_tmp%flux_surfaces(i_surf)%s(:,j) 		= surface_list%flux_surfaces(i_surf)%s(:,j)
      surface_list_tmp%flux_surfaces(i_surf)%t(:,j) 		= surface_list%flux_surfaces(i_surf)%t(:,j)
    enddo
    ! --- Copy intersections indexes
    do j=1,n_int_surf(i_surf)
      index_int_surf_tmp(i_surf,j) = index_int_surf(i_surf,j)
    enddo
    n_int_surf_tmp(i_surf) = n_int_surf(i_surf)
  enddo
  
  ! --- Copy all surfaces after added surface
  do i_surf=i_add,surface_list%n_psi
    ! --- Copy surface
    surface_list_tmp%psi_values(i_surf+1)			= surface_list%psi_values(i_surf)
    surface_list_tmp%flux_surfaces(i_surf+1)%psi		= surface_list%flux_surfaces(i_surf)%psi
    surface_list_tmp%flux_surfaces(i_surf+1)%n_pieces		= surface_list%flux_surfaces(i_surf)%n_pieces	
    surface_list_tmp%flux_surfaces(i_surf+1)%n_parts		= surface_list%flux_surfaces(i_surf)%n_parts	
    surface_list_tmp%flux_surfaces(i_surf+1)%parts_index(:)	= surface_list%flux_surfaces(i_surf)%parts_index(:)
    do j=1,surface_list%flux_surfaces(i_surf)%n_pieces
      surface_list_tmp%flux_surfaces(i_surf+1)%elm(j) 		= surface_list%flux_surfaces(i_surf)%elm(j)
      surface_list_tmp%flux_surfaces(i_surf+1)%s(:,j) 		= surface_list%flux_surfaces(i_surf)%s(:,j)
      surface_list_tmp%flux_surfaces(i_surf+1)%t(:,j) 		= surface_list%flux_surfaces(i_surf)%t(:,j)
    enddo
    ! --- Copy intersections indexes
    do j=1,n_int_surf(i_surf)
      index_int_surf_tmp(i_surf+1,j) = index_int_surf(i_surf,j)
    enddo
    n_int_surf_tmp(i_surf+1) = n_int_surf(i_surf)
  enddo
  
  ! --- Copy added surface
  surface_list_tmp%psi_values(i_add)				= surface_list_single%psi_values(1)
  surface_list_tmp%flux_surfaces(i_add)%psi		  	= surface_list_single%flux_surfaces(1)%psi
  surface_list_tmp%flux_surfaces(i_add)%n_pieces		= surface_list_single%flux_surfaces(1)%n_pieces   
  surface_list_tmp%flux_surfaces(i_add)%n_parts		  	= surface_list_single%flux_surfaces(1)%n_parts    
  surface_list_tmp%flux_surfaces(i_add)%parts_index(:)	  	= surface_list_single%flux_surfaces(1)%parts_index(:)
  do j=1,surface_list_single%flux_surfaces(1)%n_pieces
    surface_list_tmp%flux_surfaces(i_add)%elm(j)		= surface_list_single%flux_surfaces(1)%elm(j)
    surface_list_tmp%flux_surfaces(i_add)%s(:,j)		= surface_list_single%flux_surfaces(1)%s(:,j)
    surface_list_tmp%flux_surfaces(i_add)%t(:,j)		= surface_list_single%flux_surfaces(1)%t(:,j)
  enddo
  n_target = n_target + 1
  R_target(n_target) = R_out
  Z_target(n_target) = Z_out
  index_target(n_target,1) = i_add
  index_target(n_target,2) = 0 ! Hopefully we will not need the piece index later, I hope...
  index_target(n_target,3) = i_wall
  if (Z_target(n_target) .lt. Z_axis) then
    if (R_target(n_target) .lt. R_xpoint(1)) then
      index_target(n_target,4) = LowerLeft
    else
      index_target(n_target,4) = LowerRight
    endif
  else
    if (R_target(n_target) .lt. R_xpoint(2)) then
      index_target(n_target,4) = UpperLeft
    else
      index_target(n_target,4) = UpperRight
    endif
  endif
  n_int_surf_tmp(i_add) = 1
  index_int_surf_tmp(i_add,1) = n_target
  
  ! --- Reallocate surface data
  deallocate(surface_list%psi_values)
  deallocate(surface_list%flux_surfaces)
  allocate(surface_list%psi_values   (surface_list_tmp%n_psi))
  allocate(surface_list%flux_surfaces(surface_list_tmp%n_psi))
  
  ! --- Now copy back into surface_list to send back
  do i_surf=1,surface_list_tmp%n_psi
    ! --- Copy surface
    surface_list%psi_values(i_surf)				= surface_list_tmp%psi_values(i_surf)
    surface_list%flux_surfaces(i_surf)%psi			= surface_list_tmp%flux_surfaces(i_surf)%psi
    surface_list%flux_surfaces(i_surf)%n_pieces			= surface_list_tmp%flux_surfaces(i_surf)%n_pieces       
    surface_list%flux_surfaces(i_surf)%n_parts			= surface_list_tmp%flux_surfaces(i_surf)%n_parts        
    surface_list%flux_surfaces(i_surf)%parts_index(:)		= surface_list_tmp%flux_surfaces(i_surf)%parts_index(:)
    do j=1,surface_list_tmp%flux_surfaces(i_surf)%n_pieces
      surface_list%flux_surfaces(i_surf)%elm(j)			= surface_list_tmp%flux_surfaces(i_surf)%elm(j)
      surface_list%flux_surfaces(i_surf)%s(:,j)			= surface_list_tmp%flux_surfaces(i_surf)%s(:,j)
      surface_list%flux_surfaces(i_surf)%t(:,j)			= surface_list_tmp%flux_surfaces(i_surf)%t(:,j)
    enddo
    ! --- Copy intersections indexes
    do j=1,n_int_surf_tmp(i_surf)
      i_int = index_int_surf_tmp(i_surf,j)
      index_int_surf(i_surf,j) = i_int
      index_target(i_int,1) = i_surf
    enddo
    n_int_surf(i_surf) = n_int_surf_tmp(i_surf)
  enddo
  surface_list%n_psi = surface_list_tmp%n_psi
  
  deallocate(surface_list_tmp%psi_values)
  deallocate(surface_list_tmp%flux_surfaces)
  deallocate(surface_list_single%psi_values)
  deallocate(surface_list_single%flux_surfaces)
  
  ! --- The grid values need to be updated after we added flux surfaces
  if (location .eq. private)		n_private = n_private + 1
  if (location .eq. upper_private)	n_up_priv = n_up_priv + 1
  if (location .eq. SOL)		n_open    = n_open    + 1
  if (location .eq. sandwich)		n_open    = n_open    + 1
  if (location .eq. outer)		n_outer   = n_outer   + 1
  if (location .eq. inner)		n_inner   = n_inner   + 1
  n_grids(1) = n_flux	
  n_grids(3) = n_open	; n_grids(4) = n_outer  ; n_grids(5) = n_inner
  n_grids(6) = n_private; n_grids(7) = n_up_priv
  n_grids(8) = n_leg	; n_grids(9) = n_up_leg 
  
  return
end subroutine add_flux_surface
  
