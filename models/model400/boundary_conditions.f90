module mod_boundary_conditions
contains
  !*******************************************************************************
  !* Subroutine: boundary_condition                                              *
  !*******************************************************************************
  !*                                                                             *
  !* Add boundary condition on the matrix.                                       *
  !*                                                                             *
  !* Parameters:                                                                 *
  !*   my_id        - Identifier of the node in MPI_COMM_WORLD                   *
  !*   node_list    - List of nodes                                              *
  !*   element_list - List of all elements                                       *
  !*   local_elms   - List of local elements                                     *
  !*   n_local_elms - Number of local elements                                   *
  !*   index_min    - Minimal index of local elements (not with murge assembly)  *
  !*   index_max    - Maximal index of local elements (not with murge assembly)  *
  !*   xpoint2      -                                                            *
  !*   xcase2       -                                                            *
  !*   psi_axis     -                                                            *
  !*   psi_bnd      -                                                            *
  !*   Z_xpoint     -                                                            *
  !*   gmres        - boolean indicating if we are using GMRES method            *
  !*   solve_only   - Indicate if we want to perform only solve                  *
  !*                                                                             *
  !* Authors:                                                                    *
  !*   Xavier Lacoste - xavier.lacoste@inria.fr                                  *
  !*                                                                             *
  !*******************************************************************************
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list,           &
       local_elms, n_local_elms, index_min, index_max, rhs_loc, &
       xpoint2, xcase2,                                         &
       R_axis, Z_axis, psi_axis,                                &
       psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,                 &
       gmres, solve_only )

    use data_structure
    use global_distributed_matrix
    use phys_module, only: F0, GAMMA, freeboundary, tokamak_device,                  &
                           RMP_on, psi_RMP_cos, dpsi_RMP_cos_dR, dpsi_RMP_cos_dZ,    &
                           psi_RMP_sin, dpsi_RMP_sin_dR, dpsi_RMP_sin_dZ,            &
			   t_now, lambda, tset, RMP_start_time, tstep, RMP_har_cos, RMP_har_sin
    USE murge_module, ONLY : MURGE_ASSEMBLYBEGIN_WRAPPER => MURGE_ASSEMBLYBEGIN,     &
         use_murge, use_murge_element, murge_id, murge_global_n, MURGE_ASSEMBLY_OVW, &
         MURGE_ASSEMBLY_FOOL, murge_sym, murge_id_prod, murge_global_n_prod,         &
         MURGE_SUCCESS, murge_add_one_entry
    use mpi_mod

    implicit none
    !include 'mpif.h'

    ! --- Routine parameters
    integer,                   intent(in)    :: my_id
    type (type_node_list),     intent(in)    :: node_list
    type (type_element_list),  intent(in)    :: element_list
    type (type_bnd_node_list), intent(in)    :: bnd_node_list
    integer,                   intent(in)    :: local_elms(*)
    integer,                   intent(in)    :: n_local_elms
    integer,                   intent(in)    :: index_min
    integer,                   intent(in)    :: index_max
    real*8,                    intent(inout) :: rhs_loc(*)
    logical,                   intent(in)    :: xpoint2
    integer,                   intent(in)    :: xcase2
    real*8,                    intent(in)    :: R_axis
    real*8,                    intent(in)    :: Z_axis
    real*8,                    intent(in)    :: psi_axis
    real*8,                    intent(in)    :: psi_bnd
    real*8,                    intent(in)    :: R_xpoint(2)
    real*8,                    intent(in)    :: Z_xpoint(2)
    real*8,                    intent(in)    :: psi_xpoint(2)
    logical,                   intent(in)    :: gmres
    logical,                   intent(in)    :: solve_only

    ! Internal parameters
    real*8  :: bigR,  dBigR_ds, zbig
    real*8  :: R_s, R_t
    real*8  :: Z_s, Z_t
    real*8  :: xjac
    real*8  :: Vpar0, dVpar0_ds
    real*8  :: Ti0,   dTi0_ds
    real*8  :: Te0,   dTe0_ds
    real*8  :: ps0_s, ps0_t, ps0_x, ps0_y, grad_psi, Btot
    real*8  :: u0_s,  u0_t,  u0_x,  u0_y
    real*8  :: direction
    real*8  :: mach1, dmach1, d2mach1_dTi, d2mach1_dTe, mach_u, dmach_u
    integer :: i, in, iv, inode, k, kk, jj, iside, index_tmp
    integer :: index_large_i, index_node, index_node2, ielm
    integer :: ijA_position,ijA_position2, ilarge2, kv, kTi, kTe, ku, ilarge_vv, ilarge_vTi, ilarge_vTe, ilarge_vus
    integer :: ilarge_vsvs, ilarge_vsTi, ilarge_vsTe, ilarge_vsTis, ilarge_vsTes
    integer :: loop_nbr, loop, cnt, cnt_prod
    integer :: ierr
    logical :: is_local, only_count
    logical :: on_private, on_inner, on_inner_or_private
    logical :: apply_on_psi, apply_dirichlet(2), apply_RMP

    ! RMP parameters
    real*8, allocatable	:: psi_RMP_cos1(:),dpsi_RMP_cos_dR1(:),dpsi_RMP_cos_dZ1(:)
    real*8, allocatable	:: psi_RMP_sin1(:),dpsi_RMP_sin_dR1(:),dpsi_RMP_sin_dZ1(:)
    real*8  		:: Rnode, dRnode_ds, dRnode_dt
    real*8  		:: Znode, dZnode_ds, dZnode_dt
    real*8  		:: psi_node
    real*8  		:: establish_RMP
    real*8  		:: delta_psi_rmp
    real*8  		:: delta_psi_rmp_dR, delta_psi_rmp_dZ
    real*8  		:: delta_psi_rmp_ds, delta_psi_rmp_dt
    real*8  		:: delta_psi_rmp_tmp
    real*8  		:: psi_test, sigmo_fonc
    integer 		:: ilarge_vp, ilarge_vp2
    integer 		:: kp, j, err, itest 
    
    ! --- Some pre-defined integers
    ku = 2
    kv = 7
    kTi = 6
    kTe = 8
    
    ! --- Retrieve RMP profiles
    if (RMP_on .and. (n_tor .ge. 3)) then
      call tr_allocate(psi_RMP_cos1,    1, bnd_node_list%n_bnd_nodes,"psi_RMP_cos1",    CAT_UNKNOWN)
      call tr_allocate(dpsi_RMP_cos_dR1,1, bnd_node_list%n_bnd_nodes,"dpsi_RMP_cos_dR1",CAT_UNKNOWN)
      call tr_allocate(dpsi_RMP_cos_dZ1,1, bnd_node_list%n_bnd_nodes,"dpsi_RMP_cos_dZ1",CAT_UNKNOWN)
      call tr_allocate(psi_RMP_sin1,    1, bnd_node_list%n_bnd_nodes,"psi_RMP_sin1",    CAT_UNKNOWN)
      call tr_allocate(dpsi_RMP_sin_dR1,1, bnd_node_list%n_bnd_nodes,"dpsi_RMP_sin_dR1",CAT_UNKNOWN)
      call tr_allocate(dpsi_RMP_sin_dZ1,1, bnd_node_list%n_bnd_nodes,"dpsi_RMP_sin_dZ1",CAT_UNKNOWN)

      do i = 1, node_list%n_nodes
        if (node_list%node(i)%boundary .ne.0) then
          if (node_list%node(i)%boundary_index == 1 ) then
            itest = RMP_har_cos
            if (n_tor .eq. 1) itest = 1
            psi_test = node_list%node(i)%values(itest,1,1)
            if (my_id == 0) write (*,*) 'psi_bnd at previous time step', psi_test
          endif
        endif
      enddo
      
      if (abs(psi_test) .le. abs(psi_RMP_cos(1))) then
        !establish_RMP = (1.d-3)*tstep
        sigmo_fonc =  ( 1. + exp(-lambda*( t_now - RMP_start_time - tset )) )**(-1) &
                    - ( 1. + exp(-lambda*( 0. - tset                     )) )**(-1) 
        establish_RMP = (lambda * sigmo_fonc * (1 - sigmo_fonc) + 1.e-6) * tstep 
      else
         ! Other possibility (simpler) : if ( (t_now - RMP_start_time) .ge. 2.2*tset ) then establish_RMP =0.0
         establish_RMP = 0.0
      endif
    
      do j=1, bnd_node_list%n_bnd_nodes  
        psi_RMP_cos1(j)     = psi_RMP_cos(j)	 * establish_RMP
        dpsi_RMP_cos_dR1(j) = dpsi_RMP_cos_dR(j) * establish_RMP
        dpsi_RMP_cos_dZ1(j) = dpsi_RMP_cos_dZ(j) * establish_RMP
        psi_RMP_sin1(j)     = psi_RMP_sin(j)	 * establish_RMP
        dpsi_RMP_sin_dR1(j) = dpsi_RMP_sin_dR(j) * establish_RMP
        dpsi_RMP_sin_dZ1(j) = dpsi_RMP_sin_dZ(j) * establish_RMP
      end do

      if (my_id == 0) then
        write (*,'(A,3f15.7)') 'psi_RMP_cos1(1) and derivatives after multiplication in boundary conditions :',&
	                       psi_RMP_cos1(1), dpsi_RMP_cos_dR1(1), dpsi_RMP_cos_dZ1(1)
        write (*,'(A,f15.7)')  'establish_RMP', establish_RMP
      endif

    endif
    ! --- Retrieve RMP profiles (END)
    
    ! --- Loop twice if we are using Murge
    zbig = 1.d10
    if (use_murge .and. use_murge_element) then
      ! when we use murge assembly we first count entries then we had them.
      loop_nbr   = 2
      cnt	 = 0
      cnt_prod   = 0
      only_count = .true.
    else
      ! No need to do 2 loops when we build irn_glob, jcn_glob, A_glob.
      loop_nbr   = 1
      only_count = .false.
    end if

    ! --- Start boundary conditions
    do loop = 1, loop_nbr

#ifdef USE_MURGE
      if (loop == 2) then
        only_count = .false.
        write (*,*) my_id, ":: Murge Boundary Assembly phase :: ", cnt, " entries"
        if (.not. solve_only) CALL MURGE_ASSEMBLYBEGIN(murge_id,      murge_global_n,      cnt,      MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_FOOL, murge_sym, ierr)
        if (gmres)            CALL MURGE_ASSEMBLYBEGIN(murge_id_prod, murge_global_n_prod, cnt_prod, MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_FOOL, murge_sym, ierr)
      endif
#endif

      do i=1, n_local_elms

        ielm = local_elms(i)

        do iv=1, n_vertex_max

          inode = element_list%element(ielm)%vertex(iv)

          if (node_list%node(inode)%boundary .ne. 0) then

            ! --- Define (R,Z) coords and psi on node
	    Rnode     = node_list%node(inode)%x(1,1) 
            dRnode_ds = node_list%node(inode)%x(2,1) 
            dRnode_dt = node_list%node(inode)%x(3,1) 
            Znode     = node_list%node(inode)%x(1,2) 
            dZnode_ds = node_list%node(inode)%x(2,2) 
            dZnode_dt = node_list%node(inode)%x(3,2) 
	    psi_node  = node_list%node(inode)%values(1,1,1)
                    
            
	    do in=1, n_tor

              do k=1, n_var

		! ------------------------------------------------------------------------------------------------
		! ------------------------ First deal with Dirichlet boundary conditions -------------------------
		! ------------------------------------------------------------------------------------------------
		
		apply_dirichlet(1) = .false.
		apply_dirichlet(2) = .false.
                
		!------------------------------------ the open field lines (in case of x-point grid)
                if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then
		  ! Determine if we need to apply condition on psi
		  apply_on_psi = .false.
                  if (k .eq. 1) then
                    if 			      (in .eq. 1)             apply_on_psi = .true.
                    if ( (.not. RMP_on) .and. (in .ge. 2 )          ) apply_on_psi = .true.
                    if ( (RMP_on)       .and. (in .lt. RMP_har_cos) ) apply_on_psi = .true.
                    if ( (RMP_on)       .and. (in .gt. RMP_har_sin) ) apply_on_psi = .true.
		  endif
                      
		  ! Apply conditions ?
                  if (				&
                           apply_on_psi		& 
                      .or. (k .eq. 2)		&
                      .or. (k .eq. 3)		&
                      .or. (k .eq. 4)		&
                      ) apply_dirichlet(1) = .true.
		endif

                !------------------------------------ wall aligned with fluxsurface : wall (in case of x-point grid)
                if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

		  ! Determine if we are on the private or the inner boundary
		  on_private		= .false.
		  on_inner		= .false.
		  on_inner_or_private	= .false.
                  if ((xcase2 .ne. 3) .and. (psi_node .lt. psi_bnd)) then
		    on_inner_or_private = .true.
		    on_private          = .true.
		  endif
                  if  (xcase2 .eq. 3) then
                    if ( (Znode .lt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) .and. (psi_node .lt. psi_xpoint(1)) )                    on_private = .true.
                    if ( (Znode .gt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) .and. (psi_node .lt. psi_xpoint(2)) )                    on_private = .true.
                    if ( (Rnode .lt. (R_xpoint(1)+R_xpoint(2))/2.d0) .and. (psi_node .gt. max(psi_xpoint(2),psi_xpoint(2))) ) on_inner   = .true.
		  endif
		  if (on_private .or. on_inner) on_inner_or_private = .true.
		  
		  ! Determine if we need to apply condition on psi
		  apply_on_psi = .false.
                  if (k .eq. 1) then
                    if ( (freeboundary) .and. (in .eq. 1) ) apply_on_psi = .true.
                    if (.not. freeboundary) then
                      if                        (in .eq. 1)             apply_on_psi = .true.
                      if ( (.not. RMP_on) .and. (in .ge. 2)           ) apply_on_psi = .true.
                      if ( (RMP_on)	  .and. (in .lt. RMP_har_cos) ) apply_on_psi = .true.
                      if ( (RMP_on)	  .and. (in .gt. RMP_har_sin) ) apply_on_psi = .true.
		    endif
		  endif
		  
		  ! Apply conditions ?
                  if (  						    		&
                            (apply_on_psi)	 					&
                      .or.  (k .eq. 2)	 						&
                      .or.  (k .eq. 3)	 						&
                      .or.  (k .eq. 4)	 						&
                      .or.( (k .eq. 6) .and. (on_inner_or_private) )			& 
                      .or.  (k .eq. 7)	 						&
                      .or.( (k .eq. 8) .and. (on_inner_or_private) )			& 
                      ) apply_dirichlet(2) = .true.
		endif

                !------------------------------------ Apply Dirichlet boundary condition where we have chosen
		do jj=1,2
                
		  if ( apply_dirichlet(jj) ) then
		    do kk=1,2
		      
		      iside = kk
		      if ( (kk .eq. 2) .and. (jj .eq. 2) ) iside = 3
		      index_node = node_list%node(inode)%index(iside)
                      
		      if (use_murge .and. use_murge_element) then
                	call vertex_is_local(index_node, is_local)
                	if (is_local) then
                	  call murge_add_one_entry( index_node, k, in, index_node, k, in, &
                				    zbig, solve_only, gmres, cnt, cnt_prod, only_count)
                	end if
                      else
                	if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                	  call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                	  index_large_i = n_tor * n_var * (index_node - 1)
                	  ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                	  irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                	  jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                	  A_glob(ilarge2)   =  zbig
                	endif
                      end if
                    enddo
		  endif
                
		enddo
		    
		
		! ------------------------------------------------------------------------------------------------
		! ------------------------ Then deal with Mach-1 boundary conditions -----------------------------
		! ------------------------------------------------------------------------------------------------
		
                !------------------------------------ the open field lines (in case of x-point grid)
                if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then

                  if (k .eq. 7) then

                    ! --- Position of value and its derivative
		    index_node  = node_list%node(inode)%index(1)
                    index_node2 = node_list%node(inode)%index(2)

                    ! --- Define variables
		    Ti0       = node_list%node(inode)%values(1,1,6)
                    dTi0_ds   = node_list%node(inode)%values(1,2,6)
                    Te0       = node_list%node(inode)%values(1,1,8)
                    dTe0_ds   = node_list%node(inode)%values(1,2,8)
                    Vpar0     = node_list%node(inode)%values(1,1,k)
                    dVpar0_ds = node_list%node(inode)%values(1,2,k)
                    BigR      = node_list%node(inode)%x(1,1)
                    dBigR_ds  = node_list%node(inode)%x(2,1)

                    ps0_s     = node_list%node(inode)%values(1,2,1)
                    ps0_t     = node_list%node(inode)%values(1,3,1)

                    U0_s      = node_list%node(inode)%values(1,2,2)
                    U0_t      = node_list%node(inode)%values(1,3,2)

                    R_s       = node_list%node(inode)%x(2,1)
                    R_t       = node_list%node(inode)%x(3,1)
                    Z_s       = node_list%node(inode)%x(2,2)
                    Z_t       = node_list%node(inode)%x(3,2)

                    xjac  =  R_s*Z_t - R_t*Z_s
                    ps0_x = (	Z_t * ps0_s - Z_s * ps0_t ) / xjac
                    ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac
                    grad_psi = sqrt(ps0_x**2 + ps0_y**2)
                    Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR

                    u0_x = (   Z_t * u0_s - Z_s * u0_t ) / xjac
                    u0_y = ( - R_t * u0_s + R_s * u0_t ) / xjac

                    ! --- Define direction of velocity on target
                    if (tokamak_device(1:4) .eq. 'MAST') then
                      if ( (Rnode .gt. (R_xpoint(1)+R_xpoint(2))/2.d0) ) then
                        direction = 1.d0
                      else
                        direction = -1.d0
                      endif
                    else
                      direction = + ps0_x / abs(ps0_x)  	 ! temporary solution for lower x-point only
                    endif
                    if (xcase2 .eq. 2) direction = -direction
                    if ( (xcase2 .eq. 3) .and. (Znode .gt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) ) direction = -direction

                    ! --- Define mach1 velocity and its derivatives
                    mach1       = - zbig / Btot * direction                     * sqrt(GAMMA*(Ti0 + Te0))	     
                    dmach1      = - zbig / Btot * direction * 0.5d0  * GAMMA    / sqrt(GAMMA*(Ti0 + Te0))	     
                    d2mach1_dTi = + zbig / Btot * direction * 0.25d0 * GAMMA**2 / (GAMMA*(Ti0 + Te0))**(3/2) * dTi0_ds
                    d2mach1_dTe = + zbig / Btot * direction * 0.25d0 * GAMMA**2 / (GAMMA*(Ti0 + Te0))**(3/2) * dTe0_ds
                    mach_u      = - zbig * U0_s * BigR**2 / ps0_s
                    dmach_u     = - zbig        * BigR**2 / ps0_s
                    
		    ! --- Fill up matrix for main element (RHS)
		    if (use_murge .and. use_murge_element) then
                      call vertex_is_local(index_node, is_local)
                      if (is_local) then
                        call murge_add_one_entry(index_node,kv,in,index_node,  kv,  in, zbig,	 solve_only,gmres,cnt,cnt_prod,only_count)
			call murge_add_one_entry(index_node,kv,in,index_node,  kTi, in, dmach1,  solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node,kv,in,index_node,  kTe, in, dmach1,  solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node,kv,in,index_node2, ku,  in, dmach_u, solve_only,gmres,cnt,cnt_prod,only_count)
                        index_tmp = n_tor*n_var*(index_node-1) + (kv-1)*n_tor + in
                        if (.not. only_count) then
                          RHS_loc(index_tmp) = - Zbig*Vpar0 - mach_u - mach1
                        end if
                      end if
                    else
                      if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                        call locate_irn_jcn(index_node,index_node, index_min,index_max,ijA_position)
                        call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position2)

                        index_large_i = n_tor * n_var * (index_node - 1)

                        ilarge_vv  = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kv -1)*n_tor + in
                        ilarge_vTi = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kTi-1)*n_tor + in
                        ilarge_vTe = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kTe-1)*n_tor + in
                        ilarge_vus = ijA_position2 - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (ku -1)*n_tor + in

                        irn_glob(ilarge_vv)  =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                        jcn_glob(ilarge_vv)  =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                        A_glob(ilarge_vv)    =  zbig

                        irn_glob(ilarge_vTi) =  n_tor * n_var * (index_node-1) + (kv -1)*n_tor + in
                        jcn_glob(ilarge_vTi) =  n_tor * n_var * (index_node-1) + (kTi-1)*n_tor + in
                        A_glob(ilarge_vTi)   =  dmach1

                        irn_glob(ilarge_vTe) =  n_tor * n_var * (index_node-1) + (kv -1)*n_tor + in
                        jcn_glob(ilarge_vTe) =  n_tor * n_var * (index_node-1) + (kTe-1)*n_tor + in
                        A_glob(ilarge_vTe)   =  dmach1

                        irn_glob(ilarge_vus) =  n_tor * n_var * (index_node -1) + (kv-1)*n_tor + in
                        jcn_glob(ilarge_vus) =  n_tor * n_var * (index_node2-1) + (ku-1)*n_tor + in
                        A_glob(ilarge_vus)   =  dmach_u

                        index_tmp = n_tor*n_var*(index_node-1) + (kv-1)*n_tor + in
                        if (in .eq. 1) then
                          RHS_loc(index_tmp) = - Zbig*Vpar0 - mach_u - mach1
                        else
                          RHS_loc(index_tmp) = 0.d0
                        endif

                      endif
                    end if

                    ! --- Fill up matrix for linearised terms (LHS)
		    if (use_murge .and. use_murge_element) then
                      call vertex_is_local(index_node2, is_local)
                      if (is_local) then
                        call murge_add_one_entry(index_node2,kv,in,index_node2, kv,  in, zbig,        solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node2,kv,in,index_node2, kTi, in, dmach1,       solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node2,kv,in,index_node2, kTe, in, dmach1,       solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node2,kv,in,index_node,  kTi, in, d2mach1_dTi, solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node2,kv,in,index_node,  kTe, in, d2mach1_dTe, solve_only,gmres,cnt,cnt_prod,only_count)
                        index_tmp = n_tor*n_var*(index_node2-1) + (kv-1)*n_tor + in
                        if (.not. only_count) RHS_loc(index_tmp) = - Zbig * dVpar0_ds - dmach1 * (dTi0_ds + dTe0_ds)
                      end if
                    else
                      if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then

                        call locate_irn_jcn(index_node2,index_node, index_min,index_max,ijA_position)
                        call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2)

                        index_large_i = n_tor * n_var * (index_node2 - 1)

                        ilarge_vsvs  = ijA_position2 - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kv -1)*n_tor + in
                        ilarge_vsTis = ijA_position2 - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kTi-1)*n_tor + in
                        ilarge_vsTes = ijA_position2 - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kTe-1)*n_tor + in
                        ilarge_vsTi  = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kTi-1)*n_tor + in
                        ilarge_vsTe  = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kTe-1)*n_tor + in

                        irn_glob(ilarge_vsvs)  =  n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in
                        jcn_glob(ilarge_vsvs)  =  n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in
                        A_glob(ilarge_vsvs)    =  zbig

                        irn_glob(ilarge_vsTis) =  n_tor * n_var * (index_node2-1) + (kv -1)*n_tor + in
                        jcn_glob(ilarge_vsTis) =  n_tor * n_var * (index_node2-1) + (kTi-1)*n_tor + in
                        A_glob(ilarge_vsTis)   =  dmach1

                        irn_glob(ilarge_vsTes) =  n_tor * n_var * (index_node2-1) + (kv -1)*n_tor + in
                        jcn_glob(ilarge_vsTes) =  n_tor * n_var * (index_node2-1) + (kTe-1)*n_tor + in
                        A_glob(ilarge_vsTes)   =  dmach1

                        irn_glob(ilarge_vsTi)  =  n_tor * n_var * (index_node2-1) + (kv -1)*n_tor + in
                        jcn_glob(ilarge_vsTi)  =  n_tor * n_var * (index_node -1) + (kTi-1)*n_tor + in
                        A_glob(ilarge_vsTi)    =  d2mach1_dTi

                        irn_glob(ilarge_vsTe)  =  n_tor * n_var * (index_node2-1) + (kv -1)*n_tor + in
                        jcn_glob(ilarge_vsTe)  =  n_tor * n_var * (index_node -1) + (kTe-1)*n_tor + in
                        A_glob(ilarge_vsTe)    =  d2mach1_dTe

                        index_tmp = n_tor*n_var*(index_node2-1) + (kv-1)*n_tor + in
                        if (in .eq. 1) then
                          RHS_loc(index_tmp) = - Zbig * dVpar0_ds - dmach1 * (dTi0_ds + dTe0_ds)
                        else
                          Rhs_loc(index_tmp) = 0.d0
                        endif

                      endif
                    end if
                  end if

                end if

	      
		! ------------------------------------------------------------------------------------------------
		! ------------------------ Then deal with RMP boundary conditions --------------------------------
		! ------------------------------------------------------------------------------------------------
		
		apply_RMP = .false.
                
                if ( (RMP_on) .and. (k.eq.1) .and. (.not. freeboundary) .and. ((in.eq.RMP_har_cos) .or. (in.eq.RMP_har_sin)) ) then
		  !------------------------------------ the open field lines (in case of x-point grid)
                  if ( (node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3) ) apply_RMP = .true.
                  !------------------------------------ wall aligned with fluxsurface : wall (in case of x-point grid)
                  if ( (node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3) ) apply_RMP = .true.
		endif
		    
                !------------------------------------ Apply RMP boundary condition where we have chosen
                if (apply_RMP) then
		  kp=1    ! variable psi
                  kv=1    ! equation for psi
                  
                  index_node = node_list%node(inode)%index(1)  ! index in RHS (or matrix A not compressed)
                  					 
                  if (in.eq.RMP_har_cos) then
                    delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index)
                    delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index)
                    delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index)
                  else 
                    delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index)
                    delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index)
                    delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index)
                  endif
                  
                  if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                    call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                    
		    index_tmp = n_tor*n_var * (index_node-1) + (kv-1)*n_tor + in
                    Rhs_loc(index_tmp) = ZBIG * delta_psi_rmp
                    
                    ilarge_vp = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kp-1)*n_tor + in
                    irn_glob(ilarge_vp) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                    jcn_glob(ilarge_vp) =  n_tor * n_var * (index_node-1) + (kp-1)*n_tor + in
                    A_glob(ilarge_vp)	= ZBIG
                  endif
                  
		  iside = 2
                  delta_psi_rmp_ds  = delta_psi_rmp_dR * dRnode_ds + delta_psi_rmp_dZ * dZnode_ds
                  delta_psi_rmp_dt  = delta_psi_rmp_dR * dRnode_dt + delta_psi_rmp_dZ * dZnode_dt
                  delta_psi_rmp_tmp = delta_psi_rmp_ds
                  if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then
		    iside = 3
                    delta_psi_rmp_tmp = delta_psi_rmp_dt
		  endif
                  
		  index_node2 = node_list%node(inode)%index(iside)
		  index_tmp   = n_tor*n_var * (index_node2-1) + (kv-1)*n_tor + in
                  Rhs_loc(index_tmp) = ZBIG * delta_psi_rmp_tmp

                  if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then			    
                    call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2)
                    ilarge_vp2  = ijA_position2  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kp-1)*n_tor + in
                    
		    irn_glob(ilarge_vp2) =  n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in
                    jcn_glob(ilarge_vp2) =  n_tor * n_var * (index_node2-1) + (kp-1)*n_tor + in
                    A_glob(ilarge_vp2)   = ZBIG
                  endif
		
		endif
		
		
		
	      enddo
            enddo
          
	  endif
        
	enddo
      enddo

#ifdef USE_MURGE
      if (loop == 2) then
        if (.not. solve_only) CALL MURGE_ASSEMBLYEND(murge_id,      ierr)
        if (gmres)	      CALL MURGE_ASSEMBLYEND(murge_id_prod, ierr)
      end if
#endif
    
    end do
    
    return
  
  end subroutine boundary_conditions
end module mod_boundary_conditions
