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

    ! --- Internal parameters
    real*8  :: zbig
    real*8  :: Vpar0, Vpar0_s
    real*8  :: Ti0,   Ti0_s
    real*8  :: Te0,   Te0_s
    real*8  :: R, R_s, R_t
    real*8  :: Z, Z_s, Z_t
    real*8  :: xjac, Btot
    real*8  :: ps0, ps0_s, ps0_t, ps0_x, ps0_y, grad_psi
    real*8  :: u0_s, u0_t, u0_x, u0_y
    real*8  :: direction
    real*8  :: mach1, dmach1, d2mach1_dTi, d2mach1_dTe, mach_u, dmach_u, dmach_rho
    integer :: i, in, iv, inode, k_var
    integer :: index_node, index_node2, ielm, index_tmp
    integer :: ijA_position,ijA_position2, ilarge2, k_Vpar, k_Ti, k_Te, k_u, ilarge_vv, ilarge_vTi, ilarge_vTe, ilarge_vus
    integer :: ilarge_vsvs, ilarge_vsTi, ilarge_vsTe, ilarge_vsTis, ilarge_vsTes
    integer :: loop_nbr, loop, cnt, cnt_prod
    integer :: ierr
    logical :: is_local, only_count
    logical :: apply_dirichlet, apply_on_psi, on_private, on_inner, on_inner_or_private

    ! --- RMP parameters
    real*8, allocatable	:: psi_RMP_cos1(:),dpsi_RMP_cos_dR1(:),dpsi_RMP_cos_dZ1(:)
    real*8, allocatable	:: psi_RMP_sin1(:),dpsi_RMP_sin_dR1(:),dpsi_RMP_sin_dZ1(:)
    real*8  		:: establish_RMP
    real*8  		:: delta_psi_rmp, delta_psi_rmp_dR, delta_psi_rmp_dZ, delta_psi_rmp_ds, delta_psi_rmp_dt, psi_test, sigmo_fonc
    integer 		:: ilarge_vp, ilarge_vp2
    integer 		:: k_psi, j, err, itest 
    
    ! --- Constants (note: zbig is used to overwrite equations from element_matrix with the ones here)
    k_psi  = 1
    k_u    = 2
    k_Ti   = 6
    k_Vpar = 7
    k_Te   = 8
    zbig   = 1.d10
    
    ! -------------------------
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
        establish_RMP = 0.5d-3 
      else
        establish_RMP = 0.0
      endif
    
      do j=1, bnd_node_list%n_bnd_nodes  
        psi_RMP_cos1(j)     =  psi_RMP_cos(j)	 * establish_RMP
        dpsi_RMP_cos_dR1(j) = dpsi_RMP_cos_dR(j) * establish_RMP
        dpsi_RMP_cos_dZ1(j) = dpsi_RMP_cos_dZ(j) * establish_RMP
        psi_RMP_sin1(j)     =  psi_RMP_sin(j)	 * establish_RMP
        dpsi_RMP_sin_dR1(j) = dpsi_RMP_sin_dR(j) * establish_RMP
        dpsi_RMP_sin_dZ1(j) = dpsi_RMP_sin_dZ(j) * establish_RMP
      end do
    endif
    ! --- Retrieve RMP profiles (END)
    ! -------------------------------
    
    ! --- when we use murge assembly we first count entries then we had them.
    if (use_murge .and. use_murge_element) then
      loop_nbr   = 2
      cnt	 = 0
      cnt_prod   = 0
      only_count = .true.
    else
    ! --- No need to do 2 loops when we build irn_glob, jcn_glob, A_glob.
      loop_nbr   = 1
      only_count = .false.
    end if

    do loop = 1, loop_nbr

#ifdef USE_MURGE
      if (loop == 2) then
        only_count = .false.
        write (*,*) my_id, ":: Murge Boundary Assembly phase :: ", cnt, " entries"
        if (.not. solve_only) then
          CALL MURGE_ASSEMBLYBEGIN( murge_id, murge_global_n, cnt,	       &
               &		    MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW,    &
               &		    MURGE_ASSEMBLY_FOOL, murge_sym, ierr)
        endif
        if (gmres) then
          CALL MURGE_ASSEMBLYBEGIN( murge_id_prod, murge_global_n_prod, cnt_prod, &
               &		    MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW,	  &
               &		    MURGE_ASSEMBLY_FOOL, murge_sym, ierr)
        endif
      endif
#endif

      ! --- Loop on each element
      do i=1, n_local_elms
        ielm = local_elms(i)

        ! --- Take each node of element
        do iv=1, n_vertex_max
          inode = element_list%element(ielm)%vertex(iv)

          ! --- We only care about boundary elements
          if (node_list%node(inode)%boundary .ne. 0) then

            ! --- Define (R,Z) coords and psi on node
            R	 = node_list%node(inode)%x(1,1)
            R_s  = node_list%node(inode)%x(2,1)
            R_t  = node_list%node(inode)%x(3,1)
            Z	 = node_list%node(inode)%x(1,2)
            Z_s  = node_list%node(inode)%x(2,2)
            Z_t  = node_list%node(inode)%x(3,2)
            xjac = R_s*Z_t - R_t*Z_s
	    ps0  = node_list%node(inode)%values(1,1,1)
            
            do in=1, n_tor

              do k_var=1, n_var

                ! --------------------------------------------------------------------------------------------------------------
                ! ------------------------------------ the targets (in case of x-point grid) -----------------------------------
                ! --------------------------------------------------------------------------------------------------------------
                if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then

                  ! ---------------------------------------------
                  ! --- Apply RMP on target (only depends on 's')
                  if (RMP_on ) then

                    if ((k_var.eq.1) .and. ((in.eq.RMP_har_cos) .or. (in.eq.RMP_har_sin)) .and. (.not. freeboundary)) then
                        		 
                      if (in.eq.RMP_har_cos) then
                         delta_psi_rmp    =  psi_RMP_cos1   (node_list%node(inode)%boundary_index)
                         delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index)
                         delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index)
                      else 
                         delta_psi_rmp    =  psi_RMP_sin1   (node_list%node(inode)%boundary_index)
                         delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index)
                         delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index)
                      endif
                      
                      delta_psi_rmp_ds = delta_psi_rmp_dR*R_s + delta_psi_rmp_dZ*Z_s

                      ! --- Condition on nodes
                      index_node = node_list%node(inode)%index(1)
                      if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                         call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                         ilarge_vp  = ijA_position  - 1 + ((k_psi-1)*n_tor + in-1) * n_var*n_tor + (k_psi-1)*n_tor + in
                         
                         irn_glob(ilarge_vp) = n_tor * n_var * (index_node-1) + (k_psi-1)*n_tor + in
                         jcn_glob(ilarge_vp) = n_tor * n_var * (index_node-1) + (k_psi-1)*n_tor + in
                         A_glob(ilarge_vp)   = ZBIG
                         
		         index_tmp = n_tor*n_var * (index_node-1) + (k_psi-1)*n_tor + in
                         Rhs_loc(index_tmp) = ZBIG * delta_psi_rmp
                      endif
                      
                      ! --- Condition between nodes (d/ds)
                      index_node2 = node_list%node(inode)%index(2)
                      if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then 			
                         call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2)
                         
                         ilarge_vp2  = ijA_position2  - 1 + ((k_psi-1)*n_tor + in-1) * n_var*n_tor + (k_psi-1)*n_tor + in

                         irn_glob(ilarge_vp2) = n_tor * n_var * (index_node2-1) + (k_psi-1)*n_tor + in
                         jcn_glob(ilarge_vp2) = n_tor * n_var * (index_node2-1) + (k_psi-1)*n_tor + in
                         A_glob(ilarge_vp2)   = ZBIG
                         
		         index_tmp = n_tor*n_var * (index_node2-1) + (k_psi-1)*n_tor + in
                         Rhs_loc(index_tmp) = ZBIG * delta_psi_rmp_ds
                      endif
                    endif
                    
                  endif
                  
                  ! -----------------------------------------------
		  ! --- Dirichlet BCs (or Neumann if commented out)
		  apply_dirichlet = .false.
		  ! --- Determine if we need to apply condition on psi (we don't want to overwrite RMPs)
		  apply_on_psi = .false.
                  if (k_var .eq. 1) then
                    if  		      (in .eq. 1)	      apply_on_psi = .true.
                    if ( (.not. RMP_on) .and. (in .ge. 2 )	    ) apply_on_psi = .true.
                    if ( (RMP_on)	.and. (in .lt. RMP_har_cos) ) apply_on_psi = .true.
                    if ( (RMP_on)	.and. (in .gt. RMP_har_sin) ) apply_on_psi = .true.
		  endif
                  
		  ! --- Apply conditions to which variables?
                  if (  			&
                           apply_on_psi 	& 
                      .or. (k_var .eq. 2)	&
                      .or. (k_var .eq. 3)	&
                      .or. (k_var .eq. 4)  	&
                      ) apply_dirichlet = .true.

		  ! --- Apply Dirichlet if required
                  if (apply_dirichlet) then

                    ! --- MURGE
		    if (use_murge .and. use_murge_element) then
                      ! --- Condition on nodes
		      index_node = node_list%node(inode)%index(1)
                      call vertex_is_local(index_node, is_local)
                      if (is_local) call murge_add_one_entry(index_node,k_var,in,index_node,k_var,in, zbig, solve_only,gmres,cnt,cnt_prod,only_count)
                      ! --- Condition between nodes (d/ds)
                      index_node = node_list%node(inode)%index(2)
                      call vertex_is_local(index_node, is_local)
                      if (is_local) call murge_add_one_entry(index_node,k_var,in,index_node,k_var,in, zbig, solve_only,gmres,cnt,cnt_prod,only_count)
                    ! --- non-MURGE
                    else
                      ! --- Condition on nodes
		      index_node = node_list%node(inode)%index(1)
                      if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                        call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                        ilarge2 = ijA_position - 1 + ((k_var-1)*n_tor + in-1) * n_var*n_tor + (k_var-1)*n_tor + in

                        irn_glob(ilarge2) = n_tor * n_var * (index_node-1) + (k_var-1)*n_tor + in
                        jcn_glob(ilarge2) = n_tor * n_var * (index_node-1) + (k_var-1)*n_tor + in
                        A_glob(ilarge2)   = zbig
                      endif
                      ! --- Condition between nodes (d/ds)
                      index_node = node_list%node(inode)%index(2)
                      if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                        call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                        ilarge2 = ijA_position - 1 + ((k_var-1)*n_tor + in-1) * n_var*n_tor + (k_var-1)*n_tor + in

                        irn_glob(ilarge2) = n_tor * n_var * (index_node-1) + (k_var-1)*n_tor + in
                        jcn_glob(ilarge2) = n_tor * n_var * (index_node-1) + (k_var-1)*n_tor + in
                        A_glob(ilarge2)   = zbig
                      endif
                    end if
                  endif

                  ! --------------
		  ! --- Mach-1 BCs
                  if (k_var .eq. 7) then

                    ! --- Define node indices
		    index_node  = node_list%node(inode)%index(1)	     ! position of value
                    index_node2 = node_list%node(inode)%index(2)	     ! position of first deriative

                    ! --- Define variables
                    ps0_s     = node_list%node(inode)%values(1,2,1)
                    ps0_t     = node_list%node(inode)%values(1,3,1)
                    ps0_x     = (  Z_t*ps0_s - Z_s*ps0_t) / xjac
                    ps0_y     = (- R_t*ps0_s + R_s*ps0_t) / xjac
                    grad_psi  = sqrt(        ps0_x**2 + ps0_y**2)
                    Btot      = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / R
                    
		    U0_s      = node_list%node(inode)%values(1,2,2)
                    U0_t      = node_list%node(inode)%values(1,3,2)
                    u0_x      = (  Z_t*u0_s - Z_s*u0_t) / xjac
                    u0_y      = (- R_t*u0_s + R_s*u0_t) / xjac
		    
		    Ti0       = node_list%node(inode)%values(1,1,6)
                    Ti0_s     = node_list%node(inode)%values(1,2,6)
                    
		    Te0       = node_list%node(inode)%values(1,1,8)
                    Te0_s     = node_list%node(inode)%values(1,2,8)
                    
		    Vpar0     = node_list%node(inode)%values(1,1,k_Vpar)
                    Vpar0_s   = node_list%node(inode)%values(1,2,k_Vpar)

                    ! --- Define direction of Vpar on target. Careful, using ps0_x/abs(ps0_x) can be treacherous.
                    if (tokamak_device(1:4) .eq. 'MAST') then
                      if ( (node_list%node(inode)%x(1,1) .gt. (R_xpoint(1)+R_xpoint(2))/2.d0) ) then
                        direction = 1.d0
                      else
                        direction = -1.d0
                      endif
                    else
                      direction = + ps0_x / abs(ps0_x)
                    endif
                    if (xcase2 .eq. 2) direction = -direction
                    if ( (xcase2 .eq. 3) .and. (node_list%node(inode)%x(1,2) .gt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) ) direction = -direction

                    ! --- Define equations before MURGE and non-MURGE fork
		    mach1       = - zbig / Btot * direction                     * sqrt(GAMMA*(Ti0 + Te0))	     
                    dmach1      = - zbig / Btot * direction * 0.5d0  * GAMMA    / sqrt(GAMMA*(Ti0 + Te0))	     
                    d2mach1_dTi = + zbig / Btot * direction * 0.25d0 * GAMMA**2 / (GAMMA*(Ti0 + Te0))**(3/2) * Ti0_s
                    d2mach1_dTe = + zbig / Btot * direction * 0.25d0 * GAMMA**2 / (GAMMA*(Ti0 + Te0))**(3/2) * Te0_s
                    mach_u      = - zbig * U0_s * R**2 / ps0_s
                    dmach_u     = - zbig        * R**2 / ps0_s
                    !mach_u      = - zbig * (U0_s + tauIC*   Pi0_s     /rho0   ) * R**2 / ps0_s
                    !dmach_u     = - zbig                                        * R**2 / ps0_s
                    !dmach_rho   = - zbig * (     - tauIC*   Pi0_s     /rho0**2) * R**2 / ps0_s &
		    !              - zbig * (     + tauIC*(Ti0_s+Ti0)/rho0   ) * R**2 / ps0_s 
                    
                    ! --- MURGE
                    if (use_murge .and. use_murge_element) then
                      ! --- Condition on nodes
		      index_tmp = n_tor*n_var*(index_node-1) + (k_Vpar-1)*n_tor + in
                      call vertex_is_local(index_node, is_local)
                      if (is_local) then
                        call murge_add_one_entry(index_node,k_Vpar,in,index_node, k_Vpar, in, zbig,    solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node,k_Vpar,in,index_node, k_Ti,   in, dmach1,  solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node,k_Vpar,in,index_node, k_Te,   in, dmach1,  solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node,k_Vpar,in,index_node2,k_u,    in, dmach_u, solve_only,gmres,cnt,cnt_prod,only_count)
                        if (.not. only_count) then
			  RHS_loc(index_tmp) = - Zbig*Vpar0 - mach_u - mach1
                        endif
                      endif
                      
		      ! --- Condition between nodes (d/ds)
                      index_tmp = n_tor*n_var*(index_node2-1) + (k_Vpar-1)*n_tor + in
                      call vertex_is_local(index_node2, is_local)
                      if (is_local) then
                        call murge_add_one_entry(index_node2,k_Vpar,in,index_node2, k_Vpar, in, zbig,        solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node2,k_Vpar,in,index_node2, k_Ti,   in, dmach1,      solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node2,k_Vpar,in,index_node2, k_Te,   in, dmach1,      solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node2,k_Vpar,in,index_node,  k_Ti,   in, d2mach1_dTi, solve_only,gmres,cnt,cnt_prod,only_count)
                        call murge_add_one_entry(index_node2,k_Vpar,in,index_node,  k_Te,   in, d2mach1_dTe, solve_only,gmres,cnt,cnt_prod,only_count)
                        if (.not. only_count) then 
                          RHS_loc(index_tmp) = - Zbig*Vpar0_s - dmach1 * (Ti0_s + Te0_s)
                        endif
                      endif
                    else
                      ! --- Condition on nodes
		      index_tmp = n_tor*n_var*(index_node-1) + (k_Vpar-1)*n_tor + in
                      if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                        call locate_irn_jcn(index_node,index_node, index_min,index_max,ijA_position)
                        call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position2)

                        ilarge_vv  = ijA_position  - 1 + ((k_Vpar-1)*n_tor + in-1) * n_var*n_tor + (k_Vpar-1)*n_tor + in
                        ilarge_vTi = ijA_position  - 1 + ((k_Vpar-1)*n_tor + in-1) * n_var*n_tor + (k_Ti  -1)*n_tor + in
                        ilarge_vTe = ijA_position  - 1 + ((k_Vpar-1)*n_tor + in-1) * n_var*n_tor + (k_Te  -1)*n_tor + in
                        ilarge_vus = ijA_position2 - 1 + ((k_Vpar-1)*n_tor + in-1) * n_var*n_tor + (k_u   -1)*n_tor + in

                        irn_glob(ilarge_vv)  = n_tor * n_var * (index_node-1) + (k_Vpar-1)*n_tor + in
                        jcn_glob(ilarge_vv)  = n_tor * n_var * (index_node-1) + (k_Vpar-1)*n_tor + in
                        A_glob(ilarge_vv)    = zbig

                        irn_glob(ilarge_vTi) = n_tor * n_var * (index_node-1) + (k_Vpar-1)*n_tor + in
                        jcn_glob(ilarge_vTi) = n_tor * n_var * (index_node-1) + (k_Ti  -1)*n_tor + in
                        A_glob(ilarge_vTi)   = dmach1

                        irn_glob(ilarge_vTe) = n_tor * n_var * (index_node-1) + (k_Vpar-1)*n_tor + in
                        jcn_glob(ilarge_vTe) = n_tor * n_var * (index_node-1) + (k_Te  -1)*n_tor + in
                        A_glob(ilarge_vTe)   = dmach1

                        irn_glob(ilarge_vus) = n_tor * n_var * (index_node -1) + (k_Vpar-1)*n_tor + in
                        jcn_glob(ilarge_vus) = n_tor * n_var * (index_node2-1) + (k_u   -1)*n_tor + in
                        A_glob(ilarge_vus)   = dmach_u

                        if (in .eq. 1) then
                          RHS_loc(index_tmp) = - Zbig*Vpar0 - mach_u - mach1
                        else
                          RHS_loc(index_tmp) = 0.d0
                        endif
                      endif
		      
                      ! --- Condition between nodes (d/ds)
                      index_tmp = n_tor*n_var*(index_node2-1) + (k_Vpar-1)*n_tor + in
                      if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then
                        call locate_irn_jcn(index_node2,index_node, index_min,index_max,ijA_position)
                        call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2)

                        ilarge_vsvs  = ijA_position2 - 1 + ((k_Vpar-1)*n_tor + in-1) * n_var*n_tor + (k_Vpar-1)*n_tor + in
                        ilarge_vsTis = ijA_position2 - 1 + ((k_Vpar-1)*n_tor + in-1) * n_var*n_tor + (k_Ti  -1)*n_tor + in
                        ilarge_vsTes = ijA_position2 - 1 + ((k_Vpar-1)*n_tor + in-1) * n_var*n_tor + (k_Te  -1)*n_tor + in
                        ilarge_vsTi  = ijA_position  - 1 + ((k_Vpar-1)*n_tor + in-1) * n_var*n_tor + (k_Ti  -1)*n_tor + in
                        ilarge_vsTe  = ijA_position  - 1 + ((k_Vpar-1)*n_tor + in-1) * n_var*n_tor + (k_Te  -1)*n_tor + in

                        irn_glob(ilarge_vsvs)  = n_tor * n_var * (index_node2-1) + (k_Vpar-1)*n_tor + in
                        jcn_glob(ilarge_vsvs)  = n_tor * n_var * (index_node2-1) + (k_Vpar-1)*n_tor + in
                        A_glob(ilarge_vsvs)    = zbig

                        irn_glob(ilarge_vsTis) = n_tor * n_var * (index_node2-1) + (k_Vpar-1)*n_tor + in
                        jcn_glob(ilarge_vsTis) = n_tor * n_var * (index_node2-1) + (k_Ti  -1)*n_tor + in
                        A_glob(ilarge_vsTis)   = dmach1

                        irn_glob(ilarge_vsTes) = n_tor * n_var * (index_node2-1) + (k_Vpar-1)*n_tor + in
                        jcn_glob(ilarge_vsTes) = n_tor * n_var * (index_node2-1) + (k_Te  -1)*n_tor + in
                        A_glob(ilarge_vsTes)   = dmach1

                        irn_glob(ilarge_vsTi)  = n_tor * n_var * (index_node2-1) + (k_Vpar-1)*n_tor + in
                        jcn_glob(ilarge_vsTi)  = n_tor * n_var * (index_node -1) + (k_Ti  -1)*n_tor + in
                        A_glob(ilarge_vsTi)    = d2mach1_dTi

                        irn_glob(ilarge_vsTe)  = n_tor * n_var * (index_node2-1) + (k_Vpar-1)*n_tor + in
                        jcn_glob(ilarge_vsTe)  = n_tor * n_var * (index_node -1) + (k_Te  -1)*n_tor + in
                        A_glob(ilarge_vsTe)    = d2mach1_dTe

                        if (in .eq. 1) then
                          RHS_loc(index_tmp) = - Zbig*Vpar0_s - dmach1 * (Ti0_s + Te0_s)
                        else
                          Rhs_loc(index_tmp) = 0.d0
                        endif
                      endif
                    endif
                    
                  endif

                endif

                
		! ----------------------------------------------------------------------------------------------------
                ! ------------------------------------ the flux-surface boundaries -----------------------------------
                ! ----------------------------------------------------------------------------------------------------
                if (   (node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3) ) then

                  ! ---------------------------------------------
                  ! --- Apply RMP on target (only depends on 't')
                  if (RMP_on ) then

                    if ((k_var.eq.1) .and. ((in.eq.RMP_har_cos) .or. (in.eq.RMP_har_sin)) .and. (.not. freeboundary)) then
                        		 
                      if (in.eq.RMP_har_cos) then
                        delta_psi_rmp    =  psi_RMP_cos1   (node_list%node(inode)%boundary_index)
                        delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index)
                        delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index)
                      else 
                        delta_psi_rmp    =  psi_RMP_sin1   (node_list%node(inode)%boundary_index)
                        delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index)
                        delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index)
                      endif

                      delta_psi_rmp_dt = delta_psi_rmp_dR * R_t + delta_psi_rmp_dZ * Z_t
                      
                      ! --- Condition on nodes
                      index_node = node_list%node(inode)%index(1)
                      if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                        call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                        ilarge_vp  = ijA_position  - 1 + ((k_psi-1)*n_tor + in-1) * n_var*n_tor + (k_psi-1)*n_tor + in
                    
                        irn_glob(ilarge_vp) = n_tor * n_var * (index_node-1) + (k_psi-1)*n_tor + in
                        jcn_glob(ilarge_vp) = n_tor * n_var * (index_node-1) + (k_psi-1)*n_tor + in
                        A_glob(ilarge_vp)   = ZBIG
                        
                        index_tmp = n_tor*n_var * (index_node-1) + (k_psi-1)*n_tor + in
                        Rhs_loc(index_tmp) = ZBIG * delta_psi_rmp
                      endif
                        		       
                      ! --- Condition between nodes (d/ds)
                      index_node2 = node_list%node(inode)%index(3)
                      if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then 			
                        call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2)
                    
                        ilarge_vp2  = ijA_position2 - 1 + ((k_psi-1)*n_tor + in-1) * n_var*n_tor + (k_psi-1)*n_tor + in
                        		
                        irn_glob(ilarge_vp2) = n_tor * n_var * (index_node2-1) + (k_psi-1)*n_tor + in
                        jcn_glob(ilarge_vp2) = n_tor * n_var * (index_node2-1) + (k_psi-1)*n_tor + in
                        A_glob(ilarge_vp2)   = ZBIG
                        
                        index_tmp = n_tor*n_var * (index_node2-1) + (k_psi-1)*n_tor + in
                        Rhs_loc(index_tmp) = ZBIG * delta_psi_rmp_dt
                      endif
                    endif
                  endif
                  
                  ! -----------------------------------------------
		  ! --- Dirichlet BCs (or Neumann if commented out)
		  apply_dirichlet = .false.
		  
		  ! --- Determine if we are on the private or the inner boundary
		  on_private		= .false.
		  on_inner		= .false.
		  on_inner_or_private	= .false.
                  if ((xcase2 .ne. 3) .and. (ps0 .lt. psi_bnd)) then
		    on_inner_or_private = .true.
		    on_private          = .true.
		  endif
                  if  (xcase2 .eq. 3) then
                    if ( (Z .lt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) .and. (ps0 .lt. psi_xpoint(1)) )                    on_private = .true.
                    if ( (Z .gt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) .and. (ps0 .lt. psi_xpoint(2)) )                    on_private = .true.
                    if ( (R .lt. (R_xpoint(1)+R_xpoint(2))/2.d0) .and. (ps0 .gt. max(psi_xpoint(2),psi_xpoint(2))) ) on_inner   = .true.
		  endif
		  if (on_private .or. on_inner) on_inner_or_private = .true.
		  
		  ! --- Determine if we need to apply condition on psi (we don't want to overwrite RMPs)
		  apply_on_psi = .false.
                  if (k_var .eq. 1) then
                    if ( (freeboundary) .and. (in .eq. 1) ) apply_on_psi = .true.
                    if (.not. freeboundary) then
                      if                        (in .eq. 1)             apply_on_psi = .true.
                      if ( (.not. RMP_on) .and. (in .ge. 2)           ) apply_on_psi = .true.
                      if ( (RMP_on)	  .and. (in .lt. RMP_har_cos) ) apply_on_psi = .true.
                      if ( (RMP_on)	  .and. (in .gt. RMP_har_sin) ) apply_on_psi = .true.
		    endif
		  endif
		  
		  ! Apply conditions to which variables and where?
                  if (  						    		&
                            (apply_on_psi)	 					&
                      .or.  (k_var .eq. 2)	 					&
                      .or.  (k_var .eq. 3)	 					&
                      .or.  (k_var .eq. 4)	 					&
                      .or.( (k_var .eq. 6) .and. (on_private) )				& 
                      .or.  (k_var .eq. 7)	 					&
                      .or.( (k_var .eq. 8) .and. (on_private) )				& 
                      ) apply_dirichlet = .true.

		  if (apply_dirichlet) then

                    ! --- MURGE
                    if (use_murge .and. use_murge_element) then
                      ! --- Condition on nodes
                      index_node = node_list%node(inode)%index(1)
                      call vertex_is_local(index_node, is_local)
                      if (is_local) call murge_add_one_entry(index_node,k_var,in,index_node,k_var,in, zbig, solve_only,gmres,cnt,cnt_prod,only_count)
                      ! --- Condition between nodes (d/ds)
		      index_node = node_list%node(inode)%index(3)
                      call vertex_is_local(index_node, is_local)
                      if (is_local) call murge_add_one_entry(index_node,k_var,in,index_node,k_var,in, zbig, solve_only,gmres,cnt,cnt_prod,only_count)
                    ! --- non-MURGE
                    else
                      ! --- Condition on nodes
                      index_node = node_list%node(inode)%index(1)
                      if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                        call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                        ilarge2 = ijA_position - 1 + ((k_var-1)*n_tor + in-1) * n_var*n_tor + (k_var-1)*n_tor + in

                        irn_glob(ilarge2) = n_tor * n_var * (index_node-1) + (k_var-1)*n_tor + in
                        jcn_glob(ilarge2) = n_tor * n_var * (index_node-1) + (k_var-1)*n_tor + in
                        A_glob(ilarge2)   = zbig
                      endif
                      ! --- Condition between nodes (d/ds)
		      index_node = node_list%node(inode)%index(3)
                      if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                        call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                        ilarge2 = ijA_position - 1 + ((k_var-1)*n_tor + in-1) * n_var*n_tor + (k_var-1)*n_tor + in

                        irn_glob(ilarge2) = n_tor * n_var * (index_node-1) + (k_var-1)*n_tor + in
                        jcn_glob(ilarge2) = n_tor * n_var * (index_node-1) + (k_var-1)*n_tor + in
                        A_glob(ilarge2)   = zbig
                      endif
                    endif
                    
                  endif

                endif

              enddo

            enddo
          endif
        enddo
      enddo
#ifdef USE_MURGE
       if (loop == 2) then
          if (.not. solve_only) then
             CALL MURGE_ASSEMBLYEND(murge_id, ierr)
          end if
          if (gmres) then
             CALL MURGE_ASSEMBLYEND(murge_id_prod, ierr)
          end if
       end if
#endif
    end do
    return
  end subroutine boundary_conditions
end module mod_boundary_conditions
