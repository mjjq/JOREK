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
subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,    & 
     &                          n_local_elms, index_min, index_max, rhs_loc, xpoint2,   &
     &                          xcase2, psi_axis, psi_bnd, Z_xpoint, psi_xpoint, gmres, solve_only )

  use data_structure
  use global_distributed_matrix
  use phys_module, only: F0, GAMMA, freeboundary, RMP_on, psi_RMP_cos, dpsi_RMP_cos_dR, dpsi_RMP_cos_dZ, &
       psi_RMP_sin, dpsi_RMP_sin_dR, dpsi_RMP_sin_dZ, t_now, lambda, tset, RMP_start_time, tstep
  USE murge_module

  implicit none
  include 'mpif.h'

!!!! WARNING: gmres already defined in phys_module!!! Hence we use phys_module, ONLY...


  ! Subroutine parameters
  INTEGER                  :: my_id
  INTEGER                  :: local_elms(*)
  INTEGER                  :: n_local_elms
  INTEGER                  :: index_min
  INTEGER                  :: index_max
  INTEGER                  :: xcase2
  TYPE (type_node_list)    :: node_list
  TYPE (type_element_list) :: element_list
  TYPE (type_bnd_node_list):: bnd_node_list
  logical                  :: xpoint2
  REAL*8                   :: psi_axis
  REAL*8                   :: psi_bnd
  REAL*8                   :: Z_xpoint(2)
  REAL*8                   :: psi_xpoint(2)
  logical                  :: gmres
  logical                  :: solve_only
  real*8                   :: rhs_loc(*)

  ! Internal parameters
  real*8  :: zbig,  T0, Vpar0, bigR, dT0_ds, dVpar0_ds, dBigR_ds, psi_1, R_1, Z_1
  real*8  :: R_s, R_t, Z_s, Z_t, ps0_s, ps0_t, ps0_x, ps0_y, direction, xjac
  real*8  :: Btot
  real*8  :: grad_psi, u0_s, u0_t, u0_x, u0_y
  integer :: i, in, iv, inode, k
  integer :: index_large_i, index_node, index_node2, ielm
  integer :: ijA_position,ijA_position2, ilarge2, kv, kT, ku, ilarge_vv, ilarge_vT, ilarge_vus
  integer :: ilarge_vsvs, ilarge_vsTs, ilarge_vsT
  integer :: loop_nbr, loop, cnt
  integer :: first_tor, last_tor, murge_ntor, ierr
  logical :: is_local


!=============== RMP ==============
  real*8, allocatable :: psi_RMP_cos1(:),dpsi_RMP_cos_dR1(:),dpsi_RMP_cos_dZ1(:)
  real*8, allocatable :: psi_RMP_sin1(:),dpsi_RMP_sin_dR1(:),dpsi_RMP_sin_dZ1(:)
  real*8  :: Rnode, dRnode_ds, Znode, dZnode_ds, dRnode_dt, dZnode_dt, establish_RMP
  real*8  :: delta_psi_rmp, delta_psi_rmp_dR, delta_psi_rmp_dZ, delta_psi_rmp_ds, delta_psi_rmp_dt, psi_test, sigmo_fonc
  integer :: ilarge_vp, ilarge_vp2
  integer :: kp, j, err, itest

  if (RMP_on .and. (n_tor .ge. 3)) then
     allocate(psi_RMP_cos1(bnd_node_list%n_bnd_nodes), dpsi_RMP_cos_dR1(bnd_node_list%n_bnd_nodes), &
          dpsi_RMP_cos_dZ1(bnd_node_list%n_bnd_nodes), psi_RMP_sin1(bnd_node_list%n_bnd_nodes), &
          dpsi_RMP_sin_dR1(bnd_node_list%n_bnd_nodes), dpsi_RMP_sin_dZ1(bnd_node_list%n_bnd_nodes))
!=============== RMP ==============

     do i = 1, node_list%n_nodes
        if (node_list%node(i)%boundary .ne.0) then
           if (node_list%node(i)%boundary_index == 1 ) then

              if (n_tor .eq. 1) then
                 itest = 1
              else 
                 itest = 2
              endif
              
              psi_test = node_list%node(i)%values(itest,1,1)
              if (my_id == 0) then
                 write (*,*) 'psi_bnd at previous time step', psi_test
              endif
           endif
        endif
     enddo
     
     if (abs(psi_test) .le. abs(psi_RMP_cos(1))) then
     !   establish_RMP = (1.d-3)*tstep
      sigmo_fonc = ( 1. + exp(-lambda*( t_now - RMP_start_time - tset )))**(-1) &
          - ( 1. + exp(-lambda*( 0. - tset )))**(-1) 
        establish_RMP = lambda*sigmo_fonc*(1-sigmo_fonc)*tstep
     else
        establish_RMP = 0.0
     endif
     ! Other possibility (simpler) : if ( (t_now - RMP_start_time) .ge. 2.2*tset ) then establish_RMP =0.0
  
     
     do j=1, bnd_node_list%n_bnd_nodes  
        psi_RMP_cos1(j)     = psi_RMP_cos(j)     * establish_RMP
        dpsi_RMP_cos_dR1(j) = dpsi_RMP_cos_dR(j) * establish_RMP
        dpsi_RMP_cos_dZ1(j) = dpsi_RMP_cos_dZ(j) * establish_RMP
        psi_RMP_sin1(j)     = psi_RMP_sin(j)     * establish_RMP
        dpsi_RMP_sin_dR1(j) = dpsi_RMP_sin_dR(j) * establish_RMP
        dpsi_RMP_sin_dZ1(j) = dpsi_RMP_sin_dZ(j) * establish_RMP
     end do

     if (my_id == 0) then

        write (*,*) 'psi_RMP_cos1(1) and derivatives after multiplication in boundary conditions'
        write (*,*) psi_RMP_cos1(1), dpsi_RMP_cos_dR1(1), dpsi_RMP_cos_dZ1(1)
        write (*,*) 'establish_RMP', establish_RMP

     endif

  endif

  zbig = 1.d10
  if (use_murge .and. use_murge_element) then
     ! when we use murge assembly we first count entries then we had them.
     loop_nbr = 2
     cnt      = 0
  else
     ! No need to do 2 loops when we build irn_glob, jcn_glob, A_glob.
     loop_nbr = 1
  end if

  if (gmres .and. use_murge .and. use_murge_element) then 
     if (murge_harmonic == 1) then
        first_tor = 1
        last_tor = 1
        murge_ntor = 1
     else
        first_tor = 2*(murge_harmonic-1)
        last_tor = 2*(murge_harmonic-1)+1
        murge_ntor = 2
     end if
  else
     first_tor = 1
     last_tor = n_tor
     murge_ntor = n_tor
  end if

  do loop = 1, loop_nbr
     if (loop == 2)  then
        write (*,*) my_id, ":: Murge Boundary Assembly phase :: ", cnt, " entries"
        if (.not. solve_only) then
           CALL MURGE_ASSEMBLYBEGIN(murge_id, cnt, MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW, &
                MURGE_ASSEMBLY_FOOL, murge_sym, ierr)
        end if
        if (gmres) then
           CALL MURGE_ASSEMBLYBEGIN(murge_id_prod, cnt, MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW, &
                MURGE_ASSEMBLY_FOOL, murge_sym, ierr)
        end if
     end if

     do i=1, n_local_elms

        ielm = local_elms(i)

        do iv=1, n_vertex_max

           inode = element_list%element(ielm)%vertex(iv)

           if (node_list%node(inode)%boundary .ne. 0) then

              do in=first_tor, last_tor

                 do k=1, n_var

                    !------------------------------------ the open field lines (in case of x-point grid)
                    if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then




!====================================== begining RMPs at boundary ======================================================
!================================== type 1 - boundary: only depends on 's'
! ======================================================================================================================
                       
                       if (RMP_on ) then

                          if ((k.eq.1) .and. ((in.eq.2) .or. (in.eq.3)) .and. (.not. freeboundary)) then
                             ! in .eq. 2 corresponds to cos(n_perturbation)
                             ! in .eq. 3 corresponds to sin(n_perturbation)
                                       
                             kp=1    ! variable psi
                             kv=1    ! equation for psi
                         
                             index_node = node_list%node(inode)%index(1)  ! index in RHS (or matrix A not compressed)
                                                
                             Rnode     = node_list%node(inode)%x(1,1) 
                             dRnode_ds = node_list%node(inode)%x(2,1) 
                             Znode     = node_list%node(inode)%x(1,2) 
                             dZnode_ds = node_list%node(inode)%x(2,2) 
                          
                             if (in.eq.2) then
                                delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index)
                                delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index)
                                delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index)

                             else 
                                delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index)
                                delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index)
                                delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index)

                             endif
                             
                             delta_psi_rmp_ds = delta_psi_rmp_dR * dRnode_ds + delta_psi_rmp_dZ * dZnode_ds

                             if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                                                  
                                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                                
                                !-------- index dans A_glob
                                ilarge_vp  = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kp-1)*n_tor + in
                                
                                Rhs_loc(n_tor*n_var * (index_node-1) + (kv-1)*n_tor + in) = ZBIG * delta_psi_rmp
                                
                                irn_glob(ilarge_vp) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                                jcn_glob(ilarge_vp) =  n_tor * n_var * (index_node-1) + (kp-1)*n_tor + in
                                A_glob(ilarge_vp)   = ZBIG
                             endif
                             
                             index_node2 = node_list%node(inode)%index(2)

                             if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then                         
                                call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2)
                                
                                ilarge_vp2  = ijA_position2  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kp-1)*n_tor + in
                                
                                Rhs_loc(n_tor*n_var * (index_node2-1) + (kv-1)*n_tor + in) = ZBIG * delta_psi_rmp_ds

                                irn_glob(ilarge_vp2) =  n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in
                                jcn_glob(ilarge_vp2) =  n_tor * n_var * (index_node2-1) + (kp-1)*n_tor + in
                                A_glob(ilarge_vp2)   = ZBIG
                             endif
                          endif

                          
                       endif !(end RMP)
!======================================= end RMPs ==================================




                       
                       if (                  &
 !                           ((k .eq. 1) .and. (.not. RMP_on) .and. ((in.eq.2) .or. (in.eq.3)))      &
                            ((k .eq. 1) .and. (.not. RMP_on) .and. ( in .ge. 2 ))      &
!                            .or. ((k .eq. 1) .and. ((in.ne.2) .and. (in.ne.3)))  &
                            .or. ((k .eq. 1) .and. (in .eq. 1))  & 
!                            .or. ((k .eq. 3) .and. (in .eq. 1))  &                           
                            .or. (k .eq. 2)  &
                            .or. (k .eq. 3)  &
                            .or. (k .eq. 4)  &
                                !.or. (k .eq. 5) &
                                !.or. (k .eq. 6) &
                                !.or. (k .eq. 7) &
                            ) then


                          index_node = node_list%node(inode)%index(1)
                          if (use_murge .and. use_murge_element) then
                             call vertex_is_local(index_node, is_local)
                             if (is_local) then
                                if (loop /= loop_nbr) then
                                   cnt = cnt + 1
                                else
                                   call murge_add_one_entry( & 
                                        & index_node, k, in, &
                                        & index_node, k, in, &
                                        & zbig, murge_ntor, solve_only, gmres)

                                end if
                             end if
                          else
                             if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                                index_large_i = n_tor * n_var * (index_node - 1)

                                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                A_glob(ilarge2)   = zbig

                             endif
                          end if

                          index_node = node_list%node(inode)%index(2)
                          if (use_murge .and. use_murge_element) then
                             call vertex_is_local(index_node, is_local)
                             if (is_local) then
                                if (loop /= loop_nbr) then
                                   cnt = cnt + 1
                                else
                                   call murge_add_one_entry( & 
                                        & index_node, k, in, &
                                        & index_node, k, in, &
                                        & zbig, murge_ntor, solve_only, gmres)

                                end if
                             end if
                          else
                             if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                                index_large_i = n_tor * n_var * (index_node - 1)

                                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                A_glob(ilarge2)    = zbig

                             endif
                          end if
                       endif


                       if (k .eq. 7) then

                          index_node  = node_list%node(inode)%index(1)             ! position of value
                          index_node2 = node_list%node(inode)%index(2)             ! position of first deriative

                          T0        = abs(node_list%node(inode)%values(1,1,6))
                          Vpar0     = node_list%node(inode)%values(1,1,k)
                          BigR      = node_list%node(inode)%x(1,1)
                          dT0_ds    = node_list%node(inode)%values(1,2,6)
                          dVpar0_ds = node_list%node(inode)%values(1,2,k)
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
                          ps0_x = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
                          ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac

                          u0_x = (   Z_t * u0_s - Z_s * u0_t ) / xjac
                          u0_y = ( - R_t * u0_s + R_s * u0_t ) / xjac

                          direction = + ps0_x / abs(ps0_x)             ! temporary solution for lower x-point only
			  if (xcase2 .eq. 2) direction = -direction
			  if ( (xcase2 .eq. 3) .and. (node_list%node(inode)%x(1,2) .gt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) ) direction = -direction

                          grad_psi = sqrt(ps0_x**2 + ps0_y**2)

                          Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR

                          if (in .eq. 1) then

                             !                write(*,'(A,3e14.6,A,e14.6)') ' Boundary : ',Vpar0, -BigR**2 * u0_s/ps0_s, direction*sqrt(GAMMA*T0)/Btot,&
                             !                                              ' error : ',Vpar0 - BigR**2 * u0_s/ps0_s - direction*sqrt(GAMMA*T0)/Btot

                          endif

                          ku = 2
                          kv = 7
                          kT = 6

                          if (use_murge .and. use_murge_element) then
                             call vertex_is_local(index_node, is_local)
                             if (is_local) then
                                if (loop /= loop_nbr) then
                                   cnt = cnt + 3
                                else
                                   call murge_add_one_entry( & 
                                        & index_node, kv, in, &
                                        & index_node, kv, in, &
                                        & zbig, murge_ntor, solve_only, gmres)
                                   call murge_add_one_entry( & 
                                        & index_node, kv, in, &
                                        & index_node, kT, in, &
                                        & - zbig / Btot * 0.5d0 * GAMMA / sqrt(GAMMA*T0) * direction, murge_ntor, solve_only, gmres)
                                   call murge_add_one_entry( & 
                                        & index_node,  kv, in, &
                                        & index_node2, ku, in, &
                                        & - zbig * BigR**2 / ps0_s, murge_ntor, solve_only, gmres)

                                   Rhs_loc(n_tor*n_var * (index_node-1) + (kv-1)*n_tor + in) = &
                                        Zbig * ( - Vpar0 + BigR**2 * U0_s /ps0_s + direction*sqrt(GAMMA*T0) / Btot)
                                end if
                             end if
                          else
                             if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                call locate_irn_jcn(index_node,index_node, index_min,index_max,ijA_position)
                                call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position2)

                                index_large_i = n_tor * n_var * (index_node - 1)


                                ilarge_vv  = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kv-1)*n_tor + in
                                ilarge_vT  = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in
                                ilarge_vus = ijA_position2 - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (ku-1)*n_tor + in

                                irn_glob(ilarge_vv) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                                jcn_glob(ilarge_vv) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                                A_glob(ilarge_vv)   =  zbig

                                irn_glob(ilarge_vT) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                                jcn_glob(ilarge_vT) =  n_tor * n_var * (index_node-1) + (kT-1)*n_tor + in
                                A_glob(ilarge_vT)   = - zbig / Btot * 0.5d0 * GAMMA / sqrt(GAMMA*T0) * direction

                                irn_glob(ilarge_vus) =  n_tor * n_var * (index_node -1) + (kv-1)*n_tor + in
                                jcn_glob(ilarge_vus) =  n_tor * n_var * (index_node2-1) + (ku-1)*n_tor + in
                                A_glob(ilarge_vus)   = - zbig * BigR**2 / ps0_s


                                if (in .eq. 1) then
                                   RHS_loc(n_tor*n_var * (index_node-1) + (kv-1)*n_tor + in) = &
                                        Zbig * ( - Vpar0 + BigR**2 * U0_s /ps0_s + direction*sqrt(GAMMA*T0) / Btot)
                                else
                                   RHS_loc(n_tor*n_var * (index_node-1) + (kv-1)*n_tor + in) = 0.d0
                                endif

                             endif
                          end if

                          index_node  = node_list%node(inode)%index(1)
                          index_node2 = node_list%node(inode)%index(2)
                          kv = 7
                          kT = 6

                          if (use_murge .and. use_murge_element) then
                             call vertex_is_local(index_node2, is_local)
                             if (is_local) then
                                if (loop /= loop_nbr) then
                                   cnt = cnt + 1
                                else
                                   call murge_add_one_entry( & 
                                        & index_node2, kv, in, &
                                        & index_node2, kv, in, &
                                        & zbig, murge_ntor, solve_only, gmres)
                                   call murge_add_one_entry( & 
                                        & index_node2, kv, in, &
                                        & index_node2, kT, in, &
                                        & - zbig / Btot * 0.5d0 * GAMMA / sqrt(GAMMA*T0) * direction, murge_ntor, solve_only, gmres)
                                   call murge_add_one_entry( & 
                                        & index_node2, kv, in, &
                                        & index_node,  kT, in, &
                                        & + zbig / Btot * 0.25d0 * GAMMA**2 / (GAMMA*T0)**(3/2) * dT0_ds * direction, murge_ntor, solve_only, gmres)

                                   Rhs_loc(n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in) = &
                                        Zbig*(-dVpar0_ds +  0.5d0 / Btot * GAMMA / sqrt(GAMMA*T0) * dT0_ds * direction)
                                end if
                             end if
                          else
                             if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then

                                call locate_irn_jcn(index_node2,index_node,index_min,index_max,ijA_position)
                                call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2)

                                index_large_i = n_tor * n_var * (index_node2 - 1)


                                ilarge_vsvs = ijA_position2 - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kv-1)*n_tor + in
                                ilarge_vsTs = ijA_position2 - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in
                                ilarge_vsT  = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in

                                irn_glob(ilarge_vsvs) =  n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in
                                jcn_glob(ilarge_vsvs) =  n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in
                                A_glob(ilarge_vsvs)   = zbig

                                irn_glob(ilarge_vsTs) =  n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in
                                jcn_glob(ilarge_vsTs) =  n_tor * n_var * (index_node2-1) + (kT-1)*n_tor + in
                                A_glob(ilarge_vsTs)   = - zbig / Btot * 0.5d0 * GAMMA / sqrt(GAMMA*T0) * direction

                                irn_glob(ilarge_vsT) =  n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in
                                jcn_glob(ilarge_vsT) =  n_tor * n_var * (index_node -1) + (kT-1)*n_tor + in
                                A_glob(ilarge_vsT)   = + zbig / Btot * 0.25d0 * GAMMA**2 / (GAMMA*T0)**(3/2) * dT0_ds * direction

                                if (in .eq. 1) then
                                   Rhs_loc(n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in) = &
                                        Zbig*(-dVpar0_ds +  0.5d0 / Btot * GAMMA / sqrt(GAMMA*T0) * dT0_ds * direction)
                                else
                                   Rhs_loc(n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in) = 0.d0
                                endif 

                             endif
                          end if
                       end if

                    end if



                    !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                    if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

!====================================== begining RMPs at boundary ======================================================
!================================== type 2 - boundary: only depends on 't'
! ======================================================================================================================
                       
                       if (RMP_on ) then

                          if ((k.eq.1) .and. ((in.eq.2) .or. (in.eq.3)) .and. (.not. freeboundary)) then
                             ! in .eq. 2 corresponds to cos(n_perturbation)
                             ! in .eq. 3 corresponds to sin(n_perturbation)
               

                                       
                             kp=1    ! variable psi
                             kv=1    ! equation for psi
                          
                             index_node = node_list%node(inode)%index(1)  ! index in RHS (or matrix A not compressed)
                                                
                             Rnode     = node_list%node(inode)%x(1,1) 
                             dRnode_dt = node_list%node(inode)%x(3,1) 
                             Znode     = node_list%node(inode)%x(1,2) 
                             dZnode_dt = node_list%node(inode)%x(3,2) 
                          
                             if (in.eq.2) then
                                delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index)
                                delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index)
                                delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index)

                                if (node_list%node(inode)%boundary_index == 1 ) then
                                   write (*,*) 'type2_bnd: my_id, psi_RMP_cos1, Rnode, Znode'
                                   write (*,*) my_id, delta_psi_rmp, Rnode, Znode
                                   write (*,*) 'delta_psi_rmp_dR, delta_psi_rmp_dZ'      
                                   write (*,*) delta_psi_rmp_dR, delta_psi_rmp_dZ
                                endif
                             else 
                                delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index)
                                delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index)
                                delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index)

                             endif

                             delta_psi_rmp_dt = delta_psi_rmp_dR * dRnode_dt + delta_psi_rmp_dZ * dZnode_dt
                             if (in.eq.2) then

                                if (node_list%node(inode)%boundary_index == 1 ) then
                                   write (*,*) 'delta_psi_rmp_dt', delta_psi_rmp_dt
                                   write (*,*) 'delta_psi_rmp_ds', delta_psi_rmp_ds
                                endif
                             endif

                             if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                          
                                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                         
                                !-------- index dans A_glob
                                ilarge_vp  = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kp-1)*n_tor + in
                                
                                Rhs_loc(n_tor*n_var * (index_node-1) + (kv-1)*n_tor + in) = ZBIG * delta_psi_rmp
                          
                                irn_glob(ilarge_vp) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                                jcn_glob(ilarge_vp) =  n_tor * n_var * (index_node-1) + (kp-1)*n_tor + in
                                A_glob(ilarge_vp)   = ZBIG

                             endif
                                                      
                             index_node2 = node_list%node(inode)%index(3)

                             if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then                         
                                call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2)
                         
                                ilarge_vp2  = ijA_position2  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kp-1)*n_tor + in
                                
                                Rhs_loc(n_tor*n_var * (index_node2-1) + (kv-1)*n_tor + in) = ZBIG * delta_psi_rmp_dt
                                                
                                irn_glob(ilarge_vp2) =  n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in
                                jcn_glob(ilarge_vp2) =  n_tor * n_var * (index_node2-1) + (kp-1)*n_tor + in
                                A_glob(ilarge_vp2)   = ZBIG

                             endif
                          endif
                       endif
                       
!======================================= end RMPs ==================================


                       if (                                                      &
                            ((freeboundary) .and. (k .eq. 1) .and. (in .eq. 1))  &               ! exclude condition on psi (freeboundary) except n=0
!                            .or. (( .not. freeboundary) .and. (k .eq. 1) .and. (.not. RMP_on) .and. ((in.eq.2) .or. (in.eq.3)))   &
                            .or. (( .not. freeboundary) .and. (k .eq. 1) .and. (.not. RMP_on) .and. ( in .ge. 2 ))   &
                            .or. (( .not. freeboundary) .and. (k .eq. 1) .and. (in.eq.1))        &
                            .or. (k .eq. 2)    &
                            .or. (k .eq. 3)    &
                            .or. (k .eq. 4)    &
                            .or. (k .eq. 5)    &
                                !.or.( (k .eq. 5) .and.   & 
			          !(    ((xcase2 .ne. 3) .and. (node_list%node(inode)%values(1,1,1) .lt. psi_bnd)) &
			          !.or. ((xcase2 .eq. 3) .and. (node_list%node(inode)%x(1,2) .lt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) &
				  !                      .and. (node_list%node(inode)%values(1,1,1) .lt. psi_xpoint(1)) )  &
			          !.or. ((xcase2 .eq. 3) .and. (node_list%node(inode)%x(1,2) .gt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) &
				  !                      .and. (node_list%node(inode)%values(1,1,1) .lt. psi_xpoint(2)) ) ))  &  ! private region only
                            .or. (k .eq. 6)    &
                                !.or.( (k .eq. 6) .and.   & 
			          !(    ((xcase2 .ne. 3) .and. (node_list%node(inode)%values(1,1,1) .lt. psi_bnd)) &
			          !.or. ((xcase2 .eq. 3) .and. (node_list%node(inode)%x(1,2) .lt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) &
				  !                      .and. (node_list%node(inode)%values(1,1,1) .lt. psi_xpoint(1)) )  &
			          !.or. ((xcase2 .eq. 3) .and. (node_list%node(inode)%x(1,2) .gt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) &
				  !                      .and. (node_list%node(inode)%values(1,1,1) .lt. psi_xpoint(2)) ) ))  &  ! private region only
                            .or. (k .eq. 7)    &
                            ) then

                          index_node = node_list%node(inode)%index(1)
                          if (use_murge .and. use_murge_element) then
                             call vertex_is_local(index_node, is_local)
                             if (is_local) then
                                if (loop /= loop_nbr) then
                                   cnt = cnt + 1
                                else
                                   call murge_add_one_entry( & 
                                        & index_node, k, in, &
                                        & index_node, k, in, &
                                        & zbig, murge_ntor, solve_only, gmres)

                                end if
                             end if
                          else
                             if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                                index_large_i = n_tor * n_var * (index_node - 1)

                                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                A_glob(ilarge2)   = zbig

                             endif
                          end if
                          index_node = node_list%node(inode)%index(3)

                          if (use_murge .and. use_murge_element) then
                             call vertex_is_local(index_node, is_local)
                             if (is_local) then
                                if (loop /= loop_nbr) then
                                   cnt = cnt + 1
                                else
                                   call murge_add_one_entry( & 
                                        & index_node, k, in, &
                                        & index_node, k, in, &
                                        & zbig, murge_ntor, solve_only, gmres)

                                end if
                             end if
                          else
                             if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                                index_large_i = n_tor * n_var * (index_node - 1)

                                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                A_glob(ilarge2)    = zbig

                             endif
                          end if
                       endif

                    endif

                 enddo

              enddo
           endif
        enddo
     enddo
     if (loop == 2) then
        if (.not. solve_only) then
           CALL MURGE_ASSEMBLYEND(murge_id, ierr)
        end if
        if (gmres) then
           CALL MURGE_ASSEMBLYEND(murge_id_prod, ierr)
        end if
     end if
  end do

  if (RMP_on) then
     if (allocated(psi_RMP_cos1)) deallocate(psi_RMP_cos1)
     if (allocated(dpsi_RMP_cos_dR1)) deallocate(dpsi_RMP_cos_dR1)
     if (allocated(dpsi_RMP_cos_dZ1)) deallocate(dpsi_RMP_cos_dZ1)
     if (allocated(psi_RMP_sin1)) deallocate(psi_RMP_sin1)
     if (allocated(dpsi_RMP_sin_dR1)) deallocate(dpsi_RMP_sin_dR1)
     if (allocated(dpsi_RMP_sin_dZ1)) deallocate(dpsi_RMP_sin_dZ1)
  endif

  return
end subroutine boundary_conditions
