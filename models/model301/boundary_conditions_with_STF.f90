subroutine boundary_conditions(my_id,node_list,element_list,local_elms,n_local_elms,index_min,index_max, &
                               xpoint2,psi_axis,psi_bnd,Z_xpoint)			    
!---------------------------------------------------------------
! add the boundary condition to the global matrix
!---------------------------------------------------------------
use data_structure
use global_distributed_matrix
use phys_module
use mpi_mod

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node_list)    :: nodes

integer :: my_id, local_elms(*), n_local_elms, index_min, index_max, index_min_loc, index_max_loc
real*8  :: zbig, psi_axis, psi_bnd, Z_xpoint, T0, Vpar0, bigR, dT0_ds, dVpar0_ds, dBigR_ds
real*8  :: R_s, R_t, Z_s, Z_t, ps0_s, ps0_t, ps0_x, ps0_y, direction, xjac
real*8  :: Vpar0_pol_R, Vpar0_pol_Z, Vpol_R, Vpol_Z, znormal_R, znormal_Z
real*8  :: Vpar0_perp, Vpol_perp, Btot, cs_fraction, ratio, STF
real*8  :: grad_s, grad_psi, u0_s, u0_t, u0_x, u0_y, R0, R0_s, T0_s, T0_t, T0_st, Vpar0_s, ps0_st
integer :: i_bnd, i, in, ife, iv, inode, inode1, inode2, knode, j, k, l, index_ij, index_kl
integer :: index_i, index_large_i, index_large_k, index_node, index_node1, index_node2, index_node3, index_node4
integer :: i_order, k_order, ic, ielm, ierr
integer :: ijA_position, ijA_position2, ijA_position3, ijA_position4,  nz_AA2, n_AA2, ilarge2, kp, kr, kv, kT, ku
integer :: ilarge_vv,   ilarge_vT,   ilarge_vp,  ilarge_vus, ilarge_vr,  ilarge_vvs, ilarge_vrs
integer :: ilarge_vsvs, ilarge_vsTs, ilarge_vsT, ilarge_vTs, ilarge_vTt, ilarge_vTst, ilarge_vps, ilarge_vpt
logical :: xpoint2

zbig = 1.d10

do i=1, n_local_elms

  ielm = local_elms(i)

  do iv=1, n_vertex_max

    inode = element_list%element(ielm)%vertex(iv)

    if (node_list%node(inode)%boundary .ne. 0) then

      do in=1, n_tor

        do k=1, n_var

!------------------------------------ the open field lines (in case of x-point grid)
!          if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then
          if (node_list%node(inode)%boundary .eq. 1) then

            if ((k .eq.  1) .or. (k .eq. 2)  .or. (k .eq. 3)  .or. &
                (k .eq.  4) .or. (k .eq. 95) .or. (k .eq. 96) .or. (k .eq.97) ) then


              index_node = node_list%node(inode)%index(1)

              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                index_large_i = n_tor * n_var * (index_node - 1)

                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                A_glob(ilarge2)   = zbig

              endif

              index_node = node_list%node(inode)%index(2)

              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                index_large_i = n_tor * n_var * (index_node - 1)

                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                A_glob(ilarge2)    = zbig

              endif

            endif
	    
            if ((k .eq. 6) .and. (in .eq.1)) then

              index_node  = node_list%node(inode)%index(1)             ! position of value
              index_node2 = node_list%node(inode)%index(2)             ! position of s-deriative
              index_node3 = node_list%node(inode)%index(3)             ! position of t-deriative
              index_node4 = node_list%node(inode)%index(4)             ! position of st-deriative
	      
	      STF = 3.d0 ! sheath transmission factor

              R0        = abs(node_list%node(inode)%values(1,1,5))
              R0_s      = abs(node_list%node(inode)%values(1,2,5))
              Vpar0     = node_list%node(inode)%values(1,1,7)
              Vpar0_s   = node_list%node(inode)%values(1,2,7)
	      
              BigR      = node_list%node(inode)%x(1,1)

              R_s       = node_list%node(inode)%x(2,1)
              R_t       = node_list%node(inode)%x(3,1)
              Z_s       = node_list%node(inode)%x(2,2)
              Z_t       = node_list%node(inode)%x(3,2)

              ps0_s     = node_list%node(inode)%values(1,2,1)
              ps0_t     = node_list%node(inode)%values(1,3,1)
              ps0_st    = node_list%node(inode)%values(1,4,1)
	      
              T0        = abs(node_list%node(inode)%values(1,1,6))	      
              T0_s      = node_list%node(inode)%values(1,2,6)
              T0_t      = node_list%node(inode)%values(1,3,6)
              T0_st     = node_list%node(inode)%values(1,4,6)
	      
              xjac  =  R_s*Z_t - R_t*Z_s
              ps0_x = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
              ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac

              direction = -1.d0

              Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR

              if (in .eq. 1) then

                write(*,'(A,3e14.6,A,8e14.6)') ' Boundary : ',BigR,STF * R0 * T0 * Vpar0 * Btot, &
                                                direction * ZK_par * (T0_s * ps0_t - T0_t * ps0_s) / (BigR * xjac * Btot)

              endif

              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node, index_min,index_max,ijA_position)
                call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position2)
                call locate_irn_jcn(index_node,index_node3,index_min,index_max,ijA_position3)

                index_large_i = n_tor * n_var * (index_node - 1)

                kp = 1
		kr = 5
                kv = 7
                kT = 6

                ilarge_vv  = ijA_position  - 1 + ((kT-1)*n_tor + in-1) * n_var*n_tor + (kv-1)*n_tor + in
                ilarge_vT  = ijA_position  - 1 + ((kT-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in
                ilarge_vr  = ijA_position  - 1 + ((kT-1)*n_tor + in-1) * n_var*n_tor + (kr-1)*n_tor + in

                ilarge_vTs = ijA_position2 - 1 + ((kT-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in
                ilarge_vTt = ijA_position3 - 1 + ((kT-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in

                ilarge_vps = ijA_position2 - 1 + ((kT-1)*n_tor + in-1) * n_var*n_tor + (kp-1)*n_tor + in
                ilarge_vpt = ijA_position3 - 1 + ((kT-1)*n_tor + in-1) * n_var*n_tor + (kp-1)*n_tor + in

                irn_glob(ilarge_vv) =  n_tor * n_var * (index_node-1) + (kT-1)*n_tor + in
                jcn_glob(ilarge_vv) =  n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in
                A_glob(ilarge_vv)   =  zbig * STF * R0 * T0 * Btot

                irn_glob(ilarge_vT) =  n_tor * n_var * (index_node-1) + (kT-1)*n_tor + in
                jcn_glob(ilarge_vT) =  n_tor * n_var * (index_node-1) + (kT-1)*n_tor + in
                A_glob(ilarge_vT)   =  zbig * STF * Vpar0 * R0 * Btot

                irn_glob(ilarge_vr) =  n_tor * n_var * (index_node -1) + (kT-1)*n_tor + in
                jcn_glob(ilarge_vr) =  n_tor * n_var * (index_node -1) + (kr-1)*n_tor + in
                A_glob(ilarge_vr)   =  zbig * STF * Vpar0 * T0 * Btot

                irn_glob(ilarge_vTs) =  n_tor * n_var * (index_node -1) + (kT-1)*n_tor + in
                jcn_glob(ilarge_vTs) =  n_tor * n_var * (index_node2-1) + (kT-1)*n_tor + in
                A_glob(ilarge_vTs)   = - zbig * direction * ZK_par * ps0_t / (BigR * xjac * Btot)

                irn_glob(ilarge_vTt) =  n_tor * n_var * (index_node -1) + (kT-1)*n_tor + in
                jcn_glob(ilarge_vTt) =  n_tor * n_var * (index_node3-1) + (kT-1)*n_tor + in
                A_glob(ilarge_vTt)   = + zbig * direction * ZK_par * ps0_s / (BigR * xjac * Btot)

!                irn_glob(ilarge_vps) =  n_tor * n_var * (index_node -1) + (kT-1)*n_tor + in
!                jcn_glob(ilarge_vps) =  n_tor * n_var * (index_node2-1) + (kp-1)*n_tor + in
!                A_glob(ilarge_vps)   = + zbig * direction * ZK_par * T0_t / (BigR * xjac * Btot)

!                irn_glob(ilarge_vpt) =  n_tor * n_var * (index_node -1) + (kT-1)*n_tor + in
!                jcn_glob(ilarge_vpt) =  n_tor * n_var * (index_node3-1) + (kp-1)*n_tor + in
!                A_glob(ilarge_vpt)   = - zbig * direction * ZK_par * T0_s / (BigR * xjac * Btot)

                RHS_glob(n_tor*n_var * (index_node-1) + (kT-1)*n_tor + in) = &
                         Zbig *( - STF * R0 * T0 * Vpar0 * Btot + &
			         direction * ZK_par * (T0_s * ps0_t - T0_t * ps0_s) / (BigR * xjac * Btot))

              endif

              if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then

                call locate_irn_jcn(index_node2,index_node, index_min,index_max,ijA_position)
                call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2)
                call locate_irn_jcn(index_node2,index_node3,index_min,index_max,ijA_position3)
                call locate_irn_jcn(index_node2,index_node4,index_min,index_max,ijA_position4)

                kp = 1
		kr = 5
                kv = 7
                kT = 6

                ilarge_vvs  = ijA_position2  - 1 + ((kT-1)*n_tor + in-1) * n_var*n_tor + (kv-1)*n_tor + in
                ilarge_vTs  = ijA_position2  - 1 + ((kT-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in
                ilarge_vrs  = ijA_position2  - 1 + ((kT-1)*n_tor + in-1) * n_var*n_tor + (kr-1)*n_tor + in

                ilarge_vTt  = ijA_position3 - 1 + ((kT-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in
                ilarge_vTst = ijA_position4 - 1 + ((kT-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in

                irn_glob(ilarge_vvs) =  n_tor * n_var * (index_node2-1) + (kT-1)*n_tor + in
                jcn_glob(ilarge_vvs) =  n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in
                A_glob(ilarge_vvs)   =  zbig * STF * R0 * T0 * Btot

                irn_glob(ilarge_vTs) =  n_tor * n_var * (index_node2-1) + (kT-1)*n_tor + in
                jcn_glob(ilarge_vTs) =  n_tor * n_var * (index_node2-1) + (kT-1)*n_tor + in
                A_glob(ilarge_vTs)   =  zbig * STF * Vpar0 * R0 * Btot

                irn_glob(ilarge_vrs) =  n_tor * n_var * (index_node2-1) + (kT-1)*n_tor + in
                jcn_glob(ilarge_vrs) =  n_tor * n_var * (index_node2-1) + (kr-1)*n_tor + in
                A_glob(ilarge_vrs)   =  zbig * STF * Vpar0 * T0 * Btot

                irn_glob(ilarge_vTs) =  n_tor * n_var * (index_node2-1) + (kT-1)*n_tor + in
                jcn_glob(ilarge_vTs) =  n_tor * n_var * (index_node2-1) + (kT-1)*n_tor + in
                A_glob(ilarge_vTs)   = - zbig * direction * ZK_par * ps0_st / (BigR * xjac * Btot)

                irn_glob(ilarge_vTst) =  n_tor * n_var * (index_node2-1) + (kT-1)*n_tor + in
                jcn_glob(ilarge_vTst) =  n_tor * n_var * (index_node4-1) + (kT-1)*n_tor + in
                A_glob(ilarge_vTst)   = + zbig * direction * ZK_par * ps0_s / (BigR * xjac * Btot)

!                irn_glob(ilarge_vps) =  n_tor * n_var * (index_node -1) + (kT-1)*n_tor + in
!                jcn_glob(ilarge_vps) =  n_tor * n_var * (index_node2-1) + (kp-1)*n_tor + in
!                A_glob(ilarge_vps)   = + zbig * direction * ZK_par * T0_t / (BigR * xjac * Btot)

!                irn_glob(ilarge_vpt) =  n_tor * n_var * (index_node -1) + (kT-1)*n_tor + in
!                jcn_glob(ilarge_vpt) =  n_tor * n_var * (index_node3-1) + (kp-1)*n_tor + in
!                A_glob(ilarge_vpt)   = - zbig * direction * ZK_par * T0_s / (BigR * xjac * Btot)

                RHS_glob(n_tor*n_var * (index_node2-1) + (kT-1)*n_tor + in) = &
                         Zbig *( - STF * ( R0_s * T0 * Vpar0 + R0 * T0_s * Vpar0 + R0 * T0 * Vpar0_s )* Btot + &
			         direction * ZK_par * (T0_s * ps0_st - T0_st * ps0_s) / (BigR * xjac * Btot))

              endif

            endif

            if (k .eq. 7) then

              index_node  = node_list%node(inode)%index(1)             ! position of value
              index_node2 = node_list%node(inode)%index(2)             ! position of first deriative

              T0        = abs(node_list%node(inode)%values(1,1,6))
              Vpar0     = node_list%node(inode)%values(1,1,7)
              BigR      = node_list%node(inode)%x(1,1)
              dT0_ds    = node_list%node(inode)%values(1,2,6)
              dVpar0_ds = node_list%node(inode)%values(1,2,7)
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

              grad_psi = sqrt(ps0_x**2 + ps0_y**2)

              Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR

              if (in .eq. 1) then

!                write(*,'(A,3e14.6,A,e14.6)') ' Boundary : ',Vpar0, -BigR**2 * u0_s/ps0_s, direction*sqrt(GAMMA*T0)/Btot,&
!                                              ' error : ',Vpar0 - BigR**2 * u0_s/ps0_s - direction*sqrt(GAMMA*T0)/Btot

              endif

              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node, index_min,index_max,ijA_position)
                call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position2)

                index_large_i = n_tor * n_var * (index_node - 1)

                ku = 2
                kv = 7
                kT = 6

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

                RHS_glob(n_tor*n_var * (index_node-1) + (kv-1)*n_tor + in) = &
                         Zbig * ( - Vpar0 + BigR**2 * U0_s /ps0_s + direction*sqrt(GAMMA*T0) / Btot)

              endif

              index_node  = node_list%node(inode)%index(1)
              index_node2 = node_list%node(inode)%index(2)

              if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then

                call locate_irn_jcn(index_node2,index_node,index_min,index_max,ijA_position)
                call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2)

                index_large_i = n_tor * n_var * (index_node2 - 1)

                kv = 7
                kT = 6

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

                RHS_glob(n_tor * n_var * (index_node2-1) + (kv-1)*n_tor + in) = &
                        Zbig*(-dVpar0_ds +  0.5d0 / Btot * GAMMA / sqrt(GAMMA*T0) * dT0_ds * direction)

              endif

            endif

          endif

!------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
          if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

!            if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
!                (k .eq. 4) .or. (k .eq. 95) .or. (k .eq. 96) .or. (k .eq. 7) ) then
!            if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
!                (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 6) .or. (k .eq. 7) ) then

            if ( (k .eq. 1) .or. (k .eq. 2)  .or. (k .eq. 3) .or. (k .eq. 4) .or. &
                 ( (k .eq. 5) .and. (node_list%node(inode)%values(1,1,1) .lt. psi_bnd) ) .or. &  ! private region only
		 
!                 ( (k .eq. 6) .and. (node_list%node(inode)%values(1,1,1) .lt. psi_bnd) ) .or. &  ! private region only

                 (k .eq. 6)  .or. & 
                 (k .eq. 7) ) then

              index_node = node_list%node(inode)%index(1)

              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                index_large_i = n_tor * n_var * (index_node - 1)

                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                A_glob(ilarge2)   = zbig

              endif

              index_node = node_list%node(inode)%index(3)

              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                index_large_i = n_tor * n_var * (index_node - 1)

                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                A_glob(ilarge2)    = zbig

              endif

            endif

          endif

        enddo

      enddo
    endif
  enddo
enddo

return
end
