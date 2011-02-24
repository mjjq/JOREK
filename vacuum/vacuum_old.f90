subroutine vacuum_old(my_id,node_list,element_list,bnd_elm_list,bnd_node_list, &
                  index_min,index_max,xpoint2,psi_axis,psi_bnd,Z_xpoint,rhs_loc)
!---------------------------------------------------------------------
! calculates the matrix contribution of the boundary integral of the
! induction and current equations using the STARWALL vacuum response
!---------------------------------------------------------------------
  
  use parameters
  use data_structure
  use gauss
  use basis_at_gaussian
  use phys_module
  use global_distributed_matrix
  use vacuum_response
  
  
  implicit none
  
  
  integer,                      intent(in)    :: my_id          ! Number of current MPI thread.
  type (type_node_list),        intent(in)    :: node_list      ! Node list of current grid.
  type (type_element_list),     intent(in)    :: element_list   ! Element list of current grid.
  type (type_bnd_element_list), intent(in)    :: bnd_elm_list  ! List of 1D boundary elements.
  type (type_bnd_node_list),    intent(in)    :: bnd_node_list  ! List of boundary nodes.
  integer,                      intent(in)    :: index_min      ! Minimum index handled by current MPI thread
  integer,                      intent(in)    :: index_max      ! Maximum index handled by current MPI thread
  logical,                      intent(in)    :: xpoint2        ! X-point equilibrium?
  real*8,                       intent(in)    :: psi_axis       ! Psi value at magnetic axis.
  real*8,                       intent(in)    :: psi_bnd        ! Psi value at plasma boundary.
  real*8,                       intent(in)    :: Z_xpoint       ! Z-position of X-point.
  real*8,                       intent(inout) :: rhs_loc(ndof_glob) ! Part of the RHS handled by current MPI thread
  
  
  real*8     :: x_g(n_gauss), x_s(n_gauss)
  real*8     :: y_g(n_gauss), y_s(n_gauss)
  
  real*8     :: eq_g(n_plane,n_var,n_gauss), eq_s(n_plane,n_var,n_gauss), eq_p(n_plane,n_var,n_gauss)
  real*8     :: delta_g(n_plane,n_var,n_gauss), delta_s(n_plane,n_var,n_gauss)
  
  integer    :: ibnd,i, j, ms, mp, kp, kj, kbnd, k, l, jdir, kdir, ldir, korder, kv, lv,  ilarge_vv, inode, inode_bnd, inode2_bnd, ktor
  integer    :: index_node, index_node2, index_node3, index_node_bnd, index_node2_bnd, index_node3_bnd
  integer    :: ilarge_pp, ilarge_jp, ijA_position
  integer    :: in, im
  integer    :: jpsi, iwall, jnode, jtor, jbas, jnode_glob, jnode2, jelem
  real*8     :: ws, BigR, PI
  real*8     :: A_glob_31
  real*8     :: theta, zeta
  real*8     :: v, ijsize, klsize
  real*8     :: ps0, ps0_s, r0, T0
  real*8     :: bpar
  integer    :: icurr
  real*8     :: value, delta
  integer    :: knode, knode3
  integer    :: ilarge_j
  real*8     :: br_coils, bz_coils, psi_coils
  real*8     :: rhs_contrib
  real*8     :: psi_val
  integer    :: knode_bnd
  real*8, allocatable :: psi_bnd_vec(:), dpsi_bnd_vec(:)
  
  write(*,*) '************************************'
  write(*,*) '*     VACUUM boundary integral     *'
  write(*,*) '************************************'
  
  PI    = 2.d0 * asin(1.d0)
  
  theta = 0.5d0; zeta = 0.d0          ! Crank-Nicholson parameter
  !theta = 1.0d0  ; zeta = 0.0d0       ! Euler scheme
  !theta = 1.5d0   ; zeta = 0.5d0      ! BDF2 (Gears) scheme
  
  call global_matrix_structure_vacuum(node_list, bnd_node_list, index_min, index_max)
  
  ! --- Determine vector of psi boundary values.
  allocate(  psi_bnd_vec(response_index(bnd_node_list%n_bnd_nodes, n_tor, 2)), &
            dpsi_bnd_vec(response_index(bnd_node_list%n_bnd_nodes, n_tor, 2)) )
  jpsi=1
  do jnode = 1, bnd_node_list%n_bnd_nodes       ! loop over nodes
    jnode_glob = bnd_node_list%bnd_node(jnode)%index_jorek
    do jtor = 2, 3                              ! loop over toroidal harmonics ###
      do jbas = 1, 2                            ! loop over basis functions
        jdir = bnd_node_list%bnd_node(jnode)%direction(jbas)
        psi_bnd_vec ( response_index(jnode,jtor,jbas) ) = node_list%node(jnode_glob)%values(jtor, jdir, jpsi)
        dpsi_bnd_vec( response_index(jnode,jtor,jbas) ) = node_list%node(jnode_glob)%deltas(jtor, jdir, jpsi)
      end do
    end do
  end do
  
  
  !###
!  write(*,*) '###', size(wall_curr), size(diagonal_yy)
!  wall_curr=1.d0
!  do i = 1, 1000
  !###
  
  ! --- Time-integrate currents in the resistive wall.
  !     (Y: wall current potentials, Psi: poloidal flux, [YY] and [YE]: response matrices)
  !
  !     Y^n+1 = Y^n - delta t/(sigma_wall * d_wall) * [YY] * Y - [YE] * delta Psi
  !
  if ( resistive_wall ) then
    
    ! ([YY] * Y part:)
    if ( maxval(diagonal_yy) * tstep * wall_resistivity / wall_thickness > 1. ) then
      write(*,*) 'Numerically unstable! tstep must be reduced by at least a factor of ',maxval(diagonal_yy) * tstep * wall_resistivity / wall_thickness
      stop
    end if
    
    do i = 1, 1!###
      if ( wall_curr_treatment == 'explicit' ) then
        wall_curr(:) = wall_curr(:) - tstep/1. * ( diagonal_yy(:) * wall_curr(:) ) * wall_resistivity / wall_thickness
      else
        wall_curr(:) = wall_curr(:) - diagonal_r(:) * wall_curr(:)
      end if
    end do
    
    if ( wall_curr_treatment == 'explicit' ) then
      do i = 1, n_wall_curr
        wall_curr(i) = - sum( matrix_ey(i, :) * psi_bnd_vec(:) ) !################################# WRONG ##########################
      end do
    else
      do i = 1, n_wall_curr
!        wall_curr(i) = - sum( matrix_s (i, :) * psi_bnd_vec(:) )
        wall_curr(i) = wall_curr(i) - sum( matrix_s (i, :) * dpsi_bnd_vec(:) )
      end do
    end if
    
    ! ([YE] * delta Psi part:)
!    jpsi=1 ! select variable Psi
!    wall_curr = 0.d0 !###
!    do jnode = 1, bnd_node_list%n_bnd_nodes       ! loop over nodes
!      jnode_glob = bnd_node_list%bnd_node(jnode)%index_jorek
!      do jtor = 2, 3                              ! loop over toroidal harmonics ###
!        do jbas = 1, 2                            ! loop over basis functions
!          
!          jdir = bnd_node_list%bnd_node(jnode)%direction(jbas)
!          
!          if ( wall_curr_treatment == 'explicit' ) then
!            wall_curr(:) = wall_curr(:) - matrix_ye( :, response_index(jnode,jtor,jbas) ) * node_list%node(jnode_glob)%deltas(jtor, jdir, jpsi)
!          else
!            !wall_curr(:) = wall_curr(:) - matrix_s ( :, response_index(jnode,jtor,jbas) ) * node_list%node(jnode_glob)%deltas(jtor, jdir, jpsi)
!            wall_curr(:) = wall_curr(:) - matrix_s ( :, response_index(jnode,jtor,jbas) ) * node_list%node(jnode_glob)%values(jtor, jdir, jpsi)
!          end if
!          
!        end do
!      end do
!    end do
  end if
  
  !###
  if ( resistive_wall ) write(96,*) sum(abs(wall_curr))
!  end do
!  stop
  !###
  
  
  !write(35+my_id,'(4ES20.12)') sum(abs(rhs_loc)), sum(abs(A_glob))
  
  ibnd_loop: do ibnd = 1, bnd_elm_list%n_bnd_elements ! (loop over boundary elements)
  
    ! --- Determine coordinate and variable values at Gaussian points on boundary elements.
    x_g  = 0.d0; x_s  = 0.d0;
    y_g  = 0.d0; y_s  = 0.d0; 
    eq_g = 0.d0; eq_s = 0.d0; eq_p = 0.d0;
    delta_g = 0.d0; delta_s = 0.d0;
  
    do i = 1 , 2  ! (loop over nodes of boundary element)
      do j = 1, 2 ! (loop over basis functions)
        
        inode  = bnd_elm_list%bnd_element(ibnd)%vertex(i)
        jdir   = bnd_elm_list%bnd_element(ibnd)%direction(i,j)
        ijsize = bnd_elm_list%bnd_element(ibnd)%size(i,j)

        x_g(:)  = x_g(:)  + node_list%node(inode)%x(jdir,1) * ijsize * H1(i,j,:)     ! x:=R
        x_s(:)  = x_s(:)  + node_list%node(inode)%x(jdir,1) * ijsize * H1_s(i,j,:)   ! dx/ds
        y_g(:)  = y_g(:)  + node_list%node(inode)%x(jdir,2) * ijsize * H1(i,j,:)     ! y:=Z
        y_s(:)  = y_s(:)  + node_list%node(inode)%x(jdir,2) * ijsize * H1_s(i,j,:)   ! dZ/ds
        
        do mp = 1, n_plane   ! (loop over toroidal planes)
          do k = 1, n_var    ! (loop over all variables)
            do in = 1, n_tor ! (loop over toroidal harmonics)

              value = node_list%node(inode)%values(in,jdir,k)
              delta = node_list%node(inode)%deltas(in,jdir,k)

              eq_g(mp,k,:)    = eq_g(mp,k,:)    + value * ijsize * H1(i,j,:)   * HZ(in,mp)    ! variable value
              eq_s(mp,k,:)    = eq_s(mp,k,:)    + value * ijsize * H1_s(i,j,:) * HZ(in,mp)    ! dvar/ds
              eq_p(mp,k,:)    = eq_p(mp,k,:)    + value * ijsize * H1(i,j,:)   * HZ_p(in,mp)  ! dvar/dphi
              delta_g(mp,k,:) = delta_g(mp,k,:) + delta * ijsize * H1(i,j,:)   * HZ(in,mp)    ! var change in last timestep
              delta_s(mp,k,:) = delta_s(mp,k,:) + delta * ijsize * H1_s(i,j,:) * HZ(in,mp)    ! ddelta/ds

            end do
          end do
        end do
        
      end do
    end do
  
    ! --- Contribution of boundary integrals to psi- and j-equations.
    kp = 1 ! select psi variable (boundary integral in psi equation)
    kj = 3 ! select current variable (boundary integral in j equation)
    ms_loop: do ms = 1, n_gauss ! (loop over Gaussian points)
      ws = wgauss(ms)
  
      mp_loop: do mp = 1, n_plane ! (loop over toroidal planes)
        ! Note: Loops over ibnd, ms, and mp implement the boundary integral as three sums.
  
        ps0   = eq_g(mp,1,ms) ! psi (pol. flux)
        ps0_s = eq_s(mp,1,ms) ! dpsi/ds
        r0    = eq_g(mp,5,ms) ! rho (density)
        T0    = eq_g(mp,6,ms) ! T (temperature)
        BigR  = x_g(ms)       ! R (major radius)
        
        i_loop: do i = 1, 2 ! (loop over nodes of boundary element)
          j_loop: do j = 1, 2 ! (loop over basis functions)
  
            jdir           = bnd_elm_list%bnd_element(ibnd)%direction(i,j)
            inode          = bnd_elm_list%bnd_element(ibnd)%vertex(i)
            inode_bnd      = bnd_elm_list%bnd_element(ibnd)%bnd_vertex(i)
            index_node     = node_list%node(inode)%index(jdir)
            index_node_bnd = 2 * (bnd_node_list%bnd_node(inode_bnd)%index_starwall-1) + j
            ijsize         = bnd_elm_list%bnd_element(ibnd)%size(i,j)
            
            index_if: if ( (index_node >= index_min) .and. (index_node <= index_max) ) then ! Is the current MPI thread in charge?
  
              im_loop: do im = 1, n_tor ! (loop over toroidal harmonics)
                ! Note: Loops over i, j, and im select a test function.
                
                ilarge_j = n_tor * n_var * (index_node-1) + n_tor * (kj-1) + im
                v = H1(i,j,ms) * ijsize * HZ(im,mp) ! (test function)
                
                k_loop: do k = 1, 2 ! (loop over nodes in element ibnd)
                
                  l_loop: do l = 1, 2 ! (loop over basis functions)
  
                    ldir            = bnd_elm_list%bnd_element(ibnd)%direction(k,l)         ! ibnd marks the boundary element (not the node number)
                    index_node2     = node_list%node(bnd_elm_list%bnd_element(ibnd)%vertex(k))%index(ldir)
                    inode2_bnd      = bnd_elm_list%bnd_element(ibnd)%bnd_vertex(k)

                    index_node2_bnd = 2 * (bnd_node_list%bnd_node(inode2_bnd)%index_starwall-1) + l !#### to be checked

                    in_loop: do in = 2, n_tor ! (loop over toroidal harmonics)
                      ! Note: Loops over k, l, and in select the boundary degree of freedom
                      !       at which the vacuum response is calculated (***i***)
                      
                      rhs_contrib = 0.d0
  
                      klsize = bnd_elm_list%bnd_element(ibnd)%size(k,l)
                      
                      if ( use_starwall ) then
                        bpar = H1(k,l,ms)   * klsize * HZ(in,mp) * sqrt(x_s(ms)**2 + y_s(ms)**2)
                      else
                        bpar = H1_s(k,l,ms) * klsize * HZ(in,mp)
                      end if
        
                      kbnd_loop: do kbnd = 1, bnd_node_list%n_bnd_nodes ! (loop over boundary nodes)
        
                        korder_loop: do korder = 1, 2 ! (loop over basis functions)
                        do ktor = 1, n_tor !###############################################################<<<<<<<<<<<<<<<<<<<
                        
                          ! Note: Loops over kbnd, and korder select the boundary degree of freedom
                          !       which contributes to the vacuum response (***j***)

                          !###TODO: Understand why there is no loop over toroidal mode number indices necessary here
                          !### -> because no coupling of harmonics is taken into account currently
                                
                          kdir            = bnd_node_list%bnd_node(kbnd)%direction(korder)
                          knode           = bnd_node_list%bnd_node(kbnd)%index_jorek
                          index_node3     = node_list%node(knode)%index(kdir)
                          index_node3_bnd = 2*(kbnd-1) + korder
                          
                          call locate_irn_jcn(index_node,index_node3,index_min,index_max,ijA_position)   ! position in the global matrix
  
                          ilarge_pp  = ijA_position  - 1 + ((kp-1)*n_tor + im-1) * n_var*n_tor + (kp-1)*n_tor + ktor  ! to be verified (in<->im)
                          ilarge_jp  = ijA_position  - 1 + ((kj-1)*n_tor + im-1) * n_var*n_tor + (kp-1)*n_tor + ktor  ! to be verified (in<->im)
                          
                          ! --- Set row and column numbers in the main matrix (as some combinations are not set yet)
                          !do kv=1,n_var
                          !  do lv=1,n_var
                          !    ilarge_vv  = ijA_position  - 1 + ((kv-1)*n_tor + im-1) * n_var*n_tor + (lv-1)*n_tor + ktor
                          !    irn_glob(ilarge_vv) = n_tor * n_var * (index_node-1)  + (kv-1)*n_tor + im
                          !    jcn_glob(ilarge_vv) = n_tor * n_var * (index_node3-1) + (lv-1)*n_tor + ktor
                          !  end do
                          !end do
                          !irn_glob(ilarge_jp) =  n_tor * n_var * (index_node-1)  + (kj-1)*n_tor + im
                          !jcn_glob(ilarge_jp) =  n_tor * n_var * (index_node3-1) + (kp-1)*n_tor + ktor
                          
                          if  ( (ktor > 1) .and. (in > 1) ) then !###
                            
                            ! --- Contribution of the vacuum response to the left hand side, i.e., to the a-matrix
                            A_glob_31 = v * ws * bpar
                        
                            if ( .not. resistive_wall ) then
                              A_glob(ilarge_jp) =  A_glob(ilarge_jp) - A_glob_31 * vac_response(response_index(inode2_bnd,in,l), response_index(kbnd,ktor,korder))
                            else
                              if ( wall_curr_treatment == 'explicit' ) then
                                A_glob(ilarge_jp) =  A_glob(ilarge_jp) - A_glob_31 * matrix_ee(response_index(inode2_bnd,in,l), response_index(kbnd,ktor,korder))
                              else
                                A_glob(ilarge_jp) =  A_glob(ilarge_jp) - A_glob_31 * matrix_t (response_index(inode2_bnd,in,l), response_index(kbnd,ktor,korder))
                              end if
                            end if
                            
                            !write(35+my_id,*) ilarge_jp, - A_glob_31 * vac_response(response_index(inode2_bnd,in,l), response_index(kbnd,ktor,korder))
                            
                            ! --- Contribution to the rhs ([EE] * Psi part)
!                            if ( .not. resistive_wall ) then
!                              rhs_contrib = rhs_contrib + vac_response(response_index(kbnd,ktor,korder), response_index(inode2_bnd,in,l)) * node_list%node(knode)%values(ktor,kdir,kp)
!                            else
!                              rhs_contrib = rhs_contrib + matrix_ee      (response_index(kbnd,ktor,korder), response_index(inode2_bnd,in,l)) * node_list%node(knode)%values(ktor,kdir,kp)
!                            end if
                            
                          end if !###
                          
                        end do !###
                        end do korder_loop
                      end do kbnd_loop
                
                      ! --- Contribution to the rhs ([EE] * Psi part)
                      if ( in /= 1 ) then !###
                      if ( .not. resistive_wall ) then
                        rhs_contrib = rhs_contrib + sum( vac_response(response_index(inode2_bnd,in,l), :) * psi_bnd_vec(:) )
                      else
                        rhs_contrib = rhs_contrib + sum( matrix_ee   (response_index(inode2_bnd,in,l), :) * psi_bnd_vec(:) )
                      end if
                      
                      ! --- Further contribution to the rhs ([EY] * Y part)
                      if ( resistive_wall ) then
                        if ( wall_curr_treatment == 'explicit' ) then
                          rhs_contrib = rhs_contrib + sum( matrix_ey(response_index(inode2_bnd,in,l),:) * wall_curr(:) )
                        else
                          rhs_contrib = rhs_contrib + sum( matrix_u (response_index(inode2_bnd,in,l),:) * wall_curr(:) )
                        end if
                      end if
                      end if !###
                      rhs_loc(ilarge_j) = rhs_loc(ilarge_j) + v * ws * bpar * rhs_contrib
                      
                      !write(35+my_id,*) ilarge_j, v * ws * bpar * rhs_contrib
                      
                    end do in_loop
                  end do l_loop
                end do k_loop

              end do im_loop
            end if index_if
          end do j_loop
        end do i_loop
        
      end do mp_loop
    end do ms_loop
  end do ibnd_loop
  
  write(*,*) 'checksums:'
  write(*,*) sum(abs(rhs_loc))
  write(*,*) sum(abs(A_glob))
  write(*,*) sum(abs(irn_glob))
  write(*,*) sum(abs(jcn_glob))
  !write(35+my_id,'(4ES20.12)') sum(abs(rhs_loc)), sum(abs(A_glob)), real(sum(abs(irn_glob))), real(sum(abs(jcn_glob)))
  
  !do i = 1, size( irn_glob)
  !  write(33+my_id,*) irn_glob(i)
  !end do
  
  write(94,*)
  write(95,*)
  write(45,*)
  
  deallocate( psi_bnd_vec, dpsi_bnd_vec )
  
  return
end subroutine vacuum_old
