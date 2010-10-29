subroutine vacuum(my_id,node_list,element_list,boundary_list,bnd_node_list, &
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
  use vacuum_response_module
  
  
  implicit none
  
  
  integer,                      intent(in)    :: my_id          ! Number of current MPI thread.
  type (type_node_list),        intent(in)    :: node_list      ! Node list of current grid.
  type (type_element_list),     intent(in)    :: element_list   ! Element list of current grid.
  type (type_bnd_element_list), intent(in)    :: boundary_list  ! List of 1D boundary elements.
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
  
  integer    :: ibnd,i, j, ms, mp, kp, kj, kbnd, k, l, jdir, kdir, ldir, korder, kv, lv,  ilarge_vv, inode, inode_bnd, inode2_bnd
  integer    :: index_node, index_node2, index_node3, index_node_bnd, index_node2_bnd, index_node3_bnd
  integer    :: ilarge_pp, ilarge_jp, ijA_position
  integer    :: in, im
  integer    :: jpsi, iwall, jnode, jtor, jbas, jindex, jnode_glob, jnode2, jelem
  real*8     :: dPdt, YE_dPdt
  real*8     :: ws, BigR, PI
  real*8     :: A_glob_31
  real*8     :: theta, zeta
  real*8     :: v, ijsize, klsize
  real*8     :: ps0, ps0_s, r0, T0
  real*8     :: psi, psi_s, rho,  T, eta_T
  integer    :: icurr
  real*8     :: value, delta
  integer    :: knode, knode3
  integer    :: ilarge_j
  real*8     :: br_coils, bz_coils, psi_coils
  real*8     :: vac_contrib
  real*8     :: psi_val
  
  
  write(*,*) '************************************'
  write(*,*) '*     VACUUM boundary integral     *'
  write(*,*) '************************************'
  
  PI    = 2.d0*asin(1.d0)
  
  theta = 0.5d0; zeta = 0.d0          ! Crank-Nicholson parameter
  !theta = 1.0d0  ; zeta = 0.0d0       ! Euler scheme
  !theta = 1.5d0   ; zeta = 0.5d0      ! BDF2 (Gears) scheme
  
  
  
  ! --- Time-integrate currents in the resistive wall.
  !     (Y are wall current potentials, Psi denotes the poloidal flux)
  !
  !     dY/dt = - 1/(sigma * d) * [YY] * Y - [YE] * dPsi/dt
  !
  if ( resistive_wall ) then
  
    ! --- First part with YY matrix
    wall_curr(:) = wall_curr(:) - tstep * ( diagonal_yy(:) * wall_curr(:) ) * wall_resistivity / wall_thickness
    
    ! --- Second part with YE matrix
    jpsi=1 ! select poloidal flux variable
    do iwall = 1, n_wall_curr                       ! loop over wall current potentials
      YE_dPdt = 0.d0
      do jnode = 1, bnd_node_list%n_bnd_nodes       ! loop over nodes
        jnode_glob = bnd_node_list%bnd_node(jnode)%index_jorek
        do jtor = 2, 3                              ! loop over toroidal harmonics
          do jbas = 1, 2                            ! loop over basis functions

            jindex = 2*(jnode-1) +   (jbas-1) + 1   ! second index in response matrix
            
            ! --- Determine jdir of current node (essentially which coordinate direction is along the boundary).
            !     (maybe replace later by a backreference of the boundary nodes to their elements)
            JEL_LOOP: do jelem = 1, boundary_list%n_bnd_elements  ! loop over boundary elements
              do jnode2 = 1, 2                                    ! loop over the two nodes of the element
                if ( boundary_list%bnd_element(jelem)%bnd_vertex(jnode2) == jnode ) then
                  jdir = boundary_list%bnd_element(jelem)%direction(jnode2,jbas)
                  exit JEL_LOOP
                else if ( (jelem == boundary_list%n_bnd_elements) .AND. (jnode2 == 2) ) then
                  write(*,*) 'ERROR: Could not determine jdir. Did not find boundary node.'
                  stop
                end if
              end do
            end do JEL_LOOP
            
            dPdt = node_list%node(jnode_glob)%deltas(jtor, jdir, jpsi)/tstep ! dPsi/dt at the current node
            
            YE_dPdt = YE_dPdt + matrix_ye(iwall, jindex, jtor) * dPdt        ! [YE] * dP/dt (sum over all nodes)
            
          end do
        end do
      end do
      wall_curr(iwall) = wall_curr(iwall) - tstep * YE_dPdt
    end do
    
  end if
  
  
  
  ibnd_loop: do ibnd = 1, boundary_list%n_bnd_elements ! (loop over boundary elements)
  
    ! --- Determine coordinate and variable values at Gaussian points on boundary elements.
    x_g  = 0.d0; x_s  = 0.d0;
    y_g  = 0.d0; y_s  = 0.d0; 
    eq_g = 0.d0; eq_s = 0.d0; eq_p = 0.d0;
    delta_g = 0.d0; delta_s = 0.d0;
  
    do i = 1 , 2 ! (loop over nodes of boundary element)
      do j = 1, 2 ! (loop over basis functions)
        do ms = 1, n_gauss ! (loop over Gaussian points)
  
          inode  = boundary_list%bnd_element(ibnd)%vertex(i)
          jdir   = boundary_list%bnd_element(ibnd)%direction(i,j)
          ijsize = boundary_list%bnd_element(ibnd)%size(i,j)
  
          x_g(ms)  = x_g(ms)  + node_list%node(inode)%x(jdir,1) * ijsize * H1(i,j,ms)     ! x:=R
          x_s(ms)  = x_s(ms)  + node_list%node(inode)%x(jdir,1) * ijsize * H1_s(i,j,ms)   ! dx/ds
          y_g(ms)  = y_g(ms)  + node_list%node(inode)%x(jdir,2) * ijsize * H1(i,j,ms)     ! y:=Z
          y_s(ms)  = y_s(ms)  + node_list%node(inode)%x(jdir,2) * ijsize * H1_s(i,j,ms)   ! dZ/ds
  
          do mp = 1, n_plane ! (loop over toroidal planes)
            do k = 1, n_var ! (loop over all variables)
              do in = 1, n_tor ! (loop over toroidal harmonics)
  
                value = node_list%node(inode)%values(in,jdir,k)
                delta = node_list%node(inode)%deltas(in,jdir,k)
  
                eq_g(mp,k,ms)    = eq_g(mp,k,ms)    + value * ijsize * H1(i,j,ms)   * HZ(in,mp)    ! variable value
                eq_s(mp,k,ms)    = eq_s(mp,k,ms)    + value * ijsize * H1_s(i,j,ms) * HZ(in,mp)    ! dvar/ds
                eq_p(mp,k,ms)    = eq_p(mp,k,ms)    + value * ijsize * H1(i,j,ms)   * HZ_p(in,mp)  ! dvar/dphi
                delta_g(mp,k,ms) = delta_g(mp,k,ms) + delta * ijsize * H1(i,j,ms)   * HZ(in,mp)    ! var change in last timestep
                delta_s(mp,k,ms) = delta_s(mp,k,ms) + delta * ijsize * H1_s(i,j,ms) * HZ(in,mp)    ! ddelta/ds
  
              end do
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
        eta_T = eta !#### * (abs(T0)/T_0)**(-1.5d0)         ! temperature dependent resistivity
        
        i_loop: do i = 1, 2 ! (loop over nodes of boundary element)
          j_loop: do j = 1, 2 ! (loop over basis functions)
  
            jdir           = boundary_list%bnd_element(ibnd)%direction(i,j)
            inode          = boundary_list%bnd_element(ibnd)%vertex(i)
            inode_bnd      = boundary_list%bnd_element(ibnd)%bnd_vertex(i)
            index_node     = node_list%node(inode)%index(jdir)
            index_node_bnd = 2 * (bnd_node_list%bnd_node(inode_bnd)%index_starwall-1) + j
            ijsize         = boundary_list%bnd_element(ibnd)%size(i,j)
            
            index_if: if ( (index_node >= index_min) .and. (index_node <= index_max) ) then ! Is the current MPI thread in charge?
  
              im_loop: do im = 1, n_tor ! (loop over toroidal harmonics)
                ! Note: Loops over i, j, and im select a test function.
                
                v = H1(i,j,ms) * boundary_list%bnd_element(ibnd)%size(i,j) * HZ(im,mp) ! (test function)
                
                ! --- Determine the vacuum contribution to the rhs (ideal: M * Psi0; resitive: M * Psi0 + N * Y + C)
                vac_contrib = 0.d0
                !   --- M * psi0 contribution (poloidal flux)
                do kbnd = 1, bnd_node_list%n_bnd_nodes ! (loop over bnd nodes)
                  knode3 = bnd_node_list%bnd_node(kbnd)%index_jorek
                  do l = 1, 2 ! (loop over basis functions)
                    index_node3_bnd = 2 * (bnd_node_list%bnd_node(kbnd)%index_starwall-1) + l
                    kdir            = bnd_node_list%bnd_node(kbnd)%direction(l)
                    do in = im, im !###
                      psi_val = node_list%node(knode3)%values(in,kdir,kp)
                      if ( resistive_wall ) then
                        vac_contrib = vac_contrib - matrix_ee(index_node3_bnd,index_node_bnd,in)       * psi_val
                      else
                        vac_contrib = vac_contrib - vacuum_response(index_node3_bnd,index_node_bnd,in) * psi_val
                      end if
                    end do
                  end do
                end do
                if ( resistive_wall ) then
                  !   --- N * Y contribution (wall currents)
                  do icurr = 1, n_wall_curr
                    vac_contrib = vac_contrib + matrix_ey(index_node_bnd,icurr,im) * wall_curr(icurr)
                  end do
                  !   --- C contribution (coils)
                  !call pfcoils(x_g(ms), y_g(ms), br_coils, bz_coils, psi_coils)
                  !vac_contrib = vac_contrib - sqrt(x_s**2 + y_s**2) * ( bz_coils * y_s(ms) + br_coils * x_s(ms) )
                end if
                
                ! --- Add the vacuum contribution to the rhs.
                ilarge_j = n_tor * n_var * (index_node-1) + n_tor * (kj-1) + im
                !### debug output
                if ( (mp==1) .and. (im==2) ) then
                  if ( (i==1).and.(j==1) ) write(95,*) ibnd, vac_contrib
                  if ( (i==2).and.(j==1) ) write(94,*) ibnd, vac_contrib
                end if
                !### end debug output
                if ( use_starwall ) then
                  vac_contrib = vac_contrib * H1(i,j,ms)   * sqrt(x_s(ms)**2 + y_s(ms)**2)
                else
                  vac_contrib = vac_contrib * H1_s(i,j,ms)
                end if
                !rhs_loc(ilarge_j) = rhs_loc(ilarge_j) + vac_contrib * ijsize * HZ(im,mp) * tstep * v * ws
  
                kbnd_loop: do kbnd = 1, bnd_node_list%n_bnd_nodes ! (loop over boundary nodes)
  
                  korder_loop: do korder = 1, 2 ! (loop over basis functions)
                    ! Note: Loops over kbnd, and korder select the boundary degree of freedom
                    !       ???????????? which one ??????????

                    !###TODO: Understand why there is no loop over toroidal mode number indices necessary here
                    !### -> because no coupling of harmonics is taken into account currently
                          
                    kdir            = bnd_node_list%bnd_node(kbnd)%direction(korder)
                    knode           = bnd_node_list%bnd_node(kbnd)%index_jorek
                    index_node3     = node_list%node(knode)%index(kdir)
                    index_node3_bnd = 2*(kbnd-1) + korder
                    
                    k_loop: do k = 1, 2 ! (loop over nodes in element ibnd)
  
                      l_loop: do l = 1, 2 ! (loop over basis functions)
  
                        ldir            = boundary_list%bnd_element(ibnd)%direction(k,l)         ! ibnd marks the boundary element (not the node number)
                        index_node2     = node_list%node(boundary_list%bnd_element(ibnd)%vertex(k))%index(ldir)
                        inode2_bnd      = boundary_list%bnd_element(ibnd)%bnd_vertex(k)

                        index_node2_bnd = 2 * (bnd_node_list%bnd_node(inode2_bnd)%index_starwall-1) + l !#### to be checked

                        call locate_irn_jcn(index_node,index_node3,index_min,index_max,ijA_position)   ! position in the global matrix
  
                        in_loop: do in = 1, n_tor ! (loop over toroidal harmonics)
                          ! Note: Loops over k, l, and in select the boundary degree of freedom
                          !       ??????????which one?????????
  
                          ilarge_pp  = ijA_position  - 1 + ((kp-1)*n_tor + im-1) * n_var*n_tor + (kp-1)*n_tor + in  ! to be verified (in<->im)
                          ilarge_jp  = ijA_position  - 1 + ((kj-1)*n_tor + im-1) * n_var*n_tor + (kp-1)*n_tor + in  ! to be verified (in<->im)
                          
                          klsize = boundary_list%bnd_element(ibnd)%size(k,l)
                          
                          psi    = H1(k,l,ms)   * klsize * HZ(in,mp)
                          rho    = psi
                          T      = psi
                          psi_s  = H1_s(k,l,ms) * klsize * HZ(in,mp)
  
                          ! --- Set row and column numbers in the main matrix (as some combinations are not set yet)
                          do kv=1,n_var
                            do lv=1,n_var
                              ilarge_vv  = ijA_position  - 1 + ((kv-1)*n_tor + im-1) * n_var*n_tor + (lv-1)*n_tor + in
                              irn_glob(ilarge_vv) = n_tor * n_var * (index_node-1)  + (kv-1)*n_tor + im
                              jcn_glob(ilarge_vv) = n_tor * n_var * (index_node3-1) + (lv-1)*n_tor + in
                            end do
                          end do
                          irn_glob(ilarge_jp) =  n_tor * n_var * (index_node-1)  + (kj-1)*n_tor + im
                          jcn_glob(ilarge_jp) =  n_tor * n_var * (index_node3-1) + (kp-1)*n_tor + in
                          
                          if  (im .eq. in) then
                            
                            ! --- Contribution of the vacuum response to the left hand side, i.e., to the a-matrix
                            if ( use_starwall ) then
                              A_glob_31 = tstep * v * ws * psi * sqrt(x_s(ms)**2 + y_s(ms)**2)
                            else
                              A_glob_31 = tstep * v * ws * psi_s
                            end if
                            
                            if ( resistive_wall ) then
                              A_glob(ilarge_jp) =  A_glob(ilarge_jp) + A_glob_31 * matrix_ee(index_node3_bnd,index_node2_bnd,im)
                            else
                              A_glob(ilarge_jp) =  A_glob(ilarge_jp) + A_glob_31 * vacuum_response(index_node3_bnd,index_node2_bnd,im)
                            end if
                            
                          end if
                          
                        end do in_loop
                      end do l_loop
                    end do k_loop

                  end do korder_loop
                end do kbnd_loop
                
              end do im_loop
            end if index_if
          end do j_loop
        end do i_loop
        
      end do mp_loop
    end do ms_loop
  end do ibnd_loop
  
  write(94,*)
  write(95,*)
    
  return
end subroutine vacuum
