subroutine vacuum(my_id,node_list,element_list,boundary_list,index_min,index_max, &
                  xpoint2,psi_axis,psi_bnd,Z_xpoint)
!---------------------------------------------------------------------
! calculates the matrix contribution of the boundary integral of the
! induction equation using the vacuum response from STARWALL
!---------------------------------------------------------------------
use parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use global_distributed_matrix
use vacuum_response_module

implicit none

type (type_node_list)     :: node_list
type (type_element_list)  :: element_list
type (type_boundary_list) :: boundary_list

real*8     :: x_g(n_gauss), x_s(n_gauss), x_ss(n_gauss)
real*8     :: y_g(n_gauss), y_s(n_gauss), y_ss(n_gauss)

real*8     :: eq_g(n_plane,n_var,n_gauss), eq_s(n_plane,n_var,n_gauss), eq_p(n_plane,n_var,n_gauss), eq_ss(n_plane,n_var,n_gauss)
real*8     :: delta_g(n_plane,n_var,n_gauss), delta_s(n_plane,n_var,n_gauss)

integer    :: index_min, index_max
integer    :: my_id, ibnd,i, j, ms, mp, kp, kbnd, k, l, jdir, kdir, ldir, imode, korder, kv, lv,  ilarge_vv, inode
integer    :: index_node, index_node2, index_node3, index_node_bnd, index_node2_bnd, index_node3_bnd, ilarge_pp, ijA_position
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, kl1, kl2, kl3, kl4, kl5, kl6, kl7
real*8     :: ws, xjac,  BigR, PI, phi, eps_cyl
real*8     :: psi_axis, psi_bnd, Z_xpoint
real*8     :: rhs_glob_1, A_glob_11
real*8     :: psi_norm, theta, zeta, gamma_sheeth

real*8     :: v, v_x, v_y, v_s, v_p, v_ss, v_xx, v_yy, v_xs, v_ys
real*8     :: ps0, ps0_s, r0, T0
real*8     :: psi, psi_s, rho,  T, eta_T
real*8     :: amat_61, amat_65, amat_66, amat_67
integer    :: itmp1, itmp2
logical    :: xpoint2

write(*,*) '************************************'
write(*,*) '*     VACUUM boundary integral     *'
write(*,*) '************************************'

PI    = 2.d0*asin(1.d0)

theta = 0.5d0; zeta = 0.d0          ! Crank-Nicholson parameter
!theta = 1.0d0  ; zeta = 0.0d0       ! Euler scheme
!theta = 1.0d0   ; zeta = 0.5d0      ! BDF2 (Gears) scheme


do ibnd = 1, boundary_list%n_boundary                     ! loop over all boundary elements

  x_g  = 0.d0; x_s  = 0.d0;  x_ss  = 0.d0;                ! values of (x,y) and derivatives on Gaussian points
  y_g  = 0.d0; y_s  = 0.d0;  y_ss  = 0.d0;
  eq_g = 0.d0; eq_s = 0.d0;  eq_ss = 0.d0; eq_p = 0.d0;

  delta_g = 0.d0; delta_s = 0.d0;

  do i=1,2                                                ! loop over two corners of each boundary element

    do j=1,2                                              ! loop over the two basis functions (H, H_s)

      do ms=1, n_gauss                                    ! loop over Gaussian points, construct coordinates and values

        inode = boundary_list%boundary(ibnd)%vertex(i)
        jdir  = boundary_list%boundary(ibnd)%direction(i,j)

        x_g(ms)  = x_g(ms)  + node_list%node(inode)%x(jdir,1) &
                            * boundary_list%boundary(ibnd)%size(i,j) * H1(i,j,ms)

        x_s(ms)  = x_s(ms)  + node_list%node(inode)%x(jdir,1) &
                            * boundary_list%boundary(ibnd)%size(i,j) * H1_s(i,j,ms)

        y_g(ms)  = y_g(ms)  + node_list%node(inode)%x(jdir,2) &
                            * boundary_list%boundary(ibnd)%size(i,j) * H1(i,j,ms)

        y_s(ms)  = y_s(ms)  + node_list%node(inode)%x(jdir,2) &
                            * boundary_list%boundary(ibnd)%size(i,j) * H1_s(i,j,ms)

        do mp=1,n_plane                                   ! loop over toroidal planes

          do k=1,n_var                                    ! loop over all variables

            do in=1,n_tor                                 ! loop over toroidal harmonics

              inode = boundary_list%boundary(ibnd)%vertex(i)

              jdir  = boundary_list%boundary(ibnd)%direction(i,j)

              eq_g(mp,k,ms)  = eq_g(mp,k,ms) + node_list%node(inode)%values(in,jdir,k) &
                                             * boundary_list%boundary(ibnd)%size(i,j)       * H1(i,j,ms)   * HZ(in,mp)

              eq_s(mp,k,ms)  = eq_s(mp,k,ms) + node_list%node(inode)%values(in,jdir,k) &
                                             * boundary_list%boundary(ibnd)%size(i,j)       * H1_s(i,j,ms) * HZ(in,mp)

              eq_p(mp,k,ms)  = eq_p(mp,k,ms) + node_list%node(inode)%values(in,jdir,k) &
                                             * boundary_list%boundary(ibnd)%size(i,j)       * H1(i,j,ms)   * HZ_p(in,mp)

              delta_g(mp,k,ms) = delta_g(mp,k,ms)  + node_list%node(inode)%deltas(in,jdir,k) &
                                             * boundary_list%boundary(ibnd)%size(i,j)       * H1(i,j,ms)   * HZ(in,mp)

	      delta_s(mp,k,ms) = delta_s(mp,k,ms)  + node_list%node(inode)%deltas(in,jdir,k) &
                                             * boundary_list%boundary(ibnd)%size(i,j)       * H1_s(i,j,ms)   * HZ(in,mp)

            enddo                                         ! end loop over toroidal harmonics
          enddo                                           ! end loop over all variables
        enddo                                             ! end loop over toroidal planes

      enddo                                               ! end loop over Gaussian points
    enddo                                                 ! end loop over the two basis function
  enddo                                                   ! end loop over two corners of each boundary element


  do ms=1, n_gauss                                        ! loop over Gaussian points

    ws = wgauss(ms)

    do mp = 1, n_plane                                    ! sum over toroidal planes (direct integration in toroidal direction)

      ps0   = eq_g(mp,1,ms)
      ps0_s = eq_s(mp,1,ms)

      r0    = eq_g(mp,5,ms)
      T0    = eq_g(mp,6,ms)

      BigR = x_g(ms)

      eta_T   = eta * (abs(T0)/T_0)**(-1.5d0)                                 ! temperature dependent resistivity

      do i=1,2                                                                  ! loop over nodes of this piece of boundary

        do j=1,2                                                                ! loop over basis functions (test functions)

          jdir = boundary_list%boundary(ibnd)%direction(i,j)

          index_node = node_list%node(boundary_list%boundary(ibnd)%vertex(i))%index(jdir)

          index_node_bnd = 2*mod(ibnd+i-2,boundary_list%n_boundary) + (j-1) + 1   ! the index in the vacuum_response matrix (DANGEROUS ASSUMPTION!)

          if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

            do im=1,n_tor                                                       ! loop over toroidal harmonics

              v =  H1(i,j,ms) * boundary_list%boundary(ibnd)%size(i,j) * HZ(im,mp)         ! test function

              do kbnd = 1, boundary_list%n_boundary                              ! kbnd really numbers the nodes of the boundary list

                do korder = 1, 2                                                 ! loop over basis_functions at node kbnd

                  kdir = boundary_list%boundary(kbnd)%direction(1,korder)         ! ibnd marks the boundary element (not the node number) i.e. TO BE CHANGED

                  index_node3 = node_list%node(boundary_list%boundary(kbnd)%vertex(1))%index(kdir)

                  index_node3_bnd = 2*(kbnd-1) + (korder-1) + 1                  ! the index in the vacuum_response matrix (DANGEROUS ASSUMPTION!)

                  do k=1,2                                                       ! loop over nodes in element ibnd (resulting from perturbation at node k_bnd)

                    do l=1,2                                                     ! loop over basis functions

                      ldir = boundary_list%boundary(ibnd)%direction(k,l)         ! ibnd marks the boundary element (not the node number)

                      index_node2 = node_list%node(boundary_list%boundary(ibnd)%vertex(k))%index(ldir)

                      index_node2_bnd = 2*mod(ibnd+k-2,boundary_list%n_boundary) + (l-1) + 1   ! the index in the vacuum_response matrix (DANGEROUS ASSUMPTION!)

                      call locate_irn_jcn(index_node,index_node3,index_min,index_max,ijA_position)   ! position in the global matrix

                      kp = 1                            ! select psi variable

                      do in = 1, n_tor                                                       ! loop over toroidal harmonics

                        ilarge_pp  = ijA_position  - 1 + ((kp-1)*n_tor + im-1) * n_var*n_tor + (kp-1)*n_tor + in  ! to be verified (in<->im)

                        psi   =  H1(k,l,ms)   * boundary_list%boundary(ibnd)%size(k,l) * HZ(in,mp)         ! test function
                        psi_s =  H1_s(k,l,ms) * boundary_list%boundary(ibnd)%size(k,l) * HZ(in,mp)         ! test function derivative (along boundary, i.e. t)

                        rho = psi
		        T   = psi

                        do kv=1,n_var
                          do lv=1,n_var
                            ilarge_vv  = ijA_position  - 1 + ((kv-1)*n_tor + im-1) * n_var*n_tor + (lv-1)*n_tor + in
                            itmp1 = irn_glob(ilarge_vv)
                            itmp2 = jcn_glob(ilarge_vv)
                            irn_glob(ilarge_vv) = n_tor * n_var * (index_node-1)  + (kv-1)*n_tor + im
                            jcn_glob(ilarge_vv) = n_tor * n_var * (index_node3-1) + (lv-1)*n_tor + in
                          enddo
                        enddo

                        if ( use_starwall ) then
                          A_glob_11 = v * eta_T * psi * theta * tstep
                        else
                         A_glob_11 = v * eta_T * psi_s * theta * tstep
                        end if

                        A_glob_11 = A_glob_11 * vacuum_response(index_node3_bnd,index_node2_bnd,im)

                        irn_glob(ilarge_pp) =  n_tor * n_var * (index_node-1)  + (kp-1)*n_tor + im
                        jcn_glob(ilarge_pp) =  n_tor * n_var * (index_node3-1) + (kp-1)*n_tor + in

                        A_glob(ilarge_pp)   =  A_glob(ilarge_pp) + ws * A_glob_11

!                        write(*,*) imode, index_node2_bnd, index_node3_bnd, vacuum_response(index_node2_bnd,index_node3_bnd,imode)
!                        write(*,'(A,6i8,e14.6)') 'indices : ',index_node,index_node2,index_node3,ilarge_pp,irn_glob(ilarge_pp),jcn_glob(ilarge_pp),A_glob(ilarge_pp)

                      enddo   ! end loop over toroidal harmonics

                    enddo     ! end loop over basis functions
                  enddo       ! end of loop over all boundary nodes

                enddo       ! end loop over nodes
              enddo         ! end of loop over all boundary elements
            enddo           ! end loop over toroidal harmonics

	  endif           ! endif selection of local indices

        enddo             ! end of loop over basis functions
      enddo               ! end of loop over nodes
    enddo                 ! end of loop over toroidal planes

  enddo                   ! end of loop over Gaussian points

enddo                     ! end of loop over all boundary elements

return
end
