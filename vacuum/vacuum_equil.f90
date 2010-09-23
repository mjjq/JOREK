subroutine vacuum_equil(node_list,element_list,boundary_list,xpoint2,psi_axis,psi_bnd,Z_xpoint)
!---------------------------------------------------------------------
! calculates the matrix contribution of the boundary integral of the
! induction equation using the vacuum response from STARWALL
!---------------------------------------------------------------------
use parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use mumps_module
use vacuum_response_module

implicit none

type (type_node_list)        :: node_list
type (type_element_list)     :: element_list
type (type_bnd_element_list) :: boundary_list

real*8     :: x_g(n_gauss), x_s(n_gauss)
real*8     :: y_g(n_gauss), y_s(n_gauss)

real*8     :: eq_g(n_gauss), eq_s(n_gauss)

integer    :: my_id, ibnd,i, j, ms,  kp, kbnd, k, l, jdir, kdir, ldir, imode, korder, kv, lv,  ilarge_vv, inode
integer    :: index_node, index_node2, index_node3, index_node_bnd, index_node2_bnd, index_node3_bnd, ilarge_pp, ijA_position
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, kl1, kl2, kl3, kl4, kl5, kl6, kl7
real*8     :: ws, xjac,  BigR, PI, phi, eps_cyl
real*8     :: psi_axis, psi_bnd, Z_xpoint
real*8     :: rhs_glob_1, A_glob_11, A_glob_11_star
real*8     :: psi_norm, theta, zeta

real*8     :: v, v_x, v_y, v_s, v_p, v_ss, v_xx, v_yy, v_xs, v_ys
real*8     :: ps0, ps0_s, r0, T0
real*8     :: psi, psi_s, rho,  T, eta_T

integer    :: itmp1, itmp2, index_min, index_max, ilarge
logical    :: xpoint2

write(*,*) '**************************************************'
write(*,*) '*     VACUUM boundary integral (equilibrium)     *'
write(*,*) '**************************************************'

PI    = 2.d0*asin(1.d0)

ilarge = mumps_par%nz                                     ! the number of items in the mumps_par coordinate storage (before the boundary conditions)

do ibnd = 1, boundary_list%n_bnd_elements                     ! loop over all boundary elements

  x_g  = 0.d0; x_s  = 0.d0;                               ! values of (x,y) and derivatives on Gaussian points
  y_g  = 0.d0; y_s  = 0.d0;

  do i=1,2                                                ! loop over two corners of each boundary element

    do j=1,2                                              ! loop over the two basis functions (H, H_s)

      do ms=1, n_gauss                                    ! loop over Gaussian points, construct coordinates and values

        inode = boundary_list%bnd_element(ibnd)%vertex(i)
        jdir  = boundary_list%bnd_element(ibnd)%direction(i,j)

        x_g(ms)  = x_g(ms)  + node_list%node(inode)%x(jdir,1) * boundary_list%bnd_element(ibnd)%size(i,j) * H1(i,j,ms)
        x_s(ms)  = x_s(ms)  + node_list%node(inode)%x(jdir,1) * boundary_list%bnd_element(ibnd)%size(i,j) * H1_s(i,j,ms)
        y_g(ms)  = y_g(ms)  + node_list%node(inode)%x(jdir,2) * boundary_list%bnd_element(ibnd)%size(i,j) * H1(i,j,ms)
        y_s(ms)  = y_s(ms)  + node_list%node(inode)%x(jdir,2) * boundary_list%bnd_element(ibnd)%size(i,j) * H1_s(i,j,ms)

      enddo                                               ! end loop over Gaussian points
    enddo                                                 ! end loop over the two basis function
  enddo                                                   ! end loop over two corners of each boundary element


  do ms=1, n_gauss                                        ! loop over Gaussian points

    ws = wgauss(ms)

    BigR = x_g(ms)

    do i=1,2                                                                 ! loop over nodes of this piece of boundary

      do j=1,2                                                               ! loop over basis functions (test functions)

        jdir = boundary_list%bnd_element(ibnd)%direction(i,j)

        index_node = node_list%node(boundary_list%bnd_element(ibnd)%vertex(i))%index(jdir)

        index_node_bnd = 2*mod(ibnd+i-2,boundary_list%n_bnd_elements) + j        ! the index in the vacuum_response matrix (DANGEROUS ASSUMPTION!)

        v =  H1(i,j,ms) * boundary_list%bnd_element(ibnd)%size(i,j)             ! test function

 !!!!       mumps_par%rhs(index_node) = mumps_par%rhs(index_node) + ws * vac_ext(ms) ! boundary integral of external vacuum field (coils)

        do kbnd = 1, boundary_list%n_bnd_elements                                ! kbnd really numbers the nodes of the boundary list

          do korder = 1, 2                                                   ! loop over basis_functions at node kbnd

            kdir = boundary_list%bnd_element(kbnd)%direction(1,korder)          ! ibnd marks the boundary element (not the node number) i.e. TO BE CHANGED

            index_node3 = node_list%node(boundary_list%bnd_element(kbnd)%vertex(1))%index(kdir)

            index_node3_bnd = 2*(kbnd-1) + (korder-1) + 1                   ! the index in the vacuum_response matrix (DANGEROUS ASSUMPTION!)

            do k=1,2                                                        ! loop over nodes in element ibnd (resulting from perturbation at node k_bnd)

              do l=1,2                                                      ! loop over basis functions

                ldir = boundary_list%bnd_element(ibnd)%direction(k,l)          ! ibnd marks the boundary element (not the node number)

                index_node2 = node_list%node(boundary_list%bnd_element(ibnd)%vertex(k))%index(ldir)

                index_node2_bnd = 2*mod(ibnd+k-2,boundary_list%n_bnd_elements) + (l-1) + 1         ! the index in the vacuum_response matrix (DANGEROUS ASSUMPTION!)

                call locate_irn_jcn(index_node,index_node3,index_min,index_max,ijA_position)   ! position in the global matrix

                ilarge_pp  = ijA_position

                psi   =  H1(k,l,ms)   * boundary_list%bnd_element(ibnd)%size(k,l)          ! test function
                psi_s =  H1_s(k,l,ms) * boundary_list%bnd_element(ibnd)%size(k,l)          ! test function derivative (along boundary, i.e. t)

                if ( use_starwall ) then
                  A_glob_11 = v * psi * sqrt(x_s(ms)**2 + y_s(ms)**2)
!                  A_glob_11 = v * psi_s
                else
                  A_glob_11 = v * psi_s
                end if

!                irn_glob(ilarge_pp) =  n_tor * n_var * (index_node-1)  + 1
!                jcn_glob(ilarge_pp) =  n_tor * n_var * (index_node3-1) + 1

!                A_glob(ilarge_pp)   =  A_glob(ilarge_pp) + ws * A_glob_11 * vacuum_response(index_node3_bnd,index_node2_bnd,im)

                mumps_par%irn(ilarge) = (index_node-1)  + 1
                mumps_par%jcn(ilarge) = (index_node3-1) + 1
                mumps_par%A(ilarge)   = ws * A_glob_11 * vacuum_response(index_node3_bnd,index_node2_bnd,1)

                ilarge = ilarge + 1


              enddo     ! end loop over basis functions (l)
            enddo       ! end of loop elemet nodes (k)

          enddo         ! end loop over order (korder)
        enddo           ! end of loop over all boundary elements (kbnd)

      enddo             ! end loop over basis functions (j)
    enddo               ! end loop over nodes in this boundary element

  enddo                 ! end of loop over Gaussian points
enddo                   ! end of loop over all boundary elements

return
end
