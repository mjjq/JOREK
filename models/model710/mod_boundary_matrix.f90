module mod_boundary_matrix

  implicit none

contains

subroutine boundary_matrix(vertex, direction, element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                psi_bnd, R_xpoint, Z_xpoint, ELM, RHS)
!---------------------------------------------------------------------
! calculates the matrix contribution of the boundaries of one element
! implements the natural boundary conditions
!---------------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(4)        ! the two nodes containing the boundary nodes

real*8, dimension (:,:), pointer  :: ELM
real*8, dimension (:)  , pointer  :: RHS

integer    :: vertex(2), direction(2), xcase2
real*8     :: psi_axis, R_axis, Z_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
logical    :: xpoint2
real*8     :: R_g(n_gauss), R_s(n_gauss), R_t(n_gauss)
real*8     :: Z_g(n_gauss), Z_s(n_gauss), Z_t(n_gauss)

real*8     :: eq_g(n_plane,n_var,n_gauss), eq_s(n_plane,n_var,n_gauss), eq_t(n_plane,n_var,n_gauss), eq_p(n_plane,n_var,n_gauss)
real*8     :: delta_g(n_plane,n_var,n_gauss), delta_s(n_plane,n_var,n_gauss), delta_t(n_plane,n_var,n_gauss)
real*8     :: Fprofile(n_gauss)

real*8     :: Qbnd(n_var), Qjac(n_var,n_var)

integer    :: i, j, j2, ms, mt, mp, k, l, l2, index_ij, index_kl, ij, kl
integer    :: in, im, ivar, kvar
integer    :: j3, direction_perp(2)
real*8     :: ws, xjac,  BigR, phi, DL, Zbig
real*8     :: R_mid, Z_mid, R_cnt, Z_cnt
real*8     :: theta, zeta

real*8     :: T, r0, rho, c_s, cs_T
real*8     :: AR0, AR0_p, AR0_s, AR0_t, AR0_R, AR0_Z
real*8     :: AZ0, AZ0_p, AZ0_s, AZ0_t, AZ0_R, AZ0_Z     
real*8     :: A30, A30_p, A30_s, A30_t, A30_R, A30_Z
real*8     :: uR0, uR0_s, uR0_t, uR0_R, uR0_Z
real*8     :: uZ0, uZ0_s, uZ0_t, uZ0_R, uZ0_Z
real*8     :: UP0, Up
real*8     :: T0, T0_p, T0_s, T0_t, T0_R, T0_Z
real*8     :: AR, AR_p, AR_s, AR_t, AR_R, AR_Z
real*8     :: AZ, AZ_p, AZ_s, AZ_t, AZ_R, AZ_Z     
real*8     :: A3, A3_p, A3_s, A3_t, A3_R, A3_Z
real*8     :: uR, uR_s, uR_t, uR_R, uR_Z
real*8     :: uZ, uZ_s, uZ_t, uZ_R, uZ_Z
real*8     :: bf, bf_s, bf_t, bf_p, bf_R, bf_Z

real*8     :: BB2, BB2_AR, BB2_AZ, BB2_A3
real*8     :: BR0, BR0_AR, BR0_AZ, BR0_A3
real*8     :: BZ0, BZ0_AR, BZ0_AZ, BZ0_A3
real*8     :: Bp0, Bp0_AR, Bp0_AZ, Bp0_A3

real*8     :: B_dot_n,B_dot_n_AR, B_dot_n_AZ, B_dot_n_A3
real*8     :: U_dot_n

real*8     :: v, v_x, v_y, v_s, v_p, v_ss, v_xx, v_yy, v_xs, v_ys
real*8     :: element_size_ij, element_size_kl, element_size_perp
real*8     :: normal(2), normal_direction(2)
real*8     :: grad_s(2), grad_t(2)

theta = time_evol_theta
zeta  = time_evol_zeta


zbig = 1.d8

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
!!!! s is the coordinate along the boundary, t is the other direction
!!!! this can be different from the element (s,t) orientations

R_g  = 0.d0; R_s  = 0.d0;  R_t = 0.d0; 
Z_g  = 0.d0; Z_s  = 0.d0;  Z_t = 0.d0;
eq_g = 0.d0; eq_s = 0.d0;  eq_t = 0.d0; eq_p = 0.d0;

delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0;

Fprofile = 0.d0


R_mid = sum(nodes(1:2)%x(1,1)) / 2.d0     ! mid point on boundary (approx.)
Z_mid = sum(nodes(1:2)%x(1,2)) / 2.d0
R_cnt = sum(nodes(1:4)%x(1,1)) / 4.d0     ! center point within element (approx.)
Z_cnt = sum(nodes(1:4)%x(1,2)) / 4.d0

normal_direction = (/R_mid - R_cnt, Z_mid - Z_cnt /) / norm2((/R_mid - R_cnt, Z_mid - Z_cnt /))

direction_perp(1) = 6 / direction(2)     ! =3 if direction(2)=2, =3 if direction(2)=3
direction_perp(2) = 4


do i=1,2        ! loop over nodes

  do j=1,2      ! loop over basis functions

    j2 = direction(j)

    element_size_ij = element%size(vertex(i),j2)

    j3 = direction_perp(j)
    element_size_perp = - element%size(vertex(i),direction_perp(1)) * 3.d0

    do ms=1, n_gauss

      R_g(ms)  = R_g(ms)  + nodes(i)%x(j2,1) * element_size_ij * H1(i,j,ms)
      R_s(ms)  = R_s(ms)  + nodes(i)%x(j2,1) * element_size_ij * H1_s(i,j,ms)
      R_t(ms)  = R_t(ms)  + nodes(i)%x(j3,1) * element_size_ij * H1(i,j,ms)   * element_size_perp

      Z_g(ms)  = Z_g(ms)  + nodes(i)%x(j2,2) * element_size_ij * H1(i,j,ms)
      Z_s(ms)  = Z_s(ms)  + nodes(i)%x(j2,2) * element_size_ij * H1_s(i,j,ms)
      Z_t(ms)  = Z_t(ms)  + nodes(i)%x(j3,2) * element_size_ij * H1(i,j,ms)   * element_size_perp

      Fprofile(ms)   = Fprofile(ms)   + nodes(i)%Fprof_eq(j2)    * element_size_ij * H1(i,j,ms)

      do mp=1,n_plane

        do k=1,n_var

          do in=1,n_tor

            eq_g(mp,k,ms)  = eq_g(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ(in,mp)
            eq_s(mp,k,ms)  = eq_s(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1_s(i,j,ms) * HZ(in,mp)

            eq_t(mp,k,ms)  = eq_t(mp,k,ms)  + nodes(i)%values(in,j3,k) * element_size_ij * H1(i,j,ms)   * HZ(in,mp) * element_size_perp

            eq_p(mp,k,ms)  = eq_p(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ_p(in,mp)

            delta_g(mp,k,ms) = delta_g(mp,k,ms) + nodes(i)%deltas(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ(in,mp)
            delta_s(mp,k,ms) = delta_s(mp,k,ms) + nodes(i)%deltas(in,j2,k) * element_size_ij * H1_s(i,j,ms) * HZ(in,mp)

          enddo
        enddo
      enddo

    enddo
  enddo
enddo

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

   ws = wgauss(ms)

   xjac = R_s(ms) * Z_t(ms) - R_t(ms) * Z_s(ms)

   grad_s = (/   Z_t(ms), - R_t(ms) /) / xjac

   grad_t = (/ - Z_s(ms),   R_s(ms) /) / xjac

   normal = dot_product(grad_t,normal_direction) * grad_t

   normal = normal / norm2(normal)

   DL   = sqrt(R_s(ms)**2 + Z_s(ms)**2)
   BigR = R_g(ms)

   Qbnd = 0.d0
   Qjac = 0.d0

   do mp = 1, n_plane

     uR0   = eq_g(mp,var_uR,ms)
     uZ0   = eq_g(mp,var_uZ,ms)

     U_dot_n = UR0 * normal(1) + UZ0 * normal(2)

     do i=1,2                ! loop over nodes

       do j=1,2              ! loop over basis functions

         j2 = direction(j)
         element_size_ij = element%size(vertex(i),j2)

         do im=1,n_tor

           index_ij = n_tor*n_var*(n_order+1)*(vertex(i)-1) + n_tor * n_var * (j2-1) + im   ! index in the ELM matrix

           v   =  H1(i,j,ms) * element_size_ij * HZ(im,mp)         ! test function

           Qbnd(var_uR) = zbig * v * U_dot_n * normal(1) 

           Qbnd(var_uZ) = zbig * v * U_dot_n * normal(2) 

           do ivar= 1,n_var
             
             ij = index_ij + (ivar-1)*n_tor
           
             RHS(ij) =  RHS(ij) + ws * Qbnd(ivar) !* BigR * DL

           enddo

           do k=1,2                                                          ! loop over nodes

             do l=1,2                                                        ! loop over basis functions

               l2 = direction(l)

               element_size_kl   = element%size(vertex(k),l2)

               element_size_perp = - element%size(vertex(k),direction_perp(1)) * 3.d0


               do in = 1, n_tor                                              ! loop over toroidal harmonics
     
                 bf   = H1(k,l,ms)   * element_size_kl * HZ(in,mp)

                 uR    = bf    ;  uZ    = bf    ;  up    = bf


                 Qjac(var_uR,var_uR) = - zbig * v * uR * normal(1) * normal(1)
                 Qjac(var_uR,var_uZ) = - zbig * v * uZ * normal(2) * normal(1)


                 Qjac(var_uZ,var_uR) = - zbig * v * uR * normal(1) * normal(2)
                 Qjac(var_uZ,var_uZ) = - zbig * v * uZ * normal(2) * normal(2)


                 index_kl = n_tor*n_var*(n_order+1)*(vertex(k)-1) + n_tor * n_var * (l2-1) + in   ! index in the ELM matrix 

                 do ivar= 1,n_var
                   
                   do kvar= 1,n_var

                     ij = index_ij + (ivar-1)*n_tor
                     kl = index_kl + (kvar-1)*n_tor

                     ELM(ij,kl) =  ELM(ij,kl) + ws * theta * Qjac(ivar,kvar) !* BigR * DL

                   enddo
                 enddo

               enddo
             enddo
           enddo

         enddo
       enddo
     enddo

   enddo
enddo

return

end subroutine boundary_matrix
end module mod_boundary_matrix


