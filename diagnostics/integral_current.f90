subroutine integral_current(node_list,element_list,psi_axis, psi_bnd, xpoint2, xcase2, z_xpoint, current)
!---------------------------------------------------------------
!
!---------------------------------------------------------------
use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
real*8     :: ps_g(n_gauss,n_gauss)

integer    :: i, j, ms, mt, iv, inode, ife, n_elements, xcase2
real*8     :: zn
real*8     ::    dn_dpsi, dn_dz                                                ! 1st order derivatives
real*8     ::    dn_dpsi2, dn_dz2, dn_dpsi_dz                                  ! 2nd order derivatives
real*8     ::    dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3                   ! 2rd order derivatives
real*8     ::    dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz, dn_dz4      ! 4th order derivatives
real*8     ::    dn_dpsi5, dn_dpsi_dz4, dn_dpsi2_dz3, dn_dpsi3_dz2, dn_dpsi4_dz! 5th order derivatives (z5 not needed)
real*8     :: zT
real*8     ::    dT_dpsi,  dT_dz                                               ! 1st order derivatives
real*8     ::    dT_dpsi2, dT_dz2, dT_dpsi_dz                                  ! 2nd order derivatives
real*8     ::    dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,  dT_dz3                   ! 2rd order derivatives
real*8     ::    dT_dpsi4, dT_dpsi_dz3, dT_dpsi2_dz2, dT_dpsi3_dz, dT_dz4      ! 4th order derivatives
real*8     ::    dT_dpsi5, dT_dpsi_dz4, dT_dpsi2_dz3, dT_dpsi3_dz2, dT_dpsi4_dz! 5th order derivatives (z5 not needed)
real*8     :: zTi
real*8     ::    dTi_dpsi,  dTi_dz                                                  ! 1st order derivatives
real*8     ::    dTi_dpsi2, dTi_dz2, dTi_dpsi_dz                                    ! 2nd order derivatives
real*8     ::    dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz,  dTi_dz3                    ! 2rd order derivatives
real*8     ::    dTi_dpsi4, dTi_dpsi_dz3, dTi_dpsi2_dz2, dTi_dpsi3_dz, dTi_dz4      ! 4th order derivatives
real*8     ::    dTi_dpsi5, dTi_dpsi_dz4, dTi_dpsi2_dz3, dTi_dpsi3_dz2, dTi_dpsi4_dz! 5th order derivatives (z5 not needed)
real*8     :: zTe
real*8     ::    dTe_dpsi,  dTe_dz                                                  ! 1st order derivatives
real*8     ::    dTe_dpsi2, dTe_dz2, dTe_dpsi_dz                                    ! 2nd order derivatives
real*8     ::    dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz,  dTe_dz3                    ! 2rd order derivatives
real*8     ::    dTe_dpsi4, dTe_dpsi_dz3, dTe_dpsi2_dz2, dTe_dpsi3_dz, dTe_dz4      ! 4th order derivatives
real*8     ::    dTe_dpsi5, dTe_dpsi_dz4, dTe_dpsi2_dz3, dTe_dpsi3_dz2, dTe_dpsi4_dz! 5th order derivatives (z5 not needed)
real*8     :: zFFprime
real*8     ::    dFF_dpsi, dFF_dz                                                ! 1st order derivatives
real*8     ::    dFF_dpsi2, dFF_dz2, dFF_dpsi_dz                                 ! 2nd order derivatives
real*8     ::    dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3                 ! 2rd order derivatives
real*8     ::    dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4   ! 4th order derivatives
real*8     :: current, xjac, BigR, Z_xpoint(2), psi_axis, psi_bnd, wst
logical    :: xpoint2

current = 0.d0

do ife =1,  element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:)   = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0;
  y_g(:,:)   = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0;
  ps_g(:,:)  = 0.d0

  do i=1,n_vertex_max
    do j=1,n_degrees
      do ms=1, n_gauss
        do mt=1, n_gauss

          x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)

          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

          ps_g(ms,mt)  = ps_g(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H(i,j,ms,mt)

        enddo
      enddo
    enddo
  enddo

!--------------------------------------------------- sum over the Gaussian integration points

  do ms=1, n_gauss

    do mt=1, n_gauss

      call density(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, ps_g(ms,mt),psi_axis,psi_bnd,zn, &
                   dn_dpsi,  dn_dz, &                                             ! 1st order derivatives
                   dn_dpsi2, dn_dz2,      dn_dpsi_dz, &                           ! 2nd order derivatives
                   dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3, &                 ! 2rd order derivatives
                   dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz,  dn_dz4, &   ! 4th order derivatives
                   dn_dpsi5, dn_dpsi_dz4, dn_dpsi2_dz3, dn_dpsi3_dz2, dn_dpsi4_dz)! 5th order derivatives (z5 not needed)

      if (with_TiTe) then
        
        call temperature_i(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, ps_g(ms,mt),psi_axis,psi_bnd, zTi, &
                           dTi_dpsi,  dTi_dz, &                                                ! 1st order derivatives
                           dTi_dpsi2, dTi_dz2,      dTi_dpsi_dz, &                             ! 2nd order derivatives
                           dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz,  dTi_dz3, &                  ! 2rd order derivatives
                           dTi_dpsi4, dTi_dpsi_dz3, dTi_dpsi2_dz2, dTi_dpsi3_dz,  dTi_dz4, &   ! 4th order derivatives
                           dTi_dpsi5, dTi_dpsi_dz4, dTi_dpsi2_dz3, dTi_dpsi3_dz2, dTi_dpsi4_dz)! 5th order derivatives (z5 not needed)
        	       
        call temperature_e(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, ps_g(ms,mt),psi_axis,psi_bnd, zTe, &
                           dTe_dpsi,  dTe_dz, &                                                ! 1st order derivatives
                           dTe_dpsi2, dTe_dz2,      dTe_dpsi_dz, &                             ! 2nd order derivatives
                           dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz,  dTe_dz3, &                  ! 2rd order derivatives
                           dTe_dpsi4, dTe_dpsi_dz3, dTe_dpsi2_dz2, dTe_dpsi3_dz,  dTe_dz4, &   ! 4th order derivatives
                           dTe_dpsi5, dTe_dpsi_dz4, dTe_dpsi2_dz3, dTe_dpsi3_dz2, dTe_dpsi4_dz)! 5th order derivatives (z5 not needed)

        zT = zTi + zTe
        dT_dpsi = dTi_dpsi + dTe_dpsi
      
      else
        
        call temperature(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, ps_g(ms,mt),psi_axis,psi_bnd, zT, &
                         dT_dpsi,  dT_dz, &                                             ! 1st order derivatives
                         dT_dpsi2, dT_dz2,      dT_dpsi_dz, &                           ! 2nd order derivatives
                         dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,  dT_dz3, &                 ! 2rd order derivatives
                         dT_dpsi4, dT_dpsi_dz3, dT_dpsi2_dz2, dT_dpsi3_dz,  dT_dz4, &   ! 4th order derivatives
                         dT_dpsi5, dT_dpsi_dz4, dT_dpsi2_dz3, dT_dpsi3_dz2, dT_dpsi4_dz)! 5th order derivatives (z5 not needed)

      endif

      call FFprime(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, ps_g(ms,mt),psi_axis,psi_bnd, zFFprime, &
                   dFF_dpsi, dFF_dz, &                                             ! 1st order derivatives
                   dFF_dpsi2, dFF_dz2, dFF_dpsi_dz, &                              ! 2nd order derivatives
                   dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3, &              ! 2rd order derivatives
                   dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4, &! 4th order derivatives
                   .true.)
                   
      wst = wgauss(ms)*wgauss(mt)

      xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
      BigR = x_g(ms,mt)
           
      current = current + (zFFprime / x_g(ms,mt) - (zn * dT_dpsi + dn_dpsi * zT) * x_g(ms,mt)) * wst * xjac

    enddo
  enddo
enddo

current = current / (4.d-7 * PI)

write(*,'(A,f8.5,A)') ' current : ',current/1.e6,' MA'
return
end
