module mod_elt_matrix
  implicit none
contains

subroutine element_matrix(element, nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid)

!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
!                                             dP/dt  = Q                       
!
! Integrands weak form (excl R * xjac) :      amat * sol_coeff = rhs_ij
!
! where in the Beam-Warming scheme     :      amat   = (1+zeta)Pjac - tstep * theta * Qjac
!                                             rhs_ij = tstep * Qvec + zeta * Pvec_prev
!
! Pjac & Qjac are the derivatives of P and Q wrt all variables (Jacobian), 
! Pvec_prev is the vector resulting from the matrix multiplication of Pjac of the previous timestep with the deltas vector
!
!
! * Model 710 consists of the full viscoresistive MHD equations ( resistive heating is however not (yet) included )
! * Some form of documentation will be provided in a forthcoming paper
!
!---------------------------------------------------------------

use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use diffusivities, only: get_dperp, get_zkperp

implicit none

include 'mpif.h'

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8, dimension (:,:), pointer  :: ELM
real*8, dimension (:)  , pointer  :: RHS
integer, intent(in) :: tid

integer    :: i, j, k, l, index_ij, index_kl, index, xcase2, inode
integer    :: in, im, ij, kl, ivar, kvar, ms, mt, mp
real*8     :: wst, xjac, xjac3, BigR, phi
real*8     :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss), source_pellet
real*8     :: minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz

real*8     :: uR0, uR0_R, uR0_Z, uR0_p, uR0_s, uR0_t 
real*8     :: uZ0, uZ0_R, uZ0_Z, uZ0_p, uZ0_s, uZ0_t 
real*8     :: up0, up0_R, up0_Z, up0_p, up0_s, up0_t
real*8     :: AR0, AR0_R, AR0_Z, AR0_p, AR0_s, AR0_t
real*8     :: AZ0, AZ0_R, AZ0_Z, AZ0_p, AZ0_s, AZ0_t
real*8     :: A30, A30_R, A30_Z, A30_p, A30_s, A30_t
real*8     :: r0, r0_R, r0_Z, r0_p, r0_s, r0_t
real*8     :: T0, T0_R, T0_Z, T0_p, T0_s, T0_t
real*8     :: p0, p0_R, p0_Z, p0_p, p0_s, p0_t
real*8     :: uR, uR_R, uR_Z, uR_p, uR_s, uR_t
real*8     :: uZ, uZ_R, uZ_Z, uZ_p, uZ_s, uZ_t
real*8     :: up, up_R, up_Z, up_p, up_s, up_t
real*8     :: AR, AR_R, AR_Z, AR_p, AR_s, AR_t
real*8     :: AZ, AZ_R, AZ_Z, AZ_p, AZ_s, AZ_t
real*8     :: A3, A3_R, A3_Z, A3_p, A3_s, A3_t
real*8     :: T, T_R, T_Z, T_p, T_s, T_t
real*8     :: r, r_R, r_Z, r_p, r_s, r_t

real*8     :: v, v_R, v_Z, v_s, v_t, v_p
real*8     :: bf, bf_R, bf_Z, bf_s, bf_t, bf_p
   
real*8     :: B0grad_T0,    B0grad_T0_AR,    B0grad_T0_AZ,    B0grad_T0_A3
real*8     :: B0grad_r0,    B0grad_r0_AR,    B0grad_r0_AZ,    B0grad_r0_A3
real*8     :: BB2,          BB2_AR,          BB2_AZ,          BB2_A3
real*8     :: BR0,          BR0_AR,          BR0_AZ,          BR0_A3,          BR0_red
real*8     :: BZ0,          BZ0_AR,          BZ0_AZ,          BZ0_A3,          BZ0_red 
real*8     :: Bp0,          Bp0_AR,          Bp0_AZ,          Bp0_A3,          Bp0_red
real*8     :: B0grad_vstar,  B0grad_vstar_AR,  B0grad_vstar_AZ,  B0grad_vstar_A3
real*8     :: u0grad_vstar, u0grad_T0,    u0grad_r0,       gradr0grad_vstar, gradT0grad_vstar, gradbfgrad_vstar, u0grad_bf, B0grad_bf

real*8     :: u0grad_uR0, u0grad_uZ0, u0grad_up0

real*8     :: divu, divu_uR, divu_uZ, divu_up, divru, divru_uR, divru_uZ, divru_up, divru_r

real*8     :: ZK_prof, D_prof, psi_norm, theta, zeta, tht

real*8     :: Fprof, psieq_R, psieq_Z

real*8     :: eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, visco_num_T, eta_num_T, eta_R, eta_Z, eta_p, Zkpar_T, dZKpar_dt

real*8     :: Qvisc_uR, Qvisc_uZ, Qvisc_up, Qvisc_T, Q_uR_primitive, Q_uZ_primitive, Q_up_primitive

logical    :: xpoint2, viscores_heating, primitive, parallel_projection, offset_current, offset_current2

real*8, dimension(3          )        :: Vpar, Vperp
real*8, dimension(n_var      )        :: rhs_ij, Pvec_prev, Qvec, Vms, TG_NUM=0.0
real*8, dimension(n_var,n_var)        :: amat, Pjac, Qjac
real*8, dimension(n_gauss,n_gauss)    :: Fprofile,Fprofile_s,Fprofile_t, psieq, psieq_s, psieq_t

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t

real*8, dimension(:,:,:,:) , pointer :: eq_g, eq_s, eq_t,eq_p
real*8, dimension(:,:,:,:) , pointer :: delta_g, delta_s, delta_t

integer*4  :: rank
integer    :: my_id, ierr
! Variables for the VMS stabilization
! -----------------------------------
real*8 :: VdotB, BigR2
real*8 :: CvR0, CvZ0, Cvp0, CvGradAR0, CvGradAZ0, CvGradA30, CvGradr0, CvGradT0
real*8 :: VbR0, VbZ0, Vbp0, VbGradAR0, VbGradAZ0, VbGradA30, VbGradr0, VbGradT0
real*8 :: CvGraduR0, CvGraduZ0, CvGradup0, VbGraduR0, VbGraduZ0, VbGradup0

real*8 :: CvGradVi, VbGradVi, DiveRMVi, DiveZMVi, DivePMVi
real*8 :: CvGradVj, VbGradVj, DiveRMVj, DiveZMVj, DivePMVj

real*8 :: TG_NUM_Eq = 0.0, VmsCoefF, VmsCoefF_T, CoefAdv=0.0

real*8, dimension(n_var,n_var)        :: QvmsAd, QvmsF

! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------
eq_g    => thread_struct(tid)%eq_g   
eq_s    => thread_struct(tid)%eq_s   
eq_t    => thread_struct(tid)%eq_t   
eq_p    => thread_struct(tid)%eq_p   
delta_g => thread_struct(tid)%delta_g
delta_s => thread_struct(tid)%delta_s
delta_t => thread_struct(tid)%delta_t


ELM       = 0.d0
RHS       = 0.d0
rhs_ij    = 0.d0
amat      = 0.d0
Pjac      = 0.d0
Qjac      = 0.d0
QvmsAd      = 0.d0
Pvec_prev = 0.d0
Qvec      = 0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Some switches

! --- Take time evolution parameters from phys_module
theta = time_evol_theta
zeta  = time_evol_zeta

viscores_heating = .false.

! (partial) primitive weak formulation (no partial intergration where it doesn't reduce the order)
!                                      (only implemented for divergence terms now)
primitive = .true.                     

! Project in the magnetic field direction instead of the toroidal direction 
! (only implemented for momentum equation, and using the primitive formulation)
parallel_projection                 = .true.                     
if( parallel_projection ) primitive = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; 
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; 
eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_p = 0.d0

psieq   = 0.d0; psieq_s = 0.d0; psieq_t = 0.d0
delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0

current_source  = 0.d0
particle_source = 0.d0
heat_source     = 0.d0
Fprofile        = 0.d0 ; Fprofile_s        = 0.d0 ;  Fprofile_t        = 0.d0

!  Very simple stabilization :
!     the same foeficient for parallel and perp
!     and also for all variables
! ----------------------------------------------
TG_NUM = 0.25*tstep
TG_NUM(var_AR) = 0.d0
TG_NUM(var_AZ) = 0.d0
TG_NUM(var_A3) = 0.d0

TG_NUM = 0.d0

call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
my_id = rank

do i=1,n_vertex_max
 do j=1,n_order+1

   do ms=1, n_gauss
     do mt=1, n_gauss

       x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
       x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
       x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)

       y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H(i,j,ms,mt)
       y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
       y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)

#ifdef fullmhd
           Fprofile(ms,mt) = Fprofile(ms,mt) + nodes(i)%Fprof_eq(j) * element%size(i,j) * H(i,j,ms,mt)
!           Fprofile_s(ms,mt) = Fprofile(ms,mt) + nodes(i)%Fprof_eq(j) * element%size(i,j) * H_s(i,j,ms,mt)
!           Fprofile_t(ms,mt) = Fprofile(ms,mt) + nodes(i)%Fprof_eq(j) * element%size(i,j) * H_t(i,j,ms,mt)

           psieq(ms,mt)    = psieq(ms,mt)    + nodes(i)%psi_eq(j)   * element%size(i,j) * H(i,j,ms,mt)
           psieq_s(ms,mt)  = psieq_s(ms,mt)  + nodes(i)%psi_eq(j)   * element%size(i,j) * H_s(i,j,ms,mt)
           psieq_t(ms,mt)  = psieq_t(ms,mt)  + nodes(i)%psi_eq(j)   * element%size(i,j) * H_t(i,j,ms,mt)
#endif


       do mp=1,n_plane

         do k=1,n_var

           do in=1,n_tor

             eq_g(mp,k,ms,mt)    = eq_g(mp,k,ms,mt)    + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ(in,mp)
             eq_s(mp,k,ms,mt)    = eq_s(mp,k,ms,mt)    + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt) * HZ(in,mp)
             eq_t(mp,k,ms,mt)    = eq_t(mp,k,ms,mt)    + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) * HZ(in,mp)
             eq_p(mp,k,ms,mt)    = eq_p(mp,k,ms,mt)    + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ_p(in,mp)

             delta_g(mp,k,ms,mt) = delta_g(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ(in,mp)
             delta_s(mp,k,ms,mt) = delta_s(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt) * HZ(in,mp)
             delta_t(mp,k,ms,mt) = delta_t(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) * HZ(in,mp)

           enddo

         enddo

       enddo

     enddo
   enddo
 enddo
enddo

do ms=1, n_gauss
  do mt=1, n_gauss
       call current(xpoint2, xcase2, x_g(ms,mt),y_g(ms,mt), Z_xpoint, psieq(ms,mt),psi_axis,psi_bnd,current_source(ms,mt))
       call sources(xpoint2, xcase2, y_g(ms,mt)           , Z_xpoint, psieq(ms,mt),psi_axis,psi_bnd,particle_source(ms,mt),heat_source(ms,mt))
  enddo
enddo


!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

 do mt=1, n_gauss

   wst = wgauss(ms)*wgauss(mt)

   xjac    = x_s(ms,mt)*y_t(ms,mt)  - x_t(ms,mt)*y_s(ms,mt)
   
   BigR    = x_g(ms,mt)
   BigR2   = BigR*BigR

   xjac3 = BigR * xjac 

   Fprof = Fprofile(ms,mt)

   do mp = 1, n_plane

     uR0   = eq_g(mp,var_uR,ms,mt)
     uR0_p = eq_p(mp,var_uR,ms,mt)
     uR0_s = eq_s(mp,var_uR,ms,mt)
     uR0_t = eq_t(mp,var_uR,ms,mt)
     uR0_R = (   y_t(ms,mt) * uR0_s  - y_s(ms,mt) * uR0_t ) / xjac
     uR0_Z = ( - x_t(ms,mt) * uR0_s  + x_s(ms,mt) * uR0_t ) / xjac

     uZ0   = eq_g(mp,var_uZ,ms,mt)
     uZ0_p = eq_p(mp,var_uZ,ms,mt)
     uZ0_s = eq_s(mp,var_uZ,ms,mt)
     uZ0_t = eq_t(mp,var_uZ,ms,mt)
     uZ0_R = (   y_t(ms,mt) * uZ0_s  - y_s(ms,mt) * uZ0_t ) / xjac
     uZ0_Z = ( - x_t(ms,mt) * uZ0_s  + x_s(ms,mt) * uZ0_t ) / xjac

!-------------------------------------------up0 is defined u0_phi : V = .. + up0 * e_phi (physical component)
     up0   = eq_g(mp,var_up,ms,mt)
     up0_p = eq_p(mp,var_up,ms,mt)
     up0_s = eq_s(mp,var_up,ms,mt)
     up0_t = eq_t(mp,var_up,ms,mt)
     up0_R = (   y_t(ms,mt) * up0_s  - y_s(ms,mt) * up0_t ) / xjac
     up0_Z = ( - x_t(ms,mt) * up0_s  + x_s(ms,mt) * up0_t ) / xjac

     AR0   = eq_g(mp,var_AR,ms,mt)
     AR0_p = eq_p(mp,var_AR,ms,mt)
     AR0_s = eq_s(mp,var_AR,ms,mt)
     AR0_t = eq_t(mp,var_AR,ms,mt)
     AR0_R = (   y_t(ms,mt) * AR0_s  - y_s(ms,mt) * AR0_t ) / xjac
     AR0_Z = ( - x_t(ms,mt) * AR0_s  + x_s(ms,mt) * AR0_t ) / xjac

     AZ0   = eq_g(mp,var_AZ,ms,mt)
     AZ0_p = eq_p(mp,var_AZ,ms,mt)
     AZ0_s = eq_s(mp,var_AZ,ms,mt)
     AZ0_t = eq_t(mp,var_AZ,ms,mt)
     AZ0_R = (   y_t(ms,mt) * AZ0_s  - y_s(ms,mt) * AZ0_t ) / xjac
     AZ0_Z = ( - x_t(ms,mt) * AZ0_s  + x_s(ms,mt) * AZ0_t ) / xjac

!-------------------------------------------A30 == psi
!-------------------------------------------A30 is defined as A0_3 : A = .. + A03 * grad(phi)
     A30   = eq_g(mp,var_A3,ms,mt)
     A30_p = eq_p(mp,var_A3,ms,mt)
     A30_s = eq_s(mp,var_A3,ms,mt)
     A30_t = eq_t(mp,var_A3,ms,mt)
     A30_R = (   y_t(ms,mt) * A30_s  - y_s(ms,mt) * A30_t ) / xjac
     A30_Z = ( - x_t(ms,mt) * A30_s  + x_s(ms,mt) * A30_t ) / xjac
    
     r0    = eq_g(mp,var_r,ms,mt)
     r0_p  = eq_p(mp,var_r,ms,mt)
     r0_s  = eq_s(mp,var_r,ms,mt)
     r0_t  = eq_t(mp,var_r,ms,mt)
     r0_R  = (   y_t(ms,mt) * r0_s  - y_s(ms,mt) * r0_t ) / xjac
     r0_Z  = ( - x_t(ms,mt) * r0_s  + x_s(ms,mt) * r0_t ) / xjac

     T0    = eq_g(mp,var_T,ms,mt)
     T0_p  = eq_p(mp,var_T,ms,mt)
     T0_s  = eq_s(mp,var_T,ms,mt)
     T0_t  = eq_t(mp,var_T,ms,mt)
     T0_R  = (   y_t(ms,mt) * T0_s  - y_s(ms,mt) * T0_t ) / xjac
     T0_Z  = ( - x_t(ms,mt) * T0_s  + x_s(ms,mt) * T0_t ) / xjac

     p0    = r0 * T0
     p0_R  = r0_R * T0 + r0 * T0_R
     p0_Z  = r0_Z * T0 + r0 * T0_Z
     p0_s  = r0_s * T0 + r0 * T0_s
     p0_t  = r0_t * T0 + r0 * T0_t
     p0_p  = r0_p * T0 + r0 * T0_p

     psi_norm = (A30 - psi_axis) / (psi_bnd - psi_axis)

     if (xpoint2) then
       if ((psi_norm .lt. 1.d0) .and. (y_g(ms,mt) .lt. Z_xpoint(1)) .and. (xcase2 .ne. 2)) then
         psi_norm = 2.d0 - psi_norm
       endif
       if ((psi_norm .lt. 1.d0) .and. (y_g(ms,mt) .gt. Z_xpoint(2)) .and. (xcase2 .ne. 1)) then
         psi_norm = 2.d0 - psi_norm
       endif
     endif

     D_prof  = get_dperp (psi_norm)
     ZK_prof = get_zkperp(psi_norm)


     ! --- Temperature dependent resistivity
     if ( eta_T_dependent ) then

       eta_T     =   eta   * (abs(T0)/T_0)**(-1.5d0)
       deta_dT   = - eta   * (1.5d0)  * abs(T0)**(-2.5d0) * T_0**(1.5d0)
       d2eta_d2T =   eta   * (3.75d0) * abs(T0)**(-3.5d0) * T_0**(1.5d0)

       if (T0 .lt. T_min) then
         eta_T     = eta * (T_min/T_0)**(-1.5d0)
         deta_dT   = 0.d0
         d2eta_d2T = 0.d0
       endif

     else
       eta_T     = eta
       deta_dT   = 0.d0
       d2eta_d2T = 0.d0
     end if

     eta_R = deta_dT * T0_R
     eta_Z = deta_dT * T0_Z
     eta_p = deta_dT * T0_p

     ! --- Temperature dependent viscosity
     if ( visco_T_dependent ) then

       visco_T   = visco * (abs(T0)/T_0)**(-1.5d0)
       dvisco_dT = - visco * (1.5d0)  * abs(T0)**(-2.5d0) * T_0**(1.5d0)

       if (T0 .lt. T_min) then
         visco_T     = visco * (T_min/T_0)**(-1.5d0)
         dvisco_dT   = 0.d0
       endif

     else
       visco_T   = visco
       dvisco_dT = 0.d0
     end if

     if ( ZKpar_T_dependent ) then

       ZKpar_T   = ZK_par * (abs(T0)/T_0)**(+2.5d0)              ! temperature dependent parallel conductivity
       dZKpar_dT = ZK_par * (2.5d0)  * abs(T0)**(+1.5d0) * T_0**(-2.5d0)

!       if (ZKpar_T .gt. ZK_par_max) then
!         ZKpar_T   = Zk_par_max
!         dZKpar_dT = 0.d0
!       endif

       if (T0 .lt. T_min) then
         ZKpar_T   = ZK_par * (T_min/T_0)**(+2.5d0)
         dZKpar_dT = 0.d0
       endif

     else
       ZKpar_T   = ZK_par                                            ! parallel conductivity
       dZKpar_dT = 0.d0
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Some auxiliary quantities (taken out of the loops over vertices, order, and toroidal harmonics)

     BR0 = ( A30_Z - AZ0_p )/ BigR
     BZ0 = ( AR0_p - A30_R )/ BigR
     Bp0 = ( AZ0_R - AR0_Z )       +   Fprof / BigR

     BB2 = Bp0**2 + BR0**2 + BZ0**2


     ! Variables used For the VMS Stabilization
     ! -----------------------------------------------------------------
     VdotB  = uR0*BR0 + up0*Bp0 + uZ0*BZ0
     
     VbR0   = VdotB*BR0/BB2
     VbZ0   = VdotB*BZ0/BB2
     Vbp0   = VdotB*Bp0/BB2

     CvR0   = uR0 - VbR0
     CvZ0   = uZ0 - VbZ0
     Cvp0   = up0 - Vbp0
     ! -----------------------------------------------------------------

     

     ! Magnitude of the current density: |curl(B)|
!     JR0  = 
!     JZ0  = 
!     Jp0  = 
!     J2  = TODO

     divu            = uR0_R + (uR0 / BigR) + uZ0_Z + (up0_p / BigR)

     B0grad_T0       = BR0 * T0_R + BZ0 * T0_Z + Bp0 * T0_p / BigR
     B0grad_r0       = BR0 * r0_R + BZ0 * r0_Z + Bp0 * r0_p / BigR
     u0grad_T0       = uR0 * T0_R + uZ0 * T0_Z + up0 * T0_p / BigR
     u0grad_r0       = uR0 * r0_R + uZ0 * r0_Z + up0 * r0_p / BigR

     ! Variables used For the VMS Stabilization
     ! -----------------------------------------------------------------
     CvGradAR0       = CvR0 * AR0_R + CvZ0 * AR0_Z + Cvp0 * AR0_p / BigR
     CvGradAZ0       = CvR0 * AZ0_R + CvZ0 * AZ0_Z + Cvp0 * AZ0_p / BigR
     CvGradA30       = CvR0 * A30_R + CvZ0 * A30_Z + Cvp0 * A30_p / BigR
     
     CvGraduR0       = CvR0 * uR0_R + CvZ0 * uR0_Z + Cvp0 * uR0_p / BigR
     CvGraduZ0       = CvR0 * uZ0_R + CvZ0 * uZ0_Z + Cvp0 * uZ0_p / BigR
     CvGradup0       = CvR0 * up0_R + CvZ0 * up0_Z + Cvp0 * up0_p / BigR

     CvGradr0        = CvR0 * r0_R  + CvZ0 * r0_Z  + Cvp0 * r0_p / BigR
     CvGradT0        = CvR0 * T0_R  + CvZ0 * T0_Z  + Cvp0 * T0_p / BigR

     VbGraduR0       = VbR0 * uR0_R + VbZ0 * uR0_Z + Vbp0 * uR0_p / BigR
     VbGraduZ0       = VbR0 * uZ0_R + VbZ0 * uZ0_Z + Vbp0 * uZ0_p / BigR
     VbGradup0       = VbR0 * up0_R + VbZ0 * up0_Z + Vbp0 * up0_p / BigR
     
     VbGradr0        = VbR0 * r0_R  + VbZ0 * r0_Z  + Vbp0 * r0_p / BigR
     VbGradT0        = VbR0 * T0_R  + VbZ0 * T0_Z  + Vbp0 * T0_p / BigR

     ! Coeficients for Acoustic waves stabilization
     ! to remove the associated stabilization put these coeficients to zero
     ! ---------------------------------------------------------------------
     VmsCoefF        = gamma * T0 + BB2/r0
     VmsCoefF_T      = gamma * T0 

     ! -----------------------------------------------------------------

     
!     if (primitive) then
       u0grad_uR0      = uR0 * uR0_R + uZ0 * uR0_Z + up0 * uR0_p / BigR
       u0grad_uZ0      = uR0 * uZ0_R + uZ0 * uZ0_Z + up0 * uZ0_p / BigR
       u0grad_up0      = uR0 * up0_R + uZ0 * up0_Z + up0 * up0_p / BigR
       divru           = r0 * divu + u0grad_r0
!     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     phi = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)


     
     do i=1,n_vertex_max

       do j=1,n_order+1


         do im=1,n_tor

           v   = H(i,j,ms,mt)   * element%size(i,j) * HZ(im,mp)
           v_s = H_s(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_t = H_t(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_p = H(i,j,ms,mt)   * element%size(i,j) * HZ_p(im,mp)

           v_R = (  y_t(ms,mt) * v_s - y_s(ms,mt) * v_t ) / xjac
           v_Z = (- x_t(ms,mt) * v_s + x_s(ms,mt) * v_t ) / xjac 


           ! Some more auxiliary quantities, involving the test functions 
           B0grad_vstar     = BR0  * v_R  + BZ0  * v_Z  + Bp0 * v_p  / BigR
           u0grad_vstar     = uR0  * v_R  + uZ0  * v_Z  + up0 * v_p  / BigR
           gradr0grad_vstar = r0_R * v_R  + r0_Z * v_Z  + (r0_p / BigR) * (v_p  / BigR)
           gradT0grad_vstar = T0_R * v_R  + T0_Z * v_Z  + (T0_p / BigR) * (v_p  / BigR)
           ! End auxiliary quantities



           ! Variables used For the VMS Stabilization associated to the test function
           ! ------------------------------------------------------------------------
           CvGradVi        = CvR0 * v_R   + CvZ0 * v_Z   + Cvp0 * v_p / BigR     
           VbGradVi        = VbR0 * v_R   + VbZ0 * v_Z   + Vbp0 * v_p / BigR
           
           DiveRMVi        = v_R + v / BigR
           DiveZMVi        = v_z
           DivePMVi        = v_p / BigR





!###################################################################################################
!#  equation 1 (R component induction equation)                                                    #
!###################################################################################################

           Pvec_prev(var_AR) =   v * delta_g(mp,var_AR,ms,mt)

           Qvec(var_AR) = - eta_T * ( v_p * BZ0 / BigR - v_Z * Bp0 )                         &
                            +   v * ( (uZ0 +  eta_Z) * Bp0 - (up0 + eta_p / BigR ) * BZ0 )

           ! adding the stabilization contribution
           ! --------------------------------------------------------------------------------------
           Vms(var_AR)  =   (   CvGradAr0 + Cvp0 * A30 / BigR2             ) * CvGradVi        &
                &         + ( - CvGradA30 + CvR0 * A30 / BigR + Cvp0 * AR0 ) * Cvp0 * v / BigR2

           Vms(var_AR)   = CvGradAr0 * CvGradVi

           Qvec(var_AR)  =  Qvec(var_AR) - TG_NUM_Eq*TG_NUM(var_AR) * CoefAdv * Vms(var_AR)  
!###################################################################################################
!#  equation 2 (Z component induction equation)                                                    #
!###################################################################################################

           Pvec_prev(var_AZ) =   v * delta_g(mp,var_AZ,ms,mt)

           Qvec(var_AZ) = - eta_T * ( - v_p * BR0 / BigR + v_R * Bp0 )                         &
                            +   v * ( (up0 + eta_p / BigR) * BR0 - (uR0 + eta_R) * Bp0 )           

           ! adding the stabilization contribution
           ! --------------------------------------------------------------------------------------
           Vms(var_AZ)  =  CvGradAZ0 * CvGradVi 

           Qvec(var_AZ)  =   Qvec(var_AZ) - TG_NUM_Eq*TG_NUM(var_AZ) * CoefAdv * Vms(var_AZ)  
!###################################################################################################
!#  equation 3 (PHI component induction equation)                                                  #
!###################################################################################################

           Pvec_prev(var_A3) =   v * delta_g(mp,var_A3,ms,mt)

           Qvec(var_A3) = - eta_T * ( BigR * v_Z * BR0 - ( 2.d0 * v + BigR * v_R ) * BZ0 + v * current_source(ms,mt) )     &
                            + BigR * v * ( (uR0 + eta_R ) * BZ0 - (uZ0 + eta_Z) * BR0 )

           ! adding the stabilization contribution
           ! --------------------------------------------------------------------------------------
           Vms(var_A3)  =   (   CvGradAr0 + Cvp0 * A30 / BigR2             ) * Cvp0 * v                   &
                &         - ( - CvGradA30 + CvR0 * A30 / BigR + Cvp0 * AR0 ) * (CvGradVi + CvR0 * v / BigR)

           
           Vms(var_A3)   =  CvGradA30 * CvGradVi
           
           Qvec(var_A3)  =  Qvec(var_A3) - TG_NUM_Eq*TG_NUM(var_A3) * CoefAdv * Vms(var_A3)  
!###################################################################################################
!#  equation 4   (R component momentum equation)                                                   #
!###################################################################################################

           Pvec_prev(var_uR) =   v * r0 * delta_g(mp,var_uR,ms,mt)   &
                               + v * delta_g(mp,var_r,ms,mt) * uR0

           if (primitive) then ! Lorentz and some viscous terms still partially integrated to avoid second derivatives

             Q_uR_primitive   =   v * ( - r0  * ( u0grad_uR0 - up0**2 / BigR ) - uR0 * divru - p0_R ) &  ! advection + pressure

                                + v * visco_T * ( - uR0 / BigR**2 - 2.d0 * up0_p / BigR**2            ) &

                                - visco_T * (v_R * UR0_R + v_Z * UR0_Z + v_p * uR0_p / BigR**2)        ! laplacian part viscous term

             
             ! adding the stabilization contribution
             ! --------------------------------------------------------------------------------------
             Vms(var_uR)    = r0 * ( CvGraduR0 * CvGradVi + VbGraduR0 * VbGradVi  ) 
           
             Q_uR_primitive =  Q_uR_primitive - TG_NUM_Eq*TG_NUM(var_uR) * CoefAdv * Vms(var_uR)  


             Qvec(var_uR) = Q_uR_primitive - ( visco_T + visco2 ) * divu * ( v_R + v / BigR )  &                   ! primitive + nonprimitive viscous part
                           
                          - ( BR0 * B0grad_vstar + Bp0**2 * v / BigR ) + ( v_R + v / BigR ) * ( 0.5d0 * BB2 )  !           + nonprimitive Lorentz force

           else

             Qvec(var_uR) =   r0 * ( uR0 * u0grad_vstar + up0**2 * v / BigR )                      &
                            -      ( BR0 * B0grad_vstar + Bp0**2 * v / BigR )                      &
                            + ( v_R + v / BigR ) * (r0 * T0 + 0.5d0 * BB2 - visco2 * divu)

             ! Viscous terms stored for the linearization later
             Qvisc_uR     = 2.d0 * up0_p * v / BigR**2 + 2.d0 * uR0_R * v_R   &
                              + (uR0_Z + uZ0_R ) * v_Z                        &
                              + (uR0_p + BigR * up0_R - up0) * v_p / BigR**2

             Qvec(var_uR) = Qvec(var_uR) - visco_T * Qvisc_uR
           endif

!###################################################################################################
!#  equation 5   (Z component momentum equation)                                                   #
!###################################################################################################

           Pvec_prev(var_uZ) =   v * r0 * delta_g(mp,var_uZ,ms,mt) &
                               + v * delta_g(mp,var_r,ms,mt) * uZ0

           if (primitive) then

             Q_uZ_primitive   =   v * ( - r0  * u0grad_uZ0 - uZ0 * divru - p0_Z )                     &

                                - visco_T * (v_R * UZ0_R + v_Z * UZ0_Z + v_p * uZ0_p / BigR**2)        ! laplacian part viscous term

             Qvec(var_uZ) = Q_uZ_primitive - ( visco_T + visco2 ) * divu * v_Z   &

                                - BZ0 * B0grad_vstar + v_Z * ( 0.5d0 * BB2 )
             
             ! adding the stabilization contribution
             ! --------------------------------------------------------------------------------------
             Vms(var_uZ)    =  r0 * ( CvGraduZ0  * CvGradVi + VbGraduZ0  * VbGradVi   )    
         
             Q_uZ_primitive  = Q_uZ_primitive  - TG_NUM_Eq*TG_NUM(var_uZ) * CoefAdv * Vms(var_uZ)  


           else
             Qvec(var_uZ) =   r0 * uZ0 * u0grad_vstar - BZ0 * B0grad_vstar                         &
                            + v_Z * (r0 * T0 + 0.5d0 * BB2 - visco2 * divu)

             ! Viscous terms stored for the linearization later
             Qvisc_uZ     =  2.d0 * uZ0_Z * v_Z                                    &
                                        + (uZ0_R + uR0_Z ) * v_R                   &
                                        + (uZ0_p + BigR * up0_Z ) * v_p / BigR**2

             Qvec(var_uZ) = Qvec(var_uZ) - visco_T * Qvisc_uZ
           endif
                          
!###################################################################################################
!#  equation 6 (Phi component momentum equation)                                                   #
!###################################################################################################


            Pvec_prev(var_up) =   v * r0 * delta_g(mp,var_up,ms,mt) &
                                + v * delta_g(mp,var_r,ms,mt) * up0   

           if (primitive) then

             Q_up_primitive   = v * ( - r0  * ( u0grad_up0 + uR0 * up0 / BigR ) - up0 * divru - p0_p / BigR )   & 

                              + v * visco_T * ( - up0 / BigR**2 + 2.d0 * uR0_p / BigR**2            )  &     

                              - visco_T * (v_R * UP0_R + v_Z * UP0_Z + v_p * uP0_p / BigR**2)        ! laplacian part viscous term

             ! adding the stabilization contribution
             ! --------------------------------------------------------------------------------------
             
             Vms(var_up)  = - r0 * ( CvGraduR0 - Cvp0 * up0 / BigR  ) * Cvp0 * v / BigR &
                  &         + r0 * ( CvGradup0 + Cvp0 * uR0 / BigR  ) * CvGradVi        &
                  &         - r0 * ( VbGraduR0 - Vbp0 * up0 / BigR  ) * Vbp0 * v / BigR &
                  &         + r0 * ( VbGradup0 + Vbp0 * uR0 / BigR  ) * VbGradVi
             
             Vms(var_up)      = r0 * ( CvGradup0 * CvGradVi  + VbGradup0 * VbGradVi )
             
             Q_up_primitive   =  Q_up_primitive   - TG_NUM_Eq * TG_NUM(var_up) * CoefAdv * Vms(var_up)
             

             if (parallel_projection) then ! A parallel projection includes R,Z, and phi components

               Pvec_prev(var_up) =   v * delta_g(mp,var_r,ms,mt) * ( BR0 * uR0 + BZ0 * uZ0 + Bp0 * up0 ) &
                                   + v * r0 * ( BR0 * delta_g(mp,var_uR,ms,mt) + BZ0 * delta_g(mp,var_uZ,ms,mt) + Bp0 * delta_g(mp,var_up,ms,mt) )

               Qvec(var_up)      = BR0 * Q_uR_primitive + BZ0 * Q_uZ_primitive + Bp0 * Q_up_primitive &  ! The Lorentz force dissapears under a parallel projection
                                   - ( visco_T + visco2 ) * divu * B0grad_vstar                            ! the non-primitive viscous terms, using div-B = 0

             else
               Qvec(var_up) = Q_up_primitive - ( visco_T + visco2 ) * divu * ( v_p / BigR )  & 
                              - ( Bp0 * B0grad_vstar - BR0 * Bp0 * v / BigR ) + ( v_p / BigR ) * ( 0.5d0 * BB2 )
             endif


           else
             Qvec(var_up) =  r0 * ( up0 * u0grad_vstar - uR0 * up0 * v / BigR )                    &
                            -      ( Bp0 * B0grad_vstar - BR0 * Bp0 * v / BigR )                   &
                            + ( v_p / BigR ) * (r0 * T0 + 0.5d0 * BB2 - visco2 * divu)


             ! Viscous terms stored for the linearization later
             Qvisc_up     =  2.d0 * ( uR0 + up0_p ) * v_p / BigR**2                   &
                             + ( 3.d0 * up0 / BigR -up0_R - uR0_p / BigR ) * v / BigR &
                             + ( -up0 / BigR + up0_R + uR0_p / BigR ) * v_R           &
                             + (BigR * up0_Z + uZ0_p) * v_Z / BigR

           endif
           
!###################################################################################################
!#  equation 7 (Density equation)                                                                  #
!###################################################################################################

           Pvec_prev(var_r) =   v * delta_g(mp,var_r,ms,mt)

           Qvec(var_r) = - D_prof * gradr0grad_vstar                        &
                         - (D_par-D_prof) * B0grad_vstar * B0grad_r0 / BB2  &
                         + v * particle_source(ms,mt)


           if (primitive) then
             Qvec(var_r) =  Qvec(var_r) - v * ( r0 * divu + uR0 * r0_R + uZ0 * r0_Z + up0 * r0_p / BigR )
           else
             Qvec(var_r) =  Qvec(var_r) + r0 * u0grad_vstar
           endif

           ! adding the stabilization contribution
           ! --------------------------------------------------------------------------------------
           Vms(var_r)  =  CvGradr0 * CvGradVi    +  VbGradr0 * VbGradVi

           Qvec(var_r) =  Qvec(var_r)  - TG_NUM_Eq*TG_NUM(var_r) * CoefAdv * Vms(var_r)

!###################################################################################################
!#  equation 8 (Pressure equation)                                                                 #
!###################################################################################################
! NB: Excludes resistive heating

           Pvec_prev(var_T) =   v * r0 * delta_g(mp,var_T,ms,mt)        &
                              + v * delta_g(mp,var_r,ms,mt) * T0


           if (primitive) then 
             Qvec(var_T) =   v * ( - r0 * u0grad_T0 - T0 * u0grad_r0 - gamma * p0 * divu ) + heat_source(ms,mt)  &
                           + (gamma-1.d0)*( - ZK_prof * gradT0grad_vstar - (ZKpar_T-ZK_prof) * B0grad_T0 * B0grad_vstar / BB2 )
           else
             Qvec(var_T) =    gamma * p0 * u0grad_vstar                                    &
             + (gamma-1.d0)*( - ZK_prof * gradT0grad_vstar - (ZKpar_T-ZK_prof) * B0grad_T0 * B0grad_vstar / BB2               & 
                              + v * ( r0 * u0grad_T0 + T0 * u0grad_r0 + heat_source(ms,mt) ) )
           endif


           if (viscores_heating ) then         ! UNTESTED, resistive heating is not (yet) incorporated
             ! Viscous terms stored for the linearization later
             ! The following represents the Newtonian part (T_N : grad u)/visco_T
             ! which can be obtained as the sum of Qvisc_uR, Qvisc_uZ, and Qvisc_up)
             ! with v replaced with uR0, uZ0, and up0 resp.

             Qvisc_T          =    2.d0 * up0_p * uR0 / BigR**2 + 2.d0 * uR0_R * uR0_R      &
                                 + (uR0_Z + uZ0_R ) * uR0_Z                                 &
                                 + (uR0_p + BigR * up0_R - up0) * uR0_p / BigR**2           &
                                 + 2.d0 * uZ0_Z * uZ0_Z                                     &
                                 + (uZ0_R + uR0_Z ) * uZ0_R                                 &
                                 + (uZ0_p + BigR * up0_Z ) * uZ0_p / BigR**2                &
                                 + 2.d0 * ( uR0 + up0_p ) * up0_p / BigR**2                 &
                                 + ( 3.d0 * up0 / BigR -up0_R - uR0_p / BigR ) * up0 / BigR &
                                 + ( -up0 / BigR + up0_R + uR0_p / BigR ) * up0_R           &
                                 + (BigR * up0_Z + uZ0_p) * up0_Z / BigR

             Qvec(var_T) = Qvec(var_T) +  v * (gamma-1.d0) * (visco2 * divu**2 + visco_T * Qvisc_T)

           endif

           ! adding the stabilization contribution
           ! --------------------------------------------------------------------------------------
           Vms(var_T)  =  r0 *( CvGradT0 * CvGradVi    +  VbGradT0 * VbGradVi )
           
           Qvec(var_T) = Qvec(var_T) - TG_NUM_Eq * TG_NUM(var_T) * CoefAdv * Vms(var_T)

!###################################################################################################
!#  equations end                                                                                  #
!###################################################################################################

           index_ij = n_tor*n_var*(n_order+1)*(i-1) + n_tor * n_var * (j-1) + im   ! index in the ELM matrix

           do ivar= 1,n_var
             ij = index_ij + (ivar-1)*n_tor

             RHS_ij(ivar) = tstep * Qvec(ivar)  + zeta * Pvec_prev(ivar)
             
             RHS(ij)      =  RHS(ij) + rhs_ij(ivar) * wst * BigR * xjac
           enddo



           do k=1,n_vertex_max

             do l=1,n_order+1

               do in = 1, n_tor

                 bf    = H(k,l,ms,mt)   * element%size(k,l) * HZ(in,mp)
                 bf_p  = H(k,l,ms,mt)   * element%size(k,l) * HZ_p(in,mp)
                 bf_s  = H_s(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 bf_t  = H_t(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)


                 bf_R = (   y_t(ms,mt) * bf_s - y_s(ms,mt) * bf_t ) / xjac
                 bf_Z = ( - x_t(ms,mt) * bf_s + x_s(ms,mt) * bf_t ) / xjac

                 uR    = bf    ;  uZ    = bf    ;  up    = bf
                 uR_R  = bf_R  ;  uZ_R  = bf_R  ;  up_R  = bf_R
                 uR_Z  = bf_Z  ;  uZ_Z  = bf_Z  ;  up_Z  = bf_Z
                 uR_p  = bf_p  ;  uZ_p  = bf_p  ;  up_p  = bf_p
                 uR_s  = bf_s  ;  uZ_s  = bf_s  ;  up_s  = bf_s
                 uR_t  = bf_t  ;  uZ_t  = bf_t  ;  up_t  = bf_t

                 AR    = bf    ;  AZ    = bf    ;  A3    = bf    ; T    = bf    ; r    = bf
                 AR_R  = bf_R  ;  AZ_R  = bf_R  ;  A3_R  = bf_R  ; T_R  = bf_R  ; r_R  = bf_R 
                 AR_Z  = bf_Z  ;  AZ_Z  = bf_Z  ;  A3_Z  = bf_Z  ; T_Z  = bf_Z  ; r_Z  = bf_Z
                 AR_p  = bf_p  ;  AZ_p  = bf_p  ;  A3_p  = bf_p  ; T_p  = bf_p  ; r_p  = bf_p
                 AR_s  = bf_s  ;  AZ_s  = bf_s  ;  A3_s  = bf_s  ; T_s  = bf_s  ; r_s  = bf_s 
                 AR_t  = bf_t  ;  AZ_t  = bf_t  ;  A3_t  = bf_t  ; T_t  = bf_T  ; r_t  = bf_T 


                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ! Some auxiliary quantities
                 BR0_AR =   0.d0        ; BR0_AZ = - AZ_p / BigR ; BR0_A3 =   A3_Z / BigR
                 BZ0_AR =   AR_p / BigR ; BZ0_AZ =   0.d0        ; BZ0_A3 = - A3_R / BigR
                 Bp0_AR = - AR_Z        ; Bp0_AZ =   AZ_R        ; Bp0_A3 =   0.d0

                 BB2_AR = 2.d0*(BR0_AR * BR0 + BZ0_AR * BZ0 + Bp0_AR * Bp0 )
                 BB2_AZ = 2.d0*(BR0_AZ * BR0 + BZ0_AZ * BZ0 + Bp0_AZ * Bp0 )
                 BB2_A3 = 2.d0*(BR0_A3 * BR0 + BZ0_A3 * BZ0 + Bp0_A3 * Bp0 )
       
                 B0grad_vstar_AR = BR0_AR * v_R  + BZ0_AR * v_Z  + Bp0_AR * v_p / BigR
                 B0grad_vstar_AZ = BR0_AZ * v_R  + BZ0_AZ * v_Z  + Bp0_AZ * v_p / BigR
                 B0grad_vstar_A3 = BR0_A3 * v_R  + BZ0_A3 * v_Z  + Bp0_A3 * v_p / BigR

                 B0grad_T0_AR    = BR0_AR * T0_R + BZ0_AR * T0_Z + Bp0_AR * T0_p / BigR
                 B0grad_T0_AZ    = BR0_AZ * T0_R + BZ0_AZ * T0_Z + Bp0_AZ * T0_p / BigR
                 B0grad_T0_A3    = BR0_A3 * T0_R + BZ0_A3 * T0_Z + Bp0_A3 * T0_p / BigR

                 B0grad_r0_AR    = BR0_AR * r0_R + BZ0_AR * r0_Z + Bp0_AR * r0_p / BigR
                 B0grad_r0_AZ    = BR0_AZ * r0_R + BZ0_AZ * r0_Z + Bp0_AZ * r0_p / BigR
                 B0grad_r0_A3    = BR0_A3 * r0_R + BZ0_A3 * r0_Z + Bp0_A3 * r0_p / BigR


                 ! The following quantities are used for respectively the variation of 
                 ! gradT0grad_vstar, u0grad_T0, and B0grad_T0 wrt T0 (also T0 -> r0)
                 gradbfgrad_vstar = bf_R * v_R + bf_Z * v_Z + (bf_p/BigR) * (v_p/BigR)
                 B0grad_bf    = BR0 * bf_R  + BZ0 * bf_Z  + Bp0 * bf_p  / BigR
                 u0grad_bf    = uR0 * bf_R  + uZ0 * bf_Z  + up0 * bf_p  / BigR

                 divu_uR = uR_R + uR / BigR ; divu_uZ = uZ_Z ; divu_up = up_p / BigR

                 divru_uR = r0 * divu_uR + uR * r0_R 
                 divru_uZ = r0 * divu_uZ + uZ * r0_Z
                 divru_up = r0 * divu_up + up * r0_p / BigR
                 divru_r  = r * divu + u0grad_bf

                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ! Variables used For the VMS Stabilization associated to the trial function
                 ! ------------------------------------------------------------------------
                 CvGradVj       = CvR0 * bf_R + CvZ0 * bf_Z + Cvp0 * bf_p / BigR
                 VbGradVj       = VbR0 * bf_R + VbZ0 * bf_Z + Vbp0 * bf_p / BigR

                 DiveRMVj        = bf_R + bf / BigR
                 DiveZMVj        = bf_z
                 DivePMVj        = bf_p / BigR

!###################################################################################################
!#  equation 1   (R component induction equation)                                                  #
!###################################################################################################

                Pjac(var_AR,var_AR) =   v * AR            

                Qjac(var_AR,var_AR) = - eta_T * ( v_p * BZ0_AR / BigR - v_Z * Bp0_AR )                &
                                      +   v * ( (uZ0 +  eta_Z) * Bp0_AR - (up0 + eta_p / BigR ) * BZ0_AR )                            

                Qjac(var_AR,var_AZ) = - eta_T * ( v_p * BZ0_AZ / BigR - v_Z * Bp0_AZ )                &
                                      +   v * ( (uZ0 +  eta_Z) * Bp0_AZ - (up0 + eta_p / BigR ) * BZ0_AZ )                                        
      
                Qjac(var_AR,var_A3) = - eta_T * ( v_p * BZ0_A3 / BigR - v_Z * Bp0_A3 )                &
                                      +   v * ( (uZ0 +  eta_Z) * Bp0_A3 - (up0 + eta_p / BigR ) * BZ0_A3 )                                            

                Qjac(var_AR,var_uR) =   0.d0

                Qjac(var_AR,var_uZ) =   uZ * Bp0 * v  

                Qjac(var_AR,var_up) = - up * BZ0 * v                      

                Qjac(var_AR,var_r ) =   0.d0

                Qjac(var_AR,var_T ) = - deta_dT * T * ( v_p * BZ0 / BigR - v_Z * Bp0 )               &
                                      + v * ( d2eta_d2T * T0_Z * T + deta_dT * T_Z ) * Bp0           &
                                      - v * ( d2eta_d2T * T0_p * T + deta_dT * T_p ) * BZ0 / BigR

                ! Stabilization
                ! -----------------------------------------------------------------------------
                QvmsAd(var_AR,var_AR) =  CvGradVj * CvGradVi + Cvp0 * AR  * Cvp0 * v / BigR2
                QvmsAd(var_AR,var_A3) =  Cvp0 * A3 / BigR2  * CvGradVi                            &
                                      + ( - CvGradVj + CvR0 * A3 / BigR ) * Cvp0 * v / BigR2
                
                Qjac(var_AR, :) = Qjac(var_AR,:) - TG_NUM(var_AR) * CoefAdv *  QvmsAd(var_AR, :)
                
!###################################################################################################
!#  equation 2   (Z component induction equation)                                                  #
!###################################################################################################

                Pjac(var_AZ,var_AZ) =   v * AZ
              
     
                Qjac(var_AZ,var_AR) = - eta_T * ( - v_p * BR0_AR / BigR + v_R * Bp0_AR )              &
                                      + v * ( (up0 +  eta_p / BigR) * BR0_AR - (uR0 + eta_R) * Bp0_AR )

                Qjac(var_AZ,var_AZ) = - eta_T * ( - v_p * BR0_AZ / BigR + v_R * Bp0_AZ )              &
                                      + v * ( (up0 +  eta_p / BigR) * BR0_AZ - (uR0 + eta_R) * Bp0_AZ )                                            
           
                Qjac(var_AZ,var_A3) = - eta_T * ( - v_p * BR0_A3 / BigR + v_R * Bp0_A3 )              &
                                      + v * ( (up0 +  eta_p / BigR) * BR0_A3 - (uR0 + eta_R) * Bp0_A3 )                                           

                Qjac(var_AZ,var_uR) = - uR * Bp0 * v

                Qjac(var_AZ,var_uZ) =   0.d0  

                Qjac(var_AZ,var_up) =   up * BR0 * v  
                
                Qjac(var_AZ,var_r) =    0.d0

                Qjac(var_AZ,var_T) = - deta_dT * T * ( - v_p * BR0 / BigR + v_R * Bp0 )               &
                                     + v * ( d2eta_d2T * T0_p * T + deta_dT * T_p ) * BR0 / BigR      &
                                     - v * ( d2eta_d2T * T0_R * T + deta_dT * T_R ) * Bp0
                
                ! Stabilization
                ! -----------------------------------------------------------------------------
                QvmsAd(var_AZ,var_AZ) = CvGradVj * CvGradVi

                Qjac(var_AZ,:) = Qjac(var_AZ,:) - TG_NUM(var_AZ) * CoefAdv * QvmsAd(var_AZ,:)
                
!###################################################################################################
!#  equation 3   (Phi component induction equation)                                                #
!###################################################################################################

                Pjac(var_A3,var_A3) =   v * A3
                   
    
                Qjac(var_A3,var_AR) = - eta_T * ( BigR * v_Z * BR0_AR - ( 2.d0 * v + BigR * v_R ) * BZ0_AR ) &
                                      + BigR * v * ( (uR0 +  eta_R ) * BZ0_AR - (uZ0 + eta_Z) * BR0_AR )                       

                Qjac(var_A3,var_AZ) = - eta_T * ( BigR * v_Z * BR0_AZ - ( 2.d0 * v + BigR * v_R ) * BZ0_AZ ) &
                                      + BigR * v * ( (uR0 +  eta_R ) * BZ0_AZ - (uZ0 + eta_Z) * BR0_AZ )                                           
           
                Qjac(var_A3,var_A3) = - eta_T * ( BigR * v_Z * BR0_A3 - ( 2.d0 * v + BigR * v_R ) * BZ0_A3 ) &
                                      + BigR * v * ( (uR0 +  eta_R ) * BZ0_A3 - (uZ0 + eta_Z) * BR0_A3 )                                             

                Qjac(var_A3,var_uR) =   uR * BZ0 * BigR * v

                Qjac(var_A3,var_uZ) = - uZ * BR0 * BigR * v

                Qjac(var_A3,var_up) =   0.d0 

                Qjac(var_A3,var_r ) =   0.d0

                Qjac(var_A3,var_T ) = - deta_dT * T * ( BigR * v_Z * BR0 - ( 2.d0 * v + BigR * v_R ) * BZ0 &
                                                        + v * current_source(ms,mt) )                      &
                                      + BigR * v * ( d2eta_d2T * T0_R * T + deta_dT * T_R ) * BZ0          &                                   
                                      - BigR * v * ( d2eta_d2T * T0_Z * T + deta_dT * T_Z ) * BR0                                       

                ! Stabilization
                ! -----------------------------------------------------------------------------
                QvmsAd(var_A3,var_A3) =   (Cvp0 * A3 / BigR2 ) * Cvp0 * v                                   &
                     &                - ( - CvGradVj + CvR0 * A3 / BigR ) * (CvGradVi + CvR0 * v / BigR)
                
                QvmsAd(var_A3,var_AR) =    CvGradVj * Cvp0 * v                                              &
                     &                 - (  Cvp0 * AR ) * (CvGradVi + CvR0 * v / BigR) 

                Qjac(var_A3,:)      = Qjac(var_A3,:) - TG_NUM(var_A3) * CoefAdv * QvmsAd(var_A3,:)

!###################################################################################################
!#  equation 4   (R component momentum equation)                                                   #
!###################################################################################################

                 Pjac(var_uR,var_uR)      =   v * r0 * uR
                 Pjac(var_uR,var_r)       =   v * r  * uR0

                 Qjac(var_uR,var_AR) = - ( BR0_AR * B0grad_vstar + BR0 * B0grad_vstar_AR + 2.d0 * Bp0 * Bp0_AR * v / BigR ) &
                                       + ( v_R + v / BigR ) * 0.5d0 * BB2_AR

                 Qjac(var_uR,var_AZ) = - ( BR0_AZ * B0grad_vstar + BR0 * B0grad_vstar_AZ + 2.d0 * Bp0 * Bp0_AZ * v / BigR ) &
                                       + ( v_R + v / BigR ) * 0.5d0 * BB2_AZ

                 Qjac(var_uR,var_A3) = - ( BR0_A3 * B0grad_vstar + BR0 * B0grad_vstar_A3 + 2.d0 * Bp0 * Bp0_A3 * v / BigR ) &
                                       + ( v_R + v / BigR ) * 0.5d0 * BB2_A3

                 if (primitive) then


                 Qjac(var_uR,var_uR) =   v * ( - r0 * ( uR * uR0_R + u0grad_bf ) - uR * divru - uR0 * divru_uR )       & 
 
                                        + v * visco_T * ( - uR / BigR**2 ) &

                                       - visco_T * (v_R * UR_R + v_Z * UR_Z + v_p * uR_p / BigR**2)  &      ! laplacian part viscous term

                                       - ( visco_T + visco2 ) * divu_uR * ( v_R + v / BigR )

                 Qjac(var_uR,var_uZ) =   v * ( - r0  * uZ * uR0_Z - uR0 * divru_uZ )       &
                                       - ( visco_T + visco2 ) * divu_uZ * ( v_R + v / BigR )

                 Qjac(var_uR,var_up) =   v * ( - r0  * ( up * uR0_p / BigR - 2.d0 * up0 * up / BigR ) - uR0 * divru_up ) &
                                       + v * visco_T * ( - 2.d0 * up_p / BigR**2 ) - ( visco + visco2 ) * divu_up * ( v_R + v / BigR )

                 Qjac(var_uR,var_r)  =   v * ( - r  * ( u0grad_uR0 - up0**2 / BigR ) - uR0 * divru_r - ( r * T0_R + r_R * T0 ) )

                 Qjac(var_uR,var_T)  =   v * ( - ( r0 * T_R + r0_R * T ) )  &

                                       + v * dvisco_dT * T * ( - uR0 / BigR**2 - 2.d0 * up0_p / BigR**2          ) &

                                       - dvisco_dT * T * (v_R * UR0_R + v_Z * UR0_Z + v_p * uR0_p / BigR**2)       &

                                       - dvisco_dT * T * divu * ( v_R + v / BigR ) 


                 else

                 Qjac(var_uR,var_uR) =   r0 * ( uR * u0grad_vstar + uR0 * uR * v_R )                                  &
                                       - ( v_R + v / BigR ) * visco2 * divu_uR                                        &
                                       - visco_T * ( 2.d0 * uR_R * v_R + uR_Z * v_Z + uR_p * v_p / BigR**2 )

                 Qjac(var_uR,var_uZ) =   r0 * uR0 * uZ * v_Z                                                          &
                                       - ( v_R + v / BigR ) * visco2 * divu_uZ                                        &
                                       - visco_T * uZ_R * v_Z
 
                 Qjac(var_uR,var_up) =   r0 * ( uR0 * up * v_p / BigR + 2.d0 * up0 * up * v / BigR )                  &
                                       - ( v_R + v / BigR ) * visco2 * divu_up                                        &
                                       - visco_T * ( 2.d0 * up_p * v / BigR**2 + (BigR * up_R - up) * v_p / BigR**2 )

                 Qjac(var_uR,var_r)  =   r * ( uR0 * u0grad_vstar + up0**2 * v / BigR )                               &
                                       + ( v_R + v / BigR ) * r * T0

                 Qjac(var_uR,var_T)  =   ( v_R + v / BigR ) * r0 * T - dvisco_dT * T * Qvisc_uR                      
              endif

              ! Stabilization
              ! -----------------------------------------------------------------------------
              QvmsAd(var_uR,var_r ) =  0.0  ! to be completed  !$BNK
              
              QvmsAd(var_uR,var_uR) =   ( r0 * CvGradVj + CvGradr0 * uR  ) *   CvGradVi         &
                   &                  + ( r0 * Cvp0 * uR / BigR          ) * ( Cvp0 * v /BigR ) &
                   &                  + ( r0 * CvGradVj + CvGradr0 * uR  ) *   CvGradVi         &
                   &                  + ( r0 * Cvp0 * uR / BigR          ) * ( Cvp0 * v /BigR )
              
              QvmsAd(var_uR,var_up) =   ( - r0 * Cvp0 * up / BigR          ) *   CvGradVi         &
                   &                  + (   r0 * CvGradVj + CvGradr0 * up  ) * ( Cvp0 * v /BigR ) &
                   &                  + ( - r0 * Vbp0 * up / BigR          ) *   VbGradVi         &
                   &                  + (   r0 * VbGradVj + VbGradr0 * up  ) * ( Vbp0 * v /BigR )
              ! acoustic  waves
              ! ------------------------------
              QvmsF(var_uR,var_r  )  =  VmsCoefF * ( u0grad_bf     + r * divu         ) * DiveRMVi
              
              QvmsF(var_uR,var_uR )  =  VmsCoefF * ( r0 * DiveRMVj + uR * r0_R        ) * DiveRMVi
              QvmsF(var_uR,var_uZ )  =  VmsCoefF * ( r0 * DiveZMVj + uZ * r0_Z        ) * DiveRMVi
              QvmsF(var_uR,var_up )  =  VmsCoefF * ( r0 * DivePMVj + up * r0_p/BigR   ) * DiveRMVi
              
              
              Qjac(var_uR, :) = Qjac(var_uR, :) - TG_NUM(var_uR) * ( CoefAdv * QvmsAd(var_uR, :) + QvmsF(var_uR, :)  )
              
!###################################################################################################
!#  equation 5   (Z component momentum equation)                                                   #
!###################################################################################################

                 Pjac(var_uZ,var_uZ)      =   v * r0 * uZ
                 Pjac(var_uZ,var_r)       =   v * r  * uZ0


                 Qjac(var_uZ,var_AR) = - BZ0_AR * B0grad_vstar - BZ0 * B0grad_vstar_AR                         &
                                       + v_Z * 0.5d0 * BB2_AR   

                 Qjac(var_uZ,var_AZ) = - BZ0_AZ * B0grad_vstar - BZ0 * B0grad_vstar_AZ                         &
                                       + v_Z * 0.5d0 * BB2_AZ

                 Qjac(var_uZ,var_A3) = - BZ0_A3 * B0grad_vstar - BZ0 * B0grad_vstar_A3                         &
                                       + v_Z * 0.5d0 * BB2_A3

                 if (primitive) then

                 Qjac(var_uZ,var_uR) =   v * ( - r0  * uR * uZ0_R - uZ0 * divru_uR )                              & 
                                       - ( visco_T + visco2 ) * divu_uR * v_Z 
  
                 Qjac(var_uZ,var_uZ) =   v * ( - r0  * ( uZ * uZ0_Z + u0grad_bf ) - uZ * divru - uZ0 * divru_uZ ) &

                                       - visco_T * (v_R * UZ_R + v_Z * UZ_Z + v_p * uZ_p / BigR**2)   &     ! laplacian part viscous term

                                       - ( visco_T + visco2 ) * divu_uZ * v_Z 

                 Qjac(var_uZ,var_up) =   v * ( - r0  * up * uZ0_p / BigR - uZ0 * divru_up )                    &
                                       - ( visco_T + visco2 ) * divu_up * v_Z 

                 Qjac(var_uZ,var_r)  =   v * ( - r  * u0grad_uZ0 - uZ0 * divru_r - ( r * T0_Z + r_Z * T0 ) )

                 Qjac(var_uZ,var_T)  =   v * ( - ( r0 * T_Z + r0_Z * T ) )                                     &

                                      - dvisco_dT * T * (v_R * UZ0_R + v_Z * UZ0_Z + v_p * uZ0_p / BigR**2)    &
                                       
                                      - dvisco_dT * T * divu * v_Z

                 else

                 Qjac(var_uZ,var_uR) =   r0 * uZ0 * uR * v_R                                                   &
                                       - v_Z* visco2 * divu_uR                                                 &
                                       - visco_T * uR_Z * v_R
  
                 Qjac(var_uZ,var_uZ) =   r0 * uZ * u0grad_vstar + r0 * uZ0 * uZ * v_Z                          &
                                       - v_Z * visco2 * divu_uZ                                                &
                                       - visco_T * ( 2.d0 * uZ_Z * v_Z + uZ_R * v_R + uZ_p * v_p / BigR**2 )
 
                 Qjac(var_uZ,var_up) =   r0 * uZ0 * up * v_p / BigR                                            &
                                       - v_Z * visco2 * divu_up                                                &
                                       - visco_T * BigR * up_Z * v_p / BigR**2

                 Qjac(var_uZ,var_r)  =   r * uZ0 * u0grad_vstar + v_Z * r * T0

                 Qjac(var_uZ,var_T)  =   v_Z * r0 * T - dvisco_dT * T * Qvisc_uZ   
                 endif                      

                 ! Stabilization
                 ! -----------------------------------------------------------------------------

                 QvmsAd(var_uZ,var_r ) =  0.0 ! to be completed  !$BNK
                 
                 QvmsAd(var_uZ,var_uZ) =   ( r0 * CvGradVj + CvGradr0 * uZ )  * CvGradVi  &
                      &                  + ( r0 * VbGradVj + VbGradr0 * uZ )  * VbGradVi
                 
                 ! acoustic  waves
                 ! ------------------------------
                 QvmsF(var_uZ,var_r  )  =  VmsCoefF * ( u0grad_bf     + r * divu         ) * DiveZMVi
                 
                 QvmsF(var_uZ,var_uR )  =  VmsCoefF * ( r0 * DiveRMVj + uR * r0_R        ) * DiveZMVi
                 QvmsF(var_uZ,var_uZ )  =  VmsCoefF * ( r0 * DiveZMVj + uZ * r0_Z        ) * DiveZMVi
                 QvmsF(var_uZ,var_up )  =  VmsCoefF * ( r0 * DivePMVj + up * r0_p/BigR   ) * DiveZMVi
                 
                 
                 Qjac(var_uZ, :) = Qjac(var_uZ, :) - TG_NUM(var_uZ) * ( CoefAdv * QvmsAd(var_uZ, :) + QvmsF(var_uZ, :) )


!###################################################################################################
!#  equation 6   (Phi component momentum equation)                                                 #
!###################################################################################################

                 if ( .not. parallel_projection ) then

                   Pjac(var_up,var_up)      =   v * r0 * up
                   Pjac(var_up,var_r)       =   v * r  * up0

                   Qjac(var_up,var_AR) = - ( Bp0_AR * B0grad_vstar + Bp0 * B0grad_vstar_AR )                    &
                                         + ( BR0_AR * Bp0 + BR0 * Bp0_AR ) * v / BigR                           &
                                         + ( v_p / BigR ) * 0.5d0 * BB2_AR

                   Qjac(var_up,var_AZ) = - ( Bp0_AZ * B0grad_vstar + Bp0 * B0grad_vstar_AZ )                    &
                                         + ( BR0_AZ * Bp0 + BR0 * Bp0_AZ ) * v / BigR                           &
                                         + ( v_p / BigR ) * 0.5d0 * BB2_AZ

                   Qjac(var_up,var_A3) = - ( Bp0_A3 * B0grad_vstar + Bp0 * B0grad_vstar_A3 )                    &
                                         + ( BR0_A3 * Bp0 + BR0 * Bp0_A3) * v / BigR                            &
                                         + ( v_p / BigR ) * 0.5d0 * BB2_A3
                 endif

                 if (primitive) then

                   Qjac(var_up,var_uR) =   v * ( - r0  * ( uR * up0_R + uR * up0 / BigR ) - up0 * divru_uR ) &
                                         + v * visco_T * (  2.d0 * uR_p / BigR**2  )

                   Qjac(var_up,var_uZ) =   v * ( - r0  * uZ * up0_Z - up0 * divru_uZ )
 
                   Qjac(var_up,var_up) =   v * ( - r0  * ( up * up0_p / BigR + u0grad_bf + uR0 * up / BigR )              &  
                                               - up * divru - up0 * divru_up )                                            &

!                                         + v * visco_T * (  up_RR + up_R / BigR + up_ZZ + up_pp / BigR**2 - up / BigR**2  )

                                         + v * visco_T * ( - up / BigR**2  ) &

                                         - visco_T * (v_R * UP_R + v_Z * UP_Z + v_p * uP_p / BigR**2)        ! laplacian part viscous term

                   Qjac(var_up,var_r)  =   v * ( - r  * ( u0grad_up0 + uR0 * up0 / BigR ) - up0 * divru_r         &
                                             - ( r * T0_p + r_p * T0 ) / BigR )

                   Qjac(var_up,var_T)  =   v * ( - ( r0 * T_p + r0_p * T ) / BigR ) &
                                         
                                         + v * dvisco_dT * T * ( - up0 / BigR**2 + 2.d0 * uR0_p / BigR**2      )  &     

                                         - dvisco_dT * T * (v_R * UP0_R + v_Z * UP0_Z + v_p * uP0_p / BigR**2)         ! laplacian part viscous term


                   ! Stabilization
                   ! -----------------------------------------------------------------------------
                   QvmsAd(var_up,var_r ) = 0.0  ! to be completed  !$BNK
                   
                   QvmsAd(var_up,var_uR) =  -  ( r0 * CvGradVj + CvGradr0 * uR ) *   Cvp0 * v / BigR              &
                        &                   +  ( r0 * Cvp0 * uR / BigR         ) * ( CvGradVi + CvR0 * v /BigR )  &
                        &                   -  ( r0 * VbGradVj + VbGradr0 * uR ) *   Vbp0 * v / BigR              &
                        &                   +  ( r0 * Vbp0 * uR / BigR         ) * ( VbGradVi + VbR0 * v /BigR )  
                   
                   QvmsAd(var_up,var_up) =  - ( - r0 * Cvp0 * up / BigR        ) *   Cvp0 * v / BigR              &
                        &                   + ( r0 * CvGradVj + up * CvGradr0  ) * ( CvGradVi + CvR0 * v /BigR )  &
                        &                   - ( - r0 * Vbp0 * up / BigR        ) *   Vbp0 * v / BigR              &
                        &                   + ( r0 * VbGradVj + up * VbGradr0  ) * ( VbGradVi + VbR0 * v /BigR )
                   
                   ! acoustic  waves
                   ! ------------------------------
                   QvmsF(var_up,var_r  )  =  VmsCoefF * ( u0grad_bf     + r * divu         ) * DivePMVi
                   
                   QvmsF(var_up,var_uR )  =  VmsCoefF * ( r0 * DiveRMVj + uR * r0_R        ) * DivePMVi
                   QvmsF(var_up,var_uZ )  =  VmsCoefF * ( r0 * DiveZMVj + uZ * r0_Z        ) * DivePMVi
                   QvmsF(var_up,var_up )  =  VmsCoefF * ( r0 * DivePMVj + up * r0_p/BigR   ) * DivePMVi
                   
                   
                   
                   Qjac(var_up, :) = Qjac(var_up, :) - TG_NUM(var_up) * ( CoefAdv * QvmsAd(var_up, :) + QvmsF(var_up, :) )

                   if (parallel_projection) then ! A parallel projection includes R,Z, and phi components

                     Pjac(var_up,var_uR)      =   v * r0 * BR0 * uR
                     Pjac(var_up,var_uZ)      =   v * r0 * BZ0 * uZ
                     Pjac(var_up,var_up)      =   v * r0 * Bp0 * up
                     Pjac(var_up,var_r)       =   v * r  * ( BR0 * uR0 + BZ0 * uZ0 + Bp0 * up0 )


                     Qjac(var_up,var_uR) =   BR0 * Qjac(var_uR,var_uR) + BZ0 * Qjac(var_uZ,var_uR) + Bp0 * Qjac(var_up,var_uR) &
                                           - ( visco_T + visco2 ) * divu_uR * B0grad_vstar                      ! nonprimitive viscous terms
                     Qjac(var_up,var_uZ) =   BR0 * Qjac(var_uR,var_uZ) + BZ0 * Qjac(var_uZ,var_uZ) + Bp0 * Qjac(var_up,var_uZ) &
                                           - ( visco_T + visco2 ) * divu_uZ * B0grad_vstar
                     Qjac(var_up,var_up) =   BR0 * Qjac(var_uR,var_up) + BZ0 * Qjac(var_uZ,var_up) + Bp0 * Qjac(var_up,var_up) &
                                           - ( visco_T + visco2 ) * divu_up * B0grad_vstar

                     Qjac(var_up,var_r)  =  BR0 * Qjac(var_uR,var_r)  + BZ0 * Qjac(var_uZ,var_r)  + Bp0 * Qjac(var_up,var_r)
                     Qjac(var_up,var_T)  =  BR0 * Qjac(var_uR,var_T)  + BZ0 * Qjac(var_uZ,var_T)  + Bp0 * Qjac(var_up,var_T) &
                                          - dvisco_dT * T * divu * B0grad_vstar 

                  else
                     
                     Qjac(var_up,var_uR) = Qjac(var_up,var_uR) - ( visco_T + visco2 ) * divu_uR * ( v_p / BigR ) ! nonprimitive viscous terms
                     Qjac(var_up,var_uZ) = Qjac(var_up,var_uZ) - ( visco_T + visco2 ) * divu_uZ * ( v_p / BigR ) 
                     Qjac(var_up,var_up) = Qjac(var_up,var_up) - ( visco_T + visco2 ) * divu_up * ( v_p / BigR )
                     
                  endif
                  
                else ! non Primitive_formulation

                   Qjac(var_up,var_uR) =   r0 * ( up0 * uR * v_R - uR * up0 * v / BigR )                        &
                                         - ( v_p / BigR ) * visco2 * divu_uR                                    &
                                         - visco_T * ( 2.d0 * uR * v_p / BigR**2                                &
                                                       - uR_p * v / BigR**2                                     &
                                                       + uR_p * v_R / BigR                                      &
                                                     )

                   Qjac(var_up,var_uZ) =   r0 * up0 * uZ * v_Z                                                  &
                                         - ( v_p / BigR ) * visco2 * divu_uZ                                    &
                                         - visco_T * uZ_p * v_Z / BigR
 
                   Qjac(var_up,var_up) =   r0 * ( up * u0grad_vstar + up0 * up * v_p / BigR - uR0 * up * v / BigR ) &
                                         - ( v_p / BigR ) * visco2 * divu_up                                        &
                                         - visco_T * ( 2.d0 * up_p * v_p / BigR**2                                  &
                                                       + ( 3.d0 * up / BigR - up_R ) * v / BigR                     &
                                                       + ( -up / BigR + up_R ) * v_R                                &
                                                       + (BigR * up_Z ) * v_Z / BigR                                &
                                                     )

                   Qjac(var_up,var_r)  =   r * ( up0 * u0grad_vstar - uR0 * up0 * v / BigR )                   &
                                         + ( v_p / BigR ) * r * T0

                   Qjac(var_up,var_T)  =   ( v_p / BigR ) * r0 * T - dvisco_dT * T * Qvisc_up

                   ! Stabilization
                   ! -----------------------------------------------------------------------------
                   ! Qjac(var_up,var_up) =  Qjac(var_up,var_up)  &
                   !     &   - TG_NUM(var_up) * r0 * ( CvGradVj * CvGradVi + VbGradVj * VbGradVi )
                   
                   write(6,*) " Stabilization to be done for primitive without parallel !!!!!!! $BNK "
                   
                endif

!###################################################################################################
!#  equation 7   (Density equation)                                                                #
!###################################################################################################


                 Pjac(var_r,var_r)      =   v * r
        
                 
                 Qjac(var_r,var_AR) = -((D_par - D_prof) / BB2) *                                    & 
                                        ( B0grad_vstar_AR * B0grad_r0  + B0grad_vstar * B0grad_r0_AR &
                                          - BB2_AR * B0grad_r0 * B0grad_vstar / BB2 )

                 Qjac(var_r,var_AZ) = -((D_par - D_prof) / BB2) *                                    & 
                                        ( B0grad_vstar_AZ * B0grad_r0  + B0grad_vstar * B0grad_r0_AZ &
                                          - BB2_AZ * B0grad_r0 * B0grad_vstar / BB2 )

                 ! psi-dependence D_prof and mass source not yet properly taken into account
                 Qjac(var_r,var_A3) = -((D_par - D_prof) / BB2) *                                    & 
                                        ( B0grad_vstar_A3 * B0grad_r0  + B0grad_vstar * B0grad_r0_A3 &
                                          - BB2_A3 * B0grad_r0 * B0grad_vstar / BB2 )

    
                 if (primitive) then

                   Qjac(var_r,var_uR) =   - v * ( r0 * divu_uR + uR * r0_R )
                 
                   Qjac(var_r,var_uZ) =   - v * ( r0 * divu_uZ + uZ * r0_Z )

                   Qjac(var_r,var_up) =   - v * ( r0 * divu_up + up * r0_p / BigR )

                   Qjac(var_r,var_r)  = - v * ( r * divu + uR0 * r_R + uZ0 * r_Z + up0 * r_p / BigR )  &
                                        - D_prof * gradbfgrad_vstar                                    &
                                        - (D_par-D_prof) * B0grad_vstar * B0grad_bf / BB2

                 else

                   Qjac(var_r,var_uR) =   r0 * uR * v_R
                 
                   Qjac(var_r,var_uZ) =   r0 * uZ * v_Z

                   Qjac(var_r,var_up) =   r0 * up * v_p / BigR

                   Qjac(var_r,var_r)  =   r * u0grad_vstar                                 &
                                        - D_prof * gradbfgrad_vstar                        &
                                        - (D_par-D_prof) * B0grad_vstar * B0grad_bf / BB2

                 endif

                 Qjac(var_r,var_T)  =   0.

                 
                 ! Stabilization
                 ! -----------------------------------------------------------------------------
                 ! ---------------
                 QvmsAd(var_r,var_r) =  CvGradVj * CvGradVi  + VbGradVi * VbGradVj
                 
                 Qjac(var_r, :)      =  Qjac(var_r, :)  - TG_NUM(var_r) *  CoefAdv * QvmsAd(var_r, :)

!###################################################################################################
!#  equation 8   (Temperature  equation)                                                           #
!###################################################################################################


                Pjac(var_T,var_T)      =   v * r0 * T
                Pjac(var_T,var_r)      =   v * r  * T0

                
                Qjac(var_T,var_AR) = - (gamma - 1.d0 ) * ((ZKpar_T - ZK_prof) / BB2) *               & 
                                       ( B0grad_vstar_AR * B0grad_T0  + B0grad_vstar * B0grad_T0_AR &
                                     -   BB2_AR * B0grad_T0 * B0grad_vstar / BB2 )
                
                Qjac(var_T,var_AZ) = - (gamma - 1.d0 ) * ((ZKpar_T - ZK_prof) / BB2) *               & 
                                       ( B0grad_vstar_AZ * B0grad_T0  + B0grad_vstar * B0grad_T0_AZ &
                                     -   BB2_AZ * B0grad_T0 * B0grad_vstar / BB2 )   
                
                Qjac(var_T,var_A3) = - (gamma - 1.d0 ) * ((ZKpar_T - ZK_prof) / BB2) *               & 
                                       ( B0grad_vstar_A3 * B0grad_T0  + B0grad_vstar * B0grad_T0_A3 &
                                     -   BB2_A3 * B0grad_T0 * B0grad_vstar / BB2 )



                if (primitive) then

                  Qjac(var_T,var_uR) = v * ( - r0 * uR * T0_R - T0 * uR * r0_R - gamma * p0 * divu_uR )

                  Qjac(var_T,var_uZ) = v * ( - r0 * uZ * T0_Z - T0 * uZ * r0_Z - gamma * p0 * divu_uZ )

                  Qjac(var_T,var_up) = v * ( - r0 * up * T0_p / BigR - T0 * up * r0_p / BigR - gamma * p0 * divu_up)

                  Qjac(var_T,var_r)  = v * ( -  r * u0grad_T0 - T0 * u0grad_bf - gamma * r * T0 * divu ) 

                  Qjac(var_T,var_T)  = v * ( -  r0 * u0grad_bf - T * u0grad_r0 - gamma * r0 * T * divu )  & 
                                       + (gamma-1.d0)*( - ZK_prof * gradbfgrad_vstar - (ZKpar_T-ZK_prof) * B0grad_bf * B0grad_vstar / BB2 ) &
                                       + (gamma-1.d0)*(                              - (dZKpar_dT * T)   * B0grad_T0 * B0grad_vstar / BB2 )
 

                else

                  Qjac(var_T,var_uR) =    gamma * p0 * uR * v_R                               &
                                        + (gamma-1.d0) * v * ( r0 * uR * T0_R + T0 * uR * r0_R )

                  Qjac(var_T,var_uZ) =    gamma * p0 * uZ * v_Z                               &
                                        + (gamma-1.d0) * v * (    r0 * uZ * T0_Z + T0 * uZ * r0_Z  )

                  Qjac(var_T,var_up) =    gamma * p0 * up * v_p / BigR                        &
                                        + (gamma-1.d0) * v * ( r0 * up * T0_p / BigR + T0 * up * r0_p / BigR )

                  Qjac(var_T,var_r)  =  gamma * r * T0 * u0grad_vstar                       &
                                       + (gamma-1.d0) * v * ( r * u0grad_T0 + T0 * u0grad_bf )  

                  Qjac(var_T,var_T)  =    gamma * r0 * T * u0grad_vstar                       & 
                                        + (gamma-1.d0)*( - ZK_prof * gradbfgrad_vstar         & 
                                                         - (ZKpar_T-ZK_prof) * B0grad_bf * B0grad_vstar / BB2 & 
                                                         -  dZKpar_dT * T    * B0grad_T0 * B0grad_vstar / BB2 & 
                                        + v * ( r0 * u0grad_bf + T * u0grad_r0 + dvisco_dT * T * Qvisc_T ) )
                endif


                if (viscores_heating ) then

                  Qjac(var_T,var_uR) = Qjac(var_T,var_uR) + v * (gamma - 1.d0) * ( visco_T *(        &
                                            2.d0 * up0_p * uR / BigR**2 + 4.d0 * uR_R * uR0_R        &
                                          + (2.d0 * uR0_Z + uZ0_R ) * uR_Z                           &
                                          + (2.d0 * uR0_p + BigR * up0_R - up0) * uR_p / BigR**2     &
                                          + uR_Z * uZ0_R  + 2.d0 * uR * up0_p / BigR**2              &
                                          - uR_p * up0 / BigR**2 + uR_p  * up0_R / BigR              &
                                                                                            )        &
                                          + visco2 * 2.d0 * divu * divu_uR       )


                  Qjac(var_T,var_uZ) = Qjac(var_T,var_uZ) + v * (gamma - 1.d0) * ( visco_T *(        &
                                          + uZ_R * uR0_Z                                             &
                                          + 4.d0 * uZ_Z * uZ0_Z                                      &
                                          + (2.d0 * uZ0_R + uR0_Z ) * uZ_R                           &
                                          + (2.d0 * uZ0_p + BigR * up0_Z ) * uZ_p / BigR**2          &
                                          + uZ_p * up0_Z / BigR                                      &           
                                                                                          )          &
                                          + visco2 * 2.d0 * divu * divu_uZ        )

                  Qjac(var_T,var_up) = Qjac(var_T,var_up) + v * (gamma - 1.d0) * ( visco_T *(        &
                                            2.d0 * up_p * uR0   / BigR**2                            &
                                         + (up_R - up) * uR0_p / BigR**2                             &
                                         + BigR * up_Z * uZ0_p / BigR**2                             &
                                         + 2.d0 * ( uR0 + 2.d0 * up0_p ) * up_p / BigR**2            &
                                         + ( 6.d0 * up0 / BigR - up0_R - uR0_p / BigR ) * up / BigR - up_R * up0 / BigR  &
                                         + ( -up0 / BigR + 2.d0 * up0_R + uR0_p / BigR ) * up_R - up * up0_R / BigR      &
                                         + (2.d0 * BigR * up0_Z + uZ0_p) * up_Z / BigR                                   &
                                                                                          )                              &
                                        + visco2 * 2.d0 * divu * divu_up          )

               endif
               
               ! Stabilization
               ! -----------------------------------------------------------------------------
               QvmsAd(var_T,var_r) =   ( T0 * CvGradVj  + r * CvGradT0 ) * CvGradVi   &
                    &                + ( T0 * VbGradVj  + r * VbGradT0 ) * VbGradVi   
               
               
               QvmsAd(var_T,var_T) =   ( r0 * CvGradVj  + T * CvGradr0 ) * CvGradVi   &
                    &                + ( r0 * VbGradVj  + T * VbGradr0 ) * VbGradVi   
               
                   
               ! acoustic  waves
               ! ------------------------------
               QvmsF(var_T,var_r  )  =  VmsCoefF_T * ( T0*gradbfgrad_vstar + gradT0grad_vstar * r )
               QvmsF(var_T,var_T  )  =  VmsCoefF_T * ( r0*gradbfgrad_vstar + gradT0grad_vstar * r )

               
               Qjac(var_T, :) =  Qjac(var_T, :) - TG_NUM(var_T) *( CoefAdv * QvmsAd(var_T, :) + QvmsF(var_T, :) )


!###################################################################################################
!#  equations end                                                                                  #
!###################################################################################################


                 index_kl = n_tor*n_var*(n_order+1)*(k-1) + n_tor * n_var * (l-1) + in   ! index in the ELM matrix 


                 do ivar= 1,n_var
                   do kvar= 1,n_var
                     ij = index_ij + (ivar-1)*n_tor
                     kl = index_kl + (kvar-1)*n_tor

                     amat(ivar,kvar) = (1.d0+zeta)*Pjac(ivar,kvar) - tstep * theta * Qjac(ivar,kvar)

                     ELM(ij,kl) =  ELM(ij,kl) +  wst * amat(ivar,kvar)*BigR * xjac
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
enddo




return
end subroutine element_matrix
end module mod_elt_matrix
