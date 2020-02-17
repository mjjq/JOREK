subroutine F_profile(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,&
                     F_prof,dF_dpsi,dF_dz, dF_dpsi2,dF_dz2,dF_dpsi_dz, &
                     FFprime_prof,dFF_dpsi,dFF_dz, dFF_dpsi2,dFF_dz2,dFF_dpsi_dz)
!-----------------------------------------------------------------------
!
! EvdP new routine to calculate F(psi) explicitly to determine eq B_tor profile
!
!-----------------------------------------------------------------------
use phys_module

implicit none

logical :: xpoint2
integer :: xcase2
real*8  :: Fconst, profF, profF1, F_prof, dF_dpsi, dF_dz, dF_dpsi2, dF_dz2, dF_dpsi_dz
real*8  :: FFprime_prof, profFF, profFF1, FF_prof, dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz
real*8  :: dF_dpsi3, psi_edge, F_edge, sqrt_edge
real*8  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd,  psi_n, psi_barrier, sig_F, sigz, delta_psi
real*8  :: atn, datn, d2atn, d3atn
real*8  :: atn_z,   datn_z,   d2atn_z
real*8  :: atn_z_u, datn_z_u, d2atn_z_u
real*8  :: d_0, d_pert, d2_pert, d3_pert
real*8  :: tanh2, cosh3, tanh2_u, cosh3_u
real*8  :: alfa,profFFp, dprofFFp_dpsi, prof_bnd
real*8  :: prof0,     dprof0_dpsi,     dprof0_dpsi2,     dprof0_dpsi3
real*8  :: poly,      dpoly_dpsi,      dpoly_dpsi2,      dpoly_dpsi3
real*8  :: pert,      dpert_dpsi,      dpert_dpsi2,      dpert_dpsi3
real*8  :: sqrt_term, dsqrt_term_dpsi, dsqrt_term_dpsi2, dsqrt_term_dpsi3

! --- There are some rules when using FF_coefs with the F-profile in Full-MHD
if (FF_1 .ne. 0.d0) then
  write(*,*)'Full-MHD Warning!!! The F-profile does not like it if FF_1 is not zero !!!'
  write(*,*)'                    if you don,t respect this rule, we cannot guarantee that your F-profile and FFprime will be consistent!'
endif
if (FF_coef(4) .gt. 0.01) then
  write(*,*)'Full-MHD Warning!!! The tanh at FF_coef(5) with width FF_coef(4) is supposed to be a cut-off at the plasma edge !!!'
  write(*,*)'                    ie. FF_coef(5) should be the edge of your plasma, and FF_coef(4) should be very small...'
  write(*,*)'                    if you don,t respect this rule, we cannot guarantee that your F-profile and FFprime will be consistent!'
endif

! --- the cutoff of the FFprime at the edge is traditionally the tanh at FF_coef(5), not at psi_n=1.0
psi_edge  = FF_coef(5)

! --- psi_norm
psi_n = (psi - psi_axis)/(psi_bnd - psi_axis)
if (xpoint2) then
  if ((psi_n .lt. 1.d0) .and. (Z .lt. Z_xpoint(1)) .and. (xcase2 .ne. 2)) then
    psi_n = 2.d0 - psi_n
  endif
  if ((psi_n .lt. 1.d0) .and. (Z .gt. Z_xpoint(2)) .and. (xcase2 .ne. 1)) then
    psi_n = 2.d0 - psi_n
  endif
endif
delta_psi = (psi_bnd - psi_axis)  ! abs(psi_bnd - psi_axis)

! --- Polynomial part
poly        = ( psi_n  + FF_coef(1)/2.d0 * psi_n**2 + FF_coef(2)/3.d0 * psi_n**3 + FF_coef(3)/4.d0 * psi_n**4 ) * delta_psi
dpoly_dpsi  = ( 1.d0   + FF_coef(1)      * psi_n    + FF_coef(2)      * psi_n**2 + FF_coef(3)      * psi_n**3 )
dpoly_dpsi2 = (          FF_coef(1)                 + FF_coef(2)*2.0  * psi_n    + FF_coef(3)*3.0  * psi_n**2 ) / delta_psi
dpoly_dpsi3 = (                                     + FF_coef(2)*2.0             + FF_coef(3)*6.0  * psi_n    ) / delta_psi**2

! --- Perturbation part
pert        = + FF_coef(6) * tanh((psi_n - FF_coef(7))/FF_coef(8))    / 2.d0
dpert_dpsi  = + FF_coef(6) / cosh((psi_n - FF_coef(7))/FF_coef(8))**2 / (2.d0 * FF_coef(8)) / delta_psi
dpert_dpsi2 = - FF_coef(6) / cosh((psi_n - FF_coef(7))/FF_coef(8))**3 / FF_coef(8)**2       / delta_psi**2 &
                           * sinh((psi_n - FF_coef(7))/FF_coef(8))
dpert_dpsi2 = + FF_coef(6) / cosh((psi_n - FF_coef(7))/FF_coef(8))**4 / FF_coef(8)**3 * 3.0 / delta_psi**3 &
                           * sinh((psi_n - FF_coef(7))/FF_coef(8))**2                             &
              - FF_coef(6) / cosh((psi_n - FF_coef(7))/FF_coef(8))**2 / FF_coef(8)**3       / delta_psi**3
dpert_dpsi3 = - FF_coef(6) / cosh((psi_n - FF_coef(7))/FF_coef(8))**5 / FF_coef(8)**4 *12.0 / delta_psi**4 &
                           * sinh((psi_n - FF_coef(7))/FF_coef(8))**3                             &
              + FF_coef(6) / cosh((psi_n - FF_coef(7))/FF_coef(8))**3 / FF_coef(8)**4 * 3.0 / delta_psi**4 &
                     * 2.0 * sinh((psi_n - FF_coef(7))/FF_coef(8))                                &
              + FF_coef(6) / cosh((psi_n - FF_coef(7))/FF_coef(8))**3 / FF_coef(8)**4 * 2.0 / delta_psi**4 &
                           * sinh((psi_n - FF_coef(7))/FF_coef(8))

! --- Value of F at the edge
sqrt_edge =   2.0 * (FF_0 - FF_1) * (psi_edge + FF_coef(1)/2.d0 * psi_edge**2 + FF_coef(2)/3.d0 * psi_edge**3 + FF_coef(3)/4.d0 * psi_edge**4) * delta_psi &
            + 2.0 * FF_coef(6) * tanh((psi_edge - FF_coef(7))/FF_coef(8)) / 2.d0 &
            + F0**2
F_edge    = sqrt_edge**0.5 

! --- sqrt part
sqrt_term        = 2.0 * (FF_0 - FF_1) * poly        + 2.0 * pert            + F0**2 
dsqrt_term_dpsi  = 2.0 * (FF_0 - FF_1) * dpoly_dpsi  + 2.0 * dpert_dpsi 
dsqrt_term_dpsi2 = 2.0 * (FF_0 - FF_1) * dpoly_dpsi2 + 2.0 * dpert_dpsi2
dsqrt_term_dpsi3 = 2.0 * (FF_0 - FF_1) * dpoly_dpsi3 + 2.0 * dpert_dpsi3

! --- Profile and derivatives
prof0           =          sqrt_term**(+0.5)
dprof0_dpsi     = + 0.5  * sqrt_term**(-0.5) * dsqrt_term_dpsi
dprof0_dpsi2    = - 0.25 * sqrt_term**(-1.5) * dsqrt_term_dpsi**2 + 0.5  * sqrt_term**(-0.5) * dsqrt_term_dpsi2
dprof0_dpsi3    = + 0.375* sqrt_term**(-2.5) * dsqrt_term_dpsi**3 &
                  - 0.25 * sqrt_term**(-1.5) * dsqrt_term_dpsi * 2.0 * dsqrt_term_dpsi2 &
                  - 0.25 * sqrt_term**(-1.5) * dsqrt_term_dpsi       * dsqrt_term_dpsi2 &
                  + 0.5  * sqrt_term**(-0.5) * dsqrt_term_dpsi3




! --- Save F-profile
F_prof    = prof0
dF_dpsi   = dprof0_dpsi
dF_dpsi2  = dprof0_dpsi2
dF_dpsi3  = dprof0_dpsi3

! --- Save FF'-profile
FFprime_prof = F_prof * dF_dpsi
dFF_dpsi     = dF_dpsi * dF_dpsi + F_prof * dF_dpsi2
dFF_dpsi2    = 3.0 * dF_dpsi * dF_dpsi2 + F_prof * dF_dpsi3
dFF_dz       = 0.d0
dFF_dz2      = 0.d0
dFF_dpsi_dz  = 0.d0

! --- Cut-off at the plasma edge (because the tanh part if not integrable...)
if (psi_n .gt. psi_edge) then
  F_prof   = F_edge
  dF_dpsi  = 0
  dF_dpsi2 = 0
  dF_dpsi3 = 0

  FFprime_prof = 0.d0
  dFF_dpsi     = 0.d0
  dFF_dpsi2    = 0.d0
  dFF_dz       = 0.d0
  dFF_dz2      = 0.d0
  dFF_dpsi_dz  = 0.d0
endif

!--------------------------------

!write(*,'(A,5e12.4)') ' F, FFprime new, old, diff : ',F0,F_profile,FFprime_profile, profFFp, FFprime_profile - profFFp
!write(*,'(A,5e12.4)') ' F, dF_dpsi, dF_dpsi2 , FFprime: ',F_profile,dF_dpsi, dF_dpsi2, FFprime_profile
!write(*,'(A,5e12.4)') ' dFFp_dpsi new, old, diff : ',dFF_dpsi, dprofFFp_dpsi, dFF_dpsi - dprofFFp_dpsi


return
end
