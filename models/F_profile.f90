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
real*8  :: prof0, dprof0_dpsi, dprof0_dpsi2
real*8  :: prof1, dprof1_dpsi, dprof1_dpsi2
real*8  :: prof2, dprof2_dpsi, dprof2_dpsi2
real*8  :: Fconst, profF, profF1, F_prof, dF_dpsi, dF_dz, dF_dpsi2, dF_dz2, dF_dpsi_dz
real*8  :: FFprime_prof, profFF, profFF1, FF_prof, dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz
real*8  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd,  psi_n, psi_barrier, sig_F, sigz
real*8  :: atn, datn, d2atn, d3atn
real*8  :: atn_z,   datn_z,   d2atn_z
real*8  :: atn_z_u, datn_z_u, d2atn_z_u, factor
real*8  :: d_0, d_pert, d2_pert, d3_pert
real*8  :: tanh2, cosh3, tanh2_u, cosh3_u
real*8  :: alfa,profFFp, dprofFFp_dpsi, prof_bnd

psi_n     = (psi - psi_axis) / (psi_bnd - psi_axis)

psi_n = max( min(psi_n, 2.), 0. )

factor = 1.d0

if (xpoint2) then
  if ((Z .lt. Z_xpoint(1)) .and. (psi_n .lt. 1.d0)) then
    psi_n  = 2.d0 - psi_n
    factor = -1.d0
  endif
endif

prof_bnd = (1.d0 + FF_coef(1)/2.d0 + FF_coef(2)/3.d0 )


prof0           = sqrt(2.d0 * (psi_bnd - psi_axis) * (-1.d0) * (FF_0 - FF_1) &
                       
                       * (psi_n  + FF_coef(1)/2.d0 * psi_n**2 + FF_coef(2)/3.d0 * psi_n**3 - prof_bnd ) &

                       + F0**2)

dprof0_dpsi     = 1.d0/prof0 * (-1.d0) * (FF_0 - FF_1) &

                      * (1.d0   + FF_coef(1) * psi_n + FF_coef(2) * psi_n**2 )
 
dprof0_dpsi2    = -1.d0 * dprof0_dpsi/prof0**2 * (-1.d0) * (FF_0 - FF_1) & 
                 
                     * (1.d0   + FF_coef(1) * psi_n + FF_coef(2) * psi_n**2 ) & 
                  
                  +1.d0/prof0 * 1.d0/(psi_bnd - psi_axis) * (-1.d0) * (FF_0 - FF_1) & 

                     * (         FF_coef(1) + FF_coef(2) * 2.d0 * psi_n )

!-- former FF'prof:
profFFp        = -(FF_0 - FF_1) * ( 1.d0 + FF_coef(1) * psi_n +        FF_coef(2) * psi_n**2 +        FF_coef(3) * psi_n**3)
dprofFFp_dpsi  = -(FF_0 - FF_1) * (        FF_coef(1)         + 2.d0 * FF_coef(2) * psi_n    + 3.d0 * FF_coef(3) * psi_n**2)/(psi_bnd - psi_axis)


!-------------------------- choose profile:

F_prof = prof0
dF_dpsi   = dprof0_dpsi
dF_dpsi2  = dprof0_dpsi2

FFprime_prof = F_prof * dF_dpsi
dFF_dpsi     = dF_dpsi * dF_dpsi + F_prof * dF_dpsi2
dFF_dz       = 0.d0
dFF_dz2      = 0.d0
dFF_dpsi_dz  = 0.d0

if (psi_n .gt. 1.d0) then
  F_prof   = F0
  dF_dpsi  = 0
  dF_dpsi2 = 0

  FFprime_prof = 0.d0
  dFF_dpsi     = 0.d0
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
