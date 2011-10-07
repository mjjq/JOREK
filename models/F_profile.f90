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
real*8  :: alfa,profFFp, dprofFFp_dpsi

sig_F       = FF_coef(4)
psi_barrier = FF_coef(5)

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

factor = 1.d0

!-------------------------- 

d_pert  = + FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**2 / (2.d0 * FF_coef(8)) / (psi_bnd - psi_axis)
d2_pert = - FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**2 / (FF_coef(8)**2)  &
        * tanh((psi_n - FF_coef(7))/FF_coef(8)) / (psi_bnd - psi_axis)**2
d3_pert = + FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**4 / (FF_coef(8)**3)  &
        * (-2.d0 + cosh(2.d0*(psi_n-FF_coef(7))/FF_coef(8)) ) / (psi_bnd - psi_axis)**3


!--------------------------- 

prof0           = sqrt( abs(2.d0 * (psi_bnd - psi_axis) *  & 
                       -1.d0 * (FF_0 - FF_1) * &
                       (psi_n  + FF_coef(1)/2.d0 * psi_n**2 & 
                               + FF_coef(2)/3.d0 * psi_n**3 ) &
                       + F0**2))

dprof0_dpsi     = 1.d0/prof0 * &
                       -1.d0 * (FF_0 - FF_1) * & 
                       (1.d0   + FF_coef(1) * psi_n  &
                               + FF_coef(2) * psi_n**2 )
 
dprof0_dpsi2    = -1.d0 * dprof0_dpsi/prof0**2 *  & 
                       -1.d0 * (FF_0 - FF_1) * & 
                       (1.d0   + FF_coef(1) * psi_n  & 
                               + FF_coef(2) * psi_n**2 ) & 
                  +1.d0/prof0 * 1.d0/(psi_bnd - psi_axis) * & 
                       -1.d0 * (FF_0 - FF_1) * & 
                       (         FF_coef(1) + FF_coef(2) * 2.d0 * psi_n )

!-- former FF'prof:
profFFp        = -(FF_0 - FF_1) * ( 1.d0 + FF_coef(1) * psi_n +        FF_coef(2) * psi_n**2 +        FF_coef(3) * psi_n**3)
dprofFFp_dpsi  = -(FF_0 - FF_1) * (        FF_coef(1)         + 2.d0 * FF_coef(2) * psi_n    + 3.d0 * FF_coef(3) * psi_n**2)/(psi_bnd - psi_axis)

!--------------------------- arctangent

atn             = (0.5d0 - 0.5d0*tanh((psi_n - psi_barrier)/sig_F))
datn            = - 1.d0/cosh((psi_n - psi_barrier)/sig_F)**2 / (2.d0 * sig_F) / (psi_bnd - psi_axis)
d2atn           =   1.d0/cosh((psi_n - psi_barrier)/sig_F)**2 / (sig_F**2)  &
                   * tanh((psi_n - psi_barrier)/sig_F) / (psi_bnd - psi_axis)**2
d3atn           = - 1.d0/cosh((psi_n - psi_barrier)/sig_F)**4 / (sig_F**3)  &
                   * (-2.d0 + cosh(2.d0*(psi_n-psi_barrier)/sig_F) ) / (psi_bnd - psi_axis)**3

prof1        = prof0        + d_pert
dprof1_dpsi  = dprof0_dpsi  + d2_pert
dprof1_dpsi2 = dprof0_dpsi2 + d3_pert

prof2        = (prof0 + d_pert) * atn
dprof2_dpsi  = (dprof0_dpsi  + d2_pert) * atn +        (prof0 + d_pert) * datn
dprof2_dpsi2 = (dprof0_dpsi2 + d3_pert) * atn + 2.d0 * (dprof0_dpsi + d2_pert) * datn + (prof0 + d_pert) * d2atn

!-------------------------- choose profile:

F_prof = prof2
dF_dpsi   = dprof2_dpsi
dF_dpsi2  = dprof2_dpsi2

FFprime_prof = F_prof * dF_dpsi
dFF_dpsi        = dF_dpsi * dF_dpsi + F_prof * dF_dpsi2
dFF_dz          = 0.d0
dFF_dz2         = 0.d0
dFF_dpsi_dz     = 0.d0

if (xpoint2) then
  sigz    = 0.1d0

  tanh2   = tanh((Z_xpoint(1)-Z)/sigz)
  cosh3   = cosh((Z_xpoint(1)-Z)/sigz)
  tanh2_u = tanh((Z-Z_xpoint(2))/sigz)
  cosh3_u = cosh((Z-Z_xpoint(2))/sigz)
    
  atn_z 	   = (0.5d0 - 0.5d0*tanh2)
  datn_z	   =  0.5d0/cosh3**2   / sigz
  d2atn_z	   =  1.0d0/cosh3**2   / sigz**2 * tanh2
  atn_z_u	   = (0.5d0 - 0.5d0*tanh2_u)
  datn_z_u	   = -0.5d0/cosh3_u**2 / sigz
  d2atn_z_u	   =  1.0d0/cosh3_u**2 / sigz**2 * tanh2_u
  
  if(xcase2 .eq. 1) then
    atn_z_u          = 1.d0
    datn_z_u         = 0.d0
    d2atn_z_u        = 0.d0
  endif
  if(xcase .eq. 2) then
    atn_z            = 1.d0
    datn_z           = 0.d0
    d2atn_z          = 0.d0
  endif
  
  FFprime_prof  =   prof1           *	 atn_z * atn_z_u
  dFF_dpsi         =   dprof1_dpsi  *	 atn_z * atn_z_u
  dFF_dpsi2        =   dprof1_dpsi2 *	 atn_z * atn_z_u  
  dFF_dz           = + prof1        * ( datn_z * atn_z_u  +	     atn_z * datn_z_u)
  dFF_dz2          = + prof1        * (d2atn_z * atn_z_u  +  2.d0 * datn_z * datn_z_u  +  atn_z * d2atn_z_u) 
  dFF_dpsi_dz      =   dprof1_dpsi  * ( datn_z * atn_z_u  +	     atn_z * datn_z_u)

endif

!-------- simple prof -------
F_prof       = prof0
dF_dpsi         = dprof0_dpsi
dF_dpsi2        = dprof0_dpsi2
dF_dz           = 0.d0
dF_dz2          = 0.d0
dF_dpsi_dz      = 0.d0

FFprime_prof = F_prof * dF_dpsi * atn
dFF_dpsi        = (dF_dpsi * dF_dpsi + F_prof * dF_dpsi2) * atn + (F_prof * dF_dpsi) * datn
dFF_dz          = 0.d0
dFF_dz2         = 0.d0
dFF_dpsi_dz     = 0.d0
!--------------------------------

!write(*,'(A,5e12.4)') ' F, FFprime new, old, diff : ',F0,F_profile,FFprime_profile, profFFp, FFprime_profile - profFFp
!write(*,'(A,5e12.4)') ' F, dF_dpsi, dF_dpsi2 , FFprime: ',F_profile,dF_dpsi, dF_dpsi2, FFprime_profile
!write(*,'(A,5e12.4)') ' dFFp_dpsi new, old, diff : ',dFF_dpsi, dprofFFp_dpsi, dFF_dpsi - dprofFFp_dpsi


return
end
