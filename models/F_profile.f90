subroutine F_profile(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,&
                     F_prof,dF_dpsi,dF_dz, dF_dpsi2,dF_dz2,dF_dpsi_dz, &
                     FFprime_prof,dFF_dpsi,dFF_dz, dFF_dpsi2,dFF_dz2,dFF_dpsi_dz)
!-----------------------------------------------------------------------
!
! EvdP new routine to calculate F(psi) explicitly to determine eq B_tor profile
!
!-----------------------------------------------------------------------
use phys_module
use vacuum, only: current_FB_fact

implicit none

logical :: xpoint2
integer :: xcase2
real*8  :: Fconst, profF, profF1, F_prof, dF_dpsi, dF_dz, dF_dpsi2, dF_dz2, dF_dpsi_dz
real*8  :: FFprime_prof, profFF, profFF1, FF_prof, dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz
real*8  :: dF_dpsi3, psi_edge, F_edge, sqrt_edge
real*8  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd,  psi_n, psi_barrier, sig_F, sigz, delta_psi
real*8  :: psi_star
real*8  :: atn, datn, d2atn, d3atn
real*8  :: Z_star, Z_star_u
real*8  :: atn_z,   datn_z,   d2atn_z
real*8  :: atn_z_u, datn_z_u, d2atn_z_u
real*8  :: d_0, d_pert, d2_pert, d3_pert
real*8  :: tanh2, cosh3, tanh2_u, cosh3_u
real*8  :: alfa,profFFp, dprofFFp_dpsi, prof_bnd
real*8  :: prof0,     dprof0_dpsi,     dprof0_dpsi2,     dprof0_dpsi3
real*8  :: poly,      dpoly_dpsi,      dpoly_dpsi2,      dpoly_dpsi3
real*8  :: pert,      dpert_dpsi,      dpert_dpsi2,      dpert_dpsi3
real*8  :: sqrt_term, dsqrt_term_dpsi, dsqrt_term_dpsi2, dsqrt_term_dpsi3
real*8  :: no_delta_psi

! --- Jorek uses -FF' as a convention, so we need to reverse the profile before integrating
real*8  :: myFF_0, myFF_1, myFF_coef(8)

myFF_0 = - FF_0
myFF_1 = - FF_1
myFF_coef(1:8) =   FF_coef(1:8)
myFF_coef(6)   = - FF_coef(6)

! --- Initialise
F_prof     = 0.d0
dF_dpsi    = 0.d0
dF_dz      = 0.d0
dF_dpsi2   = 0.d0
dF_dz2     = 0.d0
dF_dpsi_dz = 0.d0

! --- There are some rules when using FF_coefs with the F-profile in Full-MHD
if (myFF_1 .ne. 0.d0) then
  write(*,*)'Full-MHD Warning!!! The F-profile does not like it if FF_1 is not zero !!!'
  write(*,*)'                    if you don,t respect this rule, we cannot guarantee that your F-profile and FFprime will be consistent!'
endif
if (myFF_coef(4) .gt. 0.01) then
  write(*,*)'Full-MHD Warning!!! The tanh at FF_coef(5) with width FF_coef(4) is supposed to be a cut-off at the plasma edge !!!'
  write(*,*)'                    ie. FF_coef(5) should be the edge of your plasma, and FF_coef(4) should be very small...'
  write(*,*)'                    if you don,t respect this rule, we cannot guarantee that your F-profile and FFprime will be consistent!'
endif

! --- the cutoff of the FFprime at the edge is traditionally the tanh at FF_coef(5), not at psi_n=1.0
psi_edge  = myFF_coef(5)

! --- psi_norm
psi_n = (psi - psi_axis)/(psi_bnd - psi_axis)
delta_psi = (psi_bnd - psi_axis)
no_delta_psi = 1.d0
if (FF_coef(9) .eq. 1.d0) no_delta_psi = delta_psi

! --- Polynomial part
poly        = ( psi_n  + myFF_coef(1)/2.d0 * psi_n**2 + myFF_coef(2)/3.d0 * psi_n**3 + myFF_coef(3)/4.d0 * psi_n**4 ) * delta_psi
dpoly_dpsi  = ( 1.d0   + myFF_coef(1)      * psi_n    + myFF_coef(2)      * psi_n**2 + myFF_coef(3)      * psi_n**3 )
dpoly_dpsi2 = (          myFF_coef(1)                 + myFF_coef(2)*2.0  * psi_n    + myFF_coef(3)*3.0  * psi_n**2 ) / delta_psi
dpoly_dpsi3 = (                                       + myFF_coef(2)*2.0             + myFF_coef(3)*6.0  * psi_n    ) / delta_psi**2

! --- Perturbation part
pert        = + myFF_coef(6) * tanh((psi_n - myFF_coef(7))/myFF_coef(8))    / 2.d0                                 * no_delta_psi
dpert_dpsi  = + myFF_coef(6) / cosh((psi_n - myFF_coef(7))/myFF_coef(8))**2 / (2.d0 * myFF_coef(8)) / delta_psi    * no_delta_psi
dpert_dpsi2 = - myFF_coef(6) / cosh((psi_n - myFF_coef(7))/myFF_coef(8))**3 / myFF_coef(8)**2       / delta_psi**2 * no_delta_psi &
                             * sinh((psi_n - myFF_coef(7))/myFF_coef(8))
dpert_dpsi2 = + myFF_coef(6) / cosh((psi_n - myFF_coef(7))/myFF_coef(8))**4 / myFF_coef(8)**3 * 3.0 / delta_psi**3 * no_delta_psi &
                             * sinh((psi_n - myFF_coef(7))/myFF_coef(8))**2                                                       &
              - myFF_coef(6) / cosh((psi_n - myFF_coef(7))/myFF_coef(8))**2 / myFF_coef(8)**3       / delta_psi**3 * no_delta_psi
dpert_dpsi3 = - myFF_coef(6) / cosh((psi_n - myFF_coef(7))/myFF_coef(8))**5 / myFF_coef(8)**4 *12.0 / delta_psi**4 * no_delta_psi &
                             * sinh((psi_n - myFF_coef(7))/myFF_coef(8))**3                                                       &
              + myFF_coef(6) / cosh((psi_n - myFF_coef(7))/myFF_coef(8))**3 / myFF_coef(8)**4 * 3.0 / delta_psi**4 * no_delta_psi &
                       * 2.0 * sinh((psi_n - myFF_coef(7))/myFF_coef(8))                                                          &
              + myFF_coef(6) / cosh((psi_n - myFF_coef(7))/myFF_coef(8))**3 / myFF_coef(8)**4 * 2.0 / delta_psi**4 * no_delta_psi &
                             * sinh((psi_n - myFF_coef(7))/myFF_coef(8))

! --- Value of F at the edge
sqrt_edge =   2.0 * (myFF_0 - myFF_1) * (psi_edge + myFF_coef(1)/2.d0 * psi_edge**2 &
                                                  + myFF_coef(2)/3.d0 * psi_edge**3 &
                                                  + myFF_coef(3)/4.d0 * psi_edge**4 ) * delta_psi &
            + 2.0 * myFF_coef(6) * tanh((psi_edge - myFF_coef(7))/myFF_coef(8)) / 2.d0 &
            + F0**2
F_edge    = sqrt_edge**0.5 

! --- sqrt part
sqrt_term        = 2.0 * (myFF_0 - myFF_1) * poly        + 2.0 * pert            + F0**2 
dsqrt_term_dpsi  = 2.0 * (myFF_0 - myFF_1) * dpoly_dpsi  + 2.0 * dpert_dpsi 
dsqrt_term_dpsi2 = 2.0 * (myFF_0 - myFF_1) * dpoly_dpsi2 + 2.0 * dpert_dpsi2
dsqrt_term_dpsi3 = 2.0 * (myFF_0 - myFF_1) * dpoly_dpsi3 + 2.0 * dpert_dpsi3

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





! --- Cut-off at plasma edge
sig_F     = myFF_coef(4)

psi_star  = (psi_n-psi_edge)/sig_F
psi_star  = min( max( psi_star, -40.d0), 40.d0) ! avoid floating-point exceptions

tanh2   = tanh(psi_star)
cosh3   = cosh(psi_star)

atn   = (0.5d0 - 0.5d0*tanh2)
datn  = -0.5d0/cosh3**2 / sig_F
d2atn =  1.0d0/cosh3**2 / sig_F**2 * tanh2

F_prof     = F_edge + (F_prof - F_edge) * atn
dF_dpsi    = dF_dpsi                    * atn + (F_prof - F_edge) *   datn
dF_dpsi2   = dF_dpsi2                   * atn + (F_prof - F_edge) *   d2atn + 2.0 * dF_dpsi * datn

FFprime_prof = FFprime_prof * atn
dFF_dpsi     = dFF_dpsi     * atn + FFprime_prof * datn
dFF_dpsi2    = dFF_dpsi2    * atn + FFprime_prof * d2atn + 2.0 * dFF_dpsi * datn






! --- Cut-off at X-points
if ( xpoint2 ) then
  
  sigz = 0.1d0

  if (xcase2 .eq. 1) then
    atn_z_u   = 1.d0
    datn_z_u  = 0.d0
    d2atn_z_u = 0.d0
  else
    Z_star_u  = (Z-Z_xpoint(2))/sigz
    Z_star_u  = min( max( Z_star_u, -40.d0), 40.d0) ! avoid floating-point exceptions
    
    tanh2_u   = tanh(Z_star_u)
    cosh3_u   = cosh(Z_star_u)

    atn_z_u   = (0.5d0 - 0.5d0*tanh2_u)
    datn_z_u  = -0.5d0/cosh3_u**2 / sigz
    d2atn_z_u =  1.0d0/cosh3_u**2 / sigz**2 * tanh2_u
  endif
  if (xcase2 .eq. 2) then
    atn_z   = 1.d0
    datn_z  = 0.d0
    d2atn_z = 0.d0
  else
    Z_star  = (Z_xpoint(1)-Z)/sigz
    Z_star  = min( max( Z_star, -40.d0), 40.d0) ! avoid floating-point exceptions

    tanh2   = tanh(Z_star)
    cosh3   = cosh(Z_star)
      
    atn_z   = (0.5d0 - 0.5d0*tanh2)
    datn_z  =  0.5d0/cosh3**2   / sigz
    d2atn_z =  1.0d0/cosh3**2   / sigz**2 * tanh2
  endif 
  
  F_prof     = F_edge + (F_prof - F_edge) *    atn_z * atn_z_u
  dF_dpsi    = dF_dpsi                    *    atn_z * atn_z_u
  dF_dpsi2   = dF_dpsi2                   *    atn_z * atn_z_u
  dF_dz      = (F_prof - F_edge)          * ( datn_z * atn_z_u +         atn_z * datn_z_u)
  dF_dz2     = (F_prof - F_edge)          * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  +  atn_z * d2atn_z_u)
  dF_dpsi_dz = dF_dpsi                    * ( datn_z * atn_z_u +         atn_z * datn_z_u)
 
  FFprime_prof = FFprime_prof *    atn_z * atn_z_u
  dFF_dpsi     = dFF_dpsi     *    atn_z * atn_z_u
  dFF_dpsi2    = dFF_dpsi2    *    atn_z * atn_z_u
  dFF_dz       = FFprime_prof * ( datn_z * atn_z_u +         atn_z * datn_z_u)
  dFF_dz2      = FFprime_prof * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  +  atn_z * d2atn_z_u) 
  dFF_dpsi_dz  = dFF_dpsi     * ( datn_z * atn_z_u +         atn_z * datn_z_u)

endif






if (freeboundary_equil .and. num_ffprime) then            !if the ffprime profile is given in a file and freeboundary equilibrium is on,
                                                         !the full profile is multiplied by a factor in order to iterate to a given current   
   FFprime_prof  = FFprime_prof  * current_FB_fact
   dFF_dpsi      = dFF_dpsi      * current_FB_fact
   dFF_dpsi2     = dFF_dpsi2     * current_FB_fact
   dFF_dz        = dFF_dz        * current_FB_fact
   dFF_dz2       = dFF_dz2       * current_FB_fact
   dFF_dpsi_dz   = dFF_dpsi_dz   * current_FB_fact

end if

FFprime_prof = FFprime_prof + FF_1


return
end
