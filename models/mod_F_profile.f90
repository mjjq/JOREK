module mod_F_profile

contains


subroutine F_profile(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,&
                     F_prof, &
                       dF_dpsi, dF_dz, &                                          ! 1st order derivatives
                       dF_dpsi2, dF_dz2, dF_dpsi_dz, &                            ! 2nd order derivatives
                       dF_dpsi3, dF_dpsi_dz2, dF_dpsi2_dz,  dF_dz3, &             ! 2rd order derivatives
                       dF_dpsi4, dF_dpsi_dz3, dF_dpsi2_dz2, dF_dpsi3_dz, dF_dz4, &! 4th order derivatives
                     FFprime_prof, &
                       dFF_dpsi, dFF_dz, &                                             ! 1st order derivatives
                       dFF_dpsi2, dFF_dz2, dFF_dpsi_dz, &                              ! 2nd order derivatives
                       dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3, &              ! 2rd order derivatives
                       dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4)  ! 4th order derivatives
  !-----------------------------------------------------------------------
  ! Routine to calculate F(psi) analytically by integrating FFprime formula
  !-----------------------------------------------------------------------
  use phys_module
  use vacuum, only: current_FB_fact

  implicit none

  logical :: xpoint2
  integer :: xcase2
  real*8  :: Fconst, profF, profF1, profFF, profFF1
  real*8  :: F_prof   
  real*8  ::   dF_dpsi, dF_dz                                             ! 1st order derivatives
  real*8  ::   dF_dpsi2, dF_dz2, dF_dpsi_dz                               ! 2nd order derivatives
  real*8  ::   dF_dpsi3, dF_dpsi_dz2, dF_dpsi2_dz,  dF_dz3                ! 2rd order derivatives
  real*8  ::   dF_dpsi4, dF_dpsi_dz3, dF_dpsi2_dz2, dF_dpsi3_dz, dF_dz4   ! 4th order derivatives
  real*8  ::   dF_dpsi5
  real*8  :: FFprime_prof   
  real*8  ::   dFF_dpsi, dFF_dz                                                ! 1st order derivatives
  real*8  ::   dFF_dpsi2, dFF_dz2, dFF_dpsi_dz                                 ! 2nd order derivatives
  real*8  ::   dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3                 ! 2rd order derivatives
  real*8  ::   dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4   ! 4th order derivatives
  real*8  :: psi_edge, sqrt_edge, F_edge, F_constant
  real*8  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd,  psi_n, psi_barrier, sig_F, sigz, delta_psi
  real*8  :: psi_star
  real*8  :: atn, datn, d2atn, d3atn
  real*8  :: Z_star, Z_star_u
  real*8  :: atn_z,   datn_z,   d2atn_z,   d3atn_z,   d4atn_z
  real*8  :: atn_z_u, datn_z_u, d2atn_z_u, d3atn_z_u, d4atn_z_u
  real*8  :: d_0, d_pert, d2_pert, d3_pert
  real*8  :: tanh2, cosh3, sinh3, cosh4, sinh4, tanh2_u, cosh3_u, sinh3_u, cosh4_u, sinh4_u
  real*8  :: alfa,profFFp, dprofFFp_dpsi, prof_bnd
  real*8  :: prof0,     dprof0_dpsi,     dprof0_dpsi2,     dprof0_dpsi3, dprof0_dpsi4, dprof0_dpsi5
  real*8  :: poly,      dpoly_dpsi,      dpoly_dpsi2,      dpoly_dpsi3
  real*8  :: pert,      dpert_dpsi,      dpert_dpsi2,      dpert_dpsi3
  real*8  :: sqrt_term, dsqrt_term_dpsi, dsqrt_term_dpsi2, dsqrt_term_dpsi3
  real*8  :: no_delta_psi
  ! for interpolating numerical profiles
  integer :: left, right, mid
  real*8  :: aux1, aux2

  ! --- Jorek uses -FF' as a convention, so we need to reverse the profile before integrating
  real*8  :: myFF_0, myFF_1, myFF_coef(8)

  ! --- psi_norm
  psi_n = (psi - psi_axis)/(psi_bnd - psi_axis)
  delta_psi = (psi_bnd - psi_axis)
  no_delta_psi = 1.d0
  if (FF_coef(9) .eq. 1.d0) no_delta_psi = delta_psi

  ! --- Value of F at the edge should be F0 (NEEDS TO BE DECIDED!!!)
  F_edge = F0

  ! --- use numerical representation.
  
  ! --- Interpolate profile and derivatives to position psi_n by bisections.
  left  = 1
  right = num_Fprofile_len
  do
    if ( right == left + 1 ) exit
    mid = (left + right) / 2
    if ( num_Fprofile_x(mid) >= psi_n ) then
      right = mid
    else
      left = mid
    end if
  end do
  aux1 = (psi_n - num_Fprofile_x(left)) / (num_Fprofile_x(right) - num_Fprofile_x(left))
  aux2 = (1. - aux1)
  prof0        = num_Fprofile_y0(left) * aux2 + num_Fprofile_y0(right) * aux1
  dprof0_dpsi  = num_Fprofile_y1(left) * aux2 + num_Fprofile_y1(right) * aux1
  dprof0_dpsi2 = num_Fprofile_y2(left) * aux2 + num_Fprofile_y2(right) * aux1
  dprof0_dpsi3 = num_Fprofile_y3(left) * aux2 + num_Fprofile_y3(right) * aux1
  dprof0_dpsi4 = num_Fprofile_y4(left) * aux2 + num_Fprofile_y4(right) * aux1
  dprof0_dpsi5 = num_Fprofile_y5(left) * aux2 + num_Fprofile_y5(right) * aux1
    




  ! --- Save F-profile
  F_prof    = + prof0
  dF_dpsi   = + dprof0_dpsi
  dF_dpsi2  = + dprof0_dpsi2
  dF_dpsi3  = + dprof0_dpsi3
  dF_dpsi4  = + dprof0_dpsi4
  dF_dpsi5  = + dprof0_dpsi5

  ! --- Save FF'-profile
  FFprime_prof = F_prof * dF_dpsi
  dFF_dpsi     = dF_dpsi * dF_dpsi + F_prof * dF_dpsi2
  dFF_dpsi2    = 3.0 * dF_dpsi * dF_dpsi2 + F_prof * dF_dpsi3
  dFF_dpsi3    = 3.0 * dF_dpsi2 * dF_dpsi2 + 4.0 * dF_dpsi * dF_dpsi3 + F_prof * dF_dpsi4
  dFF_dpsi4    = 3.0 * dF_dpsi3 * dF_dpsi2 + 7.0 * dF_dpsi2 * dF_dpsi3 + 5.0 * dF_dpsi * dF_dpsi4 + F_prof * dF_dpsi5
  dFF_dz       = 0.d0
  dFF_dz2      = 0.d0
  dFF_dz3      = 0.d0
  dFF_dz4      = 0.d0
  dFF_dpsi_dz  = 0.d0
  dFF_dpsi_dz2 = 0.d0
  dFF_dpsi_dz3 = 0.d0
  dFF_dpsi2_dz = 0.d0
  dFF_dpsi2_dz2= 0.d0
  dFF_dpsi3_dz = 0.d0




  ! --- Cut-off at X-points
  if ( xpoint2 ) then

    sigz = 0.1d0

    if (xcase2 .eq. 1) then
      atn_z_u   = 1.d0
      datn_z_u  = 0.d0
      d2atn_z_u = 0.d0
      d3atn_z_u = 0.d0
      d4atn_z_u = 0.d0
    else
      Z_star_u  = (Z-Z_xpoint(2))/sigz
      Z_star_u  = min( max( Z_star_u, -40.d0), 40.d0) ! avoid floating-point exceptions
      
      tanh2_u   = tanh(Z_star_u)
      cosh3_u   = cosh(Z_star_u)
      sinh3_u   = sinh(Z_star_u)
      cosh4_u   = cosh(2.0*Z_star_u)
      sinh4_u   = sinh(2.0*Z_star_u)
      
      atn_z_u   = (0.5d0 - 0.5d0*tanh2_u)
      datn_z_u  = -0.5d0/cosh3_u**2 / sigz
      d2atn_z_u =  1.0d0/cosh3_u**2 / sigz**2 * tanh2_u
      d3atn_z_u = -1.0d0/cosh3_u**4 / sigz**3 * (-2.d0 + cosh4_u) 
      d4atn_z_u = -4.0d0/cosh3_u**5 / sigz**4 * (-2.d0 + cosh4_u) * sinh3_u &
                  -1.0d0/cosh3_u**4 / sigz**4 * (-2.d0 * sinh4_u) 
    endif

    if (xcase2 .eq. 2) then
      atn_z   = 1.d0
      datn_z  = 0.d0
      d2atn_z = 0.d0
      d3atn_z = 0.d0
      d4atn_z = 0.d0
    else
      Z_star  = (Z_xpoint(1)-Z)/sigz
      Z_star  = min( max( Z_star, -40.d0), 40.d0) ! avoid floating-point exceptions
      
      tanh2   = tanh(Z_star)
      cosh3   = cosh(Z_star)
      sinh3   = sinh(Z_star)
      cosh4   = cosh(2.0*Z_star)
      sinh4   = sinh(2.0*Z_star)
      
      atn_z   = (0.5d0 - 0.5d0*tanh2)
      datn_z  =  0.5d0/cosh3**2   / sigz
      d2atn_z =  1.0d0/cosh3**2   / sigz**2 * tanh2
      d3atn_z = -1.0d0/cosh3**4 / sigz**3 * (-2.d0 + cosh4) 
      d4atn_z = -4.0d0/cosh3**5 / sigz**4 * (-2.d0 + cosh4) * sinh3 &
                -1.0d0/cosh3**4 / sigz**4 * (-2.d0 * sinh4) 
    endif 

    F_prof      = F_edge + (F_prof - F_edge) *    atn_z * atn_z_u
    dF_dpsi     = dF_dpsi                    *    atn_z * atn_z_u
    dF_dpsi2    = dF_dpsi2                   *    atn_z * atn_z_u
    dF_dpsi3    = dF_dpsi3                   *    atn_z * atn_z_u
    dF_dpsi4    = dF_dpsi4                   *    atn_z * atn_z_u
    dF_dz       = (F_prof - F_edge)          * ( datn_z * atn_z_u +         atn_z * datn_z_u)
    dF_dz2      = (F_prof - F_edge)          * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  +  atn_z * d2atn_z_u)
    dF_dz3      = (F_prof - F_edge)          * (d3atn_z * atn_z_u + 4.d0 * d2atn_z * datn_z_u + 4.d0 *  datn_z * d2atn_z_u  + atn_z * d3atn_z_u)
    dF_dz4      = (F_prof - F_edge)          * (d4atn_z * atn_z_u + 4.d0 * d3atn_z * datn_z_u + 6.d0 * d2atn_z * d2atn_z_u &
                                               + 4.d0 * datn_z * d3atn_z_u + atn_z * d4atn_z_u)
    dF_dpsi_dz  = dF_dpsi                    * ( datn_z * atn_z_u +         atn_z * datn_z_u)
    dF_dpsi_dz2 = dF_dpsi                    * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  +  atn_z * d2atn_z_u)
    dF_dpsi_dz3 = dF_dpsi                    * (d3atn_z * atn_z_u + 4.d0 * d2atn_z * datn_z_u + 4.d0 *  datn_z * d2atn_z_u  + atn_z * d3atn_z_u)
    dF_dpsi2_dz = dF_dpsi2                   * ( datn_z * atn_z_u +         atn_z * datn_z_u)
    dF_dpsi2_dz2= dF_dpsi2                   * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  +  atn_z * d2atn_z_u)
    dF_dpsi3_dz = dF_dpsi3                   * ( datn_z * atn_z_u +         atn_z * datn_z_u)

    FFprime_prof = FFprime_prof *    atn_z * atn_z_u
    dFF_dpsi     = dFF_dpsi     *    atn_z * atn_z_u
    dFF_dpsi2    = dFF_dpsi2    *    atn_z * atn_z_u
    dFF_dpsi3    = dFF_dpsi3    *    atn_z * atn_z_u
    dFF_dpsi4    = dFF_dpsi4    *    atn_z * atn_z_u
    dFF_dz       = FFprime_prof * ( datn_z * atn_z_u +         atn_z * datn_z_u)
    dFF_dz2      = FFprime_prof * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  +  atn_z * d2atn_z_u) 
    dFF_dz3      = FFprime_prof * (d3atn_z * atn_z_u + 4.d0 * d2atn_z * datn_z_u + 4.d0 *  datn_z * d2atn_z_u  + atn_z * d3atn_z_u)
    dFF_dz4      = FFprime_prof * (d4atn_z * atn_z_u + 4.d0 * d3atn_z * datn_z_u + 6.d0 * d2atn_z * d2atn_z_u &
                                  + 4.d0 * datn_z * d3atn_z_u + atn_z * d4atn_z_u)
    dFF_dpsi_dz  = dFF_dpsi     * ( datn_z * atn_z_u +         atn_z * datn_z_u)
    dFF_dpsi_dz2 = dFF_dpsi     * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  +  atn_z * d2atn_z_u)
    dFF_dpsi_dz3 = dFF_dpsi     * (d3atn_z * atn_z_u + 4.d0 * d2atn_z * datn_z_u + 4.d0 *  datn_z * d2atn_z_u  + atn_z * d3atn_z_u)
    dFF_dpsi2_dz = dFF_dpsi2    * ( datn_z * atn_z_u +         atn_z * datn_z_u)
    dFF_dpsi2_dz2= dFF_dpsi2    * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  +  atn_z * d2atn_z_u)
    dFF_dpsi3_dz = dFF_dpsi3    * ( datn_z * atn_z_u +         atn_z * datn_z_u)

  endif

  if (freeboundary_equil .and. num_ffprime) then            !if the ffprime profile is given in a file and freeboundary equilibrium is on,
                                                           !the full profile is multiplied by a factor in order to iterate to a given current   
     FFprime_prof  = FFprime_prof  * current_FB_fact
     dFF_dpsi      = dFF_dpsi      * current_FB_fact
     dFF_dpsi2     = dFF_dpsi2     * current_FB_fact
     dFF_dpsi3     = dFF_dpsi3     * current_FB_fact
     dFF_dpsi4     = dFF_dpsi4     * current_FB_fact
     dFF_dz        = dFF_dz        * current_FB_fact
     dFF_dz2       = dFF_dz2       * current_FB_fact
     dFF_dz3       = dFF_dz3       * current_FB_fact
     dFF_dz4       = dFF_dz4       * current_FB_fact
     dFF_dpsi_dz   = dFF_dpsi_dz   * current_FB_fact
     dFF_dpsi_dz2  = dFF_dpsi_dz2  * current_FB_fact
     dFF_dpsi_dz3  = dFF_dpsi_dz3  * current_FB_fact
     dFF_dpsi2_dz  = dFF_dpsi2_dz  * current_FB_fact
     dFF_dpsi2_dz2 = dFF_dpsi2_dz2 * current_FB_fact
     dFF_dpsi3_dz  = dFF_dpsi3_dz  * current_FB_fact
  
  end if

  ! --- Normally not needed since we don't allow FF_1 with model710
  FFprime_prof = FFprime_prof + FF_1

  return
end subroutine F_profile





! --- This routine integrates the FFprime numerically.
subroutine integrate_F_profile()

  use phys_module, only: xpoint, F0, FF_0, FF_1, FF_coef, n_Fprofile_internal, Fprofile_internal, Fprofile_psi_max, num_ffprime

  implicit none
  
  ! --- Internal variables
  integer :: i_prof, j, j_prev, i_step, n_integral, n_integral_max, stay_safe
  real*8  :: no_delta_psi, delta_psi
  real*8  :: integral, integrand, dx, step, step_min, step_max, psi
  real*8  :: FFprime_now, FFprime_prev, FFprime_mid
  real*8  :: psi_edge, F_edge, psi_star, psi_mid, psi_prev
  real*8  :: diff, reference, quad_tol_min, quad_tol_max
  real*8  :: Z_fake, Z_xpoint_fake(2),  psi_axis_fake, psi_bnd_fake
  real*8  :: ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12,ff13,ff14,ff15
  real*8, allocatable  :: integral_large(:), psi_large(:)

  ! --- Because when we integrate the F-profile, we don't know the plasma yet...
  Z_fake        = 0.d0
  Z_xpoint_fake = 0.d0
  psi_axis_fake = 0.d0
  psi_bnd_fake  = 1.d0

  if ( (FF_coef(9) .ne. 1.d0) .and. (.not. num_ffprime) ) then
    write(*,*)'Full-MHD Warning!!! You have to use a denormalised edge perturbation for your'
    write(*,*)'                    FFprime, or use a numerical FFprime !!!'
    write(*,*)'                    ie. FF_coef(9) should be set to 1.d0, and FF_coef(6)'
    write(*,*)'                    should be denormalised as '
    write(*,*)'                    FF_coef(6) = FF_coef(6) * (psi_bnd-psi_axis)...'
    write(*,*)'                    Aborting...'
    stop
  endif

  ! --- The equilibrium psi amplitude
  delta_psi = (psi_bnd_fake - psi_axis_fake)
  no_delta_psi = 1.d0
  if (FF_coef(9) .eq. 1.d0) no_delta_psi = delta_psi

  ! --- The cutoff of the FFprime at the edge is traditionally the tanh at FF_coef(5), not at psi_n=1.0
  ! --- However, F needs to be F0 outside the plasma, so we take a point far away.
  ! --- To keep safe, for limiter plasma, the user might have non-zero FF' outside psi_n=1.0, because we
  ! --- simply don't solve this, so in this case, we use psi_n=1.0
  psi_edge  = 1.0
  if (xpoint) psi_edge  = max(1.0,FF_coef(5) + 2.0 * FF_coef(4))

  ! --- We need an idea of the amplitude of the FFprime (for the adaptive int-step)
  psi = 0.0
  call FFprime(.false.,0,Z_fake,Z_xpoint_fake,psi,     psi_axis_fake, psi_bnd_fake, FFprime_now, &
               ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12,ff13,ff14,ff15, .false.)
  call FFprime(.false.,0,Z_fake,Z_xpoint_fake,psi_edge,psi_axis_fake, psi_bnd_fake, FFprime_prev,&
               ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12,ff13,ff14,ff15, .false.)
  FFprime_now  = -FFprime_now
  FFprime_prev = -FFprime_prev
  reference = abs(FFprime_now - FFprime_prev)
  reference = max(reference, 1.d-3)

  ! --- Compute initial FFprime
  psi  = 0.d0
  call FFprime(.false.,0,Z_fake,Z_xpoint_fake,psi,psi_axis_fake, psi_bnd_fake, FFprime_prev,&
               ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12,ff13,ff14,ff15, .false.)
  FFprime_prev = -FFprime_prev
  
  ! --- We start the integral at the plasma core, where F should be close to F0
  ! --- Note: we correct the F-profile integration constant at the end, this starting point really doesn't matter...
  integral = F0**2

  ! --- The steps for integration. Assume the full profile will be as large as minimal step
  step     = 1.d-6
  step_min = 1.d-7
  step_max = 1.d-3
  n_integral_max = int(2*psi_edge/step_min)
  allocate(integral_large(n_integral_max))
  allocate(psi_large(n_integral_max))
  integral_large = 0.d0

  ! --- Tollerance to deviate from straight line
  quad_tol_min = 1.d-6
  quad_tol_max = 0.005 ! 0.5% of FFprime amplitude should be enough...

  ! --- Start integral
  psi_prev = 0.d0
  psi      = 0.d0
  i_step   = 0
  stay_safe= 0
  do while (.true.)
  
    ! --- Make sure we are not in an infinite loop...
    stay_safe = stay_safe + 1
    if (stay_safe .gt. 10*n_integral_max) then
      write(*,*)'Warning! F-profile integration stuck in infinite loop!'
      stop
    endif
    
    ! --- New step (we go from 0.0 to ~1.0)
    if (psi+step .gt. psi_edge) then
      step = psi_edge - psi
      psi  = psi_edge
    else
      psi  = psi + step
    endif
    
    ! --- Middle of step
    psi_mid = psi - 0.5*step
    
    ! --- FFprime at this psi and at mid-point
    call FFprime(.false.,0,Z_fake,Z_xpoint_fake,psi,     psi_axis_fake, psi_bnd_fake, FFprime_now, &
                 ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12,ff13,ff14,ff15, .false.)
    call FFprime(.false.,0,Z_fake,Z_xpoint_fake,psi_mid, psi_axis_fake, psi_bnd_fake, FFprime_mid, &
                 ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12,ff13,ff14,ff15, .false.)
    FFprime_now = -FFprime_now
    FFprime_mid = -FFprime_mid
    
    ! --- Check error from straight line
    if (psi .lt. psi_edge) then
      diff = abs(FFprime_mid - (FFprime_now+FFprime_prev)/2.0)
      if (diff / reference .gt. quad_tol_max) then
        if (step / 2.0 .ge. step_min) then ! this is our threshold, otherwise continue...
          step = step / 2.0
          psi = psi_prev
          cycle
        endif
      endif
      if (diff / reference .lt. quad_tol_min) then
        if (step * 2.0 .le. step_max) then ! this is our threshold, otherwise continue...
          step = step * 2.0
          psi = psi_prev
          cycle
        endif
      endif
    endif
    
    ! --- Take the contribution (factor two because FF' = 1/2(F**2)' )
    integrand = 2.0 * FFprime_mid
    dx        = step * delta_psi
    integral  = integral + integrand * dx
    
    ! --- Check n_step
    i_step = i_step + 1
    if (i_step .eq. n_integral_max) then
      write(*,*)'Warning! F-profile not fully integrated!'
      exit
    endif
    
    ! --- Save profile
    psi_large(i_step)      = psi
    integral_large(i_step) = integral
    psi_prev = psi
    FFprime_prev = FFprime_now
    
    if (psi .eq. psi_edge) exit
    
  enddo
  n_integral = i_step
  
  ! --- WARNING! THIS IS WHERE WE DECIDE THE F-profile OFFSET! TO BE DISCUSSED/DECIDED!!!
  ! --- Correct the offset to get F0 at the plasma edge (or at the core?)
  diff = F0**2 - integral_large(n_integral)
  integral_large = integral_large + diff
    
  ! --- Squared root profile (because FF' = 1/2(F**2)' )
  integral_large(1:n_integral) = integral_large(1:n_integral)**0.5
  
  ! --- With refined profile, interpolate n_Fprofile_internal on interval psi_n=[0.0,Fprofile_psi_max]
  j_prev = 2
  do i_prof = 1,n_Fprofile_internal
    psi = Fprofile_psi_max * real(i_prof-1)/real(n_Fprofile_internal-1)
    if (psi .le. psi_large(n_integral)) then
      do j=j_prev,n_integral
        if (psi_large(j) .gt. psi) then
          j_prev = j
          diff = (psi_large(j)-psi) / (psi_large(j)-psi_large(j-1))
          Fprofile_internal(i_prof) = integral_large(j) + (integral_large(j-1)-integral_large(j)) * diff
          exit
        endif
      enddo
    else ! outside plasma, F is just constant
      Fprofile_internal(i_prof) = integral_large(n_integral)
    endif
  enddo
  
  ! --- Correct the sign of F-profile depending on F0
  if (F0 .lt. 0.d0) Fprofile_internal = - Fprofile_internal
  
  deallocate(integral_large)
  deallocate(psi_large)
  
  return
  
end subroutine integrate_F_profile





! --- This routine checks that the numerically integrated F-profile is coherent with the input FFprime
subroutine check_F_profile_accuracy()

  use phys_module, only: n_Fprofile_internal, Fprofile_internal, Fprofile_psi_max, Fprofile_tolerance

  implicit none

  ! --- Internal variables
  integer :: i_prof, n_prof_core
  real*8  :: psi, psi_n(n_Fprofile_internal)
  real*8  :: iff1(n_Fprofile_internal),ff1(n_Fprofile_internal)
  real*8  :: ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12,ff13,ff14,ff15
  real*8  :: Z_fake, Z_xpoint_fake(2)
  real*8  :: psi_axis_fake, psi_bnd_fake
  real*8  :: accumulated_error, accumulated_profile, diff_average_percent

  Z_fake        = 0.d0
  Z_xpoint_fake = 0.d0
  psi_axis_fake = 0.d0
  psi_bnd_fake  = 1.d0

  accumulated_error   = 0.d0
  accumulated_profile = 0.d0
  n_prof_core         = 0
  do i_prof = 1,n_Fprofile_internal
    psi_n(i_prof) = Fprofile_psi_max * real(i_prof-1)/real(n_Fprofile_internal-1)
    if (psi_n(i_prof) .gt. 1.01) cycle ! count just the core (outside, FFprime should be zero)
    psi   = psi_n(i_prof)
    call FFprime  (.false.,0,Z_fake,Z_xpoint_fake,psi,psi_axis_fake, psi_bnd_fake,  ff1(i_prof),&
                   ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12,ff13,ff14,ff15, .false.)
    call FFprime  (.false.,0,Z_fake,Z_xpoint_fake,psi,psi_axis_fake, psi_bnd_fake, iff1(i_prof),&
                   ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12,ff13,ff14,ff15, .true.)
    accumulated_error   = accumulated_error   + abs(ff1(i_prof)-iff1(i_prof))
    accumulated_profile = accumulated_profile + abs(ff1(i_prof))
    n_prof_core = n_prof_core + 1 ! count just the core (outside, FFprime should be zero)
  enddo
  accumulated_error    = accumulated_error   / real(n_prof_core)
  accumulated_profile  = accumulated_profile / real(n_prof_core)
  diff_average_percent = 100.0*accumulated_error/accumulated_profile


  if (diff_average_percent .gt. Fprofile_tolerance) then
    write(*,*)'*** WARNING ***'
    write(*,*)'The FFprime could not be integrated accurately. Aborting...'
    write(*,*)'Error on profile: psi_norm, input-FFprime, FFprime (from integrated F-profile), error'
    do i_prof = 1,n_Fprofile_internal
      write(*,'(A,4e15.7)')'error on profile:',psi_n(i_prof),ff1(i_prof),iff1(i_prof),abs(ff1(i_prof)-iff1(i_prof))
    enddo
    write(*,'(A,1e15.7)')'averaged error (%):',diff_average_percent
    stop
  endif

  return
end subroutine check_F_profile_accuracy






end module mod_F_profile
