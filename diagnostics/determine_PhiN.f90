!> Determine the normalized toroidal flux PhiN from flux surfaces and q-profile.
!!
!! The toroidal flux is determined from the flux surfaces and the q-profile by integration using
!! \f$ q = d\Phi/d\Psi \f$.
!!
!! - Set icorr_axis and icorr_bnd such that the discontinuities in the iota profile are corrected.
!! - To determine the input parameter surface_list, the routine find_flux_surfaces needs to be
!!   called first (see routine export_nemec for an example).
!! - To determine the input parameter q, the routine determine_q_profile needs to be called first.
!!
!! WARNING: The q-profile will be modified by this routine. It is "corrected" at magnetic axis and
!! boundary using linear extrapolation.
subroutine determine_PhiN(surface_list, q, PhiN, Phi_edge)
  
  use tr_module 
  use data_structure
  
  implicit none
  
  ! --- Routine parameters
  type (type_surface_list), intent(in)    :: surface_list             !< List of flux surfaces
  real*8,                   intent(inout) :: q(surface_list%n_psi)    !< q values at flux surfaces
  real*8,                   intent(out)   :: PhiN(surface_list%n_psi) !< Norm. tor. flux at surfaces
  real*8,                   intent(out)   :: Phi_edge                 !< Tor. flux at lcms
  
  ! --- Local variables
  real*8  :: deltaPsi
  integer :: i, icorr, n_psi, icorr_axis, icorr_bnd
  
  ! --- "Correct" how many points at the axis and at the boundary?
  n_psi=surface_list%n_psi
  icorr_axis = n_psi / 35
  icorr_bnd  = n_psi / 150
  
  ! --- "Correct" the q-profile.
  icorr = icorr_axis
  do i = 1, icorr
    q(i) = q(icorr+1) + (surface_list%psi_values(i)-surface_list%psi_values(icorr+1)) *            &
      (q(icorr+1)-q(icorr+2))/(surface_list%psi_values(icorr+1)-surface_list%psi_values(icorr+2))
  end do
  icorr = n_psi-icorr_bnd
  do i = icorr, n_psi
    q(i) = q(icorr-1) + (surface_list%psi_values(i)-surface_list%psi_values(icorr-1)) *            &
      (q(icorr-1)-q(icorr-2))/(surface_list%psi_values(icorr-1)-surface_list%psi_values(icorr-2))
  end do
  
  ! --- Determine the toroidal flux PhiN.
  PhiN(:) = 0.d0
  do i = 2, n_psi
    deltaPsi = surface_list%psi_values(i) - surface_list%psi_values(i-1)
    PhiN(i) = PhiN(i-1) + 0.5d0*(q(i)+q(i-1)) * deltaPsi
  end do
  Phi_edge = PhiN(n_psi)
  PhiN(:)  = PhiN(:) / Phi_edge
  
end subroutine determine_PhiN
