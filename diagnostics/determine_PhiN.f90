!> Determine the normalized toroidal flux PhiN from flux surfaces and
!! q-profile.
!!
!! - To determine the input parameter surface_list, the routine
!!   find_flux_surfaces needs to be called first (see routine
!!   export_nemec for an example).
!! - To determine the input parameter q, the routine determine_q_profile
!!   needs to be called first.
!!
!! WARNING: The q-profile will be "corrected" at magnetic axis and
!! boundary using linear extrapolation.
subroutine determine_PhiN(surface_list, q, PhiN, Phi_edge)
  
  use tr_module 
  use data_structure
  
  implicit none
  
  ! --- Routine parameters
  type (type_surface_list), intent(in)    :: surface_list
  real*8,                   intent(inout) :: q(surface_list%n_psi)
  real*8,                   intent(out)   :: PhiN(surface_list%n_psi)
  real*8,                   intent(out)   :: Phi_edge
  
  ! --- Local variables
  real*8  :: deltaPsi
  integer :: i, ii, n_psi
  
  ! --- "Correct" q-profile.
  n_psi=surface_list%n_psi
  ii = n_psi/50
  do i = 1, ii
    q(i) = q(ii+1) + (surface_list%psi_values(i)-surface_list%psi_values(ii+1)) * &
      (q(ii+1)-q(ii+2))/(surface_list%psi_values(ii+1)-surface_list%psi_values(ii+2))
  end do
  ii = n_psi-n_psi/150
  do i = ii, n_psi
    q(i) = q(ii-1) + (surface_list%psi_values(i)-surface_list%psi_values(ii-1)) * &
      (q(ii-1)-q(ii-2))/(surface_list%psi_values(ii-1)-surface_list%psi_values(ii-2))
  end do
  
  ! --- Determine toroidal flux.
  PhiN(:) = 0.d0
  do i = 2, n_psi
    deltaPsi = surface_list%psi_values(i) - surface_list%psi_values(i-1)
    PhiN(i) = PhiN(i-1) + 0.5d0*(q(i)+q(i-1)) * deltaPsi
  end do
  Phi_edge = PhiN(n_psi)
  PhiN(:) = PhiN(:) / Phi_edge
  
end subroutine determine_PhiN
