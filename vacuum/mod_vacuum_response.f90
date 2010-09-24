module vacuum_response_module
!-----------------------------------------------------------------------
! module contains the vacuum response matrix for all toroidal harmonics
!-----------------------------------------------------------------------

integer             :: n_dof_bnd                  ! degrees of freedom on the boundary
real*8, allocatable :: vacuum_response(:,:,:)     ! the vacuum response matrix (idrive,ireponse,itor)

real*8, allocatable :: vacuum_response2(:,:,:) !for testing only
real*8, allocatable :: vacuum_response3(:,:,:) !for testing only

! --- Resistive wall only.
!
!   dY/dt = - 1/(sigma * d) * [YY] * Y - [YE] * dP/dt
!   B_par = [EY] * Y + [EE] * P
!
!     Y: wall currents
!     P: poloidal flux
!
! Note: Matrix [EE] has the same structure as the ideal response matrix.
!
real*8              :: wall_resistivity = 0.001 !#### units etc.???
real*8              :: wall_thickness   = 0.1   !#### units etc.???

integer             :: n_wall_curr                ! Number of wall current potentials.
real*8, allocatable :: wall_curr(:)               ! Wall current potentials (Y).

real*8, allocatable :: diagonal_yy(:)
real*8, allocatable :: matrix_ye(:,:,:)
real*8, allocatable :: matrix_ey(:,:,:)
real*8, allocatable :: matrix_ee(:,:,:)

end module vacuum_response_module
