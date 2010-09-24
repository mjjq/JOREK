module vacuum_response_module
!-----------------------------------------------------------------------
! module contains the vacuum response matrix for all toroidal harmonics
!-----------------------------------------------------------------------

integer             :: n_dof_bnd                  ! degrees of freedom on the boundary
real*8, allocatable :: vacuum_response(:,:,:)     ! the vacuum response matrix (idrive,ireponse,itor)

real*8, allocatable :: vacuum_response2(:,:,:) !for testing only
real*8, allocatable :: vacuum_response3(:,:,:) !for testing only

! --- Resistive wall only.
integer             :: n_wall_curr                ! Number of wall current potentials.
real*8, allocatable :: wall_curr(:)               ! Wall current potentials.
real*8, allocatable :: diagonal_yy(:)
real*8, allocatable :: matrix_ye(:,:,:)
real*8, allocatable :: matrix_ey(:,:,:)
real*8, allocatable :: matrix_ee(:,:,:)

end module vacuum_response_module
