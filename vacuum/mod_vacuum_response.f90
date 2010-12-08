module vacuum_response_module
!-----------------------------------------------------------------------
! module contains the vacuum response matrix for all toroidal harmonics
!-----------------------------------------------------------------------
  
  ! --- Degrees of freedom of the boundary nodes
  integer             :: n_dof_bnd
  
  ! --- Ideal wall only
  real*8, allocatable :: vacuum_response(:,:,:)     ! the vacuum response matrix (idrive,ireponse,itor)
  real*8, allocatable :: vacuum_response2(:,:,:)    ! (for testing only)
  real*8, allocatable :: vacuum_response3(:,:,:)    ! (for testing only)
  
  ! --- Resistive wall only
  real*8              :: wall_resistivity = 0.001   !#### units etc.???
  real*8              :: wall_thickness   = 0.1     !#### units etc.???
  integer             :: n_wall_curr                ! Number of wall current potentials.
  real*8, allocatable :: wall_curr(:)               ! Wall current potentials (Y).
  real*8, allocatable :: diagonal_yy(:)             ! YY-matrix, see explanation below
  real*8, allocatable :: matrix_ye(:,:,:)           ! YE-matrix, see explanation below
  real*8, allocatable :: matrix_ey(:,:,:)           ! EY-matrix, see explanation below
  real*8, allocatable :: matrix_ee(:,:,:)           ! EE-matrix, see explanation below

 ! --- Coil contributions
  integer             :: n_coils                    ! the number of poloidal field coils
  real*8, allocatable :: R_coils(:), Z_coils(:)     ! positions of the poloidal field coils
  real*8, allocatable :: dR_coils(:), dZ_coils(:)   ! widths/heights of the poloidal field coils
  real*8, allocatable :: I_coils(:)                 ! the coil currents 
  real*8, allocatable :: external_field(:,:)        ! the external poloidal field  (n_dof_bnd,n_coils)
  
end module vacuum_response_module

! FURTHER INFORMATION:
! --------------------
!
!   (Y are the wall current potentials and Psi denotes the poloidal flux)
!
!
!   IDEAL WALL
!   ----------
!
!     B_par = [VR] * Psi
!
!
!   RESISTIVE WALL
!   --------------
!
!     dY/dt = - 1/(sigma * d) * [YY] * Y - [YE] * dPsi/dt
!  
!     B_par = [EY] * Y + [EE] * Psi
!  
!     (Matrix [EE] has the same structure as the ideal response matrix)
