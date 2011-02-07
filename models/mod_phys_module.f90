module phys_module
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------

  use parameters
  
  implicit none

  real*8  :: tstep, tstep_in,  tstep_n(10), eta, visco, visco_par
  real*8  :: amin, fbnd(1026), fpsi(1026)
  real*8  :: R_boundary(1026), Z_boundary(1026), psi_boundary(1026)
  real*8  :: ellip, tria_u, tria_l, quad_u, quad_l
  real*8  :: xampl, xwidth, xsig, xtheta, xshift, xleft
  real*8  :: xr1, sig1, xr2, sig2
  real*8  :: F0, GAMMA, Q_bar, sigma
  real*8  :: zjz_0, zjz_1,  zj_coef(10)
  real*8  :: T_0,   T_1,    T_coef(10)
  real*8  :: Ti_0,  Ti_1,   Ti_coef(10)
  real*8  :: Te_0,  Te_1,   Te_coef(10)
  real*8  :: rho_0, rho_1,  rho_coef(10)
  real*8  :: FF_0,  FF_1,   FF_coef(10)
  real*8  :: particlesource, heatsource,heatsource_i, heatsource_e
  real*8  :: ZK_perp(10), ZK_par, ZK_i_perp(10), ZK_e_perp(10), K_i_par, K_e_par, D_perp(10), D_par
  real*8  :: D_neutral, tauIC
  real*8  :: pellet_amplitude, pellet_R, pellet_Z, pellet_phi, pellet_radius, pellet_sig, pellet_length
  real*8  :: pellet_psi, pellet_delta_psi
  real*8  :: eta_num, visco_num, visco_par_num, D_perp_num,Zk_perp_num
  real*8  :: t_start, t_now
  integer :: nstep, nstep_n(10), n_boundary
  integer :: mf, index_start, index_now, mode(n_tor), nout
  logical :: restart, regrid, import_equil, xpoint
  logical :: refinement     ! allow mesh refinement
  logical :: freeboundary   ! use free or fixed boundary?
  logical :: use_starwall   ! use the STARWALL vacuum solution? (free boundary only)
  logical :: resistive_wall ! use a resistive or ideal wall?    (free boundary only)
  real*8, allocatable   :: xtime(:), energies(:,:,:)
  character(len=3)      :: mode_type(n_tor)  ! 'cos' or 'sin'
  
  
  
  ! --- Grid parameters
  !   --- Rectangular grid
  integer :: n_R               ! Number of grid points in R-direction
  integer :: n_Z               ! Number of grid points in Z-direction
  real*8  :: R_begin, R_end    ! Extent of grid in R-direction
  real*8  :: Z_begin, Z_end    ! Extent of grid in Z-direction
  !   --- Polar grid
  integer :: n_radial          ! Number of radial grid points
  integer :: n_pol             ! Number of poloidal grid points
  real*8  :: R_geo, Z_geo      ! Center of the grid
  !   --- Flux surface grid (no X-point)
  integer :: n_flux            ! Number of radial grid points
  integer :: n_tht             ! Number of poloidal grid points
  !   --- X-point grid (uses also parameters from flux surface grid)
  integer :: n_open            ! Number of 'radial' grid points in the open flux region
  integer :: n_private         ! Number of 'radial' grid points in the private flux region
  integer :: n_leg             ! Number of 'poloidal' grid points along the divertor legs
  real*8  :: SIG_closed        ! Width with grid accumulation
  real*8  :: SIG_open          ! -"-
  real*8  :: SIG_private       ! -"-
  real*8  :: SIG_theta         ! -"-
  real*8  :: SIG_leg_0         ! -"-
  real*8  :: SIG_leg_1         ! -"-
  real*8  :: dPSI_open         ! Delta Psi grid extends into the open flux region
  real*8  :: dPSI_private      ! Delta Psi grid extends into the private flux region
  
  
  
  ! --- Numerical input profiles
  !   --- Density
  character(len=512)  :: rho_file        ! ASCII file the profile is read from.
  logical             :: num_rho         ! is set true if rho_file /= 'none'
  integer             :: num_rho_len     ! Number of points in profile
  real*8, allocatable :: num_rho_x(:)    ! Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_rho_y0(:), num_rho_y1(:), num_rho_y2(:), num_rho_y3(:) ! values and derivatives
  !   --- Temperature
  character(len=512)  :: T_file          ! ASCII file the profile is read from.
  logical             :: num_T           ! is set true if T_file /= 'none'
  integer             :: num_T_len       ! Number of points in profile
  real*8, allocatable :: num_T_x(:)      ! Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_T_y0(:), num_T_y1(:), num_T_y2(:), num_T_y3(:) ! values and derivatives
  !   --- FFprime
  character(len=512)  :: ffprime_file    ! ASCII file the profile is read from.
  logical             :: num_ffprime     ! is set true if ffprime_file /= 'none'
  integer             :: num_ffprime_len ! Number of points in profile
  real*8, allocatable :: num_ffprime_x(:)! Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_ffprime_y0(:), num_ffprime_y1(:), num_ffprime_y2(:) ! values and derivatives
  
end module phys_module
