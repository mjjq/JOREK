module phys_module

  use parameters
  
  implicit none

  real*8  :: tstep, tstep_in, eta, visco, visco_par
  real*8  :: R_begin, R_end, Z_begin, Z_end, R_geo, Z_geo, amin, fbnd(1026), fpsi(1026)
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
  integer :: nstep, n_R, n_Z, n_radial, n_pol, n_tht, n_flux, n_boundary
  integer :: n_open, n_private, n_leg
  integer :: mf, index_start, index_now, mode(n_tor), nout
  logical :: restart, regrid, import_equil, xpoint
  logical :: refinement
  logical :: freeboundary   ! use free or fixed boundary?
  logical :: use_starwall   ! use the STARWALL vacuum solution? (free boundary only)
  logical :: resistive_wall ! use a resistive or ideal wall?    (free boundary only)
  real*8, allocatable   :: xtime(:), energies(:,:,:)
  character(len=3)      :: mode_type(n_tor)  ! 'cos' or 'sin'
  ! --- Numerical input profiles
  character(len=512)  :: rho_file, T_file, ffprime_file ! Files, the profiles are read from.
  logical             :: num_rho, num_T, num_ffprime    ! Numerical representation? (e.g., num_rho is set .true., if rho_file /= 'none')
  integer             :: num_rho_len, num_T_len, num_ffprime_len
  real*8, allocatable :: num_rho_x(:), num_rho_y0(:), num_rho_y1(:), num_rho_y2(:), num_rho_y3(:)
  real*8, allocatable :: num_T_x(:), num_T_y0(:), num_T_y1(:), num_T_y2(:), num_T_y3(:)
  real*8, allocatable :: num_ffprime_x(:), num_ffprime_y0(:), num_ffprime_y1(:), num_ffprime_y2(:)
  
end module phys_module
