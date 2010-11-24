subroutine initialise_parameters(my_id)
!-----------------------------------------------------------------------
! Initialise input parameters and read the input namelist
!-----------------------------------------------------------------------
use phys_module

implicit none

! --- Input variables.
integer, intent(in) :: my_id

! --- Internal variables.
integer :: in, i
real*8  :: psi_plot(1001),zj_plot(1001),dj_plot(1001),dz_plot(1001), z_plot(1001), z, zjz, dj_dpsi, dj_dz, psi_n

! --- Namelist with input parameters.
namelist /in1/  tstep, nstep, eta, visco, visco_par,                &
                restart,  regrid,                                   &
                n_R, n_Z, n_radial, n_pol, n_tht, n_flux,           &
                n_open,n_private,n_leg,  nout,                      &
                xr1, sig1, xr2, sig2,                               &
                R_begin, R_end, Z_begin, Z_end,                     &
                R_geo, Z_geo, amin, mf, fbnd, fpsi, mode,           &
                F0,                                                 &
                zjz_0, zjz_1, zj_coef,                              &
                rho_0, rho_1, rho_coef,                             &
                T_0,   T_1,   T_coef,                               &
                Ti_0,   Ti_1,   Ti_coef,                            &
                Te_0,   Te_1,   Te_coef,                            &
                FF_0,  FF_1,  FF_coef,                              &
                K_i_par, ZK_i_perp, K_e_par, ZK_e_perp,             &
                Zk_par, ZK_perp, D_par, D_perp,                     &
                Q_bar, sigma,                                       &
                particlesource, heatsource,                         &
                heatsource_i, heatsource_e,                         &
                eta_num, visco_num,                                 &
                ellip,tria_u,tria_l,quad_u,quad_l,                  &
                xampl,xwidth,xsig,xtheta,xshift,xleft, xpoint,      &
                rho_file, T_file, ffprime_file,                     &
                freeboundary, use_starwall, resistive_wall

if (my_id .eq. 0) then

  ! --- Preset input parameters to reasonable default values.
  tstep    = 1.d0
  tstep_in = 1.d0
  nstep    = 0

  eta   = 1.d-5
  visco = 1.d-5
  visco_par = 1.d-5

  restart      = .false.
  import_equil = .false.
  regrid       = .false.
  
  freeboundary   = .false. ! use free or fixed boundary?
  use_starwall   = .false. ! use the STARWALL vacuum solution? (freeboundary only)
  resistive_wall = .false. ! use a resistive or ideal wall?    (freeboundary only)

  n_R       = 3
  n_Z       = 3
  n_radial  = 11
  n_pol     = 6

  n_flux    = 11
  n_tht     = 17

  n_open    = 5
  n_leg     = 5
  n_private = 5

  R_geo     = 10.d0
  Z_geo     = 0.d0
  amin      = 1.d0

  F0        = 10.d0
  GAMMA     = 5.d0 / 3.d0

  mf        = 8
  fbnd      = 0.d0; fbnd(1) =2.d0

  ellip  = 1.d0
  tria_u = 0.d0
  tria_l = 0.d0
  quad_l = 0.d0
  quad_u = 0.d0

  xampl  = 0.d0
  xwidth = 0.d0
  xsig   = 0.d0
  xtheta = 0.d0
  xshift = 0.d0
  xleft  = 0.d0
  xpoint = .false.

  xr1  = 9999.d0
  sig1 = 9999.d0
  xr2  = 99999.d0
  sig2 = 99999.d0

  R_begin = -0.1d0
  R_end   =  0.1d0
  Z_begin = -0.1d0
  Z_end   = 0.1d0

  ZK_i_perp(1) = 1.d-5; ZK_i_perp(2) = 0.d0; ZK_i_perp(3)= 0.d0; ZK_i_perp(4)= 99.d0; ZK_i_perp(5) = 99.d0
  K_i_par      = 1.d+6
  ZK_e_perp(1) = 1.d-5; ZK_e_perp(2) = 0.d0; ZK_e_perp(3)= 0.d0; ZK_e_perp(4)= 99.d0; ZK_e_perp(5) = 99.d0
  K_e_par      = 1.d+7

  ZK_perp(1)   = 1.d-5; ZK_perp(2) = 0.d0; ZK_perp(3)= 0.d0; ZK_perp(4)= 99.d0; ZK_perp(5) = 99.d0
  ZK_par       = 10.
  
  D_perp(1)    = 1.d-5; D_perp(2) = 0.d0; D_perp(3)= 0.d0; D_perp(4)= 99.d0; D_perp(5) = 99.d0
  D_par        = 0.d0

  eta_num   = 0.d0
  visco_num = 0.d0

  heatsource_i   = 1.e-7
  heatsource_e   = 1.e-7
  heatsource     = 1.e-7
  particlesource = 1.e-5

  zjz_0 =  0.1173d0;  T_0   =  1.d-6  ; Ti_0   =  1.d-6  ;   Te_0   =  1.d-6  ;   rho_0 =  1.d0   ;   FF_0  =  1.d0
  zjz_1 =  0.0d0   ;  T_1   =  1.d-8  ; Ti_1   =  1.d-8  ;   Te_1   =  1.d-8  ;   rho_1 =  1.d0   ;   FF_1  =  0.d0

  zj_coef     = 0.d0;  zj_coef(1)  = -1.d0
  T_coef      = 0.d0;  T_coef(1)   = -1.d0
  Ti_coef     = 0.d0;  Ti_coef(1)  = -1.d0
  Te_coef     = 0.d0;  Te_coef(1)  = -1.d0
  rho_coef    = 0.d0;  rho_coef(1) =  0.d0
  FF_coef     = 0.d0;  FF_coef(1)  = -1.d0

  t_now       = 0.d0
  t_start     = 0.d0
  index_start = 0

  nout = 9999999

  rho_file      = 'none'
  T_file        = 'none'
  ffprime_file  = 'none'

  ! --- Read input parameters from namelist.
  if (my_id .eq. 0) read(5,in1)
  
  tstep_in = tstep

  if (nstep .gt. 0) allocate(energies(n_tor,2,nstep))
  if (nstep .gt. 0) allocate(xtime(nstep))

  write(*,*) 'USING MODEL 400'
  
  ! --- Read numerical profiles for rho, T, and ff'.
  call read_num_profiles()
  
endif

! --- Broadcast input parameters from MPI thread 0 to the others.
call broadcast_phys(my_id)

! --- Broadcast numerical input profiles from MPI thread 0 to the others.
call broadcast_num_profiles(my_id)

! --- Determine the derivatives of the numerical input profiles.
call derive_num_profiles()

return
end subroutine initialise_parameters
