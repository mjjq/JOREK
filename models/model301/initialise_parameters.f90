!> Initialise input parameters and read the input namelist
subroutine initialise_parameters(my_id)

use tr_module
use phys_module

implicit none

! --- Routine parameters
integer, intent(in) :: my_id

! --- Namelist with input parameters.
namelist /in1/  tstep, nstep, eta, visco, visco_par,                &
                restart, regrid,                                    &
                n_R, n_Z, n_radial, n_pol, n_tht, n_flux,           &
                n_open, n_private, n_leg,                           &
                SIG_closed, SIG_open, SIG_private, SIG_theta,       &
                SIG_leg_0, SIG_leg_1, dPSI_open, dPSI_private,      &
                nout, xr1, sig1, xr2, sig2,                         &
                R_begin, R_end, Z_begin, Z_end,                     &
                R_geo, Z_geo, amin, mf, fbnd, fpsi, mode,           &
                R_boundary, Z_boundary, psi_boundary, n_boundary,   &
                F0,                                                 &
                zjz_0, zjz_1, zj_coef,                              &
                rho_0, rho_1, rho_coef,                             &
                T_0,   T_1,   T_coef,                               &
                FF_0,  FF_1,  FF_coef,                              &
                ZK_par, ZK_perp, D_par, D_perp,                     &
                particlesource, heatsource,                         &
                eta_num, visco_num, visco_par_num, D_perp_num,      &
                ellip,tria_u,tria_l,quad_u,quad_l,                  &
                xampl,xwidth,xsig,xtheta,xshift,xleft, xpoint,      &
                rho_file, T_file, ffprime_file, freeboundary_equil, &
                freeboundary, use_starwall, resistive_wall,         &
                produce_live_data

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
  
  n_R       = 3
  n_Z       = 3
  n_radial  = 11
  n_pol     = 6

  n_flux    = 11
  n_tht     = 17

  n_open    = 5
  n_leg     = 5
  n_private = 5
  
  SIG_closed  = 0.1d0
  SIG_open    = 0.1d0
  SIG_private = 0.1d0
  SIG_theta   = 0.03d0
  SIG_leg_0   = 0.05d0
  SIG_leg_1   = 0.2d0
  
  dPSI_open    = 0.11
  dPSI_private = 0.03
  
  R_geo     = 10.d0
  Z_geo     = 0.d0
  amin      = 1.d0

  F0        = 10.d0
  GAMMA     = 5.d0 / 3.d0

  mf        = 8
  fbnd      = 0.d0; fbnd(1) =2.d0

  R_boundary   = 0.d0
  Z_boundary   = 0.d0
  psi_boundary = 0.d0
  n_boundary   = 0

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

  ZK_perp(1) = 1.d-5; ZK_perp(2) = 0.d0; ZK_perp(3)= 0.d0; ZK_perp(4)= 99.d0; ZK_perp(5) = 99.d0
  ZK_par     = 1.d0
  D_perp(1)  = 1.d-5; D_perp(2) = 0.d0; D_perp(3)= 0.d0; D_perp(4)= 99.d0; D_perp(5) = 99.d0
  D_par      = 0.d0

  eta_num       = 0.d0
  visco_num     = 0.d0
  visco_par_num = 0.d0
  D_perp_num    = 0.d0

  heatsource     = 1.e-7
  particlesource = 1.e-5

  zjz_0 =  0.1173d0;   T_0   =  1.d-6  ;   rho_0 =  1.d0   ;   FF_0  =  1.d0
  zjz_1 =  0.0d0   ;   T_1   =  1.d-8  ;   rho_1 =  1.d0   ;   FF_1  =  0.d0

  zj_coef     = 0.d0;  zj_coef(1)  = -1.d0
  T_coef      = 0.d0;  T_coef(1)   = -1.d0
  rho_coef    = 0.d0;  rho_coef(1) =  0.d0
  FF_coef     = 0.d0;  FF_coef(1)  = -1.d0

  t_now       = 0.d0
  t_start     = 0.d0
  index_start = 0

  nout = 9999999

  rho_file      = 'none'
  T_file        = 'none'
  ffprime_file  = 'none'

  produce_live_data = .true.
  
  ! --- Read input parameters from namelist.
  if (my_id .eq. 0) read(5,in1)
  
  tstep_in = tstep

  if (sum(nstep_n) .gt. 0) then
    nstep = sum(nstep_n)
  else
    tstep_n    = 0.d0
    tstep_n(1) = tstep
    nstep_n    = 0
    nstep_n(1) = nstep
  endif
  
  if (nstep .gt. 0) call tr_allocate(energies,1,n_tor,1,2,1,nstep,"energies")
  if (nstep .gt. 0) call tr_allocate(xtime,1,nstep,"xtime")<

  ! --- Read numerical profiles for rho, T, and ff'.
  call read_num_profiles()
  
  ! --- Determine the derivatives of the numerical input profiles.
  call derive_num_profiles()
  
endif

return
end subroutine initialise_parameters
