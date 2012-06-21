!> Presets input parameters to reasonable default values.
!! 
!! The model-specific routines initialise_parameters may overwrite
!! these defaults according to the requirements of the respective
!! model.
subroutine preset_parameters
  
  use phys_module
  use mumps_module,  only: use_mumps, no_zeros_mumps
  use murge_module,  only: use_murge, use_murge_element
  use pastix_module, only: use_pastix, no_zeros_pastix, pastix_smp_only
  use wsmp_module,   only: use_wsmp
  
  implicit none
  
  time_evol_scheme = 'Crank-Nicholson'
  
  ! --- DoubleNull flag
  xcase = LOWER_XPOINT
  
  tstep    = 1.d0
  tstep_n  = 1.d0
  nstep    = 0
  nstep_n  = 0
  
  eta_T_dependent   = .true.
  visco_T_dependent = .true.

  eta   = 1.d-5
  visco = 1.d-5
  visco_par = 1.d-5
  
  central_density = 0.d0  ! the central density in units 10^20 m^-3

  restart      = .false.
  import_equil = .false.
  regrid       = .false.

  freeboundary_equil = .false. ! use free or fixed boundary equilibrium
  freeboundary       = .false. ! use free or fixed boundary?
  resistive_wall     = .false. ! use a resistive or ideal wall?    (freeboundary only)

  bc_natural_open    = .false.! use sheath (Bohm) boundary conditions
  gamma_sheath       = 4.5d0  ! sheath transmission factor (single fluid)
  density_reflection = 0.d0   ! reflection coefficient for outgoing density
  
  n_R       = 0
  n_Z       = 0

  n_radial  = 11
  n_pol     = 16

  n_flux    = 11
  n_tht     = 16

  n_open    = 5
  n_outer   = 0
  n_inner   = 0
  n_leg     = 5
  n_private = 5
  n_up_leg  = 0
  n_up_priv = 0
  
  psi_axis_init = -0.1d0
  XR_r(:)       = 999.d0
  SIG_r(:)      = 999.d0
  XR_tht(:)     = 999.d0
  SIG_tht(:)    = 999.d0

  SIG_closed  = 0.1d0
  SIG_open    = 0.1d0
  SIG_outer   = 0.1d0
  SIG_inner   = 0.1d0
  SIG_private = 0.1d0
  SIG_up_priv = 0.1d0
  SIG_theta   = 0.03d0
  SIG_leg_0   = 0.05d0
  SIG_leg_1   = 0.2d0
  SIG_up_leg_0= 0.05d0
  SIG_up_leg_1= 0.2d0
  
  dPSI_open    = 0.11
  dPSI_outer   = 0.11
  dPSI_inner   = 0.11
  dPSI_private = 0.03
  dPSI_up_priv = 0.03
  
  R_geo     = 10.d0
  Z_geo     = 0.d0
  amin      = 1.d0

  F0        = 10.d0
  GAMMA     = 5.d0 / 3.d0

  mf        = 2
  fbnd      = 0.d0; fbnd(1) =2.d0

  R_boundary   = 0.d0
  Z_boundary   = 0.d0
  psi_boundary = 0.d0
  n_boundary   = 0

  n_pfc       = 0
  Rmin_pfc    = 0.d0
  Rmax_pfc    = 0.d0
  Zmin_pfc    = 0.d0
  Zmax_pfc    = 0.d0
  current_pfc = 0.d0

  bootstrap = .false.

  ellip  = 1.d0
  tria_u = 0.d0
  tria_l = 0.d0
  quad_l = 0.d0
  quad_u = 0.d0

  xampl  = 0.d0
  xwidth = 0.d0
  xsig   = 1.d0
  xtheta = 0.d0
  xshift = 0.d0
  xleft  = 0.d0
  xpoint = .false.
  xcase  = 1

  xr1  = 9999.d0
  sig1 = 9999.d0
  xr2  = 99999.d0
  sig2 = 99999.d0

  R_begin = -0.1d0
  R_end   =  0.1d0
  Z_begin = -0.1d0
  Z_end   = 0.1d0
  
  ZK_perp(:) = 0.d0
  ZK_perp(1) = 1.d-5; ZK_perp(2) = 0.d0; ZK_perp(3)= 0.d0; ZK_perp(4)= 99.d0; ZK_perp(5) = 99.d0
  ZK_par     = 1.d0
  D_perp(:)  = 0.d0
  D_perp(1)  = 1.d-5; D_perp(2) = 0.d0; D_perp(3)= 0.d0; D_perp(4)= 99.d0; D_perp(5) = 99.d0
  D_par      = 0.d0

  eta_num       = 0.d0
  visco_num     = 0.d0
  visco_par_num = 0.d0
  D_perp_num    = 0.d0
  ZK_perp_num   = 0.d0

  heatsource          = 1.e-7
  heatsource_psin     = 1.0d0
  heatsource_sig      = 0.1d0
  particlesource      = 1.e-5
  particlesource_psin = 1.0d0
  particlesource_sig  = 0.1d0
  
  tauIC = 0.d0

  zjz_0 =  0.1173d0;   T_0   =  1.d-6  ;   rho_0 =  1.d0   ;   FF_0  =  1.d0
  zjz_1 =  0.0d0   ;   T_1   =  1.d-8  ;   rho_1 =  1.d0   ;   FF_1  =  0.d0

  zj_coef     = 0.d0;  zj_coef(1)  = -1.d0
  T_coef      = 0.d0;  T_coef(1)   = -1.d0
  rho_coef    = 0.d0;  rho_coef(1) =  0.d0
  FF_coef     = 0.d0;  FF_coef(1)  = -1.d0

  pellet_amplitude  = 0.d0
  pellet_R          = 3.8d0
  pellet_Z          = 0.0d0
  pellet_phi        = 1.57d0
  pellet_radius     = 0.08d0
  pellet_sig        = 0.02
  pellet_length     = 0.785
  pellet_psi        = 1.0d0
  pellet_delta_psi  = 999.d0
  pellet_velocity_R = 0.d0
  pellet_velocity_Z = 0.d0
  pellet_particles  = 0.d0  
  pellet_density    = 3.d8       ! pellet density (in units 10^20 m^-3)
  use_pellet        = .false.

  t_now       = 0.d0
  t_start     = 0.d0
  index_start = 0

  nout = 9999999

  rho_file      = 'none'
  rhon_file     = 'none'
  T_file        = 'none'
  Te_file       = 'none'
  Ti_file       = 'none'
  ffprime_file  = 'none'
  d_perp_file   = 'none'
  zk_perp_file  = 'none'
  R_Z_psi_bnd_file = 'none'

  produce_live_data = .true.
  
  linear_run         = .false.
  
  export_for_nemec      = .false.
#ifdef USE_HDF5
  save_diagnostics_HDF5 = .false.
  h5_diag_nbtime        = 10.d0
#endif
  
  gmres              = .true.               ! Use iterative solver
  gmres_max_iter     = 200                  ! Max number of GMRES iterations
  use_mumps          = .false.              ! Use MUMPS solver
  use_pastix         = .true.               ! Use PASTIX solver
  use_murge          = .false.              ! Use MURGE interface to PASTIX solver
  use_murge_element  = .false.              ! Build the matrix through murge, not with a CSC.
  use_wsmp           = .false.              ! Use WSMP solver (use with care, still in development!)
  
  refinement         = .false.              ! enable mesh refinement
  
  grid_to_wall       = .false.              ! extend the grid to a physical wall
  
  adaptive_time      = .false.              ! requires no_mpi for Pastix library
  
  equil              = .true.               ! compute equilibrium
  
  bench_without_plot = .false.              ! .true. for benchmark (mesuring elapsed time without plot phases) 
  no_zeros_pastix    = .false.              ! .true. to remove nonzeros in the preconditioning matrix with MUMPS
  no_zeros_mumps     = .false.              ! .true. to remove nonzeros in the preconditioning matrix with PaStiX
  
!==== RMP parameters =====
  RMP_on             = .false.              ! .true. to activate RMPs (changes boundary conditions)
  RMP_psi_cos_file   = 'none'
  RMP_psi_sin_file   = 'none'
  lambda=0.0663
  tset = 150
  output_bnd_elements = .false.  ! writes bnd nodes and elements in output files (boundary_nodes.dat and boundary_elements.dat)

  n_limiter = 0
  R_limiter = 0.d0
  Z_limiter = 0.d0
  
 !======================MB rotation profile
  V_0=0.d0   
  V_1=0.d0    
  V_coef=0.d0
  V_coef(1)=0.d0
  V_coef(4)=0.1
  V_coef(5)=1. 
!======================MB

end subroutine preset_parameters
