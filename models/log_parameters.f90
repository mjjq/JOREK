!> Write out all relevant parameters defined in mod_parameters
!! and by the namelist input file.
subroutine log_parameters(my_id)

use phys_module
use mumps_module,  only: use_mumps, no_zeros_mumps
use murge_module,  only: use_murge, use_murge_element, murge_with_starpu, murge_cuda_nbr
use pastix_module, only: use_pastix, no_zeros_pastix, pastix_smp_only, pastix_pivot
use wsmp_module,   only: use_wsmp
use vacuum,        only: wall_resistivity

implicit none

! --- Routine parameters
integer, intent(in) :: my_id !< MPI proc id

! --- Constants
character(len=512), parameter :: REAL_FMT = "(1X,A, ' = ', 10ES12.4)"
character(len=512), parameter :: INTG_FMT = "(1X,A, ' = ', 10I12)"
character(len=512), parameter :: LOGI_FMT = "(1X,A, ' = ', 10L12)"
character(len=512), parameter :: REA2_FMT = "(1X,A, ' = ', 4ES12.4, '     ...    ', 4ES12.4)"
character(len=512), parameter :: REA3_FMT = "(1X,A, ' = ', 9ES12.4, '     ...')"
character(len=512), parameter :: VARI_FMT = "(3x,I3,': ',A)"
character(len=512), parameter :: MODE_FMT = "(3x,I3,': ',A,'(',A,'*phi)')"
character(len=512), parameter :: CHAR_FMT = "(1X,A, ' = ""', A, '""')"

! --- Local variables
integer           :: ivar, itor
character(len=10) :: mode_num

if (my_id == 0) then
  
  write(*,*)
  write(*,*) '*************************************************'
  write(*,*) '*          PARAMETERS OF THE JOREK RUN          *'
  write(*,*) '*************************************************'
  write(*,*)
  write(*,*) 'PREPROCESSOR OPTIONS'
  write(*,*) '-------------------------------------------------'
  write(*,'(1x,a)',advance='no') ' USE_MUMPS           : '
#ifdef USE_MUMPS
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,*) ' JOREK_MODEL         : ', JOREK_MODEL

  write(*,'(1x,a)',advance='no') ' USE_FFTW            : '
#ifdef USE_FFTW
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif
  
  write(*,'(1x,a)',advance='no') ' USE_PASTIX          : '
#ifdef USE_PASTIX
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif
  
  write(*,'(1x,a)',advance='no') ' USE_WSMP            : '
#ifdef USE_WSMP
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif
  
  write(*,'(1x,a)',advance='no') ' USE_MURGE           : '
#ifdef USE_MURGE
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif
  
  write(*,'(1x,a)',advance='no') ' USE_HIPS            : '
#ifdef USE_HIPS
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif
  
  write(*,'(1x,a)',advance='no') ' USE_BLOCK           : '
#ifdef USE_BLOCK
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif
  
  write(*,'(1x,a)',advance='no') ' USE_HDF5            : '
#ifdef USE_HDF5 
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif
  
  write(*,'(1x,a)',advance='no') ' MEMTRACE            : '
#ifdef MEMTRACE 
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' JECCD               : '
#ifdef JECCD
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' COMPARE_ELEMENT_MATRIX : '
#ifdef COMPARE_ELEMENT_MATRIX
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,*)
  write(*,*) 'HARD-CODED PARAMETERS'
  write(*,*) '-------------------------------------------------'
  write(*,INTG_FMT) 'jorek_model           ', jorek_model
  write(*,INTG_FMT) 'n_var                 ', n_var
  do ivar = 1, n_var
    write(*,VARI_FMT) ivar, trim(variable_names(ivar))
  end do
  write(*,INTG_FMT) 'n_dim                 ', n_dim
  write(*,INTG_FMT) 'n_order               ', n_order
  write(*,INTG_FMT) 'n_tor                 ', n_tor
  write(*,MODE_FMT) 1, 'cos', '0'
  do itor = 2, n_tor
    write(mode_num,'(I4)') mode(itor)
    write(*,MODE_FMT) itor, mode_type(itor), trim(adjustl(mode_num))
  end do
  write(*,INTG_FMT) 'n_period              ', n_period
  write(*,INTG_FMT) 'n_plane               ', n_plane
  write(*,INTG_FMT) 'n_vertex_max          ', n_vertex_max
  write(*,INTG_FMT) 'n_nodes_max           ', n_nodes_max
  write(*,INTG_FMT) 'n_elements_max        ', n_elements_max
  write(*,INTG_FMT) 'n_boundary_max        ', n_boundary_max
  write(*,INTG_FMT) 'n_pieces_max          ', n_pieces_max
  write(*,INTG_FMT) 'n_degrees             ', n_degrees
  write(*,*)
  write(*,*) 'NAMELIST INPUT PARAMETERS'
  write(*,*) '-------------------------------------------------'
  write(*,CHAR_FMT) 'time_evol_scheme      ', trim(time_evol_scheme)
  write(*,INTG_FMT) 'n_tor_fft_thresh      ', n_tor_fft_thresh
  if ( n_tor .ge. n_tor_fft_thresh ) then
    write(*,'(3x,a)') '=> fft version of element_matrix is used'
  else
    write(*,'(3x,a)') '=> no-fft version of element_matrix is used'
  end if
  write(*,REAL_FMT) 'tstep                 ', tstep
  write(*,INTG_FMT) 'nstep                 ', nstep
  write(*,REAL_FMT) 'tstep_n               ', tstep_n
  write(*,INTG_FMT) 'nstep_n               ', nstep_n
  write(*,LOGI_FMT) 'eta_T_dependent       ', eta_T_dependent
  write(*,REAL_FMT) 'eta                   ', eta
  write(*,LOGI_FMT) 'visco_T_dependent     ', visco_T_dependent
  write(*,REAL_FMT) 'visco                 ', visco
  write(*,REAL_FMT) 'visco_par             ', visco_par
  write(*,LOGI_FMT) 'restart               ', restart
  write(*,INTG_FMT) 'rst_format            ', rst_format
  write(*,INTG_FMT) 'rst_hdf5              ', rst_hdf5
  write(*,LOGI_FMT) 'regrid                ', regrid
  write(*,INTG_FMT) 'n_R                   ', n_R
  write(*,INTG_FMT) 'n_Z                   ', n_Z
  write(*,INTG_FMT) 'n_radial              ', n_radial
  write(*,INTG_FMT) 'n_pol                 ', n_pol
  
  if ( n_radial > 0 ) then
    write(*,REAL_FMT) 'psi_axis_init         ', psi_axis_init
    write(*,REAL_FMT) 'xr_r                  ', xr_r(:)
    write(*,REAL_FMT) 'sig_r                 ', sig_r(:)
    write(*,REAL_FMT) 'xr_tht                ', xr_tht(:)
    write(*,REAL_FMT) 'sig_tht               ', sig_tht(:)
  end if
  
  write(*,INTG_FMT) 'n_tht                 ', n_tht
  write(*,INTG_FMT) 'n_flux                ', n_flux
  write(*,LOGI_FMT) 'xpoint                ', xpoint
  
  if ( xpoint ) then
    write(*,INTG_FMT) 'n_open                ', n_open
    write(*,INTG_FMT) 'n_private             ', n_private
    write(*,INTG_FMT) 'n_leg                 ', n_leg
    write(*,REAL_FMT) 'SIG_closed            ', SIG_closed
    write(*,REAL_FMT) 'SIG_open              ', SIG_open
    write(*,REAL_FMT) 'SIG_private           ', SIG_private
    write(*,REAL_FMT) 'SIG_theta             ', SIG_theta
    write(*,REAL_FMT) 'SIG_leg_0             ', SIG_leg_0
    write(*,REAL_FMT) 'SIG_leg_1             ', SIG_leg_1
    write(*,REAL_FMT) 'dPSI_open             ', dPSI_open
    write(*,REAL_FMT) 'dPSI_private          ', dPSI_private
    write(*,INTG_FMT) 'xcase                 ', xcase
  end if
  
  write(*,INTG_FMT) 'nout                  ', nout
  write(*,REAL_FMT) 'xr1                   ', xr1
  write(*,REAL_FMT) 'sig1                  ', sig1
  write(*,REAL_FMT) 'xr2                   ', xr2
  write(*,REAL_FMT) 'sig2                  ', sig2
  write(*,REAL_FMT) 'R_begin               ', R_begin
  write(*,REAL_FMT) 'R_end                 ', R_end
  write(*,REAL_FMT) 'Z_begin               ', Z_begin
  write(*,REAL_FMT) 'Z_end                 ', Z_end
  write(*,REAL_FMT) 'R_geo                 ', R_geo
  write(*,REAL_FMT) 'Z_geo                 ', Z_geo
  write(*,REAL_FMT) 'amin                  ', amin
  write(*,INTG_FMT) 'mf                    ', mf
  
  if ( mf >= 0 ) then
    write(*,REA3_FMT) 'fbnd                  ', fbnd(1:MIN(9,mf))
    write(*,REA3_FMT) 'fpsi                  ', fpsi(1:MIN(9,mf))
  end if
  
  write(*,REAL_FMT) 'F0                    ', F0
  write(*,REAL_FMT) 'zjz_0                 ', zjz_0
  write(*,REAL_FMT) 'zjz_1                 ', zjz_1
  write(*,REAL_FMT) 'zj_coef               ', zj_coef
  
  if ( .not. num_rho ) then
    write(*,REAL_FMT) 'rho_0                 ', rho_0
    write(*,REAL_FMT) 'rho_1                 ', rho_1
    write(*,REAL_FMT) 'rho_coef              ', rho_coef(1:5)
  else
    write(*,CHAR_FMT) 'rho_file              ', trim(rho_file)
  end if
  
  if ( num_rot ) then
    write(*,CHAR_FMT) 'rot_file              ', trim(rot_file)
  else
    write(*,REAL_FMT) 'V_0                   ', V_0
    write(*,REAL_FMT) 'V_1                   ', V_1
    write(*,REAL_FMT) 'V_coeff               ', V_coef(1:10)
  end if

  if ( (abs(V_0) .ge. 1.d-19) .or. (num_rot) ) then  
     write(*,LOGI_FMT) 'normalized_velocity_profile', normalized_velocity_profile
  endif

  if ( .not. num_T ) then
    write(*,REAL_FMT) 'T_0                   ', T_0
    write(*,REAL_FMT) 'T_1                   ', T_1
    write(*,REAL_FMT) 'T_coef                ', T_coef(1:5)
  else
    write(*,CHAR_FMT) 'T_file                ', trim(T_file)
  end if
  
  if ( jorek_model == 400 ) then
    write(*,REAL_FMT) 'Te_0                   ', Te_0
    write(*,REAL_FMT) 'Te_1                   ', Te_1
    write(*,REAL_FMT) 'Te_coef                ', Te_coef(1:5)
    write(*,REAL_FMT) 'Ti_0                   ', Ti_0
    write(*,REAL_FMT) 'Ti_1                   ', Ti_1
    write(*,REAL_FMT) 'Ti_coef                ', Ti_coef(1:5)
    if ( .not. num_zk_e_perp ) then
      write(*,REAL_FMT) 'ZK_e_perp             ', ZK_e_perp(1:6)
    else
      write(*,CHAR_FMT) 'ZK_e_perp_file        ', trim(ZK_e_perp_file)
    end if
    if ( .not. num_zk_i_perp ) then
      write(*,REAL_FMT) 'ZK_i_perp             ', ZK_i_perp(1:6)
    else
      write(*,CHAR_FMT) 'ZK_i_perp_file        ', trim(ZK_i_perp_file)
    end if
    write(*,REAL_FMT) 'heatsource_e           ', heatsource_e
    write(*,REAL_FMT) 'heatsource_i           ', heatsource_i
    write(*,REAL_FMT) 'K_e_par                ', K_e_par
    write(*,REAL_FMT) 'K_i_par                ', K_i_par
  end if
  
  if ( .not. num_ffprime ) then
    write(*,REAL_FMT) 'FF_0                  ', FF_0
    write(*,REAL_FMT) 'FF_1                  ', FF_1
    write(*,REAL_FMT) 'FF_coef               ', FF_coef(1:8)
  else
    write(*,CHAR_FMT) 'ffprime_file          ', trim(ffprime_file)
  end if
  
  write(*,REAL_FMT) 'ZK_par                ', ZK_par
  write(*,REAL_FMT) 'ZK_par_max            ', ZK_par_max
  write(*,LOGI_FMT) 'ZKpar_T_dependent     ', ZKpar_T_dependent
  if ( .not. num_zk_perp ) then
    write(*,REAL_FMT) 'ZK_perp               ', ZK_perp(1:6)
  else
    write(*,CHAR_FMT) 'ZK_perp_file          ', trim(ZK_perp_file)
  end if
  write(*,REAL_FMT) 'D_par                 ', D_par
  if ( .not. num_d_perp ) then
    write(*,REAL_FMT) 'D_perp                ', D_perp(1:6)
  else
    write(*,CHAR_FMT) 'D_perp_file           ', trim(D_perp_file)
  end if
  write(*,REAL_FMT) 'particlesource        ', particlesource
  write(*,REAL_FMT) 'particlesource_psin   ', particlesource_psin
  write(*,REAL_FMT) 'particlesource_sig    ', particlesource_sig
  write(*,REAL_FMT) 'heatsource            ', heatsource
  write(*,REAL_FMT) 'heatsource_psin       ', heatsource_psin
  write(*,REAL_FMT) 'heatsource_sig        ', heatsource_sig
  write(*,REAL_FMT) 'tauIC                 ', tauIC
  write(*,REAL_FMT) 'eta_num               ', eta_num
  write(*,REAL_FMT) 'visco_num             ', visco_num
  write(*,REAL_FMT) 'visco_par_num         ', visco_par_num
  write(*,REAL_FMT) 'D_perp_num            ', D_perp_num
  write(*,REAL_FMT) 'ZK_perp_num           ', ZK_perp_num
  write(*,REAL_FMT) 'tgnum                 ', tgnum(:)
  write(*,REAL_FMT) 'D_prof_neg            ', D_prof_neg
  write(*,REAL_FMT) 'ZK_prof_neg           ', ZK_prof_neg
  write(*,REAL_FMT) 'T_min                 ', T_min
  write(*,LOGI_FMT) 'use_pellet            ', use_pellet
  write(*,REAL_FMT) 'corr_neg_temp_coef    ', corr_neg_temp_coef(:)
  write(*,REAL_FMT) 'corr_neg_dens_coef    ', corr_neg_dens_coef(:)
  
  if ( use_pellet ) then
    write(*,REAL_FMT) 'pellet_amplitude    ', pellet_amplitude
    write(*,REAL_FMT) 'pellet_R              ', pellet_R
    write(*,REAL_FMT) 'pellet_Z              ', pellet_Z
    write(*,REAL_FMT) 'pellet_phi            ', pellet_phi
    write(*,REAL_FMT) 'pellet_radius         ', pellet_radius
    write(*,REAL_FMT) 'pellet_sig            ', pellet_sig
    write(*,REAL_FMT) 'pellet_length         ', pellet_length
    write(*,REAL_FMT) 'pellet_psi            ', pellet_psi
    write(*,REAL_FMT) 'pellet_ellipse        ', pellet_ellipse
    write(*,REAL_FMT) 'pellet_theta          ', pellet_theta
    write(*,REAL_FMT) 'pellet_delta_psi      ', pellet_delta_psi
    write(*,REAL_FMT) 'pellet_density        ', pellet_density
    write(*,REAL_FMT) 'pellet_particles      ', pellet_particles
    write(*,REAL_FMT) 'pellet_velocity_R     ', pellet_velocity_R
    write(*,REAL_FMT) 'pellet_velocity_Z     ', pellet_velocity_Z
  end if
  
  write(*,*)
  write(*,REAL_FMT) 'ellip                 ', ellip
  write(*,REAL_FMT) 'tria_u                ', tria_u
  write(*,REAL_FMT) 'tria_l                ', tria_l
  write(*,REAL_FMT) 'quad_u                ', quad_u
  write(*,REAL_FMT) 'quad_l                ', quad_l
  write(*,REAL_FMT) 'xampl                 ', xampl
  write(*,REAL_FMT) 'xwidth                ', xwidth
  write(*,REAL_FMT) 'xsig                  ', xsig
  write(*,REAL_FMT) 'xtheta                ', xtheta
  write(*,REAL_FMT) 'xshift                ', xshift
  write(*,REAL_FMT) 'xleft                 ', xleft
  write(*,CHAR_FMT) 'R_Z_psi_bnd_file      ', trim(R_Z_psi_bnd_file)
  write(*,CHAR_FMT) 'wall_file             ', trim(wall_file)
  write(*,INTG_FMT) 'n_boundary            ', n_boundary
  if ( n_boundary > 0 ) then
    write(*,REA2_FMT) 'r_boundary            ', r_boundary(1:4), r_boundary(n_boundary-3:n_boundary)
    write(*,REA2_FMT) 'z_boundary            ', z_boundary(1:4), z_boundary(n_boundary-3:n_boundary)
    write(*,REA2_FMT) 'psi_boundary          ', psi_boundary(1:4), psi_boundary(n_boundary-3:n_boundary)
  end if
  write(*,INTG_FMT) 'n_limiter             ', n_limiter
  if ( n_limiter > 0 ) then
    write(*,REA3_FMT) 'r_limiter             ', r_limiter(1:min(9,n_limiter))
    write(*,REA3_FMT) 'z_limiter             ', z_limiter(1:min(9,n_limiter))
  end if
  
  write(*,LOGI_FMT) 'freeboundary_equil    ', freeboundary_equil
  write(*,LOGI_FMT) 'freeboundary          ', freeboundary
  if ( freeboundary ) then
    write(*,LOGI_FMT) 'resistive_wall        ', resistive_wall
    if ( resistive_wall ) then
      write(*,REAL_FMT) 'wall_resistivity      ', wall_resistivity
    end if
  end if
  
  write(*,REAL_FMT) 'Q_bar                 ', Q_bar
  write(*,REAL_FMT) 'sigma                 ', sigma
  write(*,REAL_FMT) 'density_reflection    ', density_reflection
  write(*,REAL_FMT) 'central_density       ', central_density
  write(*,REAL_FMT) 'central_mass          ', central_mass
  write(*,REAL_FMT) 'gamma_sheath          ', gamma_sheath
  write(*,LOGI_FMT) 'bc_natural_open       ', bc_natural_open
  write(*,LOGI_FMT) 'produce_live_data     ', produce_live_data
  write(*,LOGI_FMT) 'export_for_nemec      ', export_for_nemec
#ifdef USE_HDF5
  write(*,LOGI_FMT) 'save_diagnostics_HDF5 ', save_diagnostics_HDF5
  write(*,REAL_FMT) 'h5_diag_nbtime        ', h5_diag_nbtime
!  write(*,LOGI_FMT) 'h5_nbsave_all         ', h5_nbsave_all
#endif
  write(*,LOGI_FMT) 'linear_run            ', linear_run
  write(*,LOGI_FMT) 'gmres                 ', gmres
  write(*,INTG_FMT) 'gmres_max_iter        ', gmres_max_iter
  write(*,REAL_FMT) 'gmres tolerance       ', gmres_tol
  write(*,INTG_FMT) 'iter_precon           ', iter_precon
  write(*,INTG_FMT) 'gmres_m               ', gmres_m
  write(*,REAL_FMT) 'gmres_4               ', gmres_4
  write(*,LOGI_FMT) 'use_mumps             ', use_mumps
  write(*,LOGI_FMT) 'use_wsmp              ', use_wsmp
  write(*,LOGI_FMT) 'use_pastix            ', use_pastix
  write(*,LOGI_FMT) 'use_murge             ', use_murge
  write(*,LOGI_FMT) 'use_murge_element     ', use_murge_element
  write(*,LOGI_FMT) 'murge_with_starpu     ', murge_with_starpu
  write(*,INTG_FMT) 'murge_cuda_nbr        ', murge_cuda_nbr
  write(*,LOGI_FMT) 'pastix_smp_only       ', pastix_smp_only
  write(*,REAL_FMT) 'pastix_pivot          ', pastix_pivot
  write(*,LOGI_FMT) 'refinement            ', refinement
  write(*,LOGI_FMT) 'grid_to_wall          ', grid_to_wall
  write(*,LOGI_FMT) 'adaptive_time         ', adaptive_time
  write(*,LOGI_FMT) 'equil                 ', equil
  write(*,LOGI_FMT) 'bench_without_plot    ', bench_without_plot
  write(*,LOGI_FMT) 'no_zeros_mumps        ', no_zeros_mumps
  write(*,LOGI_FMT) 'no_zeros_pastix       ', no_zeros_pastix
  
  write(*,INTG_FMT) 'n_pfc                 ', n_pfc
  if ( n_pfc > 0 ) then
    write(*,REA3_FMT) 'Rmin_pfc              ', Rmin_pfc(1:min(9,n_pfc))
    write(*,REA3_FMT) 'Rmax_pfc              ', Rmax_pfc(1:min(9,n_pfc))
    write(*,REA3_FMT) 'Zmin_pfc              ', Zmin_pfc(1:min(9,n_pfc))
    write(*,REA3_FMT) 'Zmax_pfc              ', Zmax_pfc(1:min(9,n_pfc))
    write(*,REA3_FMT) 'current_pfc           ', current_pfc(1:min(9,n_pfc))
  end if

  write(*,LOGI_FMT) 'RMP_on                ', RMP_on
  if (RMP_on) then
     write(*,CHAR_FMT) 'RMP_psi_cos_file      ', trim(RMP_psi_cos_file)
     write(*,CHAR_FMT) 'RMP_psi_sin_file      ', trim(RMP_psi_sin_file)
     write(*,REAL_FMT) 'lambda                ', lambda
     write(*,REAL_FMT) 'tset                  ', tset
  endif
  write(*,LOGI_FMT) 'output_bnd_elements   ', output_bnd_elements
  write(*,LOGI_FMT) 'bootstrap             ', bootstrap
  write(*,LOGI_FMT) 'NEO                   ', NEO
  if (NEO) then
    write(*,LOGI_FMT) 'num_neo_file          ', num_neo_file
    if (num_neo_file) then
      write(*,CHAR_FMT) 'neo_file              ', trim(neo_file)
    else
      write(*,REAL_FMT) 'amu_neo_const         ', amu_neo_const
      write(*,REAL_FMT) 'aki_neo_const         ', aki_neo_const        
    endif
  endif

  if(jorek_model == 306 ) then
     write(*,REAL_FMT) 'nu_jec1_fast        ',  nu_jec1_fast
     write(*,REAL_FMT) 'nu_jec2_fast        ',  nu_jec2_fast
     write(*,REAL_FMT) 'JJ_par              ',  JJ_par
     write(*,REAL_FMT) 'jec_pos1            ',  jec_pos1
     write(*,REAL_FMT) 'jec_width           ',  jec_width
     write(*,REAL_FMT) 'jecamp              ',  jecamp
  endif

#if (JOREK_MODEL == 500) || (JOREK_MODEL == 501) || (JOREK_MODEL == 555)
     write(*,REAL_FMT) 'mgi_amplitude       ',  mgi_amplitude
     write(*,REAL_FMT) 'mgi_R               ',  mgi_R
     write(*,REAL_FMT) 'mgi_Z               ',  mgi_Z
     write(*,REAL_FMT) 'mgi_phi             ',  mgi_phi
     write(*,REAL_FMT) 'mgi_radius          ',  mgi_radius
     write(*,REAL_FMT) 'mgi_sig             ',  mgi_sig
     write(*,REAL_FMT) 'mgi_deltaphi        ',  mgi_deltaphi
     write(*,REAL_FMT) 'ksi_ion             ',  ksi_ion
     write(*,LOGI_FMT) 'JET_MGI             ',  JET_MGI
     write(*,LOGI_FMT) 'ASDEX_MGI           ',  ASDEX_MGI
     write(*,LOGI_FMT) 'gas_type            ',  gas_type
     write(*,REAL_FMT) 'A_Dmv               ',  A_Dmv
     write(*,REAL_FMT) 'K_Dmv               ',  K_Dmv
     write(*,REAL_FMT) 'V_Dmv               ',  V_Dmv
     write(*,REAL_FMT) 'L_tube              ',  L_tube
     write(*,REAL_FMT) 't_mgi               ',  t_mgi
     write(*,REAL_FMT) 'delta_n_convection  ',  delta_n_convection
     write(*,REAL_FMT) 'nimp_bg             ',  nimp_bg
#endif

  write(*,*)
  
end if

end subroutine log_parameters
