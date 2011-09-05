!> Initialise input parameters and read the input namelist
subroutine initialise_parameters(my_id)

use tr_module
use phys_module
use mumps_module,  only: use_mumps, no_zeros_mumps
use murge_module,  only: use_murge, use_murge_element
use pastix_module, only: use_pastix, no_zeros_pastix, pastix_smp_only
use vacuum,        only: vacuum_preset

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
                use_mumps, use_pastix, use_murge, use_murge_element,&
                pastix_smp_only, refinement, grid_to_wall,          &
                adaptive_time, equil, bench_without_plot,           &
                no_zeros_pastix, no_zeros_mumps,                    &
                eta_T_dependent, visco_T_dependent,                 &
                heatsource_psin, heatsource_sig,                    &
                particlesource_psin, particlesource_sig,            &
                produce_live_data, gmres, gmres_max_iter,           &
                linear_run, export_for_nemec

if (my_id .eq. 0) then

  ! --- Preset input parameters to reasonable default values.
  call preset_parameters()
  call vacuum_preset(my_id, freeboundary_equil, freeboundary, use_starwall, resistive_wall)
  
  ! --- Model-specific presets
  particlesource_psin = 100.d0
  heatsource_psin     = 0.8d0
  
  ! --- Read input parameters from namelist.
  if (my_id .eq. 0) read(5,in1)
  
  if (sum(nstep_n) .gt. 0) then
    nstep = sum(nstep_n)
  else
    tstep_n    = 0.d0
    tstep_n(1) = tstep
    nstep_n    = 0
    nstep_n(1) = nstep
  endif
  
  if (nstep .gt. 0) call tr_allocate(energies,1,n_tor,1,2,1,nstep,"energies")
  if (nstep .gt. 0) call tr_allocate(xtime,1,nstep,"xtime")

  ! --- Read numerical profiles for rho, T, and ff'.
  call read_num_profiles()
  
  ! --- Determine the derivatives of the numerical input profiles.
  call derive_num_profiles()
  
endif

return
end subroutine initialise_parameters
