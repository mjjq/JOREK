!> Initialise input parameters and read the input namelist
subroutine initialise_parameters(my_id, filename)

use tr_module
use phys_module
use pellet_module
use mumps_module,  only: use_mumps, no_zeros_mumps
use murge_module,  only: use_murge, use_murge_element, murge_with_starpu, murge_cuda_nbr
use pastix_module, only: use_pastix, no_zeros_pastix, pastix_smp_only, pastix_pivot
use vacuum
use wsmp_module,   only: use_wsmp

implicit none

! --- Routine parameters
integer,                      intent(in) :: my_id
character(len=*),             intent(in) :: filename

! --- Local variables
integer :: ierr,err,i

! --- Namelist with input parameters.
namelist /in1/  tstep, nstep, tstep_n, nstep_n,                     &
                rst_hdf5,                                           &
                eta, visco, visco_par,                              &
                restart, rst_format, regrid, bootstrap,             &
                force_horizontal_Xline,                             &
                n_R, n_Z, n_radial, n_pol, n_tht, n_flux,           &
                n_open, n_private, n_leg, n_ext,                    &
                n_outer, n_inner, n_up_priv, n_up_leg,              &
                psi_axis_init, XR_r, SIG_r, XR_tht, SIG_tht,        &
                SIG_closed, SIG_open, SIG_private, SIG_theta,       &
                SIG_leg_0, SIG_leg_1, dPSI_open, dPSI_private,      &
                SIG_up_leg_0, SIG_up_leg_1, SIG_up_priv,            &
                SIG_outer, SIG_inner,                               &
                dPSI_outer, dPSI_inner, dPSI_up_priv,               &
                nout, xr1, sig1, xr2, sig2,                         &
                R_begin, R_end, Z_begin, Z_end,                     &
                R_geo, Z_geo, amin, mf, fbnd, fpsi, mode,           &
                R_Z_psi_bnd_file,                                   &
                R_boundary, Z_boundary, psi_boundary, n_boundary,   &
                n_pfc, n_tor_fft_thresh,                            &
                Rmin_pfc, Rmax_pfc, Zmin_pfc, Zmax_pfc, current_pfc,&
                tokamak_device,                                     &
                F0, gamma_sheath, density_reflection,               &
                zjz_0, zjz_1, zj_coef,                              &
                rho_0, rho_1, rho_coef,                             &
                T_0,   T_1,   T_coef,                               &
                FF_0,  FF_1,  FF_coef,                              &
                ZK_par, ZK_par_max, ZK_perp, D_par, D_perp,         &
                particlesource, heatsource, tauIC,                  &
                eta_num, visco_num, visco_par_num, D_perp_num,      &
                ZK_perp_num,                                        &
                pellet_amplitude, pellet_R, pellet_Z, pellet_phi,   &
                pellet_radius, pellet_sig, pellet_length,           &
                pellet_psi, pellet_delta_psi, pellet_density,       &
                pellet_velocity_R, pellet_velocity_Z, pellet_theta, &
                pellet_ellipse,                                     &
                central_density, central_mass,                      &
		pellet_particles, use_pellet,                       &
                ellip,tria_u,tria_l,quad_u,quad_l,                  &
                xampl,xwidth,xsig,xtheta,xshift,xleft, xpoint,      &
                xcase, D_perp_file, ZK_perp_file,                   &
                rho_file, T_file, ffprime_file, rot_file,           &
                normalized_velocity_profile,                        &
                freeboundary_equil, freeboundary,                   &
                resistive_wall,                                     &
                wall_resistivity, wall_resistivity_fact,            &
                bc_natural_open,                                    &
                use_mumps, use_pastix, use_murge, use_murge_element,&
                murge_with_starpu, murge_cuda_nbr,                  &
                use_wsmp,                                           &
                pastix_smp_only, refinement, grid_to_wall,          &
                adaptive_time, equil, bench_without_plot,           &
                no_zeros_pastix, no_zeros_mumps,                    &
                eta_T_dependent, visco_T_dependent,                 &
                zkpar_T_dependent,                                  & 
                heatsource_psin, heatsource_sig,                    &
                particlesource_psin, particlesource_sig,            &
                edgeparticlesource, edgeparticlesource_psin,        &
                edgeparticlesource_sig,                             &
                particlesource_gauss, heatsource_gauss,             &
                heatsource_gauss_psin, heatsource_gauss_sig,        &
                particlesource_gauss_psin, particlesource_gauss_sig,&
                produce_live_data, gmres, gmres_max_iter,           &
                gmres_m, gmres_4, gmres_tol, iter_precon,           &
                tgnum,  pastix_pivot,                               &
                linear_run, export_for_nemec,                       &
#ifdef USE_HDF5
                save_diagnostics_HDF5,h5_diag_nbtime,               &
#endif
                RMP_on, RMP_har_cos,RMP_har_sin,                    &
                RMP_growth_rate, RMP_ramp_up_time,                  &
                RMP_psi_cos_file, RMP_psi_sin_file,                 &
                V_0,V_1,V_coef, output_bnd_elements,                &
                wall_file,                                          &
                n_limiter, R_limiter, Z_limiter,                    &
                first_target_point, last_target_point,		    &
                NEO, neo_file, aki_neo_const, amu_neo_const,        &
                time_evol_scheme,                                   &
                corr_neg_temp_coef, corr_neg_dens_coef,             &
                D_prof_neg, ZK_prof_neg, T_min,                     &
                amix, amix_freeb, equil_accuracy,                   &
                equil_accuracy_freeb, current_ref, FB_Ip_position,  &
                FB_Ip_integral, Z_axis_ref, FB_Zaxis_position,      &
                FB_Zaxis_derivative,FB_Zaxis_integral, start_VFB,   &
                n_feedback_current, n_feedback_vertical,            &
                n_iter_freeb, n_polcoils_nml, polcoils0,                  &
                Zaxis_find_limit, PF_pert_start_time,               &
                starwall_equil_coils

if (my_id .eq. 0) then

  ! --- Preset input parameters to reasonable default values.
  call preset_parameters()
  call vacuum_preset(my_id, freeboundary_equil, freeboundary, resistive_wall)
  
  ! --- Model-specific presets
  particlesource_psin = 100.d0
  
  ! --- Read input parameters from namelist.
  if (trim(filename) .ne. "__NO_FILENAME__" ) then
     open(42, file=filename, status='old', action='read', iostat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'ERROR: COULD NOT OPEN NAMELIST FILE "', trim(filename), '".'
      stop
    end if
    read(42,in1)
    close(42)
  else

    read(5,in1)

 endif
!    write(*,*) 'Input files: T_file = ',  trim(T_file), ',  rho_file = ', trim(rho_file) 
!    write(*,*) 'ffprime_file = ', trim(ffprime_file),  ',  R_Z_psi_bnd_file = ', trim(R_Z_psi_bnd_file)
!    if (NEO) then
!       write(*,*) 'neo_file = ', trim(neo_file)
!    endif

 !==============================R_Z_psi_bnd==========================
   if ( (n_boundary.ne.0) .and. (R_Z_psi_bnd_file /= 'none') ) then
 ! --- Open the file.
    OPEN(UNIT=243, FILE=R_Z_psi_bnd_file, FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
    if ( err /= 0 ) then
      write(*,*) 'ERROR in initialise_parameters: Cannot open file '//TRIM(R_Z_psi_bnd_file)//'.'
      write(*,*) 'Assuming data is in main input file '//TRIM(filename)//'.'
    else
      write(*,'(A)') ' boundary info from R_Z_psi_bnd_file: R_boundary, Z_boundary, psi_boundary ' 
      do i=1,n_boundary
        read(243,*) R_boundary(i),Z_boundary(i),psi_boundary(i)
        write(*,*) R_boundary(i),Z_boundary(i),psi_boundary(i)  
      enddo
    endif    
    CLOSE(243)
  endif
 !=========================================
  
 !==============================Limiter==========================
   if (n_limiter.ne.0) then
 ! --- Open the file.
    OPEN(UNIT=244, FILE=wall_file, FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
    if ( err /= 0 ) then
      write(*,*) 'ERROR in initialise_parameters: Cannot open file '//TRIM(wall_file)//'.'
      write(*,*) 'Assuming data is in main input file '//TRIM(filename)//'.'
    else
      write(*,'(A)') ' wall info from wall_file: R_wall, Z_wall ' 
      do i=1,n_limiter
        read(244,*) R_limiter(i),Z_limiter(i)
        write(*,*)  R_limiter(i),Z_limiter(i)
      enddo
    endif    
    CLOSE(244)
  endif
 !=========================================
  
  if (sum(nstep_n) .gt. 0) then
    nstep = sum(nstep_n)
  else
    tstep_n    = 0.d0
    tstep_n(1) = tstep
    nstep_n    = 0
    nstep_n(1) = nstep
  endif

  if (allocated(energies)) call tr_deallocate(energies,"energies",CAT_GRID)
  if (nstep .gt. 0) call tr_allocate(energies,1,n_tor,1,2,1,nstep,"energies",CAT_GRID)

#ifdef JECCD
  if (allocated(energies2)) call tr_deallocate(energies2,"energies2",CAT_GRID)
  if (nstep .gt. 0) call tr_allocate(energies2,1,n_tor,1,2,1,nstep,"energies2",CAT_GRID)

  if (allocated(energies3)) call tr_deallocate(energies3,"energies3",CAT_GRID)
  if (nstep .gt. 0) call tr_allocate(energies3,1,n_tor,1,2,1,nstep,"energies3",CAT_GRID)
#endif

  if (allocated(xtime)) call tr_deallocate(xtime,"xtime",CAT_GRID)
  if (nstep .gt. 0) call tr_allocate(xtime,1,nstep,"xtime",CAT_GRID)

  if (allocated(xtime_pellet_R)) call tr_deallocate(xtime_pellet_R,"xtime_pellet_R",CAT_GRID)
  if (nstep .gt. 0)              call tr_allocate(xtime_pellet_R,1,nstep,"xtime_pellet_R")
  if (allocated(xtime_pellet_Z)) call tr_deallocate(xtime_pellet_Z,"xtime_pellet_Z",CAT_GRID)
  if (nstep .gt. 0)              call tr_allocate(xtime_pellet_Z,1,nstep,"xtime_pellet_Z")
  if (allocated(xtime_pellet_psi)) call tr_deallocate(xtime_pellet_psi,"xtime_pellet_psi",CAT_GRID)
  if (nstep .gt. 0)              call tr_allocate(xtime_pellet_psi,1,nstep,"xtime_pellet_psi")
  if (allocated(xtime_pellet_particles)) call tr_deallocate(xtime_pellet_particles,"xtime_pellet_particles",CAT_GRID)
  if (nstep .gt. 0)                      call tr_allocate(xtime_pellet_particles,1,nstep,"xtime_pellet_particles")
  if (allocated(xtime_phys_ablation)) call tr_deallocate(xtime_phys_ablation,"xtime_phys_ablation",CAT_GRID)
  if (nstep .gt. 0)                   call tr_allocate(xtime_phys_ablation,1,nstep,"xtime_phys_ablation")

  if (allocated(R_axis_t)) call tr_deallocate(R_axis_t,"R_axis_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(R_axis_t,1,index_start+nstep,"R_axis_t",CAT_UNKNOWN)
  
  if (allocated(Z_axis_t)) call tr_deallocate(Z_axis_t,"Z_axis_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(Z_axis_t,1,index_start+nstep,"Z_axis_t",CAT_UNKNOWN)
  
  if (allocated(psi_axis_t)) call tr_deallocate(psi_axis_t,"psi_axis_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(psi_axis_t,1,index_start+nstep,"psi_axis_t",CAT_UNKNOWN)
  
  if (allocated(current_t)) call tr_deallocate(current_t,"current_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(current_t,1,index_start+nstep,"current_t",CAT_UNKNOWN)
  
  if (allocated(beta_p_t)) call tr_deallocate(beta_p_t,"beta_p_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(beta_p_t,1,index_start+nstep,"beta_p_t",CAT_UNKNOWN)
  
  if (allocated(beta_t_t)) call tr_deallocate(beta_t_t,"beta_t_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(beta_t_t,1,index_start+nstep,"beta_t_t",CAT_UNKNOWN)
  
  if (allocated(beta_n_t)) call tr_deallocate(beta_n_t,"beta_n_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(beta_n_t,1,index_start+nstep,"beta_n_t",CAT_UNKNOWN)
  
  if (allocated(density_in_t)) call tr_deallocate(density_in_t,"density_in_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(density_in_t,1,index_start+nstep,"density_in_t",CAT_UNKNOWN)
  
  if (allocated(density_out_t)) call tr_deallocate(density_out_t,"density_out_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(density_out_t,1,index_start+nstep,"density_out_t",CAT_UNKNOWN)
  
  if (allocated(pressure_in_t)) call tr_deallocate(pressure_in_t,"pressure_in_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(pressure_in_t,1,index_start+nstep,"pressure_in_t",CAT_UNKNOWN)
  
  if (allocated(pressure_out_t)) call tr_deallocate(pressure_out_t,"pressure_out_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(pressure_out_t,1,index_start+nstep,"pressure_out_t",CAT_UNKNOWN)
  
  if (allocated(heat_src_in_t)) call tr_deallocate(heat_src_in_t,"heat_src_in_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(heat_src_in_t,1,index_start+nstep,"heat_src_in_t",CAT_UNKNOWN)
  
  if (allocated(heat_src_out_t)) call tr_deallocate(heat_src_out_t,"heat_src_out_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(heat_src_out_t,1,index_start+nstep,"heat_src_out_t",CAT_UNKNOWN)
  
  if (allocated(part_src_in_t)) call tr_deallocate(part_src_in_t,"part_src_in_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(part_src_in_t,1,index_start+nstep,"part_src_in_t",CAT_UNKNOWN)
  
  if (allocated(part_src_out_t)) call tr_deallocate(part_src_out_t,"part_src_out_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(part_src_out_t,1,index_start+nstep,"part_src_out_t",CAT_UNKNOWN)
endif

! --- Read numerical profiles for rho, T, ff', toroidal rotation and neoclassical coefficients.
call read_num_profiles(my_id)

! --- Determine the derivatives of the numerical input profiles.
call derive_num_profiles(my_id)
  
return
end subroutine initialise_parameters
