!> Initialise input parameters and read the input namelist
subroutine initialise_parameters(my_id, filename)

use tr_module
use phys_module
use pellet_module
use mumps_module,  only: use_mumps, no_zeros_mumps, use_mumps_BLR, mumps_BLR_eps, mumps_ordering
use pastix_module, only: use_pastix, no_zeros_pastix, pastix_smp_only, pastix_pivot, &
    pastix_maxthrd
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
                rst_hdf5, rst_hdf5_version, keep_current_prof,      &
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
                particlesource, heatsource, tauIC, Wdia,            &
		U_sheath, renormalise,                              &
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
                freeboundary_equil, freeboundary,  freeb_change_indices, &
                resistive_wall,                                     &
                wall_resistivity, wall_resistivity_fact,            &
                bc_natural_open,                                    &
                use_mumps, use_mumps_BLR, mumps_BLR_eps, mumps_ordering, &
                use_pastix, use_murge, use_murge_element, use_wsmp, &
                pastix_smp_only, refinement, force_central_node,    &
                grid_to_wall,                                       &
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
                RMP_on, RMP_har_cos,RMP_har_sin,                    &
                RMP_growth_rate, RMP_ramp_up_time,                  &
                RMP_psi_cos_file, RMP_psi_sin_file,                 &
                V_0,V_1,V_coef, output_bnd_elements,                &
                wall_file,                                          &
                n_limiter, R_limiter, Z_limiter,                    &
                first_target_point, last_target_point,		    &
                NEO, neo_file, aki_neo_const, amu_neo_const,        &
                time_evol_scheme, corr_neg_temp_coef,               &
                corr_neg_dens_coef, D_prof_neg, ZK_prof_neg,        &
                D_prof_neg_thresh, ZK_prof_neg_thresh, T_min,       &
                amix, amix_freeb, equil_accuracy,                   &
                equil_accuracy_freeb, current_ref, FB_Ip_position,  &
                FB_Ip_integral, Z_axis_ref, FB_Zaxis_position,      &
                FB_Zaxis_derivative,FB_Zaxis_integral, start_VFB,   &
                n_feedback_current, n_feedback_vertical,            &
                n_iter_freeb, n_pf_coils, pf_coils,                 &
                Zaxis_find_limit, PF_pert_start_time,               &
                starwall_equil_coils, freeb_equil_iterate_area,     &
                psi_offset_freeb, diag_coils, rmp_coils,            &
                voltage_coils, vert_FB_amp, find_pf_coil_currents,  &
                pastix_maxthrd

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

  if (allocated(xtime_pellet_R))         call tr_deallocate(xtime_pellet_R,"xtime_pellet_R",CAT_GRID)
  if (nstep .gt. 0)                      call tr_allocate(xtime_pellet_R,1,nstep,"xtime_pellet_R")
  if (allocated(xtime_pellet_Z))         call tr_deallocate(xtime_pellet_Z,"xtime_pellet_Z",CAT_GRID)
  if (nstep .gt. 0)                      call tr_allocate(xtime_pellet_Z,1,nstep,"xtime_pellet_Z")
  if (allocated(xtime_pellet_psi))       call tr_deallocate(xtime_pellet_psi,"xtime_pellet_psi",CAT_GRID)
  if (nstep .gt. 0)                      call tr_allocate(xtime_pellet_psi,1,nstep,"xtime_pellet_psi")
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
  
  if (allocated(xtime_phys_ablation))    call tr_deallocate(xtime_phys_ablation,"xtime_phys_ablation",CAT_GRID)
  if (nstep .gt. 0)                      call tr_allocate(xtime_phys_ablation,1,nstep,"xtime_phys_ablation")

  if (allocated(E_tot_t)) call tr_deallocate(E_tot_t,"E_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(E_tot_t,1,index_start+nstep,"E_tot_t",CAT_UNKNOWN)

  if (allocated(helicity_tot_t)) call tr_deallocate(helicity_tot_t,"helicity_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(helicity_tot_t,1,index_start+nstep,"helicity_tot_t",CAT_UNKNOWN)

  if (allocated(kin_perp_tot_t)) call tr_deallocate(kin_perp_tot_t,"kin_perp_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(kin_perp_tot_t,1,index_start+nstep,"kin_perp_tot_t",CAT_UNKNOWN)

  if (allocated(kin_par_tot_t)) call tr_deallocate(kin_par_tot_t,"kin_par_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(kin_par_tot_t,1,index_start+nstep,"kin_par_tot_t",CAT_UNKNOWN)

  if (allocated(thermal_tot_t)) call tr_deallocate(thermal_tot_t,"thermal_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(thermal_tot_t,1,index_start+nstep,"thermal_tot_t",CAT_UNKNOWN)

  if (allocated(Wmag_tot_t)) call tr_deallocate(Wmag_tot_t,"Wmag_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(Wmag_tot_t,1,index_start+nstep,"Wmag_tot_t",CAT_UNKNOWN)

  if (allocated(ohmic_tot_t)) call tr_deallocate(ohmic_tot_t,"ohmic_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(ohmic_tot_t,1,index_start+nstep,"ohmic_tot_t",CAT_UNKNOWN)

  if (allocated(Magwork_tot_t)) call tr_deallocate(Magwork_tot_t,"Magwork_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(Magwork_tot_t,1,index_start+nstep,"Magwork_tot_t",CAT_UNKNOWN)

  if (allocated(Ip_tot_t)) call tr_deallocate(Ip_tot_t,"Ip_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(Ip_tot_t,1,index_start+nstep,"Ip_tot_t",CAT_UNKNOWN)

  if (allocated(flux_Pvn_t)) call tr_deallocate(flux_Pvn_t,"flux_Pvn_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(flux_Pvn_t,1,index_start+nstep,"flux_Pvn_t",CAT_UNKNOWN)

  if (allocated(flux_qpar_t)) call tr_deallocate(flux_qpar_t,"flux_qpar_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(flux_qpar_t,1,index_start+nstep,"flux_qpar_t",CAT_UNKNOWN)

  if (allocated(flux_qperp_t)) call tr_deallocate(flux_qperp_t,"flux_qperp_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(flux_qperp_t,1,index_start+nstep,"flux_qperp_t",CAT_UNKNOWN)

  if (allocated(flux_kinpar_t)) call tr_deallocate(flux_kinpar_t,"flux_kinpar_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(flux_kinpar_t,1,index_start+nstep,"flux_kinpar_t",CAT_UNKNOWN)

  if (allocated(dE_tot_dt)) call tr_deallocate(dE_tot_dt,"dE_tot_dt",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(dE_tot_dt,1,index_start+nstep,"dE_tot_dt",CAT_UNKNOWN)

  if (allocated(dWmag_tot_dt)) call tr_deallocate(dWmag_tot_dt,"dWmag_tot_dt",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(dWmag_tot_dt,1,index_start+nstep,"dWmag_tot_dt",CAT_UNKNOWN)

  if (allocated(dthermal_tot_dt)) call tr_deallocate(dthermal_tot_dt,"dthermal_tot_dt",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(dthermal_tot_dt,1,index_start+nstep,"dthermal_tot_dt",CAT_UNKNOWN)

  if (allocated(dkinperp_tot_dt)) call tr_deallocate(dkinperp_tot_dt,"dkinperp_tot_dt",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(dkinperp_tot_dt,1,index_start+nstep,"dkinperp_tot_dt",CAT_UNKNOWN)

  if (allocated(dkinpar_tot_dt)) call tr_deallocate(dkinpar_tot_dt,"dkinpar_tot_dt",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(dkinpar_tot_dt,1,index_start+nstep,"dkinpar_tot_dt",CAT_UNKNOWN)

  if (allocated(thmwork_tot_t)) call tr_deallocate(thmwork_tot_t,"thmwork_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(thmwork_tot_t,1,index_start+nstep,"thmwork_tot_t",CAT_UNKNOWN)

  if (allocated(viscopar_dissip_tot_t)) call tr_deallocate(viscopar_dissip_tot_t,"viscopar_dissip_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(viscopar_dissip_tot_t,1,index_start+nstep,"viscopar_dissip_tot_t",CAT_UNKNOWN)

  if (allocated(viscopar_flux_t)) call tr_deallocate(viscopar_flux_t,"viscopar_flux_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(viscopar_flux_t,1,index_start+nstep,"viscopar_flux_t",CAT_UNKNOWN)

  if (allocated(li3_t)) call tr_deallocate(li3_t,"li3_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(li3_t,1,index_start+nstep,"li3_t",CAT_UNKNOWN)

  if (allocated(li3_tot_t)) call tr_deallocate(li3_tot_t,"li3_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(li3_tot_t,1,index_start+nstep,"li3_tot_t",CAT_UNKNOWN)

  if (allocated(part_src_tot_t)) call tr_deallocate(part_src_tot_t,"part_src_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(part_src_tot_t,1,index_start+nstep,"part_src_tot_t",CAT_UNKNOWN)

  if (allocated(heat_src_tot_t)) call tr_deallocate(heat_src_tot_t,"heat_src_tot_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(heat_src_tot_t,1,index_start+nstep,"heat_src_tot_t",CAT_UNKNOWN)

  if (allocated(volume_t)) call tr_deallocate(volume_t,"volume_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(volume_t,1,index_start+nstep,"volume_t",CAT_UNKNOWN)

  if (allocated(area_t)) call tr_deallocate(area_t,"area_t",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(area_t,1,index_start+nstep,"area_t",CAT_UNKNOWN)

  if (allocated(mag_ener_src_tot)) call tr_deallocate(mag_ener_src_tot,"mag_ener_src_tot",CAT_UNKNOWN)
  if (nstep .gt. 0) call tr_allocate(mag_ener_src_tot,1,index_start+nstep,"mag_ener_src_tot",CAT_UNKNOWN)

endif

! --- Read numerical profiles for rho, T, ff', toroidal rotation and neoclassical coefficients.
call read_num_profiles(my_id)

! --- Determine the derivatives of the numerical input profiles.
call derive_num_profiles(my_id)
  
return
end subroutine initialise_parameters
