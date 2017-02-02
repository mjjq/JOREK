!> Initialise input parameters and read the input namelist
subroutine initialise_parameters(my_id, filename)

use tr_module
use phys_module
use data_structure
use constants
use mpi_mod
use corr_neg
use mumps_module,  only: use_mumps, no_zeros_mumps
use murge_module,  only: use_murge, use_murge_element
use pastix_module, only: use_pastix, no_zeros_pastix, pastix_smp_only, pastix_pivot
use vacuum
use wsmp_module,   only: use_wsmp

implicit none

! --- Routine parameters
integer,                      intent(in) :: my_id
character(len=*),             intent(in) :: filename
real*8 :: vacuum_fraction, b_over_a, a_over_b

! --- Local variables
type (type_node_list)    :: node_list
type (type_element_list) :: element_list

integer :: ierr,err,i,ifail,i_elm
integer :: err_alloc=0, err_alloc_rnd=0

real*8, dimension(2) :: P, P_s, P_t, P_phi
real*8  :: R, R_s, R_t, Z, Z_s, Z_t
real*8  :: s_out,t_out,R_out,Z_out

real*8  :: n_SI, T_eV, n_corr, T_corr
real*8  :: spi_gd_angle_01, spi_gd_angle_02        !The dispersion angles for each spi
real*8  :: spi_rotation_01, spi_rotation_02        !The rotation angle from spi coordinate to real coordinate
real*8  :: spi_Vel_totref, spi_Vel_R_tmp, spi_Vel_Z_tmp, spi_Vel_RxZ_tmp
real*8  :: spi_Vel_x, spi_Vel_y, spi_Vel_z         !Spi velocity in injection coordinate
real*8  :: spi_R_inj, spi_Z_inj, spi_phi_inj       !Injection position of SPI
real*8  :: spi_R_tmp, spi_Z_tmp, spi_phi_tmp
real*8, allocatable :: rnd(:)                      !The random number array 

! --- Namelist with input parameters.
namelist /in1/  tstep, nstep, tstep_n, nstep_n,                     &
                rst_hdf5,                                           &
                eta, visco, visco_par,                              &
                restart, rst_format, regrid, bootstrap,             &
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
                R_boundary, Z_boundary, psi_boundary, n_boundary,   &
                n_pfc, tokamak_device,                              &
                Rmin_pfc, Rmax_pfc, Zmin_pfc, Zmax_pfc, current_pfc,&
                F0, gamma_sheath, density_reflection,               &
                zjz_0, zjz_1, zj_coef,                              &
                rho_0, rho_1, rho_coef,                             &
                T_0,   T_1,   T_coef,                               &
                FF_0,  FF_1,  FF_coef,                              &
                ZK_par, ZK_par_max, ZK_perp, D_par, D_perp,         &
                particlesource, heatsource, tauIC,                  &
                eta_num, visco_num, visco_par_num, D_perp_num,      &
                ZK_perp_num, Dn_perp_num, time_evol_scheme,         &
                pellet_amplitude, pellet_R, pellet_Z, pellet_phi,   &
                pellet_radius, pellet_sig, pellet_length,           &
                pellet_psi, pellet_delta_psi, pellet_density,       &
                pellet_velocity_R, pellet_velocity_Z,               &
                central_density, central_mass,                      &
                pellet_particles, use_pellet,                       &
                ellip,tria_u,tria_l,quad_u,quad_l,                  &
                xampl,xwidth,xsig,xtheta,xshift,xleft, xpoint,      &
                xcase, D_perp_file, ZK_perp_file,                   &
                rho_file, T_file, ffprime_file,                     &
                freeboundary_equil, freeboundary,                   &
                resistive_wall, wall_resistivity,                   &
                bc_natural_open,                                    &
                use_mumps, use_pastix, use_murge, use_murge_element,&
                use_wsmp, n_tor_fft_thresh,                         &
                pastix_smp_only, refinement, grid_to_wall,          &
                adaptive_time, equil, bench_without_plot,           &
                no_zeros_pastix, no_zeros_mumps,                    &
                eta_T_dependent, visco_T_dependent,                 &
                zkpar_T_dependent,                                  & 
                heatsource_psin, heatsource_sig,                    &
                particlesource_psin, particlesource_sig,            &
                produce_live_data, gmres, gmres_max_iter,           &
                gmres_m, gmres_4, gmres_tol, iter_precon,           &
                tgnum,  pastix_pivot,                               &
                linear_run, export_for_nemec,                       &
#ifdef USE_HDF5
                save_diagnostics_HDF5,h5_diag_nbtime,               &
#endif

                V_0,V_1,V_coef, output_bnd_elements,                &
                n_limiter, R_limiter, Z_limiter,                    &
                R_Z_psi_bnd_file, wall_file,time_evol_scheme,       &
                D_prof_neg, ZK_prof_neg, T_min,                     &

                D_neutral_x, D_neutral_y, D_neutral_p,              &
                mgi_sig, mgi_deltaphi, ksi_ion, RMP_on, lambda, tset,    &  
                mgi_amplitude, mgi_R, mgi_Z, mgi_phi, mgi_radius,   &
                spi_Vel_Rref,spi_Vel_Zref, using_spi, n_spi,        &
                spi_Vel_RxZref, spi_radiusref, flag_spi,            &
                ng_radius_ratio, ng_radius_min, spi_angle,          &
                spi_L_inj, K_Dmv, A_Dmv, L_tube, V_Dmv, P_Dmv,      &
                t_mgi, JET_MGI, ASDEX_MGI,                          &
                delta_n_convection, nimp_bg,                        &
                RMP_on, lambda, tset, RMP_psi_cos_file, RMP_psi_sin_file, &
                amix, amix_freeb, equil_accuracy,                   &
                equil_accuracy_freeb, current_ref, FB_Ip_position,  &
                FB_Ip_integral, Z_axis_ref, FB_Zaxis_position,      &
                FB_Zaxis_derivative,FB_Zaxis_integral, start_VFB,   &
                n_feedback_current, n_feedback_vertical,            &
                n_iter_freeb, n_coils_nml, coils0,                  &
                Zaxis_find_limit, PF_pert_start_time

 call preset_parameters()


 if (my_id .eq. 0) then

  ! --- Preset input parameters to reasonable default values.
  !call preset_parameters()

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

  ! --- Calculate normalisation factor for MGI source (related to its toroidal shape)
  mgi_tor_norm = mgi_deltaphi * PI**0.5 * ERF(PI/mgi_deltaphi)

   if (trim(R_Z_psi_bnd_file) .ne. 'none') then

  ! --- Open the file.
     OPEN(UNIT=243, FILE=R_Z_psi_bnd_file, FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
     if ( err /= 0 ) then
       write(*,*) 'ERROR in initialise_parameters: Cannot open file '//TRIM(R_Z_psi_bnd_file)//'.'
       stop
     endif
     write(*,'(A)') ' boundary info from R_Z_psi_bnd_file: R_boundary, Z_boundary, psi_boundary '

     do i=1,n_boundary
       read(243,*) R_boundary(i),Z_boundary(i),psi_boundary(i)
       write(*,*) R_boundary(i),Z_boundary(i),psi_boundary(i)
     enddo
   endif

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
  
  if (allocated(xtime)) call tr_deallocate(xtime,"xtime",CAT_GRID)
  if (nstep .gt. 0) call tr_allocate(xtime,1,nstep,"xtime",CAT_GRID)
  
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


! --- Read numerical profiles for rho, T, and ff'.
call read_num_profiles(my_id)

! --- Determine the derivatives of the numerical input profiles.
call derive_num_profiles(my_id)

! --- Set current source term to zero
eta_T_0 = 0.0

! --- Initialize the shattered pellet position
!spi_R = mgi_R
!spi_Z = mgi_Z

if (my_id /= 0 .and. using_spi == .true.) then
  write(*,*) "Unexpected using_spi flag! my_id = ", my_id
end if

if (using_spi == .true.) then

  if (allocated(pellets)) then
    deallocate(pellets)
  end if

  allocate (pellets(n_spi),stat=err_alloc)  !< Dynamically allocate memeries for pellets

  if (err_alloc /= 0) then
    write(*,*) "Error when trying to dynamically allocate memeries for pellets, reverting to non-SPI case."
    using_spi = .false. 
  else
    if (n_spi >= 1) then

      if (my_id == 0 .and. restart == .false.) then
        open(20,file="pellets_parameters.dat",status="REPLACE")
        write(20,"(A,A11)",advance="no") "# t, "
      end if

      if (allocated(rnd)) then
        deallocate(rnd)
      end if

      allocate (rnd(2*n_spi),stat=err_alloc_rnd)  !< Dynamically allocate memeries for randoms

      if (err_alloc_rnd /= 0) then
        write(*,*) "Error when trying to dynamically allocate memeries for randoms."
        using_spi = .false.
      end if


!===================Determine the rotational transform of coordinate===============
!Here, we perform the following rotational transform from the original
!coordinate R, Z, RxZ to the so-called spi coordinate x, y ,z, with the reference
!direction of spi injection being the z axis, while y axis locates within the 
!same surface as Z and z. The rotational transform from x, y, z to R, Z, RxZ is
!as the following: first, we rotate the system around x axis clockwise, facing
!the positive x direction, for spi_rotation_01 to get coordinate X', Y', Z'. 
!Then we further rotate around Y' clockwise, facing the positive Y' direction for
!spi_rotation_02 to acquire R, Z, RxZ. Hence we have:
!R   = cos(spi_rotation_02)*x - sin(spi_rotation_02)*(-sin(spi_rotation_01)*y + cos(spi_rotation_01)*z)
!Z   = cos(spi_rotation_01)*y + sin(spi_rotation_01)*z
!RxZ = sin(spi_rotation_02)*x + cos(spi_rotation_02)*(-sin(spi_rotation_01)*y + cos(spi_rotation_01)*z)

      spi_Vel_totref  = sqrt(spi_Vel_Rref**2+spi_Vel_Zref**2+spi_Vel_RxZref**2)

      spi_R_inj       = mgi_R - spi_L_inj * (spi_Vel_Rref/spi_Vel_totref)
      spi_Z_inj       = mgi_Z - spi_L_inj * (spi_Vel_Zref/spi_Vel_totref)
      spi_phi_inj     = mgi_phi - spi_L_inj * (spi_Vel_RxZref/spi_Vel_totref)/mgi_R

      spi_rotation_01 = asin(spi_Vel_Zref/spi_Vel_totref)
      if (cos(spi_rotation_01) == 0.) then
        spi_rotation_02 = 0.
      else
        spi_rotation_02 = acos(spi_Vel_RxZref/(spi_Vel_totref*cos(spi_rotation_01)))
      end if

      write(*,*) "Rotational transform: ", spi_rotation_01, spi_rotation_02

!==========================End of rotational angles==============================

      CALL random_number(rnd)

      do i=1, n_spi

        spi_gd_angle_01 = rnd(2 * i - 1) * spi_angle / 2.0
        spi_gd_angle_02 = rnd(2 * i) * 2. * PI

        spi_Vel_x       = spi_Vel_totref * sin(spi_gd_angle_01) * cos(spi_gd_angle_02)
        spi_Vel_y       = spi_Vel_totref * sin(spi_gd_angle_01) * sin(spi_gd_angle_02)
        spi_Vel_z       = spi_Vel_totref * cos(spi_gd_angle_02) 

        spi_Vel_R_tmp   = spi_Vel_x * cos(spi_rotation_02)                          &
                          - sin(spi_rotation_02) * (-sin(spi_rotation_01)*spi_Vel_y &
                          + cos(spi_rotation_01)*spi_Vel_z)

        spi_Vel_Z_tmp   = cos(spi_rotation_01) * spi_Vel_y                          &
                          + sin(spi_rotation_01) * spi_Vel_z

        spi_Vel_RxZ_tmp = spi_Vel_x * sin(spi_rotation_02)                          &
                          - cos(spi_rotation_02) * (-sin(spi_rotation_01)*spi_Vel_y &
                          + cos(spi_rotation_01)*spi_Vel_z)

        spi_R_tmp       = spi_R_inj + spi_L_inj * (spi_Vel_R_tmp/spi_Vel_totref)
        spi_Z_tmp       = spi_Z_inj + spi_L_inj * (spi_Vel_Z_tmp/spi_Vel_totref)
        spi_phi_tmp     = spi_phi_inj + spi_L_inj * (spi_Vel_RxZ_tmp/spi_Vel_totref)/mgi_R



        pellets(i)%spi_R       = spi_R_tmp
        pellets(i)%spi_Z       = spi_Z_tmp
        pellets(i)%spi_phi     = spi_phi_tmp
        pellets(i)%spi_Vel_R   = spi_Vel_R_tmp
        pellets(i)%spi_Vel_Z   = spi_Vel_Z_tmp
        pellets(i)%spi_Vel_RxZ = spi_Vel_RxZ_tmp
        pellets(i)%spi_radius  = spi_radiusref
        pellets(i)%spi_abl     = mgi_amplitude

        if (my_id == 0 .and. i < n_spi .and. restart == .false.) then
          write(20,"(A11,I3.3)",advance="no") "abl N.: ", i
        elseif (my_id == 0 .and. i == n_spi .and. restart == .false.) then
          write(20,"(A11,I3.3)") "abl N.: ", i
        end if
        
      end do

      if (my_id == 0 .and. restart == .false.) then
        close(20)
      end if

      write(*,*) "SPI initialized successfully."

      deallocate(rnd)
    else
      write(*,*) "...... Seriously!? Reverting to non-SPI case."
      using_spi = .false.
    end if
  end if

end if
  
return
end subroutine initialise_parameters
