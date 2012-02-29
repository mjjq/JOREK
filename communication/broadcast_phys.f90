subroutine Broadcast_phys(my_id)
!----------------------------------------------------------
! Broadcast all input parameters from MPI thread 0 to the others
!----------------------------------------------------------
use tr_module
use phys_module
use mumps_module,  only: use_mumps, no_zeros_mumps
use murge_module,  only: use_murge, use_murge_element
use pastix_module, only: use_pastix, no_zeros_pastix, pastix_smp_only

implicit none

! --- Routine parameters
integer, intent(in) :: my_id

! --- internal variables
include 'mpif.h'               ! MPI fortran include file
integer                :: ierr, INT_EXT, IDBL_EXT, ILOG_EXT, CHAR_EXT, position, bufsize
character, allocatable :: buffer(:)

!----------------------------------- one line would be enough if only MPI_TYPE_STRUCT would work on IXIA
!call MPI_BCAST(phys_list,1,MPI_phys,0,MPI_COMM_WORLD,ierr)

call MPI_PACK_SIZE(1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IDBL_EXT,ierr)
call MPI_PACK_SIZE(1,MPI_INTEGER,MPI_COMM_WORLD,INT_EXT,ierr)
call MPI_PACK_SIZE(1,MPI_LOGICAL,MPI_COMM_WORLD,ILOG_EXT,ierr)
call MPI_PACK_SIZE(1,MPI_CHARACTER,MPI_COMM_WORLD,CHAR_EXT,ierr)
!==========================MB: velocity +12 double precision
bufsize = ( (224+2*max_limiter) * IDBL_EXT + (25+n_tor) * INT_EXT + 32 * ILOG_EXT + 7*512+80 * CHAR_EXT )

allocate(buffer(bufsize))
call tr_register_mem(bufsize,"bcastp_buffer")

if (my_id .eq. 0) then
  position = 0

  call MPI_PACK(tstep,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 4
  call MPI_PACK(tstep_n,               10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   ! 14
  call MPI_PACK(F0,                     1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)         ! 15
  call MPI_PACK(GAMMA,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 16
  call MPI_PACK(Q_bar,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 17
  call MPI_PACK(sigma,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 18
   
  call MPI_PACK(zjz_0,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 19
  call MPI_PACK(zjz_1,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 20
  call MPI_PACK(zj_coef,               10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   ! 30

  call MPI_PACK(T_0,                    1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)        ! 31
  call MPI_PACK(T_1,                    1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)        ! 32
  call MPI_PACK(T_coef,                10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)    ! 42

  call MPI_PACK(Ti_0,                   1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)       ! 43
  call MPI_PACK(Ti_1,                   1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)       ! 44
  call MPI_PACK(Ti_coef,               10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   ! 54

  call MPI_PACK(Te_0,                   1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)       ! 55
  call MPI_PACK(Te_1,                   1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)       ! 56
  call MPI_PACK(Te_coef,               10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   ! 66

  call MPI_PACK(rho_0,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 67
  call MPI_PACK(rho_1,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 68
  call MPI_PACK(rho_coef,              10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)  ! 78

  call MPI_PACK(FF_0,                   1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)       ! 79
  call MPI_PACK(FF_1,                   1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)       ! 80
  call MPI_PACK(FF_coef,               10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   ! 90

  call MPI_PACK(heatsource,             1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   ! 91
  call MPI_PACK(heatsource_i,           1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) ! 92
  call MPI_PACK(heatsource_e,           1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) ! 93
  call MPI_PACK(particlesource,         1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) ! 94
  
  call MPI_PACK(ZK_perp,               10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)    ! 104
  call MPI_PACK(ZK_par,                 1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 105
  call MPI_PACK(ZK_i_perp,             10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)  ! 115
  call MPI_PACK(K_i_par,                1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)     ! 116
  call MPI_PACK(ZK_e_perp,             10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)  ! 126
  call MPI_PACK(K_e_par,                1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)     ! 127
  call MPI_PACK(D_perp,                10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)     ! 137
  call MPI_PACK(D_par,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)       ! 138
  call MPI_PACK(D_neutral,              1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   ! 139
  
  call MPI_PACK(eta,                    1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)        ! 1
  call MPI_PACK(visco,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 2
  call MPI_PACK(visco_par,              1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)  ! 3

  call MPI_PACK(eta_num,                1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)         ! 144
  call MPI_PACK(visco_num,              1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)       ! 145
  call MPI_PACK(visco_par_num,          1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   ! 146
  call MPI_PACK(D_perp_num,             1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 147

  call MPI_PACK(tauIC,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)              ! 140
  call MPI_PACK(gamma_sheath,           1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)       ! 141
  call MPI_PACK(density_reflection,     1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) ! 142
  call MPI_PACK(central_density,        1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)    ! 143

  call MPI_PACK(pellet_amplitude,       1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   ! 148
  call MPI_PACK(pellet_R,               1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)           ! 149
  call MPI_PACK(pellet_Z,               1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)           ! 150
  call MPI_PACK(pellet_phi,             1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)         ! 151
  call MPI_PACK(pellet_radius,          1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 152
  call MPI_PACK(pellet_sig,             1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)         ! 153
  call MPI_PACK(pellet_length,          1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)      ! 154
  call MPI_PACK(pellet_psi,             1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)         ! 155
  call MPI_PACK(pellet_delta_psi,       1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   ! 156
  call MPI_PACK(pellet_velocity_R,      1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)  ! 157
  call MPI_PACK(pellet_velocity_Z,      1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)  ! 158
  call MPI_PACK(pellet_density,         1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)     ! 159
  call MPI_PACK(pellet_particles,       1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   ! 160
  
  call MPI_PACK(particlesource_psin,    1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)! 161
  call MPI_PACK(particlesource_sig,     1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) ! 162
  call MPI_PACK(heatsource_psin,        1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)    ! 163
  call MPI_PACK(heatsource_sig,         1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)     ! 164_R

  call MPI_PACK(rhon_0,                 1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) 
  call MPI_PACK(rhon_1,                 1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) 
  call MPI_PACK(rhon_coef,             10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  
  call MPI_PACK(D_neutral_x,            1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(D_neutral_y,            1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(D_neutral_p,            1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  
  call MPI_PACK(ksi_ion,                1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  
  
  call MPI_PACK(mgi_amplitude,          1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(mgi_R,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(mgi_Z,                  1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(mgi_phi,                1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(mgi_radius,             1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(mgi_sig,                1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(mgi_length,             1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  
  call MPI_PACK(n_zero,                 1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  
  call MPI_PACK(psi_axis_init,          1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   
  call MPI_PACK(XR_r(:),                2,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   
  call MPI_PACK(SIG_r(:),               2,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   
  call MPI_PACK(XR_tht(:),              2,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   
  call MPI_PACK(SIG_tht(:),             2,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   
  
  call MPI_PACK(SIG_closed,             1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   
  call MPI_PACK(SIG_open,               1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   
  call MPI_PACK(SIG_outer,              1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)    
  call MPI_PACK(SIG_inner,              1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)    
  call MPI_PACK(SIG_private,            1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)  
  call MPI_PACK(SIG_up_priv,            1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)  
  call MPI_PACK(SIG_leg_0,              1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   
  call MPI_PACK(SIG_leg_1,              1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   
  call MPI_PACK(SIG_up_leg_0,           1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) 
  call MPI_PACK(SIG_up_leg_1,           1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) 

  call MPI_PACK(dPSI_open,              1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   
  call MPI_PACK(dPSI_outer,             1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)  
  call MPI_PACK(dPSI_inner,             1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)   
  call MPI_PACK(dPSI_private,           1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) 
  call MPI_PACK(dPSI_up_priv,           1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) 
  call MPI_PACK(R_limiter,    max_limiter,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) 
  call MPI_PACK(Z_limiter,    max_limiter,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr) 
!==============================MB velocity +12
  call MPI_PACK(V_0,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(V_1,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(V_coef,10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
!=================
  call MPI_PACK(nstep,                  1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(nstep_n,               10,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(n_flux,                 1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(n_tht,                  1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(n_open,                 1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(n_outer,                1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(n_inner,                1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(n_leg,                  1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(n_private,              1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(n_up_leg,               1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(n_up_priv,              1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(mode,               n_tor,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(index_start,            1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(index_now,              1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(gmres_max_iter,         1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(xcase,                  1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(n_limiter,              1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(eta_T_dependent,        1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(visco_T_dependent,      1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(restart,                1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(regrid,                 1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(import_equil,           1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(xpoint,                 1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(freeboundary,           1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(resistive_wall,         1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(bc_natural_open,        1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(use_pellet,             1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(time_evol_scheme,      80,MPI_CHARACTER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(rho_file,             512,MPI_CHARACTER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(T_file,               512,MPI_CHARACTER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(ffprime_file,         512,MPI_CHARACTER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(d_perp_file,          512,MPI_CHARACTER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(zk_perp_file,         512,MPI_CHARACTER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(zk_e_perp_file,       512,MPI_CHARACTER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(zk_i_perp_file,       512,MPI_CHARACTER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(num_rho,                1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(num_T,                  1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(num_ffprime,            1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(num_d_perp,             1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(num_zk_perp,            1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(num_zk_e_perp,          1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(num_zk_i_perp,          1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(export_for_nemec,       1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(linear_run,             1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(gmres,                  1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(use_mumps,              1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(use_pastix,             1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(use_murge,              1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(use_murge_element,      1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(pastix_smp_only,        1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(refinement,             1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(grid_to_wall,           1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(adaptive_time,          1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(equil,                  1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(bench_without_plot,     1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(no_zeros_pastix,        1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(no_zeros_mumps,         1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

endif

call MPI_BCAST(buffer,bufsize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

if (my_id .ne. 0) then

  position = 0

  call MPI_UNPACK(buffer,bufsize,position,tstep,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,tstep_n,               10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,F0,                     1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,GAMMA,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,Q_bar,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,sigma,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,zjz_0,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,zjz_1,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,zj_coef,               10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,T_0,                    1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,T_1,                    1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,T_coef,                10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,Ti_0,                   1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,Ti_1,                   1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,Ti_coef,               10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,Te_0,                   1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,Te_1,                   1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,Te_coef,               10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,rho_0,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,rho_1,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,rho_coef,              10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,FF_0,                   1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,FF_1,                   1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,FF_coef,               10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,heatsource,             1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,heatsource_i,           1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,heatsource_e,           1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,particlesource,         1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,ZK_perp,               10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,ZK_par,                 1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,ZK_i_perp,             10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,K_i_par,                1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,ZK_e_perp,             10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,K_e_par,                1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,D_perp,                10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,D_par,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,D_neutral,              1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,eta,                    1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,visco,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,visco_par,              1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,eta_num,                1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,visco_num,              1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,visco_par_num,          1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,D_perp_num,             1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,tauIC,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,gamma_sheath,           1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,density_reflection,     1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,central_density,        1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,pellet_amplitude,       1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pellet_R,               1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pellet_Z,               1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pellet_phi,             1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pellet_radius,          1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pellet_sig,             1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pellet_length,          1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pellet_psi,             1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pellet_delta_psi,       1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pellet_velocity_R,      1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pellet_velocity_Z,      1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pellet_density,         1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pellet_particles,       1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,particlesource_psin,    1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,particlesource_sig,     1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,heatsource_psin,        1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,heatsource_sig,         1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,rhon_0,                 1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,rhon_1,                 1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,rhon_coef,             10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,D_neutral_x,            1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,D_neutral_y,            1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,D_neutral_p,            1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,ksi_ion,                1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,mgi_amplitude,          1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,mgi_R,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,mgi_Z,                  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,mgi_phi,                1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,mgi_radius,             1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,mgi_sig,                1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,mgi_length,             1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,n_zero,                 1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)


  call MPI_UNPACK(buffer,bufsize,position,psi_axis_init,          1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,XR_r(:),                2,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,SIG_r(:),               2,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,XR_tht(:),              2,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,SIG_tht(:),             2,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,SIG_closed,		  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,SIG_open,		  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,SIG_outer,		  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,SIG_inner,		  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,SIG_private,  	  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,SIG_up_priv,  	  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,SIG_leg_0,		  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,SIG_leg_1,		  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,SIG_up_leg_0, 	  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,SIG_up_leg_1, 	  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,dPSI_open,		  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,dPSI_outer,		  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,dPSI_inner,		  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,dPSI_private, 	  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,dPSI_up_priv, 	  1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,R_limiter,    max_limiter,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,Z_limiter,    max_limiter,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

!=================MB:initial profile of parallel velocity
  call MPI_UNPACK(buffer,bufsize,position,V_0,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,V_1,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,V_coef,10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
!==========================
  call MPI_UNPACK(buffer,bufsize,position,nstep,                  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,nstep_n,               10,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  
  call MPI_UNPACK(buffer,bufsize,position,n_flux,                 1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,n_tht,		  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,n_open,		  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,n_outer,		  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,n_inner,		  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,n_leg,		  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,n_private,		  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,n_up_leg,		  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,n_up_priv,		  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  
  call MPI_UNPACK(buffer,bufsize,position,mode,               n_tor,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,index_start,            1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,index_now,              1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,gmres_max_iter,         1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,xcase,                  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,n_limiter,              1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,eta_T_dependent,        1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,visco_T_dependent,      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,restart,                1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,regrid,                 1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,import_equil,           1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,xpoint,                 1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,freeboundary,           1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,resistive_wall,         1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,bc_natural_open,        1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,use_pellet,             1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,time_evol_scheme,      80,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,rho_file,             512,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,T_file,               512,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,ffprime_file,         512,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,d_perp_file,          512,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,zk_perp_file,         512,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,zk_e_perp_file,       512,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,zk_i_perp_file,       512,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,num_rho,                1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,num_T,                  1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,num_ffprime,            1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,num_d_perp,             1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,num_zk_perp,            1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,num_zk_e_perp,          1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,num_zk_i_perp,          1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,export_for_nemec,       1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,linear_run,             1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,gmres,                  1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,use_mumps,              1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,use_pastix,             1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,use_murge,              1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,use_murge_element,      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,pastix_smp_only,        1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,refinement,             1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,grid_to_wall,           1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,adaptive_time,          1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,equil,                  1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,bench_without_plot,     1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,no_zeros_pastix,        1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,no_zeros_mumps,         1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)

endif

call tr_unregister_mem(bufsize,"bcastp_buffer")
deallocate(buffer)

return
end subroutine Broadcast_phys
