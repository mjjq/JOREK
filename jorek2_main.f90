!> JOREK 2.0 -- Solves the (reduced) MHD equations in 3D toroidal geometry.
!!
!! - solvers implemented:
!!   - MUMPS
!!   - PastiX
!!   - GMRES (+MUMPS or PastiX preconditioner)
!!
!! - required libraries :
!!   - MPI
!!   - MUMPS
!!   - PastiX
!!   - SCOTCH (metis)
!!   - FFTW
!!   - SCALAPACK (BLACS)
!!   - LAPACK, BLAS
!!   - PPPLIB
!!
!! @author Guido Huysmans (Euratom / CEA Association)
!! @date 18-7-2008
program JOREK2

  use constants
  use mumps_module
  use pastix_module
  use murge_module,        only: murge_initialization, murge_setGraph, MURGE_Clean, use_murge,     &
    use_murge_element, murge_initialised, murge_glob2loc, murge_loc2glob, murge_id,&
    murge_termination
  use wsmp_module
  use data_structure
  use phys_module
  use parameters
  use global_distributed_matrix
  use nodes_elements
  use pellet_module
  use equil_info
  use boundary,            only: boundary_from_grid
  use vacuum,              only: vacuum_preset, vacuum_init, broadcast_vacuum, wall_curr_initialized
  use vacuum_response,     only: get_vacuum_response, update_response, init_wall_currents, I_coils
  use vacuum_equilibrium,  only: import_external_fields
  use live_data,           only: init_live_data, write_live_data, finalize_live_data

! these write additional live data (global data) used when an ECCD current is applied)
#ifdef JECCD
  use live_data2,          only: init_live_data2, write_live_data2, finalize_live_data2
  use live_data3,          only: init_live_data3, write_live_data3, finalize_live_data3
#ifdef JEC2DIAG
  use live_data4,          only: init_live_data4, write_live_data4, finalize_live_data4
#endif
#endif

  use solve_mat_n
  use tr_module
  use clock_module
#ifdef USE_HDF5
  use hdf5
  use HDF5_io_module
  use out_save_module
#endif
  use mpi_mod

#if JOREK_MODEL == 500
  use mgi_module
#endif

  use, intrinsic :: iso_c_binding
  
  implicit none

#ifdef USE_FFTW
  include 'fftw3.f03'
#endif
  
#include "r3_info.h"
#include "version.h"
  
  interface

    subroutine distribute_vector(my_id,rhs,rhs_dis,again)
      real*8               :: rhs(:), rhs_dis(:)
      integer              :: my_id
      logical              :: again
    end subroutine distribute_vector

    subroutine distribute_harmonics(my_id,my_id_n,n_cpu)
      integer              :: my_id, my_id_n,n_cpu
    end subroutine distribute_harmonics

    subroutine gmres_driver(my_id,my_id_n,i_tor,n_tor,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)
      integer :: i_tor(:), my_id, my_id_n, MPI_COMM_N, MPI_COMM_MASTER
      integer :: iter_gmres, n_tor
    end subroutine gmres_driver
    
    subroutine equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint2,xcase2, nice_q)
      use data_structure
      integer(kind=4),             intent(in)    :: my_id
      integer(kind=4),             intent(in)    :: xcase2
      type (type_node_list),       intent(inout) :: node_list
      type (type_element_list),    intent(inout) :: element_list
      type (type_bnd_node_list)   ,intent(inout) :: bnd_node_list    
      type (type_bnd_element_list),intent(inout) :: bnd_elm_list    
      logical(kind=4),             intent(in)    :: xpoint2
      logical(kind=4),             intent(in)    :: nice_q
    end subroutine equilibrium
     
    SUBROUTINE construct_matrix_murge(my_id,node_list,element_list,                &
         &                            bnd_node_list, local_elms,                   &
         &                            n_local_elms, xpoint2, xcase2,               &
         &                            minRad, R_axis, Z_axis, psi_axis,            &
         &                            psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,     &
         &                            gmres, i_tor, n_cpu,                         &
         &                            mpi_comm_n, MPI_COMM_TRANS,                  &
         &                            my_id_trans, n_cpu_trans, solve_only)
      use data_structure, only : type_node,type_element,type_element_list,type_bnd_node_list,      &
        type_node_list, thread_struct
      integer(kind=4) :: n_cpu
      integer(kind=4) ,target :: n_local_elms
      integer(kind=4) ,target :: my_id
      type (type_node_list) ,target :: node_list
      type (type_element_list) ,target :: element_list
      type (type_bnd_node_list) ,target :: bnd_node_list
      integer(kind=4) ,target :: local_elms(n_local_elms)
      logical(kind=4) ,target :: xpoint2
      integer(kind=4) ,target :: xcase2
      real(kind=8) ,target :: minrad
      real(kind=8) ,target :: r_axis
      real(kind=8) ,target :: z_axis
      real(kind=8) ,target :: psi_axis
      real(kind=8) ,target :: psi_bnd
      real(kind=8) ,target :: r_xpoint(:)
      real(kind=8) ,target :: z_xpoint(:)
      real(kind=8) ,target :: psi_xpoint(:)
      logical(kind=4) ,target :: gmres
      integer(kind=4) :: i_tor(n_cpu)
      integer(kind=4) :: mpi_comm_n
      integer(kind=4) ,target :: mpi_comm_trans
      integer(kind=4) ,target :: my_id_trans
      integer(kind=4) ,target :: n_cpu_trans
      logical(kind=4) ,target :: solve_only
    end subroutine construct_matrix_murge
  end interface
  
  type (type_surface_list) :: surface_list, flux_list
  type (t_equil_state)     :: equil_state
  real*8                   :: W_mag(n_tor), W_kin(n_tor), growth_mag, growth_kin, growth_mag0, growth_kin0
#ifdef JECCD
  real*8                   :: A_tem(n_tor), A_den(n_tor), A_jen(n_tor), A_jec(n_tor),A_jec1(n_tor), A_jec2(n_tor)
#endif
  real*8                   :: psi_lim, R_lim, Z_lim
  real*8                   :: t_matrix, t_send, t_solve
  type(clcktype)           :: t_itstart, t0, t1
  real*8                   :: psi_bnd, psi_axis, R_axis, Z_axis, s_axis, t_axis, minRad
  real*8                   :: psi_xpoint(2), R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2), mindelta, maxdelta
  integer                  :: my_id, my_id_n, my_id_master
  integer                  :: istep,jstep,ierr,i,itor,inode, i_elm_axis, i_elm_xpoint(2)
  integer                  :: n_local_ELMs
  integer                  :: i_rank(n_tor), n_cpu, n_cpu_n, n_cpu_master, m_cpu, n_masters, n_cpu_trans, my_id_trans
  integer                  :: iter_gmres
  integer                  :: MPI_COMM_N, MPI_GROUP_MASTER, MPI_GROUP_WORLD, MPI_COMM_MASTER, MPI_COMM_TRANS
  integer                  :: i_find, i_elm_find(8)
  real*8                   :: Router,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
  real*8                   :: Zouter,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
  real*8                   :: s_find(8), t_find(8)
  character*8              :: label, itlabel
  character*14             :: fileout
  integer                  :: required,provided,StatInfo
  integer, allocatable     :: local_elms(:), i_tor(:), index_min(:), index_max(:)
  real*8                   :: zjz, E_min, E_max
  logical                  :: solve_only
  integer*4                :: rank, comm_size 
  real*8                   :: zn,  dn_dpsi,  dn_dz,  dn_dpsi2,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi_dz2,  dn_dpsi2_dz
  real*8                   :: zT,  dT_dpsi,  dT_dz,  dT_dpsi2,  dT_dz2,  dT_dpsi_dz,  dT_dpsi3,  dT_dpsi_dz2,  dT_dpsi2_dz
  real*8                   :: zTi, dTi_dpsi, dTi_dz, dTi_dpsi2, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz
  real*8                   :: zTe, dTe_dpsi, dTe_dz, dTe_dpsi2, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz
  real*8                   :: zFFprime, dFFprime_dpsi, dFFprime_dz, dFFprime_dpsi_dz,dFFprime_dpsi2,dFFprime_dz2
  real*8                   :: Rp, Zp, R_out,Z_out,s_out,t_out,P_s,P_t,P_st,P_ss,P_tt, psi
  real*8                   :: Rp_start, Rp_end, density_tot,density_in,density_out,pressure_tot,pressure_in,pressure_out
  real*8,allocatable       :: xp(:), yp1(:), yp2(:), yp3(:)
  integer                  :: nplot, iplot, i_elm, ifail, ivar, iter_big, n_aa, iter_prev
  logical                  :: is_local, file_exists
  integer                  :: i_elem, inode1, i_order, index_node1
  type (type_element)      :: element
  integer                  :: index_size, id_elements
  integer                  :: list_to_be_refined(n_ref_list), n_to_be_refined    
  REAL*8                   :: max_time, min_time, tsecond
  integer, allocatable     :: tab_n_local_elems(:)
  real*8                   :: t_this, sum_deltas
  integer                  :: h5_nbsave_current,h5_nbsave,h5_nbsave_previous
! =================== plot NEO coeffs ==================
  real*8                   :: amu_neo_node, aki_neo_node
  real*8,allocatable       :: mu_neo(:), ki_neo(:)
! ======================================================
#ifdef USE_FFTW
  real*8     :: in_fft(1:n_plane)
  complex*16 :: out_fft(1:n_plane)
#endif

  real*8  :: DUMMY_REAL(1:1)
  integer :: DUMMY_INT (1:1)
  !***********************************************************************
  !*                  intialisation                                      *
  !***********************************************************************
  ! --- Initialize OpenMP threads before MPI_init
  !call init_threads()
  
  ! --- Initialise MPI / threaded MPI
#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif
#ifdef STAN_FLAG
required = 0
#endif
  call MPI_Init_thread(required, provided, StatInfo)

  call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread
  
  ! --- Determine ID of each MPI proc
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  my_id = rank
  
  ! --- Determine number of MPI procs
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  n_cpu = comm_size

  if (my_id == 0) then
    write(*,*) '*************************************************'
    write(*,*) '*   3D Reduced MHD : JOREK_2.0                  *'
    write(*,*) '*************************************************'
    write(*,*) ' MPI processes       : ', n_cpu
    write(*,*) ' OpenMP threads      : ', nbthreads
    write(*,*) ' SVN revision        : ', SVN_VERSION
    111 format(2x,a,': ',a)
    write(*,111) 'compile_time        ', trim(adjustl(compile_time))
    write(*,111) 'compile_user        ', trim(adjustl(compile_user))
    write(*,111) 'compile_machine     ', trim(adjustl(compile_machine))
    write(*,111) 'compile_dir         ', trim(adjustl(compile_dir))
    write(*,111) 'compile_command     ', trim(adjustl(compile_command))
    write(*,111) 'compile_flags       ', trim(adjustl(compile_flags))
    write(*,111) 'compile_includes    ', trim(adjustl(compile_includes))
    write(*,111) 'compile_defines     ', trim(adjustl(compile_defines))
    write(*,111) 'compile_libs        ', trim(adjustl(compile_libs))
  end if

  ! --- Initialise memory tracing
  call tr_meminit(my_id, n_cpu)

  ! --- Initialise timing
  call clck_init()
  call r3_info_init ()
  
  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  ! --- Remove file STOP_NOW if it exists
  if ( my_id == 0 ) then
    open(42, file='STOP_NOW', iostat=ierr)
    if ( ierr == 0 ) close(42, status='delete')
  end if

  ! --- Preset some solver variables
  pastix_initialised = .false.
  pastix_analysed    = .false.
  murge_initialised  = .false.
  
  ! --- Preset input parameters to reasonable defaults, then read the input file.
  call initialise_and_broadcast_parameters(my_id, "__NO_FILENAME__")
  
  ! --- Initialize the vacuum part.
  call vacuum_init(my_id, freeboundary_equil, freeboundary, resistive_wall)

  ! --- MURGE with ntor=1 doesn't work up to now because i_tor is not allocated correctly
  if (n_tor == 1) then
    gmres     = .false.
    use_murge = .false. 
  end if
  
  ! --- Write out all parameters defined in mod_parameters and the namelist input file.
  call log_parameters(my_id)
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  ! --- Some checks not to waste any cpu time
  if (required .ne. provided) then
    write(*,*) 'FATAL : MPI_THREAD_MULTIPLE (provided < required)', my_id, required, provided
    call MPI_FINALIZE(IERR)
    stop
  else if ( (.not. use_mumps) .and. (.not. use_pastix) .and. (.not. use_wsmp) ) then
    write(*,*) ' FATAL : specify a valid solver'
    call MPI_FINALIZE(IERR)
    stop
  else if ( n_plane < 2*(n_tor-1) ) then
    write(*,*) ' FATAL: n_plane >= 2 * (n_tor-1) required to avoid aliasing.'
    call MPI_FINALIZE(IERR)
    stop
#ifndef USE_FFTW
  else if ( ( n_tor >= n_tor_fft_thresh ) .and. ( iand(n_plane,n_plane-1) /= 0 ) ) then
    write(*,*) ' FATAL: If n_tor >= n_tor_fft_thresh, n_plane must be a power of 2.'
    write(*,*) ' Hint: USE_FFTW removes this constraint.'
    call MPI_FINALIZE(IERR)
    stop
#endif
  else if ( gmres .and. (nstep > 0) .and. (mod(n_cpu,(n_tor-1)/2+1) /= 0) ) then
    write(*,'(A,i4,A,i4,A)') ' FATAL : need a multiple of ',(n_tor-1)/2+1,' cpus for ',            &
      (n_tor-1)/2+1,' harmonics'
    call MPI_FINALIZE(IERR)
    stop
  else if ( use_mumps ) then
#ifndef USE_MUMPS
    write(*,*) 'FATAL : use_mumps=.true. requires USE_MUMPS=1 in Makefile.inc'
    call MPI_FINALIZE(IERR)
    stop
#endif
  else if ( use_pastix ) then
#ifndef USE_PASTIX
     write(*,*) 'FATAL : use_pastix=.true. requires USE_PASTIX=1 in Makefile.inc'
     call MPI_FINALIZE(IERR)
     stop
#endif
     if ( use_murge ) then
#ifndef USE_MURGE
        write(*,*) 'FATAL : use_murge=.true. requires USE_PASTIX_MURGE=1 in Makefile.inc'
        call MPI_FINALIZE(IERR)
        stop
#endif
     endif
  else if ( use_wsmp ) then
#ifndef USE_WSMP
    write(*,*) 'FATAL : use_wsmp=.true. requires USE_WSMP=1 in Makefile.inc'
    call MPI_FINALIZE(IERR)
    stop
#endif
#ifdef USE_BLOCK
    write(*,*) 'FATAL : USE_BLOCK=1 in Makefile.inc is currently not possible with use_wsmp'
    call MPI_FINALIZE(IERR)
    stop
#endif
      if ( .not. restart ) then
      write(*,*) 'FATAL : use_wsmp is currently not supported for the equilibrium'
      call MPI_FINALIZE(IERR)
      stop
    end if
    if ( use_pastix ) then
      write(*,*) 'FATAL : you should only select one of use_wsmp or use_pastix'
      call MPI_FINALIZE(IERR)
      stop
    end if
  end if
  
  ! --- Initialize live data file which will be filled during the code run
  if ( (my_id == 0) .and. (.not. bench_without_plot) ) call init_live_data()
#ifdef JECCD
  if ( (my_id == 0) .and. (.not. bench_without_plot) ) call init_live_data2()
  if ( (my_id == 0) .and. (.not. bench_without_plot) ) call init_live_data3()
  if (my_id == 0) write(6,*) "initializing live data"
#ifdef JEC2DIAG
   if ( (my_id == 0) .and. (.not. bench_without_plot) ) call init_live_data4()
#endif
#endif
  
  ! --- Initialise ppplib plotting library
  if (my_id == 0)  call begplt('jorek2.ps')
  
  ! --- Define the basis functions at the Gaussian points
  call initialise_basis()
  
  call tr_print_memsize("InitStep")

  !***********************************************************************
  !*                  read restart file                                  *
  !***********************************************************************
  
  if ( restart .and. (my_id == 0) ) then
    
    ! --- Read the restart file (jorek_restart.rst)
    call import_restart(node_list, element_list, 'jorek_restart.rst', rst_format, ierr)
    if ( ierr /= 0 ) stop

#ifdef USE_HDF5
    if (save_diagnostics_HDF5 .and. (my_id .eq. 0) ) then
       write(*,*) ' '
       write(*,*) '*******************************************************************************'
       write(*,*) '******* Read and initialise quantites for HDF5 saving --RESTART MODE-- ********'
       ! --- Read and initialise quantites for HDF5 saving
       ! "h5_nbsave_previous" exists, so set it to "h5_nbsave_all" = number of HDF5 files
       ! that have been written in the previous run(s)
       h5_nbsave_previous = h5_nbsave_all
       write(*,*) '  h5_nbsave_previous = ',h5_nbsave_previous
       ! number of HDF5 files that have actually been written in the current run
       h5_nbsave_current  = 0
       write(*,*) '  h5_nbsave_current  = ',h5_nbsave_current
       ! consistency test required: "h5_diag_nbtime" cannot be smaller than 1 Alfven
       ! time otherwise the "modulo" below fails
       if ( h5_diag_nbtime < 1.d0 ) then
          h5_diag_nbtime = 1.d0
          write(*,*) '  -----> your "h5_diag_nbtime" value is stupid and has been set to 1.d0 '
       else
          write(*,*) '  h5_diag_nbtime     = ',h5_diag_nbtime
       endif
       ! number of HDF5 files that are going to be written in the current run
       ! if everything goes right
       t_this = tstep*nstep
       h5_nbsave = int((t_this)/h5_diag_nbtime)-1 + min(1,mod( floor(t_this),floor(h5_diag_nbtime) ))
       write(*,*) '  h5_nbsave          = ',h5_nbsave
       write(*,*) '*******************************************************************************'
    endif
#endif
    
    ! --- Write live data for previous time-steps
    if ( .not. bench_without_plot ) then
      do index_now = 1, index_start
        call write_live_data(index_now)
#ifdef JECCD
        call write_live_data2(index_now)
        call write_live_data3(index_now)
#ifdef JEC2DIAG
        call write_live_data4(index_now)
#endif
#endif
      end do
    end if
    
    ! --- Optional: Redo flux aligned grid (DOES NOT WORK CURRENTLY)
    if (regrid) then
      if (xpoint)  then
        if (xcase .ge. 2) then
	  call grid_double_xpoint(node_list, element_list)
        else
	  call grid_xpoint(node_list,element_list,n_flux,n_open,n_private,n_leg,n_tht,   &
        		   SIG_open,SIG_closed,SIG_private,SIG_theta,SIG_leg_0,SIG_leg_1,dPSI_open,dPSI_private, xcase)
        endif
      else
        call grid_flux_surface(xpoint,xcase, node_list, element_list, surface_list, n_flux, n_tht, xr1,  &
          sig1, xr2, sig2, refinement)
      end if
    end if
    
  end if
  if ( restart .and. freeboundary ) call broadcast_vacuum(my_id, resistive_wall)
  
  !***********************************************************************
  !*                  define grid / equilibrium                          *
  !***********************************************************************
  
  if_not_restart: if (.not. restart) then
    call tr_resetfile()
    element_list%n_elements      = 0
    bnd_elm_list%n_bnd_elements  = 0
    node_list%n_nodes            = 0
    
#ifdef USE_HDF5
    if (save_diagnostics_HDF5 .and. (my_id .eq. 0) ) then
       write(*,*) ' '
       write(*,*) '*******************************************************************************'
       write(*,*) '******* Read and initialise quantites for HDF5 saving --INITIAL STEP-- ********'
       ! --- Read and initialise quantites for HDF5 saving
       ! "h5_nbsave_previous" does not exist, so = 0
       h5_nbsave_previous = 0
       write(*,*) '  h5_nbsave_previous = ',h5_nbsave_previous
       ! number of HDF5 files that have actually been written in the current run
       h5_nbsave_current  = 0
       write(*,*) '  h5_nbsave_current  = ',h5_nbsave_current
       ! consistency test required: "h5_diag_nbtime" cannot be smaller than 1 Alfven
       ! time otherwise the "modulo" below fails
       if ( h5_diag_nbtime < 1.d0 ) then
          h5_diag_nbtime = 1.d0
          write(*,*) '  -----> your "h5_diag_nbtime" value is invalid and has been set to 1.d0 '
       else
          write(*,*) '  h5_diag_nbtime     = ',h5_diag_nbtime
       endif
       write(*,*) '*******************************************************************************'
    endif
#endif

    if (my_id == 0) then
      
      ! --- Define the boundary of the initial grid
      call define_boundary()
      
      if ((n_R > 0) .and. (n_Z > 0) .and. (n_radial > 0)) then
        
        call grid_bezier_square_polar(n_R, n_Z, n_radial, R_begin, R_end, Z_begin, Z_end, R_geo,   &
          Z_geo, amin, fbnd, fpsi, mf, .true., node_list, element_list)
        
      else if ((n_R > 0) .and. (n_Z > 0) ) then
        
        call grid_bezier_square(n_R, n_Z, R_begin, R_end, Z_begin, Z_end, .true., node_list,       &
          element_list)
        
      else if ((n_radial > 0) .and. (n_pol > 0) ) then
        
        call grid_polar_bezier(R_geo, Z_geo, amin, 0.d0, 0.d0, fbnd, fpsi, mf, n_radial, n_pol,    &
          node_list, element_list)
        
      else
        write(*,*) ' FATAL : no valid combination of grid-sizes specified'
        call MPI_FINALIZE(IERR)
        stop
      end if 
      
      ! --- Determine boundary information from the grid
      call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

      call tr_debug_write("JMAIN:Def_grid elt_list",element_list%n_elements)
      call tr_debug_write("JMAIN:Def_grid node_list",node_list%n_nodes)
      call tr_debug_write("JMAIN:Def_grid bnd_elt_list",bnd_elm_list%n_bnd_elements)
      
    end if
    
    ! --- Send boundary elements and nodes to other MPI procs
    call broadcast_boundary(my_id,bnd_elm_list,bnd_node_list)
    
    ! --- Fill the vacuum response matrices for freeboundary computations
    if ( freeboundary ) then
      call get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list, freeboundary_equil,  &
        resistive_wall)
      call update_response(tstep, freeboundary_equil, resistive_wall)
      call import_external_fields('coil_field.dat', my_id)
      if ( (.not. restart) .or. (.not. wall_curr_initialized) ) call init_wall_currents(my_id, resistive_wall)
    end if
    
    ! --- Plot the grid  
    if ( (my_id == 0) .and. (.not. bench_without_plot) ) then
      call plot_grid(node_list,element_list,bnd_elm_list,bnd_node_list,.true.,.false.,'initial')
    end if
    
#ifdef USE_MUMPS
    ! --- Initialize MUMPS solver (used for equilibrium)
    call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
    call MPI_GROUP_INCL(MPI_GROUP_WORLD,1,0,MPI_GROUP_MUMPS_EQUIL,ierr)
    call MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_GROUP_MUMPS_EQUIL,MPI_COMM_MUMPS_EQUIL,ierr)
    if (my_id == 0) call initialise_mumps(MPI_COMM_MUMPS_EQUIL)
#endif

    if (my_id == 0) then
      
      ! --- Compute the plasma equilibrium
      if (equil) then
        call equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint,xcase, .true.) 
        if (export_for_nemec) call export_nemec(node_list, element_list, xpoint, xcase)
      end if

      ! --- Determine a flux surface aligned grid
      if (n_flux > 1) then
        
        if (xpoint)  then

!          if (.not. grid_to_wall) then
          if (xcase .ge. 2) then
	    call grid_double_xpoint(node_list, element_list)
          else
	  
	    if (.not. grid_to_wall) then
	      call grid_xpoint(node_list,element_list,n_flux,n_open,n_private,n_leg,n_tht,   &
                               SIG_open,SIG_closed,SIG_private,SIG_theta,SIG_leg_0,SIG_leg_1,dPSI_open,dPSI_private, xcase)
	    else
!!! works only for ITER wall for the moment
 !            write(*,*) 'ITER wall started'
              call grid_xpoint_wall(node_list,element_list,n_flux,n_open,n_private,n_leg,n_tht, n_ext,  &
                                    SIG_open,SIG_closed,SIG_private,SIG_theta,SIG_leg_0,SIG_leg_1,dPSI_open,dPSI_private)
	    endif
	           
          endif
                   
          call plot_grid(node_list,element_list,bnd_elm_list,bnd_node_list,.false.,.false.,'xpoint')
          
        else
          
          call grid_flux_surface(xpoint,xcase, node_list, element_list, surface_list, n_flux, n_tht,     &
                                 xr1, sig1, xr2, sig2,refinement)
          
          call plot_grid(node_list, element_list, bnd_elm_list, bnd_node_list, .true., .false.,'fluxsurface')
          
          ! --- Refine elements (equilibrum)
          if (refinement) then
            n_to_be_refined=0
            call Refine_Elem_List(node_list, element_list, list_to_be_refined, n_to_be_refined)
            call Ref_Update_Index(element_list, node_list)
          end if
             
        end if
        
        ! --- Determine boundary information from the grid
        call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.) 
        
	call export_boundary(node_list, bnd_elm_list, bnd_node_list)
        
        ! --- Compute the plasma equilibrium
        call equilibrium(my_id, node_list, element_list, bnd_node_list, bnd_elm_list, xpoint,xcase, .false.)
        
      end if
      
      ! --- Set initial conditions for time-evolution
      call initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint,xcase)

!      call remove_centre(node_list,element_list,n_tht,67*(n_tht-1))

      ! --- Determine initial energies
      call energy(node_list,element_list,W_mag,W_kin)
      write(*,'(A,12e16.8)') ' initial energies : ', W_mag, W_kin

#ifdef JECCD
      call temp(node_list,element_list,A_tem,A_den,A_jen,A_jec,A_jec1,A_jec2)
      write(*,'(A,12e16.8)') ' initial energies2 : ',A_tem,A_den
      write(*,'(A,12e16.8)') ' initial energies3 : ',A_jen,A_jec
#ifdef JEC2DIAG
      write(*,'(A,12e16.8)') ' initial energies4 : ',A_jec1,A_jec2
#endif
#endif
    end if ! (my_id == 0)
    
#ifdef USE_MUMPS
    ! --- Clean up this instance of mumps (used for equilibrium)
    mumps_par%JOB = -2
    if (my_id == 0) call DMUMPS(mumps_par)
#endif
    if (allocated(pastix_perm_vars))  call tr_deallocate(pastix_perm_vars,"pastix_perm_vars",CAT_UNKNOWN)
    if (allocated(pastix_iperm_vars)) call tr_deallocate(pastix_iperm_vars,"pastix_iperm_vars",CAT_UNKNOWN)
  
  end if if_not_restart
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  ! --- Determine boundary information from the grid
  if ( my_id == 0 ) call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, output_bnd_elements)
  call broadcast_boundary(my_id, bnd_elm_list, bnd_node_list)
  
  ! --- Fill the vacuum response matrices for freeboundary computations
!   if ( freeboundary_equil ) call import_external_fields('coil_field.dat')
  if ( freeboundary ) then
    call get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list, freeboundary_equil,    &
      resistive_wall)
    call update_response(tstep, freeboundary_equil, resistive_wall)
    call import_external_fields('coil_field.dat', my_id)
    if ( (.not. restart) .or. (.not. wall_curr_initialized) ) call init_wall_currents(my_id, resistive_wall)
  end if
  
  call tr_print_memsize("AfterEquilibrium")

  if (RMP_on) then
     if (my_id == 0) then
        call read_RMP_profiles(bnd_node_list)
     endif

  endif
  
  ! --- Broadcast grid information and input parameters to other MPI procs
  call broadcast_elements(my_id, element_list)                ! elements
  if (RMP_on) then
     call broadcast_RMP_profiles(my_id, bnd_node_list)        ! psi_RMP profiles
  endif

  call broadcast_nodes(my_id, node_list)                      ! nodes

  call broadcast_phys(my_id)                                  ! physics parameters
  n_AA = 0  
  do inode = 1, node_list%n_nodes  
    n_AA = max(n_AA,node_list%node(inode)%index(4))  
  end do
  mumps_par%n = n_AA

   ! --- Initialize FFTW
#ifdef USE_FFTW
  call dfftw_plan_dft_r2c_1d(fftw_plan,n_plane,in_fft,out_fft,FFTW_PATIENT)
#endif

! if (RMP_on) then
!    print*, 'bnd_node_list%n_bnd_nodes', bnd_node_list%n_bnd_nodes
!    !print*, 'psi_RMP_cos after broadcast RMP3, my_id', psi_RMP_cos(3), my_id
!    print*, 'psi_RMP_cos after broadcast RMP3, my_id', psi_RMP_cos(bnd_node_list%n_bnd_nodes)
!    !print*, 'dpsi_RMP_cos_dR after broadcast RMP3, my_id', dpsi_RMP_cos_dR(3), my_id
!    print*, 'dpsi_RMP_cos_dR after broadcast RMP3, my_id', dpsi_RMP_cos_dR(bnd_node_list%n_bnd_nodes), my_id
!    !print*, 'dpsi_RMP_cos_dZ after broadcast RMP3, my_id', dpsi_RMP_cos_dZ(3), my_id
!    print*, 'dpsi_RMP_cos_dZ after broadcast RMP3, my_id', dpsi_RMP_cos_dZ(bnd_node_list%n_bnd_nodes), my_id
!    !print*, 'psi_RMP_sin after broadcast RMP3, my_id', psi_RMP_sin(3), my_id
!    print*, 'psi_RMP_sin after broadcast RMP3, my_id', psi_RMP_sin(bnd_node_list%n_bnd_nodes), my_id
!    !print*, 'dpsi_RMP_sin_dR after broadcast RMP3, my_id', dpsi_RMP_sin_dR(3), my_id
!    print*, 'dpsi_RMP_sin_dR after broadcast RMP3, my_id', dpsi_RMP_sin_dR(bnd_node_list%n_bnd_nodes), my_id
!    !print*, 'dpsi_RMP_sin_dZ after broadcast RMP3, my_id', dpsi_RMP_sin_dZ(3), my_id
!    print*, 'dpsi_RMP_sin_dZ after broadcast RMP3, my_id', dpsi_RMP_sin_dZ(bnd_node_list%n_bnd_nodes), my_id
! endif
! 
  call tr_debug_write("JMAIN:End_init elt_list",element_list%n_elements)
  call tr_debug_write("JMAIN:End_init bnd_elt_list",bnd_elm_list%n_bnd_elements)
  call tr_debug_write("JMAIN:End_init node_list",node_list%n_nodes)
  call tr_debug_write("JMAIN:End_init nAA",n_AA)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  call update_equil_state(node_list, element_list, bnd_elm_list, xpoint, xcase, equil_state)
  if ( my_id == 0 ) then
    call print_equil_state(equil_state, .true.)
    call save_special_points(equil_state, 'special_equilibrium_points.dat', .false., ierr)
  end if

  !***********************************************************************
  !*                 end of initilisation/equilibrium                    *
  !***********************************************************************
  
  t_now     = t_start      ! t_now: current time in the simulation
  psi_bnd   = 0.d0
  
  if (nstep > 0) then

    !### THINGS LIKE THIS SHOULD BE REPLACED BY update_equil_state in the future:
    psi_bnd = 0.d0
    if (xpoint) then
      call find_xpoint(my_id,node_list, element_list, psi_xpoint, R_xpoint, Z_xpoint,             &
        i_elm_xpoint, s_xpoint, t_xpoint, xcase, ifail)
      psi_bnd  = psi_xpoint(1)
      if( (xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
        psi_bnd = psi_xpoint(2)
      endif
    else
      call find_limiter(my_id, node_list, element_list, bnd_elm_list, psi_lim, R_lim, Z_lim)
      psi_bnd = psi_lim
    end if
    !###
    
    !*******************************************************
    !*      create groups /communicators		   *
    !* MPI_COMM_N      : group for each harmonic	   *
    !* MPI_COMM_TRANS  : Transversal communicator	   *
    !*   (ie : all first proc of MPI_COMM_N, all second,   *
    !*         all third...)				   *
    !* MPI_COMM_MASTER : group of masters of each harmonic *
    !*  		 (i.e id=0 from each MPI_COMM_N)   *
    !*******************************************************
    if (gmres) then

       N_masters = (n_tor+1)/2
       if (MOD(n_cpu, N_masters) == 0) then
    	  M_cpu = n_cpu / (N_masters)
       else
    	  M_cpu = (n_cpu - MOD(n_cpu, N_masters))/N_masters +1
       end if

       call tr_allocate(i_tor,1,n_cpu,"i_tor",CAT_UNKNOWN)
       
       do i = 1, n_cpu 
    	  i_tor(i) =  MOD(i-1, M_cpu)+1
       end do
       call MPI_COMM_SPLIT(MPI_COMM_WORLD,i_tor(my_id+1),my_id,MPI_COMM_TRANS,ierr)

       do i=1,n_cpu
    	  i_tor(i) = ((i-1) - MOD(i-1, M_cpu))/ M_cpu  + 1
       enddo

       call MPI_COMM_SPLIT(MPI_COMM_WORLD,i_tor(my_id+1),my_id,MPI_COMM_N,ierr)
       
       do i=1,N_masters
    	  i_rank(i) = (i-1) * M_cpu
       enddo
 
       call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
       call MPI_GROUP_INCL(MPI_GROUP_WORLD,N_masters,i_rank,MPI_GROUP_MASTER,ierr)

       call MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_GROUP_MASTER,MPI_COMM_MASTER,ierr)

       call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)		     ! the id of each cpu
       call MPI_COMM_SIZE(MPI_COMM_N, n_cpu_n, ierr)		     ! the number of cpus
       call MPI_COMM_RANK(MPI_COMM_TRANS, my_id_trans, ierr)	     ! the id of each cpu
       call MPI_COMM_SIZE(MPI_COMM_TRANS, n_cpu_trans, ierr)	     ! the number of cpus
       ! TODO : MPI_COMM_MASTER = MPI_COMM_TRANS
       if (my_id_n .eq. 0) then
    	  call MPI_COMM_RANK(MPI_COMM_MASTER, my_id_master, ierr)     ! the id of each cpu
    	  call MPI_COMM_SIZE(MPI_COMM_MASTER, n_cpu_master, ierr)     ! the number of cpus
       endif
    else
       my_id_n = my_id
       MPI_COMM_N = MPI_COMM_WORLD
    endif

    !***********************************************************************
    !*  	  distribute nodes and elements over cpu's		   *
    !***********************************************************************
    if ( use_pastix .and. use_murge .and. use_murge_element .and. gmres ) then
       index_size  = n_cpu_n
       id_elements = my_id_n
    else
       index_size  = n_cpu
       id_elements = my_id
    endif

    call tr_allocate(local_elms,1,element_list%n_elements,"local_elms",CAT_FEM)
    call tr_allocate(index_min,1,index_size,"index_min",CAT_FEM)
    call tr_allocate(index_max,1,index_size,"index_max",CAT_FEM)
    if ( .not. (use_pastix .and. use_murge .and. use_murge_element .and. gmres) ) then
       call tr_allocate(local_index_start,1,n_cpu,"local_index_start",CAT_FEM)
       call tr_allocate(local_index_end,1,n_cpu,"local_index_end",CAT_FEM)
    end if
    !
    ! Construct index_min, index_max and local_elems
    !
    call distribute_nodes_elements(id_elements,index_size,node_list,element_list,local_elms,	  &
    	 n_local_elms,ndof_glob,index_min,index_max)

    node_list%n_dof = ndof_glob
    if ( .not. (use_pastix .and. use_murge  .and. use_murge_element .and. gmres) ) then
       local_index_start = index_min
       local_index_end   = index_max
    end if
    ! Build ijA_index, ijA_size and irn_jcn

    ! TODO : ne pas appeler avec MURGE si pas utile
    call global_matrix_structure(my_id_n,node_List,element_list,bnd_elm_list, freeboundary,&
    	 local_elms,n_local_elms,index_min(id_elements+1),index_max(id_elements+1))

    if ( use_pastix .and. use_murge .and. use_murge_element ) then
       
       write (*,*) "--- Murge initilisation ---"

       ! TODO : deplacer dans un subroutine, dans mod_murge.f90

       ! --- Murge initialisation and graph definition edge by edge
       if (use_murge_element) call murge_initialization(gmres, my_id, MPI_COMM_N, i_tor)
       ! --- Build the graph
       call murge_setgraph(gmres, mumps_par%n, local_elms, n_local_elms,      &
            &              element_list, node_list, n_aa, my_id, my_id_trans, &
            &              n_cpu_trans, MPI_COMM_N, MPI_COMM_TRANS)

    END IF
    if (use_mumps) then
       if (.not. gmres) then
    	  call initialise_mumps(MPI_COMM_WORLD)    ! start MUMPS sparse matrix solver all cpus
       else
    	  call initialise_mumps(MPI_COMM_N)	   ! start MUMPS sparse matrix solver on local groups
       endif
    endif

 endif ! (nstep >0)
  
  ! --- Export a restart file before the first timestep
  if ( (my_id == 0) .and. (.not. restart) ) then
    call export_restart(node_list,element_list,'jorek00000.rst')
  end if

  !***********************************************************************
  !***********************************************************************
  !*                          time stepping                              *
  !***********************************************************************
  !***********************************************************************
  
  if (nstep > 0) call update_deltas(my_id, node_list) ! create list of delta values in local_matrix module

  iter_gmres  = 999
  iter_big    = gmres_max_iter
  iter_prev   = 0

  call tr_print_memsize("BeforeTimeStepping")
  call r3_info_print (-2, -2, 'INITIALIZATION')    ! timing
  
  index_now = index_start  ! index_now: Index of current timestep

  jstep_loop: do jstep = 1, 10 ! Go through the different values of the tstep_n and nstep_n arrays
  istep_loop: do istep = 1, nstep_n(jstep)
    call clck_time(t_itstart)
    t0 = t_itstart

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call flushc !flush the output stream
    call tr_debug_write("JMAIN:Index_now",index_now)

    index_now = index_now + 1
    
    tstep = tstep_n(jstep)
    
    if ( freeboundary ) call update_response(tstep, freeboundary_equil, resistive_wall)
    
    if ( my_id == 0 ) then
      write(*,*) '******************************************************'
      write(*,'(A17,3i7,f14.5,A)') ' *   time step : ',jstep,istep,index_now,tstep,'  *'
      write(*,*) '******************************************************'
    end if

    ! --- Initialise the buffers needed by OpenMP threads. The values of n_tor, 
    ! --- n_plane, n_var have to remain the same until the end of the program.
    call new_thread_buffers()

    call find_axis(99,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

    psi_bnd = 0.d0
    if (xpoint) then
      call find_xpoint(99,node_list, element_list, psi_xpoint, R_xpoint, Z_xpoint,             &
        i_elm_xpoint, s_xpoint, t_xpoint, xcase, ifail)
      psi_bnd  = psi_xpoint(1)
      if( (xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
        psi_bnd = psi_xpoint(2)
      endif
    else
      call find_limiter(99, node_list, element_list, bnd_elm_list, psi_lim, R_lim, Z_lim)
      psi_bnd = psi_lim
    end if
    
    call update_equil_state(node_list, element_list, bnd_elm_list, xpoint, xcase, equil_state)
    if ( my_id == 0 ) call print_equil_state(equil_state, .false.)
    psi_bnd = equil_state%psi_bnd
    
    if (bootstrap) then
      flux_list%n_psi = 1
      call tr_allocate(flux_list%psi_values,1,flux_list%n_psi,"flux_list%psi_values",CAT_GRID)
      flux_list%psi_values(1) = psi_bnd
      call find_flux_surfaces(xpoint,xcase,node_list,element_list,flux_list)
      call find_theta_surface(node_list, element_list, flux_list, 1, 0.0, R_axis, Z_axis,i_elm_find,s_find,t_find,i_find)
      call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
        	     Router,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,  &
        	     Zouter,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      call tr_deallocate(flux_list%psi_values,"flux_list%psi_values",CAT_GRID)
      minRad = Router - R_axis
    else
      minRad = 0.0
    endif
    
    call tr_debug_write("JMAIN:Find_axis_R",R_axis)
    call tr_debug_write("JMAIN:Find_axis_Z",Z_axis)
    call tr_debug_write("JMAIN:Find_axis_T",T_axis)
    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
!    if (my_id .eq. 0) then
!       write(*,FMT_TIMING)  my_id, '# Elapsed time init_time_step :',tsecond
!    end if

    ! Build the matrix 
    call clck_time(t0)
    if (gmres) then
       solve_only = .false.
       if ((gmres) .and. (istep .gt. 1)) then
    	  solve_only = .true.
    	  if (iter_gmres+iter_prev .gt. 2*iter_precon) then			   ! redo preconditioner
    	     solve_only = .false.
    	  endif
       endif
    endif
    
    if (use_pellet) then	    ! calculating the pellet_volume (total_pellet_volume)
      pellet_volume = 3.1415926 * pellet_radius**2 * 2.d0 * 3.1415926535 * pellet_R
      call Integrals_3D(my_id, node_list,element_list,density_tot,density_in,density_out,pressure_tot,pressure_in,pressure_out)
    endif
    call tr_debug_write("JMAIN:Debconstruct_n_elms",n_local_elms)
    
    ! --- construct the matrix from elemental matrices
    if ( use_pastix .and. use_murge .and. use_murge_element ) then

       call construct_matrix_murge(my_id, node_list, element_list, bnd_node_list, local_elms,      &
    	 n_local_ELms, xpoint, xcase, minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint,         &
         Z_xpoint, psi_xpoint, gmres, i_tor, n_cpu, mpi_comm_n, mpi_comm_trans, my_id_trans,       &
         n_cpu_trans, solve_only)
    else

       call construct_matrix(my_id, local_elms, n_local_ELms, index_min(my_id+1),                  &
         index_max(my_id+1), xpoint, xcase, minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint,   &
         Z_xpoint, psi_xpoint)
    endif
    

    call clck_time_barrier(t1)
    if (my_id .eq. 0) then
       call clck_ldiff(t0,t1,tsecond)
      write(*,FMT_TIMING) my_id, '# Elapsed time construct_matrix :',tsecond
    endif     

    ! --- Free the buffers needed by OpenMP threads (ELM-RHS etc.)
    call del_thread_buffers()

    if (.not. gmres) then
       if (use_mumps) then

    	  call solve_mumps_all(my_id)

       else

    	  ! Recuperer la solution
    	  if (use_murge) then
    	     call solve_murge_all(n_cpu,my_id,index_min(my_id+1),index_max(my_id+1), i_tor, gmres, my_id_n, mpi_comm_n, mpi_comm_master)
    	  else
    	     call solve_pastix_all(n_cpu,my_id,index_min(my_id+1),index_max(my_id+1))
    	  endif

       endif

    else
       call clck_time(t0)
       if (.not. solve_only) then
    	  ! with murge elementary assembly harmonic distribution is already done.
    	  IF ( .not. ( use_pastix .and. use_murge .and. use_murge_element ) ) THEN
    	     call distribute_harmonics(my_id,my_id_n,n_cpu)
    	  ELSE
    	     call distribute_vector(my_id,rhs_glob,mumps_par%rhs,.false.)	       
    	  END IF
       else
          call distribute_vector(my_id,rhs_glob,mumps_par%rhs,.true.)	       
       endif
       call clck_time_barrier(t1)
       call clck_ldiff(t0,t1,tsecond)
       if (my_id .eq. 0) then
          write(*,FMT_TIMING) my_id, '# Elapsed time distribute :',tsecond
       end if

       call clck_time(t0)
       if (use_murge .and. use_murge_element) then
    	  call solve_murge_all(n_cpu,my_id,index_min(my_id_n+1),index_max(my_id_n+1), i_tor, gmres, my_id_n, mpi_comm_n, mpi_comm_master)
       else
    	  call solve_matrix_n(my_id,i_tor,MPI_COMM_N,MPI_COMM_MASTER,solve_only)    ! factorise preconditioning matrices
       end if
       call clck_time_barrier(t1)
       call clck_ldiff(t0,t1,tsecond)
       if (my_id .eq. 0) then
          write(*,FMT_TIMING) my_id, '# Elapsed time first solve :',tsecond
       end if
    endif

    call clck_time_barrier(t0)
    if (gmres) then
      iter_prev = iter_gmres
      iter_gmres = gmres_max_iter
      call gmres_driver(my_id,my_id_n,i_tor, n_tor,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)
    endif
    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0) then
       write(*,FMT_TIMING)  my_id, '# Elapsed time gmres/solve :',tsecond
    end if

    call clck_time_barrier(t0)
    if ( (gmres .and. (iter_gmres .lt. iter_big)) .or. (.not.gmres) ) then

       if (use_pellet) then
         pellet_volume = total_pellet_volume
         call update_pellet(my_id,node_list,element_list)

!         if (my_id == 0) then
!           xtime_pellet_R(index_now)         = pellet_R
!           xtime_pellet_Z(index_now)         = pellet_Z
!           xtime_pellet_particles(index_now) = pellet_particles
!           xtime_phys_ablation(index_now)    = phys_ablation
!          endif
       endif

#if JOREK_MODEL == 500
       call update_mgi(my_id,node_list,element_list)
#endif

       call update_values(my_id,element_list,node_list,deltas)         ! add solution to node values
       call update_deltas(my_id,node_list)
 
       !***********************************************************************
       !*                          output saving for diagnostics              *
       !*                                                                     *
       !*  ===> set boolean "save_diagnostics_HDF5" to "true" if wanted       *
       !*       in the input file, to "false" if not                          *
       !*  ===> the diagnostics are saved every "h5_diag_nbtime" Alfven times *
       !***********************************************************************
#ifdef USE_HDF5
       !*   0D-1D and 2D diagnostics saving in HDF5 format  *
       if (save_diagnostics_HDF5 .and. (my_id .eq. 0) ) then
          ! the number of HDF5 files that:
          !   - have been written in the previous runs: "h5_nbsave_previous"
          !   - have been written so far: "h5_nbsave_current"
          !   - should be written if everything goes right: "h5_nbsave"
          ! are computed above, around line 288
          if ( ( (   (mod( floor(t_now),floor(h5_diag_nbtime) ).eq.0).or.(t_now.eq.t_this)   ) &
               .and. (h5_nbsave_current .le. h5_nbsave) ) &
               .or. (h5_nbsave_previous .eq. 0)         ) then
             write(*,*) ' '
             write(*,*) '*******************************************************************************'
             write(*,*) '*     BEGIN --- writing the HDF5 diagnostics                                  *'
             write(*,*) '*******************************************************************************'
             ! compute quantities in (R,Z) and (psi,theta) coordinates
!             call HDF5_compute_R_Z_psi_th()
             write(*,*) ' ===> writing the basic parameters.............................................'
             call HDF5_basics_save(index_now,t_now)
             write(*,*) ' ===> writing the n_tor profiles...............................................'
             call HDF5_ntor_profiles_save(index_now)
             write(*,*) ' ===> writing the radial (psi) profiles........................................'
             !call HDF5_radial_profiles_save(index_now)
             write(*,*) '*******************************************************************************'
             write(*,*) '*     END --- writing the HDF5 diagnostics                                    *'
             write(*,*) '*******************************************************************************'
             write(*,*) ' '
          endif
          h5_nbsave_current = h5_nbsave_current + 1
          ! this quantity is now saved in the restart file and becomes the new h5_nbsave_previous
          h5_nbsave_all     = h5_nbsave_previous + h5_nbsave_current
       endif
#endif

          t_now = t_now + tstep

       else
          if ( my_id == 0 ) then
             write(*,*)
             write(*,'(a,i6.6,a)') '>>>>> NO CONVERGENCE AFTER ', iter_gmres, ' ITERATIONS. ABORTING <<<<<'
             write(*,*)
          end if
          index_now = index_now - 1 ! Undo the time step
          exit jstep_loop
       end if

    !-------------------------------------------------------- adapt time step (in progress...)
    mindelta = minval(deltas); maxdelta = maxval(deltas);

    if (gmres .and. adaptive_time) then        ! experimental
       if (iter_gmres .ge. iter_big) then
    	  tstep = tstep /2.d0
    	  write(*,*) my_id,' REDUCTION TIMESTEP : ',tstep
       elseif (max(abs(mindelta),abs(maxdelta)) .gt. 0.05) then
    	  !	 tstep = tstep /2.d0
    	  !	 iter_gmres = 99999
    	  !	 write(*,*) my_id,' REDUCTION TIMESTEP : ',tstep
       elseif (max(abs(mindelta),abs(maxdelta)) .lt. 0.001) then
    	  !	 tstep = tstep * 2.d0
    	  !	 iter_gmres = 99999
    	  !	 write(*,*) my_id,' INCREASE TIMESTEP : ',tstep
       endif
    endif

    !--------------------------------------------------------- energies
    if ( (my_id == 0) .and. (.not. bench_without_plot) ) then
       call energy(node_list,element_list,W_mag,W_kin)

       xtime(index_now) = t_now
       energies(1:n_tor,1,index_now) = W_mag(1:n_tor)
       energies(1:n_tor,2,index_now) = W_kin(1:n_tor)

#ifdef JECCD
       call temp(node_list,element_list,A_tem,A_den,A_jen,A_jec,A_jec1,A_jec2)
       write(*,'(A,12e16.8)') ' current energies2 : ',A_tem,A_den
       write(*,'(A,12e16.8)') ' current energies3 : ',A_jen,A_jec
#ifdef JEC2DIAG
       write(*,'(A,12e16.8)') ' current energies4 : ',A_jec1,A_jec2
#endif

       energies2(1:n_tor,1,index_now) = A_tem(1:n_tor)
       energies2(1:n_tor,2,index_now) = A_den(1:n_tor)

       energies3(1:n_tor,1,index_now) = A_jen(1:n_tor)
       energies3(1:n_tor,2,index_now) = A_jec(1:n_tor)

#ifdef JEC2DIAG
       energies4(1:n_tor,1,index_now) = A_jec1(1:n_tor)
       energies4(1:n_tor,2,index_now) = A_jec2(1:n_tor)
#endif

       write(6,*) ' exiting current energies '
#endif
       
       ! --- Output some information about the current timestep
       130 format(1x,a,i5.5,a,es10.3,a)
       131 format(1x,a,2(2(es10.2,' ...',es10.2,',')))
       132 format(1x,'-------------------------------------------------------------------')
       133 format(1x,a,2(es10.2,' at ',i10,','))
       write(*,*)
       write(*,132)
       write(*,130) 'After step ', istep, ' (t_now=', t_now, '):'
       write(*,132)
       write(*,133) 'min,max deltas  =', mindelta, minloc(deltas), maxdelta, maxloc(deltas)
       write(*,131) 'W_mag,_kin      =', W_mag(1), W_mag(n_tor), W_kin(1), W_kin(n_tor)
       Growth_mag  = 0.d0; Growth_kin  = 0.d0; Growth_mag0 = 0.d0; Growth_kin0 = 0.d0
       if (index_now > index_start+1) then
    	 Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,index_now)/energies(n_tor,1,index_now-1)))/ tstep
    	 Growth_kin  = 0.5d0*log(abs(energies(n_tor,2,index_now)/energies(n_tor,2,index_now-1)))/ tstep
    	 Growth_mag0 = 0.5d0*log(abs(energies(1,1,index_now)/energies(1,1,index_now-1)))/ tstep
    	 Growth_kin0 = 0.5d0*log(abs(energies(1,2,index_now)/energies(1,2,index_now-1)))/ tstep
    	 write(*,131) 'Growth_mag,_kin =', Growth_mag0, Growth_mag, Growth_kin0, Growth_kin
       endif
       write(*,132)
       write(*,*)

       ! --- Output energies and growth_rates to text files during the code run
       if ( .not. bench_without_plot ) call write_live_data(index_now)
#ifdef JECCD
       if ( .not. bench_without_plot ) call write_live_data2(index_now)
       if ( .not. bench_without_plot ) call write_live_data3(index_now)
#ifdef JEC2DIAG
       if ( .not. bench_without_plot ) call write_live_data4(index_now)
#endif
#endif
    endif

    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0) then
       write(*,FMT_TIMING)  my_id, '# Diagnostics :',tsecond
    end if
    !---------------------------------------------------------timing
    if ( istep == 1 ) then
       call r3_info_print (-3, -2, 'ITERATION	 1')
    else
       call r3_info_print (istep, -2, 'ITERATION')
    endif
    write(itlabel,'(I8)') istep
    call tr_print_memsize("AfterIter"//itlabel)
    
    ! --- Write a restart file every nout timesteps
    if ( (my_id == 0) .and. (mod(index_now,nout) == 0) ) then
      write(fileout,'(A5,i5.5,A4)') 'jorek',index_now,'.rst'
      call export_restart(node_list,element_list,fileout)
    endif
    
    ! --- Exit the code if a file "STOP_NOW" exists in the run directory.
    inquire(file='STOP_NOW', exist=file_exists)
    if ( file_exists ) then
      if ( my_id == 0 ) then
    	write(*,*)
    	write(*,*) '>>>>> FOUND FILE STOP_NOW: EXITING THE CODE <<<<<'
    	write(*,*)
      end if
      exit jstep_loop
    end if
    
    ! --- Exit the code if NaNs are detected.
    if ( allocated(deltas) ) then
      sum_deltas = sum(deltas)
      if ( sum_deltas /= sum_deltas ) then
    	write(*,*)
    	write(*,*) '>>>>> NaNs DETECTED: EXITING THE CODE <<<<<'
    	write(*,*)
        exit jstep_loop
      end if
    end if
    
    call clck_time_barrier(t1)
    call clck_ldiff(t_itstart,t1,tsecond)
    if (my_id .eq. 0) then
       write(*,FMT_TIMING)  my_id, '# Elapsed time ITERATION :',tsecond
    end if

  enddo istep_loop
  enddo jstep_loop
  
  !***********************************************************************
  !*                         cleanup  (solvers)                          *
  !***********************************************************************

  if (nstep .gt.0) then
    if (use_mumps) then
#ifdef USE_MUMPS
      mumps_par%JOB = -2                            ! clean up this instance of mumps
      call DMUMPS(mumps_par)
#endif
    elseif (use_pastix) then
       if ( use_murge ) then 
    	 call murge_termination(gmres)
       else
    	  pastix_iparm(2)     = 7			! Clean-up
    	  pastix_iparm(3)     = 7

    	  if (.not. gmres) then

    	     call pastix_fortran(pastix_data,MPI_COMM_WORLD,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
    		  pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

    	  elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0))  ) then

            call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,&
                  DUMMY_INT, DUMMY_INT, DUMMY_REAL, &
                  pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
   
       	  endif
          
       end if

    elseif (use_wsmp) then

#ifdef USE_WSMP
      call PWGSMP__deallocate()
#endif

    endif
  endif
  
  ! --- Close open files
  if ( (my_id == 0) .and. (.not. bench_without_plot) ) call finalize_live_data()
#ifdef JECCD
  if ( (my_id == 0) .and. (.not. bench_without_plot) ) call finalize_live_data2()
  if ( (my_id == 0) .and. (.not. bench_without_plot) ) call finalize_live_data3()
#ifdef JEC2DIAG
  if ( (my_id == 0) .and. (.not. bench_without_plot) ) call finalize_live_data4()
#endif
#endif

  !***********************************************************************
  !*                          plots etc.                                 *
  !***********************************************************************

  if (my_id .eq. 0)  then

    call export_restart(node_list,element_list,'jorek_restart.rst')

    if (.not. bench_without_plot) then
       
       do ivar=1,n_var
    	  call plot_solution(node_list,element_list,ivar,-1,1,variable_names(ivar))
       enddo

       do i=1,n_tor,2
          write(label,'(A4,i3,A1)') '(n =',((i-1)/2)*n_period,')'

    	  do ivar=1,n_var
          if ((ivar .ne. 3) .and. (ivar .ne. 4)) then
             call plot_solution(node_list,element_list,ivar,i,1,variable_names(ivar)//label)
          endif
    	  enddo

       enddo
    endif

    if (index_now .gt. 1) then

       E_min =  1.d20
       E_max = -1.d20
       E_max = max(E_max,maxval(energies(1,2,1:index_now)))
       E_min = min(E_min,minval(energies(1,2,1:index_now)))
       do i=2,n_tor
    	  E_max = max(E_max,maxval(energies(i,1,1:index_now)))
    	  E_min = min(E_min,minval(energies(i,1,1:index_now)))
    	  E_max = max(E_max,maxval(energies(i,2,1:index_now)))
    	  E_min = min(E_min,minval(energies(i,2,1:index_now)))
       enddo

       call nframe(1,1,2,xtime(1),xtime(index_now),E_min,E_max,'energies',7,'time',4,' ',1)

       do i=1,n_tor
	 if (mod(i,2) .eq. 0) then
	   call lincol(mod(i/2,10))
	 else
	   call lincol(mod((i-1)/2,10))
	 endif
    	  call lplot(1,1,2,xtime(1:index_now),energies(i,1,1:index_now),-index_now,1,'Magnetic Energie',16,'time',4,'Emag',4)
    	  call lincol(4)
    	  if (n_tor .eq. 3) call lincol(2)
    	  call lplot(1,1,2,xtime(1:index_now),energies(i,2,1:index_now),-index_now,1,'Kinetic Energie',15,'time',4,'Ekin',4)
       enddo
       call lincol(3)
       call lplot(1,1,2,xtime(1:index_now),energies(1,2,1:index_now),-index_now,1,'Kinetic Energie',15,'time',4,'Ekin',4)
       call lincol(0)
    endif

!---------------------------------------------- plot equilibrium current profile (to be removed)
    call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis, ifail)

    nplot = 501
    call tr_allocate(xp,1,nplot,"xp",CAT_GRID)
    call tr_allocate(yp1,1,nplot,"yp1",CAT_GRID)
    call tr_allocate(yp2,1,nplot,"yp2",CAT_GRID)
    call tr_allocate(yp3,1,nplot,"yp3",CAT_GRID)
! ---- plot neoclassical coefficients -----
    if (NEO) then
       call tr_allocate(mu_neo,1,nplot,"mu_neo",CAT_GRID)
       call tr_allocate(ki_neo,1,nplot,"ki_neo",CAT_GRID)
    endif
    iplot = 0

    if (xpoint) then
       call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
    else
       psi_bnd = psi_xpoint(1)
       if( (xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
    	 psi_bnd = psi_xpoint(2)
       endif
    endif

    Rp_start = R_axis - amin*2.d0
    Rp_end   = R_axis + amin*2.d0

    Zp = Z_axis

    do i=1,nplot

       Rp =  Rp_start + float(i-1)/float(nplot-1) * (Rp_end - Rp_start)

       call find_RZ(node_list,element_list,Rp,Zp,R_out,Z_out,i_elm,s_out,t_out,ifail)

       if (ifail .eq. 0) then

    	  call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi,P_s,P_t,P_st,P_ss,P_tt)

    	  call density(    xpoint,xcase, Zp, Z_xpoint, psi,psi_axis,psi_bnd,	       &
    	       zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz)
    	  if (jorek_model .eq. 400) then	     
    	    call temperature_i(xpoint,xcase, Zp, Z_xpoint, psi,psi_axis,psi_bnd, &
    			     zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2,dTi_dpsi2_dz)			   
    	    call temperature_e(xpoint,xcase, Zp, Z_xpoint, psi,psi_axis,psi_bnd, &
    	     zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz)	     
      zT = zTi + zTe
    	    dT_dpsi = dTi_dpsi + dTe_dpsi	    
    	  else
      call temperature(xpoint,xcase, Zp, Z_xpoint, psi,psi_axis,psi_bnd, &
    		   zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)
    	  endif
    	  call FFprime(    xpoint,xcase, Zp, Z_xpoint, psi,psi_axis,psi_bnd,	       &
    	       zFFprime,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2,dFFprime_dpsi_dz)

       if (NEO) then
          if (num_neo_file) then
             call neo_coef (xpoint, xcase, Zp, Z_xpoint, psi, psi_axis,psi_bnd, &
                  amu_neo_node, aki_neo_node)
          endif
       endif

    	  zjz	= (zFFprime - Rp*Rp * (zn * dT_dpsi + dn_dpsi * zT)) / Rp

    	  iplot = iplot + 1

    	  xp(iplot)  = Rp
    	  yp1(iplot) = zFFprime / Rp
    	  yp2(iplot) = zjz
    	  yp3(iplot) = - Rp*Rp * (zn * dT_dpsi + dn_dpsi * zT) / Rp

    	  !	 write(*,'(A,8e16.8)') ' profiles : ',xp(iplot),psi,psi_axis,psi_bnd,yp2(iplot),yp1(iplot),yp3(iplot)
           if (NEO) then
              if ( num_neo_file) then
                 mu_neo(iplot) = amu_neo_node
                 ki_neo(iplot) = aki_neo_node
                 write(*,'(A,8e16.8)') ' profiles : ',xp(iplot),psi,psi_axis,psi_xpoint,mu_neo(iplot),ki_neo(iplot)
              endif
           endif

        endif

    enddo

    call lplot6(1,1,xp,yp2,iplot,' ')
    call lincol(1)
    call lplot6(1,1,xp,yp1,-iplot,' ')
    call lincol(2)
    call lplot6(1,1,xp,yp3,-iplot,' ')
    call lincol(0)
    if (NEO) then
       if ( num_neo_file) then
          call lplot6(1,1,xp,mu_neo,iplot,' ')
          call lincol(1)
          call lplot6(1,1,xp,ki_neo,iplot,' ')
          call lincol(0)
       end if
    endif
    call finplt 					 ! close plot file

!  cll export_POV(node_list,element_list,3,1)	       ! export to POVray native bezier patch format

    call export_helena(node_list,element_list,bnd_elm_list)

    if (allocated(energies))  call tr_deallocate(energies,"energies",CAT_UNKNOWN)
    if (allocated(xtime))     call tr_deallocate(xtime,"xtime",CAT_UNKNOWN)

#ifdef JECCD
    if (allocated(energies2)) call tr_deallocate(energies2,"energies2",CAT_UNKNOWN)
    if (allocated(energies3)) call tr_deallocate(energies3,"energies3",CAT_UNKNOWN)
#ifdef JEC2DIAG
    if (allocated(energies4)) call tr_deallocate(energies4,"energies4",CAT_UNKNOWN)
#endif
#endif
  endif
  
#ifdef USE_FFTW
  call dfftw_destroy_plan(fftw_plan)
#endif

  call r3_info_summary ()                                ! timing
  call MPI_FINALIZE(IERR)                                ! clean up MPI

end program JOREK2


! --- The following comment block defines the start page of the Doxygen code documentation ---
!
!> \mainpage
!!
!! \section JOREK About JOREK
!!
!! JOREK solves the (reduced) MHD equations in 3D toroidal geometry using a discretization with
!! Bezier finite elements in the poloidal plane and a Fourier expansion in toroidal direction.
!!
!! \section DOCU About this documentation
!!
!! This documentation covers the code structure, i.e.,
!! - <a href="dirs.html">directories</a>,
!! - <a href="files.html">source files</a>,
!! - <a href="annotated.html">Fortran modules</a>,
!! - <a href="globals_func.html">subroutines and functions</a>, and
!! - routine parameters.
!!
!! <img src="JOREK_DOC.jpeg" style="border:1px solid black;float:right"/>
!! Background information on the following topics is covered in a seperate documentation
!! available in the docu/ folder of the JOREK repository as LaTeX source or
!! <a href="JOREK_DOC.pdf">online as PDF</a>:
!! - Implemented equations
!! - Numerical methods
!! - Using the version control system Subversion
!! - Compiling and running the code
!! - Bezier elements
!! - Free boundary extension
!! - ...
!! 
!! \section Repository Browse Subversion Repository online
!! 
!! The <a href="http://subversion.apache.org/">Subversion</a> repository of JOREK can
!! be <a href="https://gforge.inria.fr/scm/viewvc.php?view=rev&root=aster">browsed online</a>
!! after logging in if you already have an account for it. On the Server,
!! <a href="http://www.viewvc.org/">ViewVC</a> is installed for this purpose.
!!
!! \section Doxygen Code documentation with Doxygen
!!
!! This documentation is generated directly from the source code using
!! <a href="http://www.doxygen.org">Doxygen</a>.
!! For this to work properly, you should follow some
!! simple rules when writing comments in the code.
!!
!! - An example for the proper documentation of a subroutine:                               \n
!!                                                                                          \n
!! <code>
!! !\> Brief documentation for the test routine                                             \n
!! !!                                                                                       \n
!! !! More details on the functionality of the                                              \n
!! !! test routines can be added like this.                                                 \n
!! !!                                                                                       \n
!! subroutine test(i,r)                                                                     \n
!!                                                                                          \n
!!   ! --- Routine parameters                                                               \n
!!   integer, intent(in)  :: i !< Information on parameter i                                \n
!!   real*8,  intent(out) :: r !< Information regarding parameter r                         \n
!!                                                                                          \n
!!   ! --- Local variables                                                                  \n
!!   integer :: k   ! Information on k (Local variables are not documented by Doxygen)      \n
!!                                                                                          \n
!!   ...                                                                                    \n
!!                                                                                          \n
!! end subroutine test                                                                      \n
!! </code>
!!
!! - The documentation of a module works in the same way:                                   \n
!!                                                                                          \n
!! <code>
!! !\> Brief documentation for the test module                                              \n
!! !!                                                                                       \n
!! !! More details on the functionality of the                                              \n
!! !! test module can be added like this.                                                   \n
!! !!                                                                                       \n
!! module some_test_module                                                                  \n
!!                                                                                          \n
!!  implicit none                                                                           \n
!!                                                                                          \n
!!  !\> \@name Rectangular Grid                                                             \n
!!  integer :: n_R               !< Number of grid points in R-direction                    \n
!!  integer :: n_Z               !< Number of grid points in Z-direction                    \n
!!                                                                                          \n
!!  !\> \@name Polar Grid                                                                   \n
!!  !! Parameters defining a non flux-aligned polar grid in the poloidal plane.             \n
!!  integer :: n_radial          !< Number of radial grid points                            \n
!!  integer :: n_pol             !< Number of poloidal grid points                          \n
!!                                                                                          \n
!!  contains                                                                                \n
!!                                                                                          \n
!!  !\> Description for routine                                                             \n
!!  subroutine test()                                                                       \n
!!  ...                                                                                     \n
!! </code>
!! Note, that the \@name command allows to create groups of variables. An example for this
!! is the module ::phys_module.
!!
!! - HTML markup is possible, e.g., \<b\>some text\</b\> will display in bold face:
!!   <b>some text</b>
!!
!! - For further information, refer to the
!!   <a href="http://www.stack.nl/~dimitri/doxygen/manual.html">Doxygen manual</a>.
!!
