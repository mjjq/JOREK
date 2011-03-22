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
  
  use mumps_module
  use pastix_module
  use murge_module,        only: murge_initialization, murge_setGraph, MURGE_Clean, use_murge,     &
    use_murge_element, murge_initialised, murge_harmonic, murge_glob2loc, murge_loc2glob, murge_id
  use data_structure
  use phys_module
  use global_distributed_matrix
  use nodes_elements
  use boundary,            only: boundary_from_grid
  use vacuum_response,     only: vacuum_preset, vacuum_init, get_vacuum_response
  use vacuum_equilibrium,  only: import_external_fields
  use live_data,           only: init_live_data, write_live_data, finalize_live_data
  
  implicit none
  
  include 'mpif.h'
#include "r3_info.h"
  
  interface
    subroutine gmres_driver(my_id,my_id_n,i_tor,n_tor,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)
      integer :: i_tor(:), my_id, my_id_n, MPI_COMM_N, MPI_COMM_MASTER
      integer :: iter_gmres, n_tor
    end subroutine gmres_driver
    
    subroutine equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint2)
      use data_structure
      integer(kind=4),             intent(in)    :: my_id
      type (type_node_list),       intent(inout) :: node_list
      type (type_element_list),    intent(inout) :: element_list
      type (type_bnd_node_list)   ,intent(inout) :: bnd_node_list    
      type (type_bnd_element_list),intent(inout) :: bnd_elm_list    
      logical(kind=4),             intent(in)    :: xpoint2
    end subroutine equilibrium
     
    subroutine update_rhs_n(my_id,my_id_n,i_tor,mpi_comm_master)
      integer(kind=4) :: my_id
      integer(kind=4) :: my_id_n
      integer(kind=4) :: i_tor(:)
      integer(kind=4) :: mpi_comm_master
    end subroutine update_rhs_n
    
    subroutine construct_matrix_murge(my_id,node_list,                  &
      element_list,local_elms,n_local_elms,xpoint2,psi_axis,psi_bnd,    &
      z_xpoint,gmres,i_tor,n_cpu,mpi_comm_n,mpi_comm_trans,my_id_trans, &
      n_cpu_trans,solve_only)
      use data_structure, only : type_node, type_element,               &
        type_element_list, type_node_list
      integer(kind=4) :: n_cpu
      integer(kind=4), target :: n_local_elms
      integer(kind=4), target :: my_id
      type (type_node_list), target :: node_list
      type (type_element_list), target :: element_list
      integer(kind=4), target :: local_elms(n_local_elms)
      logical(kind=4), target :: xpoint2
      real(kind=8), target :: psi_axis
      real(kind=8), target :: psi_bnd
      real(kind=8), target :: z_xpoint
      logical(kind=4), target :: gmres
      integer(kind=4) :: i_tor(n_cpu)
      integer(kind=4) :: mpi_comm_n
      integer(kind=4), target :: mpi_comm_trans
      integer(kind=4), target :: my_id_trans
      integer(kind=4), target :: n_cpu_trans
      logical(kind=4), target :: solve_only
    end subroutine construct_matrix_murge
  end interface
  
  type (type_surface_list) :: surface_list
  real*8                   :: W_mag(n_tor), W_kin(n_tor), growth_mag, growth_kin, growth_mag0, growth_kin0
  real*8                   :: t_matrix_0, t_matrix_1, PI
  real*8                   :: t_send_0, t_send_1, t_solve_0, t_solve_1, t_solve_2
  real*8                   :: psi_bnd, psi_axis, R_axis, Z_axis, s_axis, t_axis
  real*8                   :: psi_xpoint, R_xpoint, Z_xpoint, s_xpoint, t_xpoint, mindelta, maxdelta
  integer                  :: my_id, my_id_n, my_id_master
  integer                  :: istep,jstep,ierr,i,itor,inode, i_elm_axis, i_elm_xpoint
  integer                  :: n_local_ELMs, index_total
  integer                  :: i_rank(n_tor), n_cpu, n_cpu_n, n_cpu_master, m_cpu, n_masters, n_cpu_trans, my_id_trans
  integer                  :: iter_gmres
  integer                  :: MPI_COMM_N, MPI_GROUP_MASTER, MPI_GROUP_WORLD, MPI_COMM_MASTER, MPI_COMM_TRANS
  character*8              :: label
  character*14             :: fileout
  integer                  :: required,provided,StatInfo
  integer, allocatable     :: local_elms(:), i_tor(:), index_min(:), index_max(:)
  real*8                   :: zjz, E_min, E_max
  logical                  :: gmres, solve_only, adaptive_time
  integer*4                :: rank, comm_size 
  real*8 :: zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz
  real*8 :: zT, dT_dpsi, dT_dz, dT_dpsi2, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz
  real*8 :: zFFprime, dFFprime_dpsi, dFFprime_dz, dFFprime_dpsi_dz,dFFprime_dpsi2,dFFprime_dz2
  real*8 :: Rp, Zp, R_out,Z_out,s_out,t_out,P_s,P_t,P_st,P_ss,P_tt, psi
  real*8 :: Rp_start, Rp_end
  real*8,allocatable       :: xp(:), yp1(:), yp2(:), yp3(:)
  integer                  :: nplot, iplot, i_elm, ifail, ivar, iter_big, iter_precon, n_aa
  logical                  :: is_local
  integer                  :: i_elem, inode1, i_order, index_node1
  type (type_element)      :: element
  integer                  :: index_size, id_elements
  integer                  :: list_to_be_refined(n_ref_list), n_to_be_refined    
  logical                  :: bench_without_plot
  integer                  :: t0,t1,nb_periodes_max,nb_periodes_sec, nb_periods
  character(len=20), parameter :: FMT_TIMING = "(I2,A70,F7.2)"
  
  !***********************************************************************
  !*                  intialisation                                      *
  !***********************************************************************
  
  ! --- Initialise MPI / threaded MPI
  !call MPI_INIT(IERR)
  required=MPI_THREAD_MULTIPLE
  call MPI_Init_thread(required,provided,StatInfo)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)      ! id of each MPI proc
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr) ! number of MPI procs
  my_id = rank
  n_cpu = comm_size
  
  if (my_id == 0) then
    write(*,*) '****************************************'
    write(*,*) '*   3D Reduced MHD : JOREK_2.0         *'
    write(*,*) '****************************************'
    write(*,*) ' number of cpus : ',n_cpu
  end if
  
  ! --- Initialise timing
  call system_clock(count_rate=nb_periodes_sec, count_max=nb_periodes_max)
  call r3_info_init ()
  
  ! --- Select solver
  gmres              = .true.             ! .true. for gmres, .false. for direct solver
  use_mumps          = .false.            ! Use MUMPS solver
  use_pastix         = (.not. use_mumps)  ! Use PASTIX solver
  use_murge          = .false.            ! Use MURGE interface to PASTIX solver
  use_murge_element  = .false.
  pastix_initialised = .false.
  pastix_analysed    = .false.
  murge_initialised  = .false.
  pastix_smp_only    = .false.            ! Implies that each MPI group resides within one node!
  if (n_tor == 1) gmres = .false.
  
  refinement         = .false.            ! Enable mesh refinement?
  
  adaptive_time      = .false.            ! Requires no_mpi for Pastix library
  
  ! --- Flag from HSLT
  bench_without_plot              = .false.    ! .true. for benchmark (mesuring elapsed time without plot phases) 
  use_matrix_whitout_zeros_pastix = .false.    ! .true. to remove nonzeros in the preconditioning matrix with MUMPS
  use_matrix_whitout_zeros_mumps  = .false.    ! .true. to remove nonzeros in the preconditioning matrix with PaStiX
  
  ! --- Preset input parameters to reasonable defaults, then read the input file.
  call vacuum_preset(my_id, freeboundary_equil, freeboundary, use_starwall, resistive_wall)
  call initialise_parameters(my_id)
  call vacuum_init(my_id, freeboundary_equil, freeboundary, use_starwall, resistive_wall)
  
  ! --- Fill the arrays mode (toroidal mode number n) and mode_type (cos or sin).
  do itor=1, n_tor
    mode(itor)        = int(itor / 2) * n_period
    if ( (itor==1) .or. (mod(itor,2)==0) ) then
      mode_type(itor) = 'cos'
    else
      mode_type(itor) = 'sin'
    end if
  end do
  
  ! --- Write out all parameters defined in mod_parameters and by the input file.
  call log_parameters(my_id)
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  ! --- Some checks not to waste any cpu time
  if (required.ne.provided) then
    write(*,*) 'FATAL : MPI_THREAD_MULTIPLE (provided < required)', my_id, required, provided
    call MPI_FINALIZE(IERR)
    stop
  end if
  
  if ( (.not. use_mumps) .and. (.not. use_pastix) ) then
    write(*,*) ' FATAL : specify a valid solver'
    call MPI_FINALIZE(IERR)
    stop
  end if
  
  if ((n_plane < n_tor + 1) .and. (n_tor > 1)) then
    write(*,*) ' FATAL : n_plane too small ',n_plane,n_tor
    call MPI_FINALIZE(IERR)
    stop
  end if
  
  if ( gmres .and. (nstep > 0) .and. (mod(n_cpu,(n_tor-1)/2+1) /= 0) ) then
    write(*,'(A,i4,A,i4,A)') ' FATAL : need a multiple of ',(n_tor-1)/2+1,' cpus for ',            &
      (n_tor-1)/2+1,' harmonics'
    call MPI_FINALIZE(IERR)
    stop
  end if
  
  ! --- Open files which will be filled during the code run
  if ( (my_id == 0) .and. (.not. bench_without_plot) ) call init_live_data()
  
  ! --- Initialise ppplib plotting library
  if (my_id == 0)  call begplt('jorek2.ps')
  
  ! --- Define the basis functions at the Gaussian points
  call initialise_basis()
  
  ! --- Read the restart file to continue a previous JOREK run (if restart is .true.)
  if ( restart .and. (my_id == 0) ) then
    
    ! --- Read the restart file (jorek_restart.rst)
    call import_restart(node_list, element_list)
    tstep = tstep_in
    
    ! --- Output energies and growth_rates to text files (values from restart file)
    do index_now = 1, index_start
      call write_live_data(index_now)
    end do
    
    ! --- Optional: Redo flux aligned grid (DOES NOT WORK CURRENTLY)
    if (regrid) then
      if (xpoint)  then
        call grid_xpoint(node_list, element_list, n_flux, n_open, n_private, n_leg, n_tht,         &
          SIG_open, SIG_closed, SIG_private, SIG_theta, SIG_leg_0, SIG_leg_1, dPSI_open,           &
          dPSI_private)
      else
        call grid_flux_surface(xpoint, node_list, element_list, surface_list, n_flux, n_tht, xr1,  &
          sig1, xr2, sig2)
      end if
    end if
    
  end if
  
  !***********************************************************************
  !*                  define grid / equilibrium                          *
  !***********************************************************************
  
  if_not_restart: if (.not. restart) then
    
    element_list%n_elements      = 0
    bnd_elm_list%n_bnd_elements  = 0
    node_list%n_nodes            = 0
    
    if (my_id == 0) then
      
      call define_boundary()
      
      if ((n_R > 0) .and. (n_Z > 0) .and. (n_radial > 0)) then
        
        call grid_bezier_square_polar(n_R, n_Z, n_radial, R_begin, R_end, Z_begin, Z_end, amin,    &
          fbnd, fpsi, mf, .true., node_list, element_list)
        
      else if ((n_R > 0) .and. (n_Z > 0) ) then
        
        call grid_bezier_square(n_R, n_Z, R_begin, R_end, Z_begin, Z_end, .true., node_list,       &
          element_list)
        
      else if ((n_radial > 0) .and. (n_pol > 0) ) then
        
        call grid_polar_bezier(R_geo, Z_geo, amin, 0.d0, fbnd, fpsi, mf, n_radial, n_pol,          &
          node_list, element_list)
        
      else
        write(*,*) ' FATAL : no valid combination of grid-sizes specified'
        call MPI_FINALIZE(IERR)
        stop
      end if 
      
      ! --- Determine boundary information from the grid
      call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list)
      
    end if
    
    ! --- Send boundary elements to different MPI threads
    call broadcast_boundary(my_id,bnd_elm_list,bnd_node_list)
    
    ! --- Fill the vacuum response matrices for freeboundary computations
    if ( freeboundary_equil ) call import_external_fields()
    if ( freeboundary ) call get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list,    &
      freeboundary_equil, use_starwall, resistive_wall)
    
    ! --- Plot the grid  
    if ( (my_id == 0) .and. (.not. bench_without_plot) ) then
      call plot_grid(node_list,element_list,bnd_elm_list,bnd_node_list,.true.,.false.)
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
      call equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint) 
      
      ! --- Determine a flux surface aligned grid
      if (n_flux > 1) then
        
        if (xpoint)  then
          
          call grid_xpoint(node_list, element_list, n_flux, n_open, n_private, n_leg, n_tht,       &
            SIG_open, SIG_closed, SIG_private, SIG_theta, SIG_leg_0, SIG_leg_1, dPSI_open,         &
            dPSI_private)
          
          call plot_grid(node_list,element_list,bnd_elm_list,bnd_node_list,.false.,.false.)
          
        else
          
          call grid_flux_surface(xpoint, node_list, element_list, surface_list, n_flux, n_tht,     &
            xr1, sig1, xr2, sig2)
          
          call plot_grid(node_list, element_list, bnd_elm_list, bnd_node_list, .true., .false.)  
          
          ! --- Refine elements (equilibrum)
          if (refinement) then
            n_to_be_refined=0
            call Refine_Elem_List(node_list, element_list, list_to_be_refined, n_to_be_refined)
            call Ref_Update_Index(element_list, node_list)
          end if
             
        end if
        
        ! --- Determine boundary information from the grid         
        call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list) 
        
        ! --- Compute the plasma equilibrium
        call equilibrium(my_id, node_list, element_list, bnd_node_list, bnd_elm_list, xpoint)
        
      end if
      
      ! --- Set initial conditions for time-evolution
      call initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint)
      
      ! --- Determine initial energies
      call energy(node_list,element_list,W_mag,W_kin)
      write(*,'(A,12e16.8)') ' initial energies : ', W_mag, W_kin
      
    end if
    
#ifdef USE_MUMPS
    ! --- Clean up this instance of mumps (used for equilibrium)
    mumps_par%JOB = -2
    if (my_id == 0) call DMUMPS(mumps_par)
#endif
    if (allocated(pastix_perm_vars))  deallocate(pastix_perm_vars)
    if (allocated(pastix_iperm_vars)) deallocate(pastix_iperm_vars)
    
  end if if_not_restart
  
  ! --- Broadcast grid information and input parameters to other MPI procs
  call broadcast_elements(my_id, element_list)                ! elements
  call broadcast_boundary(my_id, bnd_elm_list, bnd_node_list) ! boundary elements
  call broadcast_nodes(my_id, node_list)                      ! nodes
  call broadcast_phys(my_id)                                  ! physics parameters
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  !***********************************************************************
  !*                 end of initilisation/equilibrium                    *
  !***********************************************************************
  
  t_now     = t_start      ! t_now: current time in the simulation
  psi_bnd   = 0.d0
  
  if (nstep > 0) then

     call find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

     if (xpoint) then
        call find_xpoint(node_list, element_list, psi_xpoint, R_xpoint, Z_xpoint, i_elm_xpoint,    &
          s_xpoint, t_xpoint, ifail)
        psi_bnd = psi_xpoint
     else
        psi_bnd = 0.d0
     end if


     !*******************************************************
     !*      create groups /communicators                   *
     !* MPI_COMM_N      : group for each harmonic           *
     !* MPI_COMM_TRANS  : Transversal communicator          *
     !*   (ie : all first proc of MPI_COMM_N, all second,   *
     !*         all third...)                               *
     !* MPI_COMM_MASTER : group of masters of each harmonic *
     !*                   (i.e id=0 from each MPI_COMM_N)   *
     !*******************************************************
     murge_harmonic = 1
     if (gmres) then

        N_masters = (n_tor+1)/2
        if (MOD(n_cpu, N_masters) == 0) then
           M_cpu = n_cpu / (N_masters)
        else
           M_cpu = (n_cpu - MOD(n_cpu, N_masters))/N_masters +1
        end if

        allocate(i_tor(n_cpu))
        
        do i = 1, n_cpu 
           i_tor(i) =  MOD(i-1, M_cpu)+1
        end do
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,i_tor(my_id+1),my_id,MPI_COMM_TRANS,ierr)

        do i=1,n_cpu
           i_tor(i) = ((i-1) - MOD(i-1, M_cpu))/ M_cpu  + 1
        enddo
        murge_harmonic = i_tor(my_id+1)

        call MPI_COMM_SPLIT(MPI_COMM_WORLD,i_tor(my_id+1),my_id,MPI_COMM_N,ierr)
        
        do i=1,N_masters
           i_rank(i) = (i-1) * M_cpu
        enddo
 

        call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
        call MPI_GROUP_INCL(MPI_GROUP_WORLD,N_masters,i_rank,MPI_GROUP_MASTER,ierr)

        call MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_GROUP_MASTER,MPI_COMM_MASTER,ierr)

        call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)                 ! the id of each cpu
        call MPI_COMM_SIZE(MPI_COMM_N, n_cpu_n, ierr)                 ! the number of cpus
        call MPI_COMM_RANK(MPI_COMM_TRANS, my_id_trans, ierr)         ! the id of each cpu
        call MPI_COMM_SIZE(MPI_COMM_TRANS, n_cpu_trans, ierr)         ! the number of cpus
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
     !*            distribute nodes and elements over cpu's                 *
     !***********************************************************************
     if ( use_pastix .and. use_murge  .and. use_murge_element .and. gmres ) then
        index_size  = n_cpu_n
        id_elements = my_id_n
     else
        index_size  = n_cpu
        id_elements = my_id
     endif

     allocate(local_elms(element_list%n_elements))
     allocate(index_min(index_size),index_max(index_size))
     if ( .not. (use_pastix .and. use_murge .and. use_murge_element .and. gmres) ) then
        allocate(local_index_start(n_cpu),local_index_end(n_cpu))
     end if
     !
     ! Construct index_min, index_max and local_elems
     !
     call distribute_nodes_elements(id_elements,index_size,node_list,element_list,local_elms,      &
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
        !   TODO : Avoid doubles
        call murge_setgraph(gmres, mumps_par%n, local_elms, n_local_elms, &
          element_list, node_list, n_aa, my_id)
        call system_clock(count=t0)
        ! Build local_elms from loc2glob
        n_local_elms = 0

        DO i_elem = 1, element_list%n_elements

           element = element_list%element(i_elem)
           DO i=1,n_vertex_max

              inode1 = element%vertex(i)

              DO i_order = 1, n_order+1

                 index_node1 = node_list%node(inode1)%index(i_order)

                 call vertex_is_local(index_node1, is_local)
                 IF (is_local) THEN      
                    n_local_elms = n_local_elms + 1
                    EXIT
                 END IF
              END DO
              IF (is_local) THEN      
                 n_local_elms = n_local_elms + 1
                 EXIT
              END IF
           END DO
           IF (.not. is_local) then
              print *, my_id, ":", i_elem
           END IF
        END DO

        IF (ALLOCATED(local_elms)) DEALLOCATE(local_elms)
        ! Build local_elms from loc2glob
        ALLOCATE(local_elms(n_local_elms))

        n_local_elms = 0
        DO i_elem = 1, element_list%n_elements

           element = element_list%element(i_elem)
           L_I: DO i=1,n_vertex_max

              inode1         = element%vertex(i)

              DO i_order = 1, n_order+1

                 index_node1 = node_list%node(inode1)%index(i_order)

                 !index_large_i = n_tor * n_var * (index_node1 - 1)

                 !call vertex_is_local(index_node1*n_tor * n_var, loc2glob, local_n, is_local)
                 call vertex_is_local(index_node1, is_local)
                 IF (is_local) THEN      
                    n_local_elms = n_local_elms + 1
                    local_elms(n_local_elms) = i_elem
                    exit L_I
                 END IF
              END DO
           END DO L_I
        END DO
        call system_clock(count=t1)   
        nb_periods = t1-t0
        if (t1<t0) nb_periods = nb_periods + nb_periodes_max
        write(*,FMT_TIMING) my_id, ' system_clock elapsed time computing new local element list ',REAL(nb_periods)/nb_periodes_sec

        index_total = -1
        do inode=1, node_list%n_nodes
           index_total = max(index_total,maxval(node_list%node(inode)%index))
        enddo

        ndof_glob  = index_total * n_tor * n_var

        node_list%n_dof = ndof_glob
     END IF
     if (use_mumps) then
        if (.not. gmres) then
           call initialise_mumps(MPI_COMM_WORLD)    ! start MUMPS sparse matrix solver all cpus
        else
           call initialise_mumps(MPI_COMM_N)        ! start MUMPS sparse matrix solver on local groups
        endif
     endif

  endif
  
  ! --- Export a restart file before the first timestep
  if ( (my_id == 0) .and. (nstep > 0) .and. (.not. restart) ) then
    call export_restart(node_list,element_list,'jorek00000.rst')
  end if

  !***********************************************************************
  !***********************************************************************
  !*                          time stepping                              *
  !***********************************************************************
  !***********************************************************************

  if (nstep > 0) call update_deltas(my_id, node_list) ! create list of delta values in local_matrix module

  iter_gmres  = 999
  iter_big    = 200
  iter_precon = 22

  call r3_info_print (-2, -2, 'INITIALIZATION')    ! timing
  
  index_now = index_start  ! index_now: Index of current timestep

  do jstep = 1, 10 ! Go through the different values of the tstep_n and nstep_n arrays
  do istep = 1, nstep_n(jstep)

     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call flushc !flush the output stream

     index_now = index_now + 1
     
     tstep = tstep_n(jstep)
     
     if ( my_id == 0 ) then
       write(*,*) '******************************************************'
       write(*,'(A17,3i7,f14.5,A)') ' *   time step : ',jstep,istep,index_now,tstep,'  *'
       write(*,*) '******************************************************'
     end if

     call find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

     if (xpoint) then
       call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,ifail)
       psi_bnd = psi_xpoint
     else
       psi_bnd = 0.d0
     endif

     call cpu_time(t_matrix_0)

     ! Build the matrix 
     
     call system_clock(count=t0)
     if (gmres) then
        solve_only = .false.
        if ((gmres) .and. (istep .gt. 1)) then
           solve_only = .true.
           if (iter_gmres .gt. iter_precon) then                        ! redo preconditioner
              solve_only = .false.
           endif
        endif
     endif
     IF ( use_pastix .and. use_murge .and. use_murge_element ) THEN
        call construct_matrix_murge(my_id, node_list, element_list, local_elms, &
             n_local_ELms,  xpoint, &
             psi_axis, psi_bnd, Z_xpoint, gmres, i_tor, n_cpu, mpi_comm_n, &
             mpi_comm_trans, my_id_trans, n_cpu_trans, solve_only)        ! construct the matrix from elemental matrices
     ELSE
        call construct_matrix(my_id, local_elms, &
             n_local_ELms, index_min(my_id+1),index_max(my_id+1), &
             xpoint,psi_axis,psi_bnd,Z_xpoint)        ! construct the matrix from elemental matrices
     END IF

     call system_clock(count=t1)   
     nb_periods = t1-t0
     if (t1<t0) nb_periods = nb_periods + nb_periodes_max
     write(*,FMT_TIMING) my_id, ' system_clock elapsed time in construct_matrix ',REAL(nb_periods)/nb_periodes_sec
     call cpu_time(t_matrix_1)

     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     if (my_id .eq. 0) write(*,'(i3,A,f8.3)') my_id,' matrix  : ',t_matrix_1-t_matrix_0

     call cpu_time(t_solve_0)

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

        if (.not. solve_only) then

           call cpu_time(t_send_0)
           ! with murge elementary assembly harmonic distribution is already done.
           IF ( .not. ( use_pastix .and. use_murge .and. use_murge_element ) ) THEN
              call distribute_harmonics(my_id,my_id_n,n_cpu)
           ELSE
              call distribute_vector(my_id,rhs_glob,mumps_par%rhs)              
           END IF
           call cpu_time(t_send_1)

           if (my_id .eq. 0) write(*,'(i3,A,f8.3)') my_id,' distribute  : ',t_send_1-t_send_0

           call MPI_Barrier(MPI_COMM_WORLD,ierr)
        endif

        if (.not. gmres) call update_rhs_n(my_id,my_id_n, i_tor, MPI_COMM_MASTER)      ! correct the RHS with the previous solution (deltas)
        if (use_murge .and. use_murge_element) then
           call solve_murge_all(n_cpu,my_id,index_min(my_id_n+1),index_max(my_id_n+1), i_tor, gmres, my_id_n, mpi_comm_n, mpi_comm_master)
        else
           call solve_matrix_n(my_id,i_tor,MPI_COMM_N,MPI_COMM_MASTER,solve_only)    ! factorise preconditioning matrices
        end if
     endif

     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     call cpu_time(t_solve_1)

     if (gmres) call gmres_driver(my_id,my_id_n,i_tor, n_tor,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)

     call cpu_time(t_solve_2)

     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     if (my_id .eq. 0) write(*,'(i3,A,f8.3)') my_id,' solve : ',t_solve_1-t_solve_0
     if (my_id .eq. 0) write(*,'(i3,A,f8.3)') my_id,' gmres : ',t_solve_2-t_solve_1

     if ( (gmres .and. (iter_gmres .lt. iter_big)) .or. (.not.gmres) ) then

        call update_values(my_id,element_list,node_list,deltas)         ! add solution to node values
        call update_deltas(my_id,node_list)

        t_now = t_now + tstep

     else
        write(*,*) ' TIME STEP SKIPPED !', iter_gmres
	exit
     endif

     !-------------------------------------------------------- adapt time step (in progress...)
     mindelta = minval(deltas); maxdelta = maxval(deltas);
     if (my_id .eq. 0) write(*,'(A,2e16.8,2i12)') ' min/max deltas : ',mindelta,maxdelta,minloc(deltas),maxloc(deltas)

     if (gmres .and. adaptive_time) then        ! experimental
        if (iter_gmres .ge. iter_big) then
           tstep = tstep /2.d0
           write(*,*) my_id,' REDUCTION TIMESTEP : ',tstep
        elseif (max(abs(mindelta),abs(maxdelta)) .gt. 0.05) then
           !      tstep = tstep /2.d0
           !      iter_gmres = 99999
           !      write(*,*) my_id,' REDUCTION TIMESTEP : ',tstep
        elseif (max(abs(mindelta),abs(maxdelta)) .lt. 0.001) then
           !      tstep = tstep * 2.d0
           !      iter_gmres = 99999
           !      write(*,*) my_id,' INCREASE TIMESTEP : ',tstep
        endif
     endif

     !--------------------------------------------------------- energies
     if ( (my_id == 0) .and. (.not. bench_without_plot) ) then
        call energy(node_list,element_list,W_mag,W_kin)

        xtime(index_now) = t_now
        energies(1:n_tor,1,index_now) = W_mag(1:n_tor)
        energies(1:n_tor,2,index_now) = W_kin(1:n_tor)

        Growth_mag  = 0.d0; Growth_kin  = 0.d0; Growth_mag0 = 0.d0; Growth_kin0 = 0.d0

        if (index_start+istep .gt. 1) then
           Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,index_now)/energies(n_tor,1,index_now-1)))/ tstep
           Growth_kin  = 0.5d0*log(abs(energies(n_tor,2,index_now)/energies(n_tor,2,index_now-1)))/ tstep
           Growth_mag0 = 0.5d0*log(abs(energies(1,1,index_now)/energies(1,1,index_now-1)))/ tstep
           Growth_kin0 = 0.5d0*log(abs(energies(1,2,index_now)/energies(1,2,index_now-1)))/ tstep
        endif

        write(*,'(i5,12e14.6)') istep,t_now,W_mag(1),W_kin(1),W_mag(n_tor),W_kin(n_tor),Growth_kin0,Growth_kin

        ! --- Output energies and growth_rates to text files during the code run
        call write_live_data(index_now)
        
     endif

     !---------------------------------------------------------timing
     if ( istep == 1 ) then
        call r3_info_print (-3, -2, 'ITERATION    1')
     else
        call r3_info_print (istep, -2, 'ITERATION')
     endif
     
     ! --- Write a restart file every nout timesteps
     if ( (my_id == 0) .and. (mod(index_now,nout) == 0) ) then
       write(fileout,'(A5,i5.5,A4)') 'jorek',index_now,'.rst'
       call export_restart(node_list,element_list,fileout)
     endif

  enddo                                              ! end of time stepping
  enddo
  
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


           IF ( use_murge_element ) THEN
              IF (ALLOCATED(murge_glob2loc)) DEALLOCATE(murge_glob2loc)
              IF (ALLOCATED(murge_loc2glob)) DEALLOCATE(murge_loc2glob)
           END IF
           CALL MURGE_Clean(murge_id, ierr)
        else
           pastix_iparm(2)     = 7                       ! Clean-up
           pastix_iparm(3)     = 7

           if (.not. gmres) then

              call pastix_fortran(pastix_data,MPI_COMM_WORLD,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                   pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

           elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0))  ) then

              call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                   pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
           endif


        end if
     endif
  endif
  
  ! --- Close open files
  if ( (my_id == 0) .and. (.not. bench_without_plot) ) call finalize_live_data()

  !***********************************************************************
  !*                          plots etc.                                 *
  !***********************************************************************
  
  if (my_id .eq. 0)  then

     call export_restart(node_list,element_list,'jorek_restart.rst')

     if (.not. bench_without_plot) then
        PI = 2.d0*asin(1.d0)

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
           call lincol(1)
           call lplot(1,1,2,xtime(1:index_now),energies(i,1,1:index_now),-index_now,1,'Magnetic Energie',16,'time',4,'Emag',4)
           call lincol(2)
           call lplot(1,1,2,xtime(1:index_now),energies(i,2,1:index_now),-index_now,1,'Kinetic Energie',15,'time',4,'Ekin',4)
        enddo
        call lincol(3)
        call lplot(1,1,2,xtime(1:index_now),energies(1,2,1:index_now),-index_now,1,'Kinetic Energie',15,'time',4,'Ekin',4)
        call lincol(0)
     endif

!---------------------------------------------- plot equilibrium current profile (to be removed)
     call find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis, ifail)

     nplot = 501
     allocate(xp(nplot),yp1(nplot),yp2(nplot),yp3(nplot))
     iplot = 0

     if (xpoint) then
        call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,ifail)
     else
        psi_xpoint = psi_bnd
     endif

     Rp_start = R_axis - amin*2.d0
     Rp_end   = R_axis + amin*2.d0

     Zp = Z_axis

     do i=1,nplot

        Rp =  Rp_start + float(i-1)/float(nplot-1) * (Rp_end - Rp_start)

        call find_RZ(node_list,element_list,Rp,Zp,R_out,Z_out,i_elm,s_out,t_out,ifail)

        if (ifail .eq. 0) then

           call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi,P_s,P_t,P_st,P_ss,P_tt)

           call density(    xpoint, Zp, Z_xpoint, psi,psi_axis,psi_xpoint,           &
                zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz)
           call temperature(xpoint, Zp, Z_xpoint, psi,psi_axis,psi_xpoint,           &
                zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)
           call FFprime(    xpoint, Zp, Z_xpoint, psi,psi_axis,psi_xpoint,           &
                zFFprime,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2,dFFprime_dpsi_dz)

           zjz   = (zFFprime - Rp*Rp * (zn * dT_dpsi + dn_dpsi * zT)) / Rp

           iplot = iplot + 1

           xp(iplot)  = Rp
           yp1(iplot) = zFFprime / Rp
           yp2(iplot) = zjz
           yp3(iplot) = - Rp*Rp * (zn * dT_dpsi + dn_dpsi * zT) / Rp

           !      write(*,'(A,8e16.8)') ' profiles : ',xp(iplot),psi,psi_axis,psi_xpoint,yp2(iplot),yp1(iplot),yp3(iplot)

        endif

     enddo

     call lplot6(1,1,xp,yp2,iplot,' ')
     call lincol(1)
     call lplot6(1,1,xp,yp1,-iplot,' ')
     call lincol(2)
     call lplot6(1,1,xp,yp3,-iplot,' ')
     call lincol(0)
     call finplt                                          ! close plot file

!  call export_POV(node_list,element_list,3,1)          ! export to POVray native bezier patch format

     call export_helena(node_list,element_list,bnd_elm_list)

     if (allocated(energies))  deallocate(energies)
     if (allocated(xtime))     deallocate(xtime)
  endif
  
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
!! <ul>
!! <li>
!!   <a href="http://www.ipp.mpg.de/~mhoelzl/jorek/dirs.html">directories</a>,
!! </li><li>
!!   <a href="http://www.ipp.mpg.de/~mhoelzl/jorek/files.html">source files</a>,
!! </li><li>
!!   <a href="http://www.ipp.mpg.de/~mhoelzl/jorek/namespaces.html">Fortran modules</a>,
!! </li><li>
!!   <a href="http://www.ipp.mpg.de/~mhoelzl/jorek/globals_func.html">subroutines and functions</a>, and
!! </li><li>
!!   routine parameters.
!! </li>
!! </ul>
!!
!! Background information on the implemented equations,
!! the numerical methods, compiling and running the code, ..., is covered in a LaTeX
!! documentation that is available in the docu/ folder of the JOREK repository as LaTeX
!! source and PDF.
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
