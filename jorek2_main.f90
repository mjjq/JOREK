!***************************************************************************************************
!*                                         JOREK 2.0                                               *
!***************************************************************************************************
!*   program solves the (reduced) MHD equations in 3D toroidal geometry                            *
!*                                                                                                 *
!*   solvers implemented:                                                                          *
!*     - MUMPS                                                                                     *
!*     - PastiX                                                                                    *
!*     - GMRES (+MUMPS or PastiX preconditioner)                                                   *
!*                                                                                                 *
!*   required libraries :                                                                          *
!*     - MPI                                                                                       *
!*     - MUMPS                                                                                     *
!*     - PastiX                                                                                    *
!*     - SCOTCH (metis)                                                                            *
!*     - FFTW                                                                                      *
!*     - SCALAPACK (BLACS)                                                                         *
!*     - LAPACK, BLAS                                                                              *
!*     - PPPLIB                                                                                    *
!*                                                                                                 *
!*  Author : Guido Huysmans (Euratom / CEA Association)                                            *
!*  Date   : 18-7-2008                                                                             *
!***************************************************************************************************
program JOREK2

  use mumps_module
  use pastix_module
  use murge_module
  use data_structure
  use phys_module
  use global_distributed_matrix
  use nodes_elements
  use vacuum_response_module

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

  end interface
  
  type (type_surface_list) :: surface_list
  logical                  :: grid_changed
  real*8                   :: W_mag(n_tor), W_kin(n_tor), growth_mag, growth_kin, growth_mag0, growth_kin0
  real*8                   :: t_matrix_0, t_matrix_1, PI
  real*8                   :: t_send_0, t_send_1, t_solve_0, t_solve_1, t_solve_2
  real*8                   :: psi_bnd, psi_axis, R_axis, Z_axis, s_axis, t_axis
  real*8                   :: psi_xpoint, R_xpoint, Z_xpoint, s_xpoint, t_xpoint, mindelta, maxdelta
  integer                  :: my_id, my_id_n, my_id_master, my_id_eq, n_cpu_eq
  integer                  :: istep,ierr,i,j,k,itor,inode, i_elm_axis, i_elm_xpoint
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
  integer                  :: nplot, iplot, i_elm, ifail, ivar, iter_big, iter_precon
  logical                  :: is_local
  integer                  :: i_elem, inode1, i_order, k_order, index_node1, index_node2, knode
  type (type_element)      :: element
  integer                  :: index_size, id_elements

  integer                  ::  list_to_be_refined(n_ref_list), n_to_be_refined    
  
  logical                  :: bench_without_plot
  integer                  :: t0,t1,nb_periodes_max,nb_periodes_sec, nb_periods
  character(len=20), parameter :: FMT_TIMING = "(I2,A70,F7.2)"

  !***********************************************************************
  !*                  intialisation                                      *
  !***********************************************************************

  ! --- Initialise MPI
  !call MPI_INIT(IERR)                                     ! initialise MPI
  required=MPI_THREAD_MULTIPLE
  call MPI_Init_thread(required,provided,StatInfo)         ! initialise threaded MPI (openMPI)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)           ! the id of each cpu
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)      ! the number of cpus
  my_id = rank
  n_cpu = comm_size
 
  if (my_id .eq. 0) then
    write(*,*) '****************************************'
    write(*,*) '*   3D Reduced MHD : JOREK_2.0         *'
    write(*,*) '****************************************'
    write(*,*) ' number of cpus : ',n_cpu
  endif

  ! --- Initialise timing
  call system_clock(count_rate=nb_periodes_sec,count_max=nb_periodes_max) ! elapsed time
  call r3_info_init ()                                     ! timing

  ! --- Select solver
  gmres              = .true.                              ! .true. for gmres, .false. for direct
  use_mumps          = .false.
  use_pastix         = (.not. use_mumps)
  use_murge          = .false.
  use_murge_element  = .false.
  pastix_initialised = .false.
  pastix_analysed    = .false.
  murge_initialised  = .false.
  pastix_smp_only    = .false.         ! implies that each MPI group resides within one node!    
  
  refinement    = .false.              ! enable mesh refinement

  adaptive_time = .false.              ! requires no_mpi for Pastix library

  if (n_tor .eq. 1) gmres  = .false.

  ! --- Flag from HSLT
  bench_without_plot              = .false.    ! .true. for benchmark (mesuring elapsed time without plot phases) 
  use_matrix_whitout_zeros_pastix = .false.    ! .true. to remove nonzeros in the preconditioning matrix with MUMPS
  use_matrix_whitout_zeros_mumps  = .false.    ! .true. to remove nonzeros in the preconditioning matrix with PaStiX

  ! --- Preset input parameters to reasonable defaults; then read the input file.
  call initialise_parameters(my_id)

  ! --- Fill the arrays mode (toroidal mode number n) and mode_type (cos or sin).
  do itor=1, n_tor
    mode(itor)      = + int(itor / 2) * n_period
    if ( (itor==1) .or. (mod(itor,2)==0) ) then
      mode_type(itor) = 'cos'
    else
      mode_type(itor) = 'sin'
    end if
    !write(*,*) ' toroidal mode numbers : ',itor,mode(itor)
  enddo
  
  ! --- Write out all parameters defined in mod_parameters and by the input file.
  call log_parameters(my_id)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Some checks not to waste any cpu time
  if (required.ne.provided) then
    write(*,*) 'FATAL : MPI_THREAD_MULTIPLE (provided is smaller than required)',my_id,required,provided
    call MPI_FINALIZE(IERR)
    stop
  endif

  if ( (.not. use_mumps) .and. (.not. use_pastix) ) then
    write(*,*) ' FATAL : specify a valid solver'
    call MPI_FINALIZE(IERR)                                ! clean up MPI
    stop
  end if

  if ((n_plane .lt. n_tor + 1) .and. (n_tor .gt. 1)) then
    write(*,*) ' FATAL : n_plane too small ',n_plane,n_tor
    call MPI_FINALIZE(IERR)                                ! clean up MPI
    stop
  end if

  if ((gmres) .and. (nstep .gt. 0)) then
    if (n_cpu .lt. (n_tor-1)/2+1) then
      write(*,'(A,i4,A,i4,A)') ' FATAL : need at least',(n_tor-1)/2+1,' cpus for ',(n_tor-1)/2+1,' harmonics'
  !    call MPI_FINALIZE(IERR)                             ! clean up MPI
  !    stop
    end if
    if (mod(n_cpu,(n_tor-1)/2+1) .ne. 0) then
      write(*,'(A,i4,A,i4,A)') ' FATAL : need a multiple of ',(n_tor-1)/2+1,' cpus for ',(n_tor-1)/2+1,' harmonics'
      call MPI_FINALIZE(IERR)                              ! clean up MPI
      stop
    end if
  end if
  
  ! --- Initialise ppplib plotting library
  if (my_id .eq. 0)  call begplt('jorek2.ps')

  ! --- Define the basis functions at the Gaussian points
  call initialise_basis()

  ! --- Read the restart file to continue a previous JOREK run (if restart is .true.)
  if ( restart .and. (my_id .eq. 0) ) then

    call import_restart(node_list,element_list)    ! read restart file
    tstep = tstep_in

    ! --- Optional: redo fluxsurface grid (DOES NOT WORK CURRENTLY)
    if (regrid) then
      if (xpoint)  then
        call grid_xpoint(node_list,element_list,n_flux,n_open,n_private,n_leg,n_tht,            &
          SIG_open,SIG_closed,SIG_private,SIG_theta,SIG_leg_0,SIG_leg_1,dPSI_open,dPSI_private)
      else
        call grid_flux_surface(xpoint,node_list,element_list,surface_list,n_flux,n_tht,xr1,sig1,xr2,sig2)
      end if
    end if
    
  end if

  !***********************************************************************
  !*                  define grid / equilibrium                          *
  !***********************************************************************

  if (.not. restart) then

    element_list%n_elements      = 0
    bnd_elm_list%n_bnd_elements  = 0
    node_list%n_nodes            = 0

    if (my_id .eq. 0) then

      call define_boundary

      if ((n_R .gt. 0) .and. (n_Z .gt. 0) .and. (n_radial .gt.0)) then

        call grid_bezier_square_polar(n_R,n_Z,n_radial,R_begin,R_end,Z_begin,Z_end,amin,fbnd,fpsi,mf,.true.,node_list,element_list)

      elseif ((n_R .gt. 0) .and. (n_Z .gt. 0) ) then

        call grid_bezier_square(n_R,n_Z,R_begin,R_end,Z_begin,Z_end,.true.,node_list,element_list)

      elseif ((n_radial .gt. 0) .and. (n_pol .gt. 0) ) then

        call grid_polar_bezier(R_geo,Z_geo,amin,0.d0,fbnd,fpsi,mf,n_radial,n_pol,node_list,element_list)

      else
        write(*,*) ' FATAL : no valid combination of grid-sizes specified'
      endif 

      call boundary_from_grid(node_list,element_list,bnd_node_list,bnd_elm_list)                                                ! Determine boundary information from the grid
        
      if ( freeboundary ) then         
        call get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list,.true.)  ! Fill the vacuum response matrix/matrices
      endif
        
      if (.not. bench_without_plot) call plot_grid(node_list,element_list,bnd_elm_list,bnd_node_list,.true.,.false.)    ! plot the grid
            
    endif
    
#ifdef USE_MUMPS
    call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
    call MPI_GROUP_INCL(MPI_GROUP_WORLD,1,0,MPI_GROUP_MUMPS_EQUIL,ierr)
    call MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_GROUP_MUMPS_EQUIL,MPI_COMM_MUMPS_EQUIL,ierr)

    if (my_id .eq. 0) call initialise_mumps(MPI_COMM_MUMPS_EQUIL)                   ! start MUMPS sparse matrix solver
#endif

    if (my_id .eq. 0) then

      call equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint) 

      if (n_flux .gt. 1) then                                                          ! flux surface grid
        
        if (xpoint)  then

          call grid_xpoint(node_list,element_list,n_flux,n_open,n_private,n_leg,n_tht,   &
                           SIG_open,SIG_closed,SIG_private,SIG_theta,SIG_leg_0,SIG_leg_1,dPSI_open,dPSI_private)

          call plot_grid(node_list,element_list,bnd_elm_list,bnd_node_list,.false.,.false.)
          
        else

          call grid_flux_surface(xpoint,node_list,element_list,surface_list,n_flux,n_tht,xr1,sig1,xr2,sig2)
          
          call plot_grid(node_list,element_list,bnd_elm_list,bnd_node_list,.true.,.false.)  
          
          !****************************************************************************
          !                        Refining elements (equilibrum)                     *
          !****************************************************************************

          if (refinement) then
  
            n_to_be_refined=0
            call Refine_Elem_List(node_list, element_list,list_to_be_refined,n_to_be_refined)
            call Ref_Update_Index( element_list,node_list)
          
          endif
             
        endif
                                               ! Determine boundary information from the grid         
        call boundary_from_grid(node_list,element_list,bnd_node_list,bnd_elm_list) 

        call equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint)

      endif
  
      call initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint) ! initial conditions
 
      call energy(node_list,element_list,W_mag,W_kin)
      write(*,'(A,12e16.8)') ' initial energies : ', W_mag, W_kin
  
    endif

#ifdef USE_MUMPS
    mumps_par%JOB = -2                                       ! clean up this instance of mumps
    if (my_id .eq. 0) call DMUMPS(mumps_par)
#endif
    if (allocated(pastix_perm_vars))  deallocate(pastix_perm_vars)
    if (allocated(pastix_iperm_vars)) deallocate(pastix_iperm_vars)

  endif

  call broadcast_elements(my_id,element_list)                ! sending all elements
  call broadcast_boundary(my_id,bnd_elm_list,bnd_node_list)  ! sending boundary elements
  call broadcast_nodes(my_id,node_list)                      ! sending all nodes
  call broadcast_phys(my_id)                                 ! sending the physics parameters

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  !************************************************************************
  !*                        vacuum initialisation                         *
  !************************************************************************
    
  ! --- Fill the vacuum response matrix/matrices
!  if ( freeboundary ) call get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list)

  !***********************************************************************
  !*                 end of initilisation/equilibrium                    *
  !***********************************************************************

  t_now         = t_start
  index_now     = index_start

  if (nstep .gt. 0) then

     grid_changed  = .true.
     call find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

     if (xpoint) then
        call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,ifail)
        psi_bnd = psi_xpoint
     else
        psi_bnd = 0.d0
     endif


     !*******************************************************
     !*      create groups /communicators                   *
     !* MPI_COMM_N      : group for each harmonic           *
     !* MPI_COMM_TRANS  : Transversal communicator          *
     !*   (ie : all first proc of MPI_COMM_N, all second,   *
     !*         all third...)                               *
     !* MPI_COMM_MASTER : group of masters of each harmonic *
     !*                   (i.e id=0 from each MPI_COMM_N)   *
     !*******************************************************

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
     if ( use_pastix .and. use_murge  .and. use_murge_element .and. gmres ) then
     else
        allocate(local_index_start(n_cpu),local_index_end(n_cpu))
     endif
     !
     ! Construct index_min, index_max and local_elems
     !
     call distribute_nodes_elements(id_elements,index_size,node_list,element_list,local_elms,n_local_elms, &
          ndof_glob,index_min,index_max)

     node_list%n_dof = ndof_glob
     if ( use_pastix .and. use_murge  .and. use_murge_element .and. gmres ) then
     else
        local_index_start = index_min
        local_index_end   = index_max
     end IF
     ! Build ijA_index, ijA_size and irn_jcn
     
     call global_matrix_structure(my_id_n,node_List,element_list,bnd_elm_list, freeboundary,&
          local_elms,n_local_elms,index_min(id_elements+1),index_max(id_elements+1))

     if ( use_pastix .and. use_murge .and. use_murge_element ) then

        write (*,*) "--- Murge initilisation ---"

        !
        ! Murge initialisation and 
        ! graph definition edge by edge
        !
        IF (use_murge_element .and. .NOT. murge_initialised) THEN
           !
           ! Init murge
           !
           if (gmres) then
              CALL MURGE_Initialize(N_masters, ierr)
              murge_id = i_tor(my_id+1)-1
              CALL MURGE_SetCommunicator(murge_id, MPI_COMM_N, ierr)
           else
              CALL MURGE_Initialize(1, ierr)
              murge_id = 0;
           end if

           CALL MURGE_GetSolver(murge_solver, ierr)
           IF (murge_solver == MURGE_SOLVER_PASTIX) THEN
              CALL MURGE_SetDefaultOptions(murge_id, 0, ierr)
              write (*,*) "--- Murge_setdefaultoptions ---"
              CALL MURGE_SetOptionINT(murge_id, IPARM_VERBOSE,             API_VERBOSE_YES, ierr)
              CALL MURGE_SetOptionINT(murge_id, IPARM_MATRIX_VERIFICATION, API_YES,         ierr)
              ! refinement : max number of iterations
              CALL MURGE_SetOptionINT(murge_id, IPARM_ITERMAX,             murge_iter,      ierr) 
              ! degrees of freedom per node (not correct)
              if (gmres) then
                 if ( i_tor(my_id+1) == 1 ) then
                    murge_ndof = n_var
                 else
                    murge_ndof = 2*n_var
                 end if
              CALL MURGE_SetOptionINT(murge_id, IPARM_DOF_COST,     2*n_var,      ierr) 
              else
                 murge_ndof = n_tor*n_var
              end if
              CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_DOF,          murge_ndof,      ierr) 

              CALL MURGE_SetOptionINT(murge_id, IPARM_THREAD_NBR,          murge_nthrd,     ierr)
              CALL MURGE_SetOptionINT(murge_id, IPARM_LEVEL_OF_FILL,       murge_iluk,      ierr)
              CALL MURGE_SetOptionINT(murge_id, IPARM_INCOMPLETE,          murge_ricar,     ierr)
              CALL MURGE_SetOptionINT(murge_id, IPARM_AMALGAMATION_LEVEL,  murge_amalg,     ierr)
              CALL MURGE_SetOptionINT(murge_id, IPARM_MATRIX_VERIFICATION, API_YES,         ierr)

              CALL MURGE_SetOptionREAL(murge_id, DPARM_EPSILON_MAGN_CTRL,  murge_pivot,     ierr)

           ENDIF


           CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_SYM,         murge_sym,   ierr)
           CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_BASEVAL,     1,   ierr)

           CALL MURGE_SetOptionREAL(murge_id, MURGE_RPARAM_EPSILON_ERROR, murge_epsilon, ierr)

           murge_initialised = .TRUE.

        ENDIF

        !
        ! Build the graph
        !
        ! TODO: Avoid doubles
        !
        !if (my_id == 0) then
        CALL MURGE_GRAPHBEGIN(murge_id, mumps_par%n, &
             n_local_elms*(n_order+1)*(n_order+1)*n_vertex_max*n_vertex_max, ierr)
        IF (ierr /= MURGE_SUCCESS) THEN
           write (*,*) "ERROR in MURGE_GRAPHBEGIN"
           STOP
        END IF

        call system_clock(count=t0)

        DO i_elem = 1, n_local_elms

           element = element_list%element(local_elms(i_elem))
           DO i=1,n_vertex_max

              inode1         = element%vertex(i)

              DO i_order = 1, n_order+1

                 index_node1 = node_list%node(inode1)%index(i_order)



                 ! Build nodes Matrices
                 DO k=1,n_vertex_max

                    knode         = element%vertex(k)

                    DO k_order = 1, n_order+1

                       index_node2 = node_list%node(knode)%index(k_order)


                       CALL MURGE_GRAPHEDGE(murge_id,  &
                            index_node1,         &
                            index_node2,         &
                            ierr)
                       IF (ierr /= MURGE_SUCCESS) THEN
                          write (*,*) "N", mumps_par%n, &
                               "I", index_node1, &
                               "J", index_node2
                          STOP
                       END IF


                    END DO
                 END DO
              END DO
           END DO
        END DO

        call system_clock(count=t1)
        nb_periods = t1-t0
        if (t1<t0) nb_periods = nb_periods + nb_periodes_max
        write(*,FMT_TIMING) my_id, ' system_clock elapsed time entering graph ',REAL(nb_periods)/nb_periodes_sec

        call system_clock(count=t0)
        CALL MURGE_GRAPHEND(murge_id, ierr)
        IF (ierr /= MURGE_SUCCESS) THEN
           write (*,*) "ERROR in MURGE_GRAPHEND"
           STOP
        END IF
        call system_clock(count=t1)
        nb_periods = t1-t0
        if (t1<t0) nb_periods = nb_periods + nb_periodes_max
        write(*,FMT_TIMING) my_id, ' system_clock elapsed time in MURGE_GRAPHEND ',REAL(nb_periods)/nb_periodes_sec


        call system_clock(count=t0)
        CALL MURGE_GETLOCALNODENBR(murge_id, murge_local_n, ierr)
        IF (ierr /= MURGE_SUCCESS) THEN
           write (*,*) "ERROR in MURGE_GETLOCALNODENBR"
           STOP
        END IF
        call system_clock(count=t1)
        nb_periods = t1-t0
        if (t1<t0) nb_periods = nb_periods + nb_periodes_max
        write(*,FMT_TIMING) my_id, ' system_clock elapsed time in MURGE_GETLOCALNODENBR ',REAL(nb_periods)/nb_periodes_sec

        
        ALLOCATE(murge_loc2glob(murge_local_n))
        murge_global_n = mumps_par%n
        write (*,*) "Local number of nodes", murge_local_n, "global",  mumps_par%n
        call system_clock(count=t0)
        CALL MURGE_GETLOCALNODELIST(murge_id, murge_loc2glob, ierr)
        if (allocated(murge_glob2loc)) deallocate(murge_glob2loc)

        call system_clock(count=t1)
        nb_periods = t1-t0
        if (t1<t0) nb_periods = nb_periods + nb_periodes_max
        write(*,FMT_TIMING) my_id, ' system_clock elapsed time in MURGE_GETLOCALNODELIST ',REAL(nb_periods)/nb_periodes_sec
        IF (ierr /= MURGE_SUCCESS) THEN
           write (*,*) "ERROR in MURGE_GETLOCALNODELIST"
           STOP
        END IF

        call system_clock(count=t0)
        ! Build local_elms from loc2glob
        n_local_elms = 0

        DO i_elem = 1, element_list%n_elements

           element = element_list%element(i_elem)
           DO i=1,n_vertex_max

              inode1         = element%vertex(i)

              DO i_order = 1, n_order+1

                 index_node1 = node_list%node(inode1)%index(i_order)

                 call vertex_is_local(index_node1, is_local)
                 IF (is_local) THEN      
                    n_local_elms = n_local_elms + 1
                    GOTO 10
                 END IF
              END DO
           END DO
10         continue
      END DO 
        IF (ALLOCATED(local_elms)) DEALLOCATE(local_elms)
        ! Build local_elms from loc2glob
        ALLOCATE(local_elms(n_local_elms))

        n_local_elms = 0
        DO i_elem = 1, element_list%n_elements

           element = element_list%element(i_elem)
           DO i=1,n_vertex_max

              inode1         = element%vertex(i)

              DO i_order = 1, n_order+1

                 index_node1 = node_list%node(inode1)%index(i_order)

                 !index_large_i = n_tor * n_var * (index_node1 - 1)

                 !call vertex_is_local(index_node1*n_tor * n_var, loc2glob, local_n, is_local)
                 call vertex_is_local(index_node1, is_local)
                 IF (is_local) THEN      
                    n_local_elms = n_local_elms + 1
                    local_elms(n_local_elms) = i_elem
                    GOTO 20
                 END IF
              END DO
           END DO
20         continue
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


  !***********************************************************************
  !***********************************************************************
  !*                          time stepping                              *
  !***********************************************************************
  !***********************************************************************


  if (nstep .gt. 0) call update_deltas(my_id,node_list)        ! create list of delta values in local_matrix module

  iter_gmres  = 999
  iter_big    = 200
  iter_precon = 22

  call r3_info_print (-2, -2, 'INITIALIZATION')    ! timing

  do istep = 1, nstep

     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     index_now = index_start + istep

     if (my_id .eq. 0)  write(*,*) '********************************************'
     if (my_id .eq. 0)  write(*,'(A17,2i7,f8.2,A)') ' *   time step : ',istep,index_now,tstep,'     *'
     if (my_id .eq. 0)  write(*,*) '********************************************'

     !  call find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

     !  if (xpoint) then
     !    call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,ifail)
     !    psi_bnd = psi_xpoint
     !  else
     !    psi_bnd = 0.d0
     !  endif

     call cpu_time(t_matrix_0)

     ! Build the matrix 
     
     call system_clock(count=t0)
     IF ( use_pastix .and. use_murge .and. use_murge_element ) THEN
        call construct_matrix_murge(my_id, node_list, element_list, local_elms, &
             n_local_ELms,  xpoint, &
             psi_axis, psi_bnd, Z_xpoint, gmres, i_tor, n_cpu, mpi_comm_n, &
             mpi_comm_trans, my_id_trans, n_cpu_trans)        ! construct the matrix from elemental matrices
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

        solve_only = .false.
        if ((gmres) .and. (istep .gt. 1)) then
           solve_only = .true.
           if (iter_gmres .gt. iter_precon) then                        ! redo preconditioner
              solve_only = .false.
           endif
        endif

        if (.not. solve_only) then

           call cpu_time(t_send_0)
           IF ( .not. ( use_pastix .and. use_murge .and. use_murge_element ) ) THEN
              call distribute_harmonics(my_id,my_id_n,n_cpu)
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

     if (gmres) then

        call gmres_driver(my_id,my_id_n,i_tor, n_tor,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)     ! gmres solution

     endif

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
     endif

     if ( freeboundary .and. (.not. resistive_wall) ) call boundary_check()

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
     if (my_id .eq. 0 .and. .not. bench_without_plot)  then
        call energy(node_list,element_list,W_mag,W_kin)

        xtime(index_start+istep) = t_now
        energies(1:n_tor,1,index_start+istep) = W_mag(1:n_tor)
        energies(1:n_tor,2,index_start+istep) = W_kin(1:n_tor)

        Growth_mag  = 0.d0; Growth_kin  = 0.d0; Growth_mag0 = 0.d0; Growth_kin0 = 0.d0

        if (index_start+istep .gt. 1) then
           Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,index_start+istep)/energies(n_tor,1,index_start+istep-1)))/ tstep
           Growth_kin  = 0.5d0*log(abs(energies(n_tor,2,index_start+istep)/energies(n_tor,2,index_start+istep-1)))/ tstep
           Growth_mag0 = 0.5d0*log(abs(energies(1,1,index_start+istep)/energies(1,1,index_start+istep-1)))/ tstep
           Growth_kin0 = 0.5d0*log(abs(energies(1,2,index_start+istep)/energies(1,2,index_start+istep-1)))/ tstep
        endif

        write(*,'(i5,12e14.6)') istep,t_now,W_mag(1),W_kin(1),W_mag(n_tor),W_kin(n_tor),Growth_kin0,Growth_kin

        !### FOR TESTING: OUTPUT SOME ADDITIONAL INFORMATION
        !open(42, file='times_and_steps.dat', status='REPLACE',action='WRITE')
        !do j = 1, index_now
        !  write(42,*) j, xtime(j)
        !end do
        !close(42)
        !
        !open(42, file='energies.dat', status='REPLACE',action='WRITE')
        !do i = 1, n_tor
        !  do j = 1, index_now
        !    write(42,*) xtime(j), energies(i,1:2,j)
        !  end do
        !  write(42,*)
        !  write(42,*)
        !end do
        !close(42)
        !
        !open(42, file='growth_rates.dat', status='REPLACE',action='WRITE')
        !do i = 1, n_tor
        !  do j = 2, index_now
        !   write(42,*) ( xtime(j) + xtime(j-1) ) / 2., 0.5 * ( LOG(energies(i,1:2,j)) - LOG(energies(i,1:2,j-1)) ) / ( xtime(j) - xtime(j-1) )
        !  end do
        !  write(42,*)
        !  write(42,*)
        !end do
        !close(42)
        !### END FOR TESTING

        !  print*,"enrj",abs(energies(1,2,index_start+istep-1))
     endif

     !---------------------------------------------------------timing
     if (istep.eq.1) then
        call r3_info_print (-3, -2, 'ITERATION    1')
     else
        call r3_info_print (istep, -2, 'ITERATION')
     endif

     if (my_id .eq. 0) then
        if (mod(index_now,nout) .eq. 0) then
           write(fileout,'(A5,i5.5,A4)') 'jorek',index_now,'.rst'
           call export_restart(node_list,element_list,fileout)
        endif
     endif

  enddo                                              ! end of time stepping

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

endif
call MPI_FINALIZE(IERR)                                ! clean up MPI
  stop
if (my_id.eq.0) then
     if (allocated(energies))  deallocate(energies)
     if (allocated(xtime))     deallocate(xtime)

  endif
  call r3_info_summary ()                                ! timing
  call MPI_FINALIZE(IERR)                                ! clean up MPI

end program JOREK2
