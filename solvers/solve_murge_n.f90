subroutine solve_murge_n(my_id,i_tor,MPI_COMM_N,MPI_COMM_MASTER,solve_only)
  !---------------------------------------------------------------------
  ! subroutine solves the system of equation for each harmonic
  ! using mumps with centralised matrix on the group mpi_group_n (mpi_comm_n)
  !---------------------------------------------------------------------
  use parameters

  use mumps_module
  use pastix_module
  USE murge_module

  use global_distributed_matrix
  implicit none
  include 'mpif.h'

  integer :: i, my_id, i_tor(*), i_reduced, j_reduced, n_i, n_j, index, index1, index2
  integer :: MPI_COMM_N, MPI_COMM_MASTER, my_id_n, n_cpu_n, ierr, my_id_master, n_cpu_master
  real*8  :: t_analysis_0, t_analysis_1, t_fact_0, t_fact_1, t_solv_0, t_solv_1
  real*8, allocatable :: RHS_tmp(:)
  logical :: solve_only

  integer, external :: omp_get_num_threads, omp_get_thread_num

  if (my_id .eq. 0) then
     write(*,*) my_id,'*********************************'
     write(*,*) my_id,'*      solve local matrix  (n)  *'
     write(*,*) my_id,'*********************************'

     if (use_mumps)  write(*,*) my_id,'*       using solver MUMPS      *'
     if (use_pastix) write(*,*) my_id,'*       using solver PastiX     *'

     write(*,*) my_id,'*********************************'
  endif

  call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)     ! the id of each cpu
  call MPI_COMM_SIZE(MPI_COMM_N, n_cpu_n, ierr)     ! the number of cpus

  if (my_id_n .eq. 0) then
     call MPI_COMM_RANK(MPI_COMM_MASTER, my_id_master, ierr)     ! the id of each cpu
     call MPI_COMM_SIZE(MPI_COMM_MASTER, n_cpu_master, ierr)     ! the number of cpus
  endif


  if (.not. solve_only) then


     if  ((.NOT. use_murge_element) .AND. (.not. murge_initialised))  then

        if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

           !          if (pastix_smp_only) pastix_nthrd = n_cpu_n                ! use the size of the MPIgroup for the number of threads

           !$omp parallel default(none) shared(pastix_nthrd)    
           !pastix_nthrd = omp_get_num_threads()
           !$omp end parallel
           !CALL coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)

           CALL MURGE_Initialize(1, ierr)
           id = 0;
           IF (solver == MURGE_SOLVER_PASTIX) THEN
              CALL MURGE_SetDefaultOptions(id, 0, ierr)
              CALL MURGE_SetOptionINT(id, IPARM_VERBOSE,             API_VERBOSE_YES,  ierr)
              CALL MURGE_SetOptionINT(id, IPARM_MATRIX_VERIFICATION, API_YES,         ierr)
              CALL MURGE_SetOptionINT(id, IPARM_ITERMAX,             murge_iter,      ierr) ! refinement : max number of iterations
              !    CALL MURGE_SetOptionINT(id, MURGE_IPARAM_DOF     , n_tor * n_var,  ierr)   
              ! degrees of freedom per node (not correct)
              CALL MURGE_SetOptionINT(id, IPARM_THREAD_NBR,          murge_nthrd,     ierr)
              CALL MURGE_SetOptionINT(id, IPARM_LEVEL_OF_FILL,       murge_iluk,      ierr)
              CALL MURGE_SetOptionINT(id, IPARM_INCOMPLETE,          murge_ricar,     ierr)
              CALL MURGE_SetOptionINT(id, IPARM_AMALGAMATION_LEVEL,  murge_amalg,     ierr)
              CALL MURGE_SetOptionINT(id, IPARM_MATRIX_VERIFICATION, API_YES, ierr)

              CALL MURGE_SetOptionREAL(id, DPARM_EPSILON_MAGN_CTRL,    murge_pivot,   ierr)

           ENDIF


           CALL MURGE_SetOptionINT(id, MURGE_IPARAM_SYM,         murge_sym,   ierr)
           CALL MURGE_SetOptionINT(id, MURGE_IPARAM_BASEVAL,     1,   ierr)

           CALL MURGE_SetOptionREAL(id, MURGE_RPARAM_EPSILON_ERROR, murge_epsilon, ierr)

           pastix_initialised = .TRUE.

        endif

     endif

     if ((.not. use_murge_element).and. (.not. pastix_analysed)) then

!        if (my_id_n .eq. 0) then
!
!           if (allocated(sparskit_work)) deallocate(sparskit_work)
!           allocate(sparskit_work(mumps_par%N + 1))
!
!           call coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)
!
!           deallocate(sparskit_work)
!
!        endif
!
        !$omp parallel default(none) shared(pastix_nthrd)    
        pastix_nthrd = omp_get_num_threads()
        !$omp end parallel

        write(*,'(i5,A,i5)') my_id,' PastiX n_threads : ',pastix_nthrd 
!        if (.not. pastix_smp_only) then 
!           call MPI_BCAST(mumps_par%n,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
!           call MPI_BCAST(mumps_par%nz,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
!
!           if (my_id_n .gt. 0) then
!              if (associated(mumps_par%irn)) deallocate(mumps_par%irn)
!              if (associated(mumps_par%jcn)) deallocate(mumps_par%jcn)
!              if (associated(mumps_par%A))   deallocate(mumps_par%A)
!              if (associated(mumps_par%rhs)) deallocate(mumps_par%rhs)
!              allocate(mumps_par%irn(mumps_par%nz),mumps_par%jcn(mumps_par%nz),mumps_par%a(mumps_par%nz),mumps_par%rhs(mumps_par%n))
!           endif
!
!           call MPI_BCAST(mumps_par%IRN,mumps_par%nz,MPI_INTEGER,0,MPI_COMM_N,ierr)
!           call MPI_BCAST(mumps_par%JCN,mumps_par%nz,MPI_INTEGER,0,MPI_COMM_N,ierr)
!           call MPI_BCAST(mumps_par%A,  mumps_par%nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
!           call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
!
!
!        endif
!
        !if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

           CALL CPU_TIME(t_analysis_0)

           WRITE(*,*) '***********************************'
           WRITE(*,*) '* analyse Murge                   *'
           WRITE(*,*) '***********************************'
           ! this processor enters the A(myfirstrow:mylastrow, *) 
           ! part of the matrix non-zero pattern
           !CALL MURGE_GraphGlobalIJV(id, n_glob, nz_glob, IRN_glob, JCN_glob, -1, ierr);
           !CALL MURGE_MatrixGlobalIJV(id, n_glob, nz_glob, IRN_glob, JCN_glob, A_glob, -1, MURGE_ASSEMBLY_OVW, murge_sym, ierr);
           write (*,*) 'colptr[n]-1 ', mumps_par%jcn(mumps_par%n+1)-1
           CALL MURGE_GraphGlobalIJV(id, mumps_par%n,mumps_par%nz, mumps_par%jcn, mumps_par%irn, 0, ierr)
           write (*,*) "ierr", ierr
           CALL MURGE_MatrixGlobalIJV(id, mumps_par%n, mumps_par%nz, mumps_par%jcn, mumps_par%irn, mumps_par%A, 0, MURGE_ASSEMBLY_OVW, murge_sym, ierr)
           CALL CPU_TIME(t_analysis_1)

           pastix_analysed = .TRUE.

        !endif

     endif ! use_pastix, analysis

  endif   ! (.not., solve_only)


  if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

     call cpu_time(t_fact_0)

     CALL MURGE_SetGlobalRhs(id, mumps_par%rhs, -1,MURGE_ASSEMBLY_OVW , ierr)
     CALL MURGE_GetGlobalSolution(id, mumps_par%rhs, -1, ierr)
     call cpu_time(t_fact_1)

     if (my_id_n .eq.0)   write(*,'(i3,A,f8.3)')  my_id,' PastiX, fact      : ',t_fact_1-t_fact_0

  endif


  if (my_id_n .eq. 0) then

     if (allocated(deltas)) deallocate(deltas)
     allocate(deltas(ndof_glob))
     deltas = 0.d0

     allocate(rhs_tmp(ndof_glob))

     rhs_tmp = 0.d0

     if (my_id .eq. 0 ) then

        rhs_tmp(1:ndof_glob:n_tor) = mumps_par%rhs(1:mumps_par%n)

     else

        rhs_tmp(2*i_tor(my_id+1)-2:ndof_glob:n_tor) = mumps_par%rhs(1:mumps_par%n:2)
        rhs_tmp(2*i_tor(my_id+1)-1:ndof_glob:n_tor) = mumps_par%rhs(2:mumps_par%n:2)

     endif

     call MPI_AllReduce(RHS_tmp,deltas,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_MASTER,ierr)

     deallocate(rhs_tmp)

  endif

  return
end subroutine solve_murge_n
