module solve_mat_n

contains
  subroutine pastix_bind_threads(my_id)

    use pastix_module

    integer, intent(in) :: my_id
    integer*4, dimension(1:pastix_nthrd) :: thread_map
    integer*4 k, iplace
    
!    do k = 1, pastix_nthrd/2
!      thread_map(1+2*(k-1)) = my_id * (pastix_nthrd/2) + k-1
!      thread_map(2+2*(k-1)) = my_id * (pastix_nthrd/2) + k-1
!    end do
    do k = 1, pastix_nthrd
      thread_map(k) = mod(my_id * pastix_nthrd,68) + k-1
    end do
    call pastix_fortran_bindthreads(pastix_data, pastix_nthrd, thread_map(1:))

  end subroutine pastix_bind_threads

  !> Solves the system of equation for each harmonic using mumps, pastix, or wsmp
  subroutine solve_matrix_n(my_id,i_tor,MPI_COMM_N,MPI_COMM_MASTER,solve_only)

    use tr_module
    use mod_parameters
    use mumps_module
    use murge_module, only: API_NO, IPARM_MODIFY_PARAMETER, MURGE_EPSILON,     &
         &                  MURGE_RPARAM_EPSILON_ERROR, MURGE_IPARAM_BASEVAL,  &
         &                  MURGE_SYM, MURGE_IPARAM_SYM, MURGE_PIVOT,          &
         &                  DPARM_EPSILON_MAGN_CTRL, IPARM_PID, MURGE_AMALG,   &
         &                  IPARM_AMALGAMATION_LEVEL, MURGE_RICAR,             &
         &                  IPARM_INCOMPLETE, MURGE_ILUK, IPARM_LEVEL_OF_FILL, &
         &                  MURGE_ILUK, IPARM_LEVEL_OF_FILL, MURGE_NTHRD,      &
         &                  IPARM_THREAD_NBR, MURGE_ITER, IPARM_ITERMAX,       &
         &                  API_YES, IPARM_MATRIX_VERIFICATION,                &
         &                  API_VERBOSE_NO, IPARM_VERBOSE,                     &
         &                  MURGE_SOLVER_PASTIX, MURGE_SOLVER, murge_id,       &
         &                  use_murge, MURGE_SUCCESS, DPARM_MEM_MAX,           &
         &                  API_TASK_SOLVE, API_TASK_NUMFACT,                  &
         &                  MURGE_ASSEMBLY_OVW, API_TASK_ANALYSE,              &
         &                  API_TASK_ORDERING, IPARM_DOF_NBR,                  &
         &                  DPARM_EPSILON_REFINEMENT, IPARM_SYM,               &
         &                  IPARM_RHS_MAKING, IPARM_FACTORIZATION,             &
#ifdef WORLDWAR2
         &                  API_THREAD_MULTIPLE, API_THREAD_FUNNELED,          &
#endif
         &                  IPARM_THREAD_COMM_MODE,                            &
         &                  IPARM_END_TASK, API_TASK_INIT, IPARM_START_TASK,   &
         &                  IPARM_BINDTHRD, API_BIND_TAB
#ifdef USE_MURGE
    use murge_module, only : MURGE_MatrixReset
    USE murge_module, only : MURGE_Initialize
    USE murge_module, only : MURGE_GetSolver
    USE murge_module, only : MURGE_SetDefaultOptions
    USE murge_module, only : MURGE_SetOptionINT
    USE murge_module, only : MURGE_SetOptionREAL
    USE murge_module, only : MURGE_SetCommunicator
    USE murge_module, only : MURGE_GraphGlobalCSC
    USE murge_module, only : MURGE_MatrixGlobalCSC
    USE murge_module, only : MURGE_SetGlobalRhs
    USE murge_module, only : MURGE_GetGlobalSolution
#endif
    use wsmp_module
    use pastix_module
    use global_distributed_matrix
    use mpi_mod 
    use mod_clock
    use phys_module, only : index_now
    use mod_coicsr
    implicit none

#include "r3_info.h"

    integer, intent(in) :: my_id
    integer, dimension(:), intent(in) :: i_tor(:)
    integer, intent(in) :: MPI_COMM_N, MPI_COMM_MASTER
    logical, intent(in) :: solve_only

    integer :: i, j, k, my_id_n, n_cpu_n, ierr, my_id_master, n_cpu_master
    integer :: i_reduced, j_reduced, n_i, n_j, index, index1, index2
    type(clcktype) :: t_itstart, t0, t1, t2, t3
    real*8  :: tsecond
    real*8, allocatable :: RHS_tmp(:)
    integer, external :: omp_get_num_threads, omp_get_thread_num
    !Split broadcast
    character*8 :: type
    INTEGER :: increment
    real*8  :: DUMMY_REAL(1:1)
    integer :: DUMMY_INT (1:1)
    CHARACTER(LEN=128) :: fname

    !+increment because of difference between murge.inc (0 based) and pastix_fortran.h (1 based)
    !in old version of PaStiX
    increment=1-IPARM_MODIFY_PARAMETER

    call r3_info_begin (r3_info_index_0, 'solve_matrix_n')                  ! timing
    call tr_print_memsize("BeforeSolveN")
    call tr_debug_writei("smn_A_mumps_par%n",mumps_par%n)

    if (my_id .eq. 0) then
      write(*,*) my_id,'*********************************'
      write(*,*) my_id,'*      solve local matrix  (n)  *'
      write(*,*) my_id,'*********************************'

      if (use_mumps)  write(*,*) my_id,'*       using solver MUMPS      *'
      if (use_pastix) write(*,*) my_id,'*       using solver PastiX     *'
      if (use_wsmp)   write(*,*) my_id,'*       using solver WSMP       *'

      write(*,*) my_id,'*********************************'
    endif

    call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)     ! the id of each cpu
    call MPI_COMM_SIZE(MPI_COMM_N, n_cpu_n, ierr)     ! the number of cpus

    if (my_id_n .eq. 0) then
      call MPI_COMM_RANK(MPI_COMM_MASTER, my_id_master, ierr)     ! the id of each cpu
      call MPI_COMM_SIZE(MPI_COMM_MASTER, n_cpu_master, ierr)     ! the number of cpus
    endif


    if (.not. solve_only) then

      !---------------------------------------- column scaling 
      if (my_id_n .eq. 0) then

        if (allocated(column_scaling))  call tr_deallocate(column_scaling,"column_scaling",CAT_DMATRIX)
        call tr_allocate(column_scaling,1,mumps_par%N,"column_scaling",CAT_DMATRIX)

        column_scaling = 1.d-20
        do k=1,mumps_par%nz
          j = mumps_par%jcn(k)
          column_scaling(j) = min(max(column_scaling(j),abs(mumps_par%A(k))),1d20)
        enddo
        if (my_id .eq. 0) then
           write(fname,'(A,I6.6)')"column_scaling",index_now
           call tr_vdump(fname,column_scaling,mumps_par%N)
        end if
        !CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        write(*,'(2i4,A,2e12.4)') my_id,my_id_n,' COLUMN SCALING : ',minval(column_scaling),maxval(column_scaling)
        do k=1,mumps_par%nz
          j = mumps_par%jcn(k)
          mumps_par%A(k) = mumps_par%A(k) / column_scaling(j)
        enddo
      endif

      if (my_id_n .eq. 0) then                          ! elapsed time analysis start
         call MPI_Barrier(MPI_COMM_MASTER,ierr)
         call clck_time(t0)
      endif


      if (use_mumps) then
#ifdef USE_MUMPS

        mumps_par%JOB = 1                                 ! Analysis, only needed when grid has changed

        mumps_par%icntl(7)  = 4                            ! reorderign option (7:automatic, 3:Scotch, 4:PORD, 5:METIS)
        mumps_par%icntl(8)  = 7                            ! row and column scaling
        mumps_par%icntl(14) = 30                           ! MAXS
        mumps_par%icntl(18) = 0

        call DMUMPS(mumps_par)

#endif
      else ! .not. use_mumps --> use_pastix or use_murge or use_wsmp

        if (my_id_n .eq. 0) then           

          if (my_id_n .eq. 0) then                
            call MPI_Barrier(MPI_COMM_MASTER,ierr)
            call clck_time(t2)
          endif

#ifdef USE_BLOCK
          !---------------------------- reduce IRN,JCN to make use of blocksize ntor*nvar
          !                             temporary solution before using blocks everywhere

          block_size  = n_var
          if (my_id .ne. 0) block_size = 2*n_var

          block_size2 = block_size**2
          !---------------------------- reduce IRN,JCN to make use of blocksize ntor*nvar
          n_block   = mumps_par%n  / block_size
          nnz_block = mumps_par%nz / block_size2

          do i=1,nnz_block  
            mumps_par%irn(i) = (mumps_par%irn((i-1)*block_size2+1) - 1) / block_size + 1 
            mumps_par%jcn(i) = (mumps_par%jcn((i-1)*block_size2+1) - 1) / block_size + 1 
          enddo

          if (allocated(sparskit_work)) deallocate(sparskit_work)
          allocate(sparskit_work(n_block+1))
          call coicsr2(n_block,nnz_block,mumps_par%A,mumps_par%IRN(1:nnz_block),mumps_par%JCN(1:nnz_block),block_size,sparskit_work)

          ! WARNING:  USE_BLOCK does not (yet) work with WSMP!!!
          if (use_wsmp) then
#ifdef USE_WSMP
            call PWGSMP__allocate(n_block, nnz_block, my_id_n)
            call PWGSMP__initialize_matrix(n_block, nnz_block,                                     &
              mumps_par%a, mumps_par%jcn, mumps_par%irn, my_id_n )
#endif
          endif

#else
          if (allocated(sparskit_work)) deallocate(sparskit_work)
          allocate(sparskit_work(mumps_par%N + 1))

          call coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)

          if (use_wsmp) then
#ifdef USE_WSMP
            call PWGSMP__allocate(mumps_par%N, mumps_par%NZ, my_id_n)
            call PWGSMP__initialize_matrix(mumps_par%N, mumps_par%NZ,                              &
              mumps_par%a, mumps_par%jcn, mumps_par%irn, my_id_n )
#endif
          endif
#endif

          if (allocated(sparskit_work)) deallocate(sparskit_work)

          if (my_id_n .eq. 0) then
            call MPI_Barrier(MPI_COMM_MASTER,ierr) 
            call clck_time(t3)
            call clck_ldiff(t2,t3,tsecond)
            write(*,FMT_TIMING) my_id, '### Elapsed time coicsr :', tsecond
          endif

        else  ! (my_id_n > 0) below
#ifdef USE_WSMP
          if (use_wsmp) call PWGSMP__allocate(0, 0, my_id_n)
#endif
        endif ! end (my_id_n .eq. 0)



        ! --- Dstribute data to the MPI "slave" tasks (>0)
        !     (When using WSMP, this is *not necessary* in 0-master mode!)
        if ((.not. use_wsmp).and.(.not. pastix_smp_only)) then

          !$omp parallel default(none) shared(pastix_nthrd)    
          !$omp master
          pastix_nthrd = omp_get_num_threads()/2
          !$omp end master
          !$omp end parallel
          
          if (my_id .eq. 0) write(*,'(I5,A,i5)') my_id,' Second PastiX n_threads : ',pastix_nthrd, "OMP", omp_get_num_threads() 

          call MPI_BCAST(mumps_par%n,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
          call MPI_BCAST(mumps_par%nz,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
#ifdef USE_BLOCK
          call MPI_BCAST(block_size,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
          call MPI_BCAST(n_block,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
          call MPI_BCAST(nnz_block,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
#endif
          if (my_id_n .gt. 0) then
            if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%irn",CAT_DMATRIX)
            if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%jcn",CAT_DMATRIX)
            if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)
            if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs",CAT_DMATRIX)
            call tr_allocatep(mumps_par%irn,1,mumps_par%nz,"mumps_par%irn",CAT_DMATRIX)
            call tr_allocatep(mumps_par%jcn,1,mumps_par%nz,"mumps_par%jcn",CAT_DMATRIX)
            call tr_allocatep(mumps_par%a,1,mumps_par%nz,"mumps_par%a",CAT_DMATRIX)
            call tr_allocatep(mumps_par%rhs,1,mumps_par%n,"mumps_par%rhs",CAT_DMATRIX)
          endif

          ! Split MPI_BCAST if MPI buffer beyond 2Go
          type='intIRN'
          call split_broadcast(type,MPI_COMM_N)
          type='intJCN'
          call split_broadcast(type,MPI_COMM_N)
          type='double'
          call split_broadcast(type,MPI_COMM_N)

        endif

        if  (.not. pastix_initialised)  then

          if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

            if (use_murge) then

#ifdef USE_MURGE
              CALL MURGE_Initialize(n_tor, ierr)
              if (ierr /= MURGE_SUCCESS) then 
                write (*,*) "ERROR in MURGE_Initialize"; 
                STOP
              end if
              murge_id = 0
              write (*,*) "murge_id : ", murge_id
              CALL MURGE_GetSolver(murge_solver, ierr)
              IF (murge_solver == MURGE_SOLVER_PASTIX) THEN
                CALL MURGE_SetDefaultOptions(murge_id, 0, ierr)
                CALL MURGE_SetOptionINT(murge_id, IPARM_VERBOSE,             API_VERBOSE_NO,  ierr)
                CALL MURGE_SetOptionINT(murge_id, IPARM_MATRIX_VERIFICATION, API_YES,         ierr)
                CALL MURGE_SetOptionINT(murge_id, IPARM_ITERMAX,             murge_iter,      ierr) ! refinement : max number of iterations
                !    CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_DOF     , n_tor * n_var,  ierr)   ! degrees of freedom per node (not correct)
                CALL MURGE_SetOptionINT(murge_id, IPARM_THREAD_NBR,          murge_nthrd,     ierr)
                CALL MURGE_SetOptionINT(murge_id, IPARM_LEVEL_OF_FILL,       murge_iluk,      ierr)
                CALL MURGE_SetOptionINT(murge_id, IPARM_INCOMPLETE,          murge_ricar,     ierr)
                CALL MURGE_SetOptionINT(murge_id, IPARM_AMALGAMATION_LEVEL,  murge_amalg,     ierr)
                CALL MURGE_SetOptionINT(murge_id, IPARM_MATRIX_VERIFICATION, API_YES, ierr)
                CALL MURGE_SetOptionINT(murge_id, IPARM_PID,                 murge_id, ierr);

                CALL MURGE_SetOptionREAL(murge_id, DPARM_EPSILON_MAGN_CTRL,    murge_pivot,   ierr)

              ENDIF


              CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_SYM,         murge_sym,   ierr)
              CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_BASEVAL,     mumps_par%jcn(1),   ierr)

              CALL MURGE_SetOptionREAL(murge_id, MURGE_RPARAM_EPSILON_ERROR, murge_epsilon, ierr)
              CALL MURGE_SetCommunicator(murge_id, MPI_COMM_N, ierr)
              write (*,*) murge_id, "MPI_COMM_N", MPI_COMM_N
#else
              print *, "Binary built without murge"
              call abort()
#endif
            elseif (use_pastix) then

              !          if (pastix_smp_only) pastix_nthrd = n_cpu_n                ! use the size of the MPIgroup for the number of threads

              !$omp parallel default(none) shared(pastix_nthrd)    
              !$omp master
              pastix_nthrd = omp_get_num_threads()/2
              !$omp end master
              !$omp end parallel
          if (my_id .eq. 0) write(*,'(I5,A,i5)') my_id,' First PastiX n_threads : ',pastix_nthrd, "OMP", omp_get_num_threads() 
              pastix_iparm(IPARM_MODIFY_PARAMETER+increment) = API_NO         ! insert default values
              pastix_iparm(IPARM_START_TASK+increment)       = API_TASK_INIT  ! initializse
              pastix_iparm(IPARM_END_TASK+increment)         = API_TASK_INIT
!              pastix_iparm(IPARM_BINDTHRD+increment)         = API_NO
              if (.not. pastix_smp_only) call MPI_BCAST(mumps_par%n,1,MPI_INTEGER,0,MPI_COMM_N,ierr)

#ifdef USE_BLOCK
              call tr_allocate(pastix_perm_vars,1,n_block,"pastix_perm_vars",CAT_UNKNOWN)
              call tr_allocate(pastix_iperm_vars,1,n_block,"pastix_iperm_vars",CAT_UNKNOWN)

              call pastix_fortran(pastix_data,MPI_COMM_N,n_block,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#else
              call tr_allocate(pastix_perm_vars,1 ,mumps_par%n,"pastix_perm_vars",CAT_UNKNOWN)
              call tr_allocate(pastix_iperm_vars,1,mumps_par%n,"pastix_iperm_vars",CAT_UNKNOWN)

              call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#endif

              pastix_iparm(IPARM_VERBOSE+increment)            = pastix_verb              
              pastix_iparm(IPARM_ITERMAX+increment)            = pastix_iter                ! refinement : max number of iterations

              pastix_iparm(IPARM_FACTORIZATION+increment)      = pastix_facto
              pastix_iparm(IPARM_THREAD_NBR+increment)         = pastix_nthrd               ! number of threads
              pastix_iparm(IPARM_RHS_MAKING+increment)         = pastix_rhs                 ! right hand side (0 : use RHS)

              pastix_iparm(IPARM_SYM+increment)                = pastix_sym

              pastix_iparm(IPARM_INCOMPLETE+increment)         = pastix_ricar
              pastix_iparm(IPARM_LEVEL_OF_FILL+increment)      = pastix_iluk
              pastix_iparm(IPARM_AMALGAMATION_LEVEL+increment) = pastix_amalg

#ifdef WORLDWAR2
#ifdef FUNNELED
              pastix_iparm(IPARM_THREAD_COMM_MODE+increment)  = API_THREAD_FUNNELED
#else
              pastix_iparm(IPARM_THREAD_COMM_MODE+increment)  = API_THREAD_MULTIPLE
#endif
#endif

              pastix_dparm(DPARM_EPSILON_REFINEMENT+increment) = pastix_epsilon             ! error level refinement
              pastix_dparm(DPARM_EPSILON_MAGN_CTRL+increment)  = pastix_pivot               ! pivot threshold
#ifdef USE_BLOCK
              pastix_iparm(IPARM_DOF_NBR+increment)            = block_size                 ! block size
#else
              pastix_iparm(IPARM_DOF_NBR+increment)            = 1
#endif
            else if (use_wsmp) then
#ifdef USE_WSMP
              call PWGSMP__initialize_solver(my_id_n, MPI_COMM_N)
#endif
            end if
            pastix_initialised = .true.

          endif

        endif !.not. pastix_initialised



        if (.not. pastix_analysed) then

          if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

            if (use_murge) then

              WRITE(*,*) '***********************************'
              WRITE(*,*) '* analyse Murge                   *'
              WRITE(*,*) '***********************************'
#ifdef USE_MURGE
              CALL MURGE_GraphGlobalCSC(murge_id, mumps_par%n, mumps_par%jcn, mumps_par%irn, -1, ierr)
#else
              print *, "Binary built without murge"
              call abort()
#endif
              if (ierr /= MURGE_SUCCESS) then 
                write (*,*) "ERROR in MURGE_GraphGlobalCSC"; 
                STOP
              end if

            else if (use_pastix) then

              pastix_iparm(IPARM_START_TASK+increment) = API_TASK_ORDERING
              pastix_iparm(IPARM_END_TASK+increment)   = API_TASK_ANALYSE
!              pastix_iparm(IPARM_BINDTHRD+increment)   = API_NO
#ifdef USE_BLOCK
              call pastix_fortran(pastix_data,MPI_COMM_N, n_block, &
                mumps_par%jcn(1:n_block+1), mumps_par%irn(1:nnz_block), mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#else
              call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#endif
            else if (use_wsmp) then
              ! do nothing
            endif

            pastix_analysed = .true.
          endif

        endif ! .not. pastix_analysed

      endif   ! (else, use_mumps)


      if (my_id_n .eq.0) then                            ! elapsed time analysis end
         call MPI_Barrier(MPI_COMM_MASTER,ierr)
         call clck_time(t1)
         call clck_ldiff(t0,t1,tsecond)
         write(*, FMT_TIMING) my_id,' ## Elapsed time, analysis :',tsecond
         call clck_time(t0)                              ! elapsed time facto start 
      endif

      if (use_mumps) then
#ifdef USE_MUMPS

        mumps_par%JOB = 2                                   ! factorisation

        call DMUMPS(mumps_par)

        if (my_id_n .eq.0)   write(*,'(i3,A,i8)')    my_id,' MUMPS, mem       : ',mumps_par%info(16)
#endif
      elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

        if (use_murge) then
          WRITE(*,*) '***********************************'
          WRITE(*,*) '* Matrix Murge                    *'
          WRITE(*,*) '***********************************'
#ifdef USE_MURGE
          CALL MURGE_MatrixGlobalCSC(murge_id, mumps_par%n, mumps_par%jcn, mumps_par%irn, mumps_par%A, &
            -1, MURGE_ASSEMBLY_OVW, murge_sym, ierr)
#else
          print *, "Binary built without murge"
          call abort()
#endif
          if (ierr /= MURGE_SUCCESS) then 
            write (*,*) "ERROR in MURGE_MatrixGlobalCSC"; 
            STOP
          end if

          call tr_deallocatep(mumps_par%A,"special:mumps_par%A",CAT_DMATRIX)
          call tr_deallocatep(mumps_par%irn,"special:mumps_par%irn",CAT_DMATRIX)
          call tr_deallocatep(mumps_par%jcn,"special:mumps_par%jcn",CAT_DMATRIX)

        else if (use_pastix) then

          pastix_iparm(IPARM_START_TASK+increment) = API_TASK_NUMFACT
          pastix_iparm(IPARM_END_TASK+increment)   = API_TASK_NUMFACT
          pastix_iparm(IPARM_BINDTHRD+increment)   = API_BIND_TAB
#ifdef USE_BLOCK
          call pastix_bind_threads(my_id)
          call pastix_fortran(pastix_data,MPI_COMM_N, n_block, &
            mumps_par%jcn, mumps_par%irn, mumps_par%A, &
            pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#else	   
          call pastix_bind_threads(my_id)
          call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
            pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#endif
        else if (use_wsmp) then
#ifdef USE_WSMP
          call PWGSMP__LU_factorization(my_id_n)
#endif
        end if

      endif
      
      if (my_id_n .eq.0) then                            ! elapsed time facto end
         call MPI_Barrier(MPI_COMM_MASTER,ierr)
         call clck_time(t1)
         call clck_ldiff(t0,t1,tsecond)
         write(*, FMT_TIMING) my_id,' ## Elapsed time, facto :',tsecond
         call clck_time(t0)
      end if
   endif
   call tr_debug_writei("smn_B_mumps_par%n",mumps_par%n)


   if (my_id_n .eq. 0) then                          ! elapsed time solve start
      call MPI_Barrier(MPI_COMM_MASTER,ierr)
      call clck_time(t0)
   endif
   if (use_mumps) then
#ifdef USE_MUMPS
      mumps_par%JOB = 3                                   ! Solve

      call DMUMPS(mumps_par)

#endif
    elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

      if (use_murge) then
#ifdef USE_MURGE
        call tr_locvnorms("smn_rhs",mumps_par%rhs,mumps_par%n)
        CALL MURGE_SetGlobalRhs(murge_id, mumps_par%rhs, 0,MURGE_ASSEMBLY_OVW , ierr)
        if (ierr /= MURGE_SUCCESS) then 
          write (*,*) "ERROR in MURGE_SetGlobalRhs"; 
          STOP
        end if
        CALL MURGE_GetGlobalSolution(murge_id, mumps_par%rhs, 0, ierr)
        if (ierr /= MURGE_SUCCESS) then 
          write (*,*) "ERROR in MURGE_GetGlobalSolution"; 
          STOP
        end if
#else
        print *, "Binary built without murge"
        call abort()
#endif
      else if (use_pastix) then

        pastix_iparm(IPARM_START_TASK+increment) = API_TASK_SOLVE
        pastix_iparm(IPARM_END_TASK+increment)   = pastix_endsolve
!        pastix_iparm(IPARM_BINDTHRD+increment)   = API_NO
        if (.not. pastix_smp_only) then
           call tr_debug_writei("smn_C_mumps_par%n",mumps_par%n)
           call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
        end if
        if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%irn",CAT_DMATRIX)
        if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%jcn",CAT_DMATRIX)
        if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)

        call tr_locvnorms("smn_rhs",mumps_par%rhs,mumps_par%n)

#ifdef USE_BLOCK
        call pastix_fortran(pastix_data,MPI_COMM_N, n_block,                &
!             mumps_par%jcn,mumps_par%irn,mumps_par%A, &
             DUMMY_INT,DUMMY_INT,DUMMY_REAL, &
             pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#else
        call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,&
!             mumps_par%jcn,mumps_par%irn,mumps_par%A, &
          DUMMY_INT,DUMMY_INT,DUMMY_REAL, &
          pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#endif
      else if (use_wsmp) then
#ifdef USE_WSMP
        call PWGSMP__back_substitution(mumps_par%rhs, my_id_n)
#endif
     end if
   endif

    if (my_id_n .eq.0) then                            ! elapsed time solve end
       call MPI_Barrier(MPI_COMM_MASTER,ierr)
       call clck_time(t1)
       call clck_ldiff(t0,t1,tsecond)
       write(*, FMT_TIMING) my_id,' ## Elapsed time, solve :',tsecond
       call clck_time(t0)
    end if


    if (my_id_n .eq. 0) then

      !------------------------------------------ undo column scaling
      do k=1,mumps_par%n
        mumps_par%rhs(k) =  mumps_par%rhs(k) / column_scaling(k)
      enddo

      if (allocated(deltas)) call tr_deallocate(deltas,"deltas",CAT_PRECOND)
      call tr_allocate(deltas,1,ndof_glob,"deltas",CAT_PRECOND)
      deltas = 0.d0

      call tr_allocate(rhs_tmp,1,ndof_glob,"rhs_tmp",CAT_PRECOND)

      rhs_tmp = 0.d0

      if (my_id .eq. 0 ) then
        !        rhs_tmp(1:ndof_glob:n_tor) = mumps_par%rhs(1:mumps_par%n)
        do i=0, mumps_par%n-1
          rhs_tmp(1+i*n_tor)=mumps_par%rhs(1+i)
        end do
      else
        !        rhs_tmp(2*i_tor(my_id+1)-2:ndof_glob:n_tor) = mumps_par%rhs(1:mumps_par%n:2)
        !        rhs_tmp(2*i_tor(my_id+1)-1:ndof_glob:n_tor) = mumps_par%rhs(2:mumps_par%n:2)
        do i=0, mumps_par%n/2-1
          rhs_tmp(2*i_tor(my_id+1)-2+i*n_tor) = mumps_par%rhs(1+i*2)
          rhs_tmp(2*i_tor(my_id+1)-1+i*n_tor) = mumps_par%rhs(2+i*2)
        end do

      endif

      call MPI_AllReduce(RHS_tmp,deltas,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_MASTER,ierr)
      call tr_deallocate(rhs_tmp,"rhs_tmp",CAT_PRECOND)

      call tr_locvnorms("smn_res",mumps_par%rhs,mumps_par%n)
      call tr_locvnorms("smn_delta",deltas,ndof_glob)
    endif
    if (.not. use_murge) then
      call tr_set_precondmem(pastix_dparm(DPARM_MEM_MAX+increment)) 
    end if
    call tr_print_memsize("AfterSolveN")
    call r3_info_end (r3_info_index_0)         ! timing
    return
  end subroutine solve_matrix_n
end module solve_mat_n
