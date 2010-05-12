subroutine solve_matrix_n(my_id,i_tor,MPI_COMM_N,MPI_COMM_MASTER,solve_only)
  !---------------------------------------------------------------------
  ! subroutine solves the system of equation for each harmonic
  ! using mumps with centralised matrix on the group mpi_group_n (mpi_comm_n)
  !---------------------------------------------------------------------
  use parameters

  use mumps_module
  use murge_module
  use pastix_module

  use global_distributed_matrix
  implicit none
  include 'mpif.h'
#include "r3_info.h"

  integer :: i, my_id, i_tor(*), i_reduced, j_reduced, n_i, n_j, index, index1, index2
  integer :: MPI_COMM_N, MPI_COMM_MASTER, my_id_n, n_cpu_n, ierr, my_id_master, n_cpu_master
  real*8  :: t_analysis_0, t_analysis_1, t_fact_0, t_fact_1, t_solv_0, t_solv_1
  real*8, allocatable :: RHS_tmp(:)
  logical :: solve_only
  integer::time_ini_1,time_ini_0,time_facto_1,time_facto_0,nb_periods,nb_periodes_max,nb_periodes_sec
  integer, external :: omp_get_num_threads, omp_get_thread_num
  !Split broadcast
  character*8 :: type
  !Matrix without zeros
  integer                 :: nz2,k,kk
  integer,allocatable     :: irn2(:),jcn2(:),tmp(:)
  real*8,allocatable      :: A2(:)

  call system_clock(count_rate=nb_periodes_sec,count_max=nb_periodes_max) ! elapsed time
  call r3_info_begin (r3_info_index_0, 'solve_matrix_n')                  ! timing

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

     if (use_mumps) then

        !----------------Pre-processing to remove nonzeros in A (MUMPS)
        if (use_matrix_whitout_zeros_mumps) then
           ALLOCATE(tmp(mumps_par%NZ))
           tmp(:)=0
           !
           nz2=0
           DO i=1,mumps_par%NZ
              IF(mumps_par%A(i).NE.0.d0)THEN
                 nz2=nz2+1
                 tmp(nz2)=i
              ENDIF
           ENDDO
           WRITE(*,*) '% zeros',my_id,mumps_par%NZ,nz2,REAL(nz2)/REAL(mumps_par%NZ)*100
           !
           ALLOCATE(A2(nz2))
           ALLOCATE(irn2(nz2))
           ALLOCATE(jcn2(nz2))
           !
           A2(:) = 0.d0
           irn2(:) = 0
           jcn2(:) = 0
           !
           DO i=1,nz2
              A2(i)   = mumps_par%A(tmp(i))
              irn2(i) = mumps_par%irn(tmp(i))
              jcn2(i) = mumps_par%jcn(tmp(i))        
           ENDDO
           !
           DEALLOCATE(tmp)  
           !
           IF (ASSOCIATED(mumps_par%A))   DEALLOCATE(mumps_par%A)
           IF (ASSOCIATED(mumps_par%irn)) DEALLOCATE(mumps_par%irn)
           IF (ASSOCIATED(mumps_par%jcn)) DEALLOCATE(mumps_par%jcn)
           !
           mumps_par%NZ = nz2
           !
           ALLOCATE(mumps_par%A(mumps_par%NZ))
           ALLOCATE(mumps_par%irn(mumps_par%NZ))
           ALLOCATE(mumps_par%jcn(mumps_par%NZ))
           !
           mumps_par%A(:)   = A2(:)
           mumps_par%irn(:) = irn2(:)
           mumps_par%jcn(:) = jcn2(:)
           !
           DEALLOCATE(A2)
           DEALLOCATE(irn2)
           DEALLOCATE(jcn2)
           !
        endif
        !----------------End pre-processing (MUMPS)

        if (my_id_n .eq. 0) then                          ! elapsed time analysis start
           call MPI_Barrier(MPI_COMM_MASTER,ierr)
           call system_clock(count=time_ini_0)
        endif
        mumps_par%JOB = 1                                 ! Analysis, only needed when grid has changed
        call cpu_time(t_analysis_0)

        mumps_par%icntl(7)  = 4                            ! reorderign option (7:automatic, 3:Scotch, 4:PORD, 5:METIS)
        mumps_par%icntl(8)  = 7                            ! row and column scaling
        mumps_par%icntl(14) = 30                           ! MAXS
        mumps_par%icntl(18) = 0

        call DMUMPS(mumps_par)

        call cpu_time(t_analysis_1)

        if (my_id_n .eq.0) write(*, '(i3,A,f8.3)') my_id,' MUMPS, analysis  : ',t_analysis_1-t_analysis_0
        if (my_id_n .eq.0) then                            ! elapsed time analysis end
           call MPI_Barrier(MPI_COMM_MASTER,ierr)
           call system_clock(count=time_ini_1)
           nb_periods = time_ini_1-time_ini_0
           if (time_ini_1<time_ini_0) nb_periods = nb_periods + nb_periodes_max
           write(*,*) 'system_clock elapsed time analysis',REAL(nb_periods)/nb_periodes_sec
        endif

     else ! .not. use_mumps --> use_pastix or use_murge

        if (my_id_n .eq. 0) then

           !----------------Pre-processing to remove nonzeros in A (PASTIX)
           if (use_matrix_whitout_zeros_pastix) then
              ALLOCATE(tmp(mumps_par%NZ))
              tmp(:)=0
              !
              nz2=0
              DO i=1,mumps_par%NZ
                 IF(mumps_par%A(i).NE.0.d0)THEN
                    nz2=nz2+1
                    tmp(nz2)=i
                 ENDIF
              ENDDO
              WRITE(*,*) '% zeros',my_id,mumps_par%NZ,nz2,REAL(nz2)/REAL(mumps_par%NZ)*100
              !
              ALLOCATE(A2(nz2))
              ALLOCATE(irn2(nz2))
              ALLOCATE(jcn2(nz2))
              !
              A2(:) = 0.d0
              irn2(:) = 0
              jcn2(:) = 0
              !
              DO i=1,nz2
                 A2(i)   = mumps_par%A(tmp(i))
                 irn2(i) = mumps_par%irn(tmp(i))
                 jcn2(i) = mumps_par%jcn(tmp(i))        
              ENDDO
              !
              DEALLOCATE(tmp)  
              !
              IF (ASSOCIATED(mumps_par%A))   DEALLOCATE(mumps_par%A)
              IF (ASSOCIATED(mumps_par%irn)) DEALLOCATE(mumps_par%irn)
              IF (ASSOCIATED(mumps_par%jcn)) DEALLOCATE(mumps_par%jcn)
              !
              kk = mumps_par%NZ  !tmp, just for write(*,*) below
              mumps_par%NZ = nz2
              !
              ALLOCATE(mumps_par%A(mumps_par%NZ))
              ALLOCATE(mumps_par%irn(mumps_par%NZ))
              ALLOCATE(mumps_par%jcn(mumps_par%NZ))
              !
              mumps_par%A(:)   = A2(:)
              mumps_par%irn(:) = irn2(:)
              mumps_par%jcn(:) = jcn2(:)
              !
              DEALLOCATE(A2)
              DEALLOCATE(irn2)
              DEALLOCATE(jcn2)
              !
           endif
           !----------------End pre-processing (PASTIX)

           if (allocated(sparskit_work)) deallocate(sparskit_work)
           allocate(sparskit_work(mumps_par%N + 1))

           call coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)

           deallocate(sparskit_work)

        endif


        if (.not. pastix_smp_only) then

           !$omp parallel default(none) shared(pastix_nthrd)    
           pastix_nthrd = omp_get_num_threads()
           !$omp end parallel

           write(*,'(i5,A,i5)') my_id,' PastiX n_threads : ',pastix_nthrd 

           call MPI_BCAST(mumps_par%n,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
           call MPI_BCAST(mumps_par%nz,1,MPI_INTEGER,0,MPI_COMM_N,ierr)

           if (my_id_n .gt. 0) then
              if (associated(mumps_par%irn)) deallocate(mumps_par%irn)
              if (associated(mumps_par%jcn)) deallocate(mumps_par%jcn)
              if (associated(mumps_par%A))   deallocate(mumps_par%A)
              if (associated(mumps_par%rhs)) deallocate(mumps_par%rhs)
              allocate(mumps_par%irn(mumps_par%nz),mumps_par%jcn(mumps_par%nz),mumps_par%a(mumps_par%nz),mumps_par%rhs(mumps_par%n))
           endif

!!$           call MPI_BCAST(mumps_par%IRN,mumps_par%nz,MPI_INTEGER,0,MPI_COMM_N,ierr)
!!$           call MPI_BCAST(mumps_par%JCN,mumps_par%nz,MPI_INTEGER,0,MPI_COMM_N,ierr)
!!$           call MPI_BCAST(mumps_par%A,mumps_par%nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)

           ! Split MPI_BCAST if MPI buffer beyond 2Go
           type='intIRN'
           call split_brodcast(type,MPI_COMM_N)
           type='intJCN'
           call split_brodcast(type,MPI_COMM_N)
           type='double'
           call split_brodcast(type,MPI_COMM_N)

           !----------------PaStiX need an input matrix with symmetric structure   
           IF (use_matrix_whitout_zeros_pastix) then
              CALL pastix_fortran_checkmatrix(pastix_data,MPI_COMM_WORLD,pastix_verb, &
                   pastix_sym,1,mumps_par%N,mumps_par%JCN,mumps_par%IRN,mumps_par%A,-1,1)

              k = mumps_par%JCN(mumps_par%N+1)-1
              IF(k/=nz2)THEN
                 WRITE(*,*) 'New nnz to symmetrize',my_id,k,REAL(k)/REAL(kk)*100
                 DEALLOCATE(mumps_par%IRN)
                 DEALLOCATE(mumps_par%A)
                 ALLOCATE(mumps_par%A(k))
                 ALLOCATE(mumps_par%IRN(k))
                 CALL pastix_fortran_checkmatrix_end(pastix_data,pastix_verb,mumps_par%IRN,mumps_par%A,1)
              ENDIF
           ENDIF
           !----------------End symmetric

        endif

        if  (.not. pastix_initialised)  then

           if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

              if (use_murge) then

                 CALL MURGE_Initialize(n_tor, ierr)
                 if (ierr /= MURGE_SUCCESS) then 
                    write (*,*) "ERROR in MURGE_Initialize"; 
                    STOP
                 end if
                 id = i_tor(my_id+1) -1
                 write (*,*) "id : ", id
                 CALL MURGE_GetSolver(solver, ierr)
                 IF (solver == MURGE_SOLVER_PASTIX) THEN
                    CALL MURGE_SetDefaultOptions(id, 0, ierr)
                    CALL MURGE_SetOptionINT(id, IPARM_VERBOSE,             API_VERBOSE_NO,  ierr)
                    CALL MURGE_SetOptionINT(id, IPARM_MATRIX_VERIFICATION, API_YES,         ierr)
                    CALL MURGE_SetOptionINT(id, IPARM_ITERMAX,             murge_iter,      ierr) ! refinement : max number of iterations
                    !    CALL MURGE_SetOptionINT(id, MURGE_IPARAM_DOF     , n_tor * n_var,  ierr)   ! degrees of freedom per node (not correct)
                    CALL MURGE_SetOptionINT(id, IPARM_THREAD_NBR,          murge_nthrd,     ierr)
                    CALL MURGE_SetOptionINT(id, IPARM_LEVEL_OF_FILL,       murge_iluk,      ierr)
                    CALL MURGE_SetOptionINT(id, IPARM_INCOMPLETE,          murge_ricar,     ierr)
                    CALL MURGE_SetOptionINT(id, IPARM_AMALGAMATION_LEVEL,  murge_amalg,     ierr)
                    CALL MURGE_SetOptionINT(id, IPARM_MATRIX_VERIFICATION, API_YES, ierr)
                    CALL MURGE_SetOptionINT(id, IPARM_PID,                 id, ierr);

                    CALL MURGE_SetOptionREAL(id, DPARM_EPSILON_MAGN_CTRL,    murge_pivot,   ierr)

                 ENDIF


                 CALL MURGE_SetOptionINT(id, MURGE_IPARAM_SYM,         murge_sym,   ierr)
                 CALL MURGE_SetOptionINT(id, MURGE_IPARAM_BASEVAL,     mumps_par%jcn(1),   ierr)

                 CALL MURGE_SetOptionREAL(id, MURGE_RPARAM_EPSILON_ERROR, murge_epsilon, ierr)
                 CALL MURGE_SetCommunicator(id, MPI_COMM_N, ierr)
                 write (*,*) id, "MPI_COMM_N", MPI_COMM_N
              else

                 !          if (pastix_smp_only) pastix_nthrd = n_cpu_n                ! use the size of the MPIgroup for the number of threads

                 !$omp parallel default(none) shared(pastix_nthrd)    
                 pastix_nthrd = omp_get_num_threads()
                 !$omp end parallel

                 pastix_iparm(1)     = 0                                     ! insert default values
                 pastix_iparm(2)     = 0                                     ! initializse
                 pastix_iparm(3)     = 0

                 if (.not. pastix_smp_only) call MPI_BCAST(mumps_par%n,1,MPI_INTEGER,0,MPI_COMM_N,ierr)

                 allocate(pastix_perm_vars(mumps_par%n),pastix_iperm_vars(mumps_par%n))

                 call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
                 pastix_iparm(63) =    i_tor(my_id+1) -1 
		 pastix_iparm(4) = pastix_verb              
              end if
              pastix_initialised = .true.

           endif

        endif !.not. pastix_initialised



        if (.not. pastix_analysed) then

           if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

              if (use_murge) then
                 CALL CPU_TIME(t_analysis_0)

                 WRITE(*,*) '***********************************'
                 WRITE(*,*) '* analyse Murge                   *'
                 WRITE(*,*) '***********************************'
                 CALL MURGE_GraphGlobalCSC(id, mumps_par%n, mumps_par%jcn, mumps_par%irn, -1, ierr)
                 if (ierr /= MURGE_SUCCESS) then 
                    write (*,*) "ERROR in MURGE_GraphGlobalCSC"; 
                    STOP
                 end if
                 CALL CPU_TIME(t_analysis_1)

                 IF (my_id .EQ. 0)  WRITE(*,'(A,f8.3)') ' PASTIX, analysis  : ',t_analysis_1-t_analysis_0

              else
                 if (my_id_n .eq. 0) then                     ! elapsed time analysis start
                    call MPI_Barrier(MPI_COMM_MASTER,ierr)
                    call system_clock(count=time_ini_0)
                 endif
                 pastix_iparm(2) = 1
                 pastix_iparm(3) = 3
                 pastix_iparm(6) = pastix_iter                 ! refinement : max number of iterations

                 pastix_iparm(31) = pastix_facto
                 pastix_iparm(35) = pastix_nthrd               !   numthreads   ! number of threads
                 pastix_iparm(39) = pastix_rhs                 ! right hand side (0 : use RHS)

                 pastix_iparm(41) = pastix_sym

                 pastix_iparm(42) = pastix_ricar
                 pastix_iparm(37) = pastix_iluk
                 pastix_iparm(14) = pastix_amalg

                 pastix_dparm(6)  = pastix_epsilon             ! error level refinement
                 pastix_dparm(11) = pastix_pivot               ! pivot threshold?
                 pastix_iparm(63) =    i_tor(my_id+1) -1                                
                 call cpu_time(t_analysis_0)

                 call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

                 call cpu_time(t_analysis_1)

                 if (my_id_n .eq.0) write(*, '(i3,A,f8.3)') my_id,' PASTIX, analysis  : ',t_analysis_1-t_analysis_0
                 if (my_id_n .eq. 0) then
                    call MPI_Barrier(MPI_COMM_MASTER,ierr)     ! elapsed time analysis end
                    call system_clock(count=time_ini_1)
                    nb_periods = time_ini_1-time_ini_0
                    if (time_ini_1<time_ini_0) nb_periods = nb_periods + nb_periodes_max
                    write(*,*) 'system_clock elapsed time analysis',REAL(nb_periods)/nb_periodes_sec
                 endif

              endif

              pastix_analysed = .true.
           endif

        endif ! .not. pastix_analysed

     endif   ! (else, use_mumps)


     if (my_id_n .eq. 0) then                              ! elapsed time factorisation start
        call MPI_Barrier(MPI_COMM_MASTER,ierr)
        call system_clock(count=time_facto_0)
     endif

     if (use_mumps) then

        call cpu_time(t_fact_0)

        mumps_par%JOB = 2                                   ! factorisation

        call DMUMPS(mumps_par)

        call cpu_time(t_fact_1)

        if (my_id_n .eq.0)   write(*,'(i3,A,f8.3)')  my_id,' MUMPS, fact      : ',t_fact_1-t_fact_0
        if (my_id_n .eq.0)   write(*,'(i3,A,i8)')    my_id,' MUMPS, mem       : ',mumps_par%info(16)

     elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

        if (use_murge) then
           WRITE(*,*) '***********************************'
           WRITE(*,*) '* Matrix Murge                    *'
           WRITE(*,*) '***********************************'
           CALL MURGE_MatrixGlobalCSC(id, mumps_par%n, mumps_par%jcn, mumps_par%irn, mumps_par%A, &
                -1, MURGE_ASSEMBLY_OVW, murge_sym, ierr)
           if (ierr /= MURGE_SUCCESS) then 
              write (*,*) "ERROR in MURGE_MatrixGlobalCSC"; 
              STOP
           end if
        else
           call cpu_time(t_fact_0)

           pastix_iparm(2) = 4
           pastix_iparm(3) = 4
           pastix_iparm(6) = pastix_iter          ! refinement : max number of iterations

           pastix_iparm(31) = pastix_facto
           pastix_iparm(35) = pastix_nthrd         ! numthreads   ! number of threads
           pastix_iparm(39) = pastix_rhs         ! right hand side (0 : use RHS)
           pastix_iparm(41) = pastix_sym

           pastix_iparm(42) = pastix_ricar
           pastix_iparm(37) = pastix_iluk
           pastix_iparm(14) = pastix_amalg

           pastix_dparm(6)  = pastix_epsilon    ! error level refinement
           pastix_dparm(11) = pastix_pivot    ! pivot threshold?
           pastix_iparm(63) =    i_tor(my_id+1) -1                          
           call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

           call cpu_time(t_fact_1)

           if (my_id_n .eq.0)   write(*,'(i3,A,f8.3)')  my_id,' PastiX, fact      : ',t_fact_1-t_fact_0
        end if

     endif

  endif

  if (my_id_n .eq. 0) then                          ! elapsed time factorisation end
     call MPI_Barrier(MPI_COMM_MASTER,ierr)
     call system_clock(count=time_facto_1)
     nb_periods = time_facto_1-time_facto_0
     if (time_facto_1<time_facto_0) nb_periods = nb_periods + nb_periodes_max
     write(*,*) 'system_clock elapsed time factorization',REAL(nb_periods)/nb_periodes_sec
  endif


  if (use_mumps) then

     mumps_par%JOB = 3                                   ! Solve

     call cpu_time(t_solv_0)

     call DMUMPS(mumps_par)

     call cpu_time(t_solv_1)

     if (my_id_n .eq.0)   write(*,'(i3,A,f8.3)')  my_id,' MUMPS, solv      : ',t_solv_1-t_solv_0

  elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

     if (use_murge) then
        if (.not. pastix_smp_only) call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)

        CALL MURGE_SetGlobalRhs(id, mumps_par%rhs, -1,MURGE_ASSEMBLY_OVW , ierr)
        if (ierr /= MURGE_SUCCESS) then 
           write (*,*) "ERROR in MURGE_SetGlobalRhs"; 
           STOP
        end if
        CALL MURGE_GetGlobalSolution(id, mumps_par%rhs, -1, ierr)
        if (ierr /= MURGE_SUCCESS) then 
           write (*,*) "ERROR in MURGE_GetGlobalSolution"; 
           STOP
        end if

     else
        call cpu_time(t_solv_0)

        pastix_iparm(2) = 5
        pastix_iparm(3) = pastix_endsolve
        pastix_iparm(6) = pastix_iter          ! refinement : max number of iterations

        pastix_iparm(31) = pastix_facto
        pastix_iparm(35) = pastix_nthrd ! numthreads   ! number of threads
        pastix_iparm(39) = pastix_rhs            ! right hand side (0 : use RHS)
        pastix_iparm(41) = pastix_sym

        pastix_iparm(42) = pastix_ricar
        pastix_iparm(37) = pastix_iluk
        pastix_iparm(14) = pastix_amalg

        pastix_dparm(6)  = pastix_epsilon    ! error level refinement
        pastix_dparm(11) = pastix_pivot    ! pivot threshold?
        pastix_iparm(63) =    i_tor(my_id+1) -1                       
        if (.not. pastix_smp_only) call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)

        call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
             pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

        call cpu_time(t_solv_1)

        if (my_id_n .eq.0)   write(*,'(i3,A,f8.3)')  my_id,' PastiX, solv      : ',t_solv_1-t_solv_0
     end if
     CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
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
  call r3_info_end (r3_info_index_0)         ! timing
  return
end subroutine solve_matrix_n
