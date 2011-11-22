module solve_mat_n
contains
subroutine solve_matrix_n(my_id,i_tor,MPI_COMM_N,MPI_COMM_MASTER,solve_only)
  !---------------------------------------------------------------------
  ! subroutine solves the system of equation for each harmonic
  ! using mumps with centralised matrix on the group mpi_group_n (mpi_comm_n)
  !---------------------------------------------------------------------
  use tr_module 
  use parameters

  use mumps_module
  use murge_module
  use pastix_module

  use global_distributed_matrix
  implicit none
  include 'mpif.h'
#include "r3_info.h"
  
  integer, intent(in) :: my_id
  integer, dimension(:), intent(in) :: i_tor(:)
  integer, intent(in) :: MPI_COMM_N, MPI_COMM_MASTER
  logical, intent(in) :: solve_only

  integer :: i, my_id_n, n_cpu_n, ierr, my_id_master, n_cpu_master
  integer :: i_reduced, j_reduced, n_i, n_j, index, index1, index2
  real*8  :: t_analysis_0, t_analysis_1, t_fact_0, t_fact_1, t_solv_0, t_solv_1
  real*8, allocatable :: RHS_tmp(:)
  integer :: t1, t0, time_ini_1,time_ini_0,time_facto_1,time_facto_0, time_solve_0
  integer :: time_solve_1,nb_periods,nb_periodes_max,nb_periodes_sec
  integer, external :: omp_get_num_threads, omp_get_thread_num
  !Split broadcast
  character*8 :: type
  !Matrix without zeros
  integer                 :: nz2,k,kk
  integer,allocatable     :: irn2(:),jcn2(:),tmp(:)
  real*8,allocatable      :: A2(:)
  
  real*8,  pointer :: dummy_real(:)
  integer, pointer :: dummy_int(:)
  dummy_real => NULL()
  dummy_int  => NULL()

  call system_clock(count_rate=nb_periodes_sec,count_max=nb_periodes_max) ! elapsed time
  call r3_info_begin (r3_info_index_0, 'solve_matrix_n')                  ! timing
  call tr_print_memsize("BeforeSolveN")
     
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
#ifdef USE_MUMPS
        !----------------Pre-processing to remove nonzeros in A (MUMPS)
        if (no_zeros_mumps) then
           call tr_allocate(tmp,1,mumps_par%NZ,"tmp",CAT_PRECOND)
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
           call tr_allocate(A2,1,nz2,"A2",CAT_PRECOND)
           call tr_allocate(irn2,1,nz2,"irn2",CAT_PRECOND)
           call tr_allocate(jcn2,1,nz2,"jcn2",CAT_PRECOND)
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
           call tr_deallocate(tmp,"tmp",CAT_PRECOND)  
           !
           IF (ASSOCIATED(mumps_par%A))   call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)
           IF (ASSOCIATED(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%irn",CAT_DMATRIX)
           IF (ASSOCIATED(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%jcn",CAT_DMATRIX)
           !
           mumps_par%NZ = nz2
           !
           call tr_allocatep(mumps_par%A,1,mumps_par%NZ,"mumps_par%A",CAT_DMATRIX)
           call tr_allocatep(mumps_par%irn,1,mumps_par%NZ,"mumps_par%irn",CAT_DMATRIX)
           call tr_allocatep(mumps_par%jcn,1,mumps_par%NZ,"mumps_par%jcn",CAT_DMATRIX)
           !
           mumps_par%A(:)   = A2(:)
           mumps_par%irn(:) = irn2(:)
           mumps_par%jcn(:) = jcn2(:)
           !
           call tr_deallocate(a2,"a2",CAT_PRECOND)
           call tr_deallocate(irn2,"irn2",CAT_PRECOND)
           call tr_deallocate(jcn2,"jcn2",CAT_PRECOND)
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
#endif
     else ! .not. use_mumps --> use_pastix or use_murge

        if (my_id_n .eq. 0) then

           !----------------Pre-processing to remove nonzeros in A (PASTIX)
           if (no_zeros_pastix) then
              call tr_allocate(tmp,1,mumps_par%NZ,"tmp",CAT_PRECOND)
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
              call tr_allocate(A2,1,nz2,"A2",CAT_PRECOND)
              call tr_allocate(irn2,1,nz2,"irn2",CAT_PRECOND)
              call tr_allocate(jcn2,1,nz2,"jcn2",CAT_PRECOND)
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
              call tr_deallocate(tmp,"tmp",CAT_PRECOND)  
              !
              IF (ASSOCIATED(mumps_par%A))   call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)
              IF (ASSOCIATED(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%irn",CAT_DMATRIX)
              IF (ASSOCIATED(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%jcn",CAT_DMATRIX)
              !
              kk = mumps_par%NZ  !tmp, just for write(*,*) below
              mumps_par%NZ = nz2
              !
              call tr_allocatep(mumps_par%A,1,mumps_par%NZ,"mumps_par%A",CAT_DMATRIX)
              call tr_allocatep(mumps_par%irn,1,mumps_par%NZ,"mumps_par%irn",CAT_DMATRIX)
              call tr_allocatep(mumps_par%jcn,1,mumps_par%NZ,"mumps_par%jcn",CAT_DMATRIX)
              !
              mumps_par%A(:)   = A2(:)
              mumps_par%irn(:) = irn2(:)
              mumps_par%jcn(:) = jcn2(:)
              !
              call tr_deallocate(a2,"a2",CAT_PRECOND)
              call tr_deallocate(irn2,"irn2",CAT_PRECOND)
              call tr_deallocate(jcn2,"jcn2",CAT_PRECOND)
              !
           endif
           !----------------End pre-processing (PASTIX)

           
           
           if (my_id_n .eq. 0) then                     ! elapsed time analysis start
              call MPI_Barrier(MPI_COMM_MASTER,ierr)
              call system_clock(count=t0)
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

#else
           if (allocated(sparskit_work)) deallocate(sparskit_work)
           allocate(sparskit_work(mumps_par%N + 1))

           call coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)

#endif

           if (my_id_n .eq. 0) then
              call MPI_Barrier(MPI_COMM_MASTER,ierr)     ! elapsed time analysis end
              call system_clock(count=t1)
              nb_periods = t1-t0
              if (t1<t0) nb_periods = nb_periods + nb_periodes_max
              write(*,*) 'system_clock elapsed time coicsr',REAL(nb_periods)/nb_periodes_sec
           endif
           deallocate(sparskit_work)

        endif


        if (.not. pastix_smp_only) then

           !$omp parallel default(none) shared(pastix_nthrd)    
           pastix_nthrd = omp_get_num_threads()
           !$omp end parallel

           if (my_id .eq. 0) write(*,'(i5,A,i5)') my_id,' PastiX n_threads : ',pastix_nthrd 

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
           IF (no_zeros_pastix) then
              CALL pastix_fortran_checkmatrix(pastix_data,MPI_COMM_WORLD,pastix_verb, &
                   pastix_sym,1,mumps_par%N,mumps_par%JCN,mumps_par%IRN,mumps_par%A,-1,1)

              k = mumps_par%JCN(mumps_par%N+1)-1
              IF(k/=nz2)THEN
                 WRITE(*,*) 'New nnz to symmetrize',my_id,k,REAL(k)/REAL(kk)*100
                 call tr_deallocatep(mumps_par%IRN,"mumps_par%IRN",CAT_DMATRIX)
                 call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)
                 call tr_allocatep(mumps_par%A,1,k,"mumps_par%A",CAT_DMATRIX)
                 call tr_allocatep(mumps_par%IRN,1,k,"mumps_par%IRN",CAT_DMATRIX)
                 CALL pastix_fortran_checkmatrix_end(pastix_data,pastix_verb,mumps_par%IRN,mumps_par%A,1)
              ENDIF
           ENDIF
           !----------------End symmetric

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
              else

                 !          if (pastix_smp_only) pastix_nthrd = n_cpu_n                ! use the size of the MPIgroup for the number of threads

                 !$omp parallel default(none) shared(pastix_nthrd)    
                 pastix_nthrd = omp_get_num_threads()
                 !$omp end parallel
                 !+1 because there is a difference between murge.inc (0 based) and pastix_fortran.h (1 based) 
                 pastix_iparm(IPARM_MODIFY_PARAMETER+1) = API_NO         ! insert default values
                 pastix_iparm(IPARM_START_TASK+1)       = API_TASK_INIT  ! initializse
                 pastix_iparm(IPARM_END_TASK+1)         = API_TASK_INIT

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

		 pastix_iparm(IPARM_VERBOSE+1)            = pastix_verb              
                 pastix_iparm(IPARM_ITERMAX+1)            = pastix_iter                ! refinement : max number of iterations

                 pastix_iparm(IPARM_FACTORIZATION+1)      = pastix_facto
                 pastix_iparm(IPARM_THREAD_NBR+1)         = pastix_nthrd               ! number of threads
                 pastix_iparm(IPARM_RHS_MAKING+1)         = pastix_rhs                 ! right hand side (0 : use RHS)

                 pastix_iparm(IPARM_SYM+1)                = pastix_sym

                 pastix_iparm(IPARM_INCOMPLETE+1)         = pastix_ricar
                 pastix_iparm(IPARM_LEVEL_OF_FILL+1)      = pastix_iluk
                 pastix_iparm(IPARM_AMALGAMATION_LEVEL+1) = pastix_amalg

                 pastix_dparm(DPARM_EPSILON_REFINEMENT+1) = pastix_epsilon             ! error level refinement
                 pastix_dparm(DPARM_EPSILON_MAGN_CTRL+1)  = pastix_pivot               ! pivot threshold
#ifdef USE_BLOCK
                 pastix_iparm(IPARM_DOF_NBR+1)            = block_size                 ! block size
#else
                 pastix_iparm(IPARM_DOF_NBR+1)            = 1
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

              else
                 if (my_id_n .eq. 0) then                     ! elapsed time analysis start
                    call MPI_Barrier(MPI_COMM_MASTER,ierr)
                    call system_clock(count=time_ini_0)
                 endif
                 pastix_iparm(IPARM_START_TASK+1) = API_TASK_ORDERING
                 pastix_iparm(IPARM_END_TASK+1)   = API_TASK_ANALYSE
                 call cpu_time(t_analysis_0)

#ifdef USE_BLOCK
  
                 call pastix_fortran(pastix_data,MPI_COMM_N, n_block, &
                      mumps_par%jcn(1:n_block+1), mumps_par%irn(1:nnz_block), mumps_par%A, &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#else

                 call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#endif
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



     if (use_mumps) then
#ifdef USE_MUMPS
        call cpu_time(t_fact_0)

        mumps_par%JOB = 2                                   ! factorisation

        call DMUMPS(mumps_par)

        call cpu_time(t_fact_1)

        if (my_id_n .eq.0)   write(*,'(i3,A,f8.3)')  my_id,' MUMPS, fact      : ',t_fact_1-t_fact_0
        if (my_id_n .eq.0)   write(*,'(i3,A,i8)')    my_id,' MUMPS, mem       : ',mumps_par%info(16)
#endif
     elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

        CALL CPU_TIME(t_analysis_0)
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

           CALL CPU_TIME(t_analysis_1)
                 
           IF (my_id_n .EQ. 0)  WRITE(*,'(i3,A,f8.3)') my_id, ' MURGE_MatrixGlobalCSC  : ',t_analysis_1-t_analysis_0
        else
           call cpu_time(t_fact_0)

           pastix_iparm(IPARM_START_TASK+1) = API_TASK_NUMFACT
           pastix_iparm(IPARM_END_TASK+1)   = API_TASK_NUMFACT
           if (my_id_n .eq. 0) then                              ! elapsed time factorisation start
              call MPI_Barrier(MPI_COMM_MASTER,ierr)
              call system_clock(count=time_facto_0)
           endif
	   
#ifdef USE_BLOCK
           call pastix_fortran(pastix_data,MPI_COMM_N, n_block, &
                      mumps_par%jcn(1:n_block+1), mumps_par%irn(1:nnz_block), mumps_par%A, &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#else	   
           call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#endif
           call cpu_time(t_fact_1)

           if (my_id_n .eq.0)   write(*,'(i3,A,f8.3)')  my_id,' PastiX, fact      : ',t_fact_1-t_fact_0
           
           if (my_id_n .eq. 0) then                          ! elapsed time factorisation end
              call MPI_Barrier(MPI_COMM_MASTER,ierr)
              call system_clock(count=time_facto_1)
              nb_periods = time_facto_1-time_facto_0
              if (time_facto_1<time_facto_0) nb_periods = nb_periods + nb_periodes_max
              write(*,*) 'system_clock elapsed time factorization',REAL(nb_periods)/nb_periodes_sec
           endif
        end if

     endif

  endif


  if (my_id_n .eq. 0) then                              ! elapsed time solve start
     call MPI_Barrier(MPI_COMM_MASTER,ierr)
     call system_clock(count=time_solve_0)
  endif

  if (use_mumps) then
#ifdef USE_MUMPS
     mumps_par%JOB = 3                                   ! Solve

     call cpu_time(t_solv_0)

     call DMUMPS(mumps_par)

     call cpu_time(t_solv_1)

     if (my_id_n .eq.0)   write(*,'(i3,A,f8.3)')  my_id,' MUMPS, solv      : ',t_solv_1-t_solv_0
#endif
  elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

     if (use_murge) then
#ifdef USE_MURGE
        CALL MURGE_SetGlobalRhs(murge_id, mumps_par%rhs, 0,MURGE_ASSEMBLY_OVW , ierr)
        if (ierr /= MURGE_SUCCESS) then 
           write (*,*) "ERROR in MURGE_SetGlobalRhs"; 
           STOP
        end if
        call cpu_time(t_fact_0)
        CALL MURGE_GetGlobalSolution(murge_id, mumps_par%rhs, 0, ierr)
        call cpu_time(t_fact_1)
        if (ierr /= MURGE_SUCCESS) then 
           write (*,*) "ERROR in MURGE_GetGlobalSolution"; 
           STOP
        end if
        if (my_id_n .eq.0)   write(*,'(i3,A,f8.3)')  my_id,' MURGE_GetGlobalSolution    : ',t_fact_1-t_fact_0

#else
        print *, "Binary built without murge"
        call abort()
#endif
     else
        call cpu_time(t_solv_0)

        pastix_iparm(IPARM_START_TASK+1) = API_TASK_SOLVE
        pastix_iparm(IPARM_END_TASK+1)   = pastix_endsolve
        if (.not. pastix_smp_only) call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
        if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%irn",CAT_DMATRIX)
        if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%jcn",CAT_DMATRIX)
        if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)
#ifdef USE_BLOCK
        call pastix_fortran(pastix_data,MPI_COMM_N, n_block,                                         &
                    NULL(),NULL(),NULL(), &
                    pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#else
        call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,NULL(),NULL(),NULL(), &
             pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#endif
        call cpu_time(t_solv_1)

        if (my_id_n .eq.0)  write(*,'(i3,A,f8.3)')  my_id,' PastiX, solv      : ',t_solv_1-t_solv_0
     end if
  endif


  CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (my_id_n .eq. 0) then                          ! elapsed time factorisation end
     call MPI_Barrier(MPI_COMM_MASTER,ierr)
     call system_clock(count=time_solve_1)
     nb_periods = time_solve_1-time_solve_0
     if (time_solve_1<time_solve_0) nb_periods = nb_periods + nb_periodes_max
     write(*,*) 'system_clock elapsed time resolution',REAL(nb_periods)/nb_periodes_sec
  endif

  if (my_id_n .eq. 0) then

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

  endif
  if (.not. use_murge) then
     call tr_set_precondmem(pastix_dparm(DPARM_MEM_MAX+1)) 
  end if
  call tr_print_memsize("AfterSolveN")
  call r3_info_end (r3_info_index_0)         ! timing
  return
end subroutine solve_matrix_n
end module solve_mat_n
