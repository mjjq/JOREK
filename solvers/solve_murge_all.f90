!*******************************************************************************
!* Subroutine: solve_murge_all                                                 *
!*******************************************************************************
!*                                                                             *
!* Solves AX=B using Murge interface (PaStiX).                                 *
!*                                                                             *
!* Parameters:                                                                 *
!*   n_cpu           - Number of processes in MPI_COMM_WORLD                   *
!*   my_id           - Identifier of the node in MPI_COMM_WORLD                *
!*   index_min       - Minimal index of local nodes                            *
!*   index_max       - Maximal index of local nodes                            *
!*   i_tor           - Tor number                                              *
!*   gmres           - Solve method (.true. for gmres, .false for 'direct')    *
!*   my_id_n         - Identifier of the node in solver communicator.          *
!*   mpi_comm_master - masters MPI communicator.                               *
!*                                                                             *
!* Authors:                                                                    *
!*   Xavier Lacoste - xavier.lacoste@inria.fr                                  *
!*                                                                             *
!*******************************************************************************
SUBROUTINE solve_murge_all(n_cpu,my_id,index_min,index_max, i_tor,  gmres, & 
     my_id_n, mpi_comm_n, mpi_comm_master)
  !---------------------------------------------------------------------
  ! subroutine solves the complete system of equation using pastix with
  ! distributed matrix on the main group mpi_comm_world
  !---------------------------------------------------------------------
  USE tr_module 
  USE parameters
  USE mumps_module
  USE pastix_module
  USE murge_module
  USE global_distributed_matrix
  IMPLICIT NONE
  INCLUDE 'mpif.h'
#include "r3_info.h"
  ! Subroutine parameters:
  INTEGER                  :: n_cpu
  INTEGER                  :: my_id
  INTEGER                  :: index_min, index_max       ! global index_min, index_max for this cpu
  INTEGER                  :: i_tor(n_cpu)
  LOGICAL                  :: gmres
  INTEGER                  :: my_id_n
  INTEGER                  :: mpi_comm_master, MPI_COMM_N

  ! local variables:
  REAL*8,ALLOCATABLE       :: column_local(:)
  REAL*8                   :: t_analysis_0, t_analysis_1, t_fact_0, t_fact_1, t_comm_0, t_comm_1
  REAL*8                   :: t_scale_0, t_scale_1
  INTEGER                  :: i, k, j, ierr, m_loc
  INTEGER,ALLOCATABLE      :: counts(:), displacements(:)
  REAL*8, ALLOCATABLE      :: rhs_tmp(:)
  integer, external :: omp_get_num_threads, omp_get_thread_num
  integer                  :: t0,t1,nb_periodes_max,nb_periodes_sec, nb_periods
  CHARACTER(LEN=20), PARAMETER :: FMT_TIMING = "(I2,A70,F7.2)"

#ifdef USE_MURGE
  WRITE(*,*) my_id,'*********************************'
  WRITE(*,*) my_id,'*  solve global matrix (PastiX) *'
  WRITE(*,*) my_id,'*********************************'

  call system_clock(count_rate=nb_periodes_sec,count_max=nb_periodes_max) ! elapsed time
  call r3_info_begin (r3_info_index_0, 'solve_matrix_n')                  ! timing
  if (use_murge_element) then
     m_loc = murge_local_n * n_tor * n_var
  else
     m_loc = (index_max - index_min + 1) * n_tor * n_var
  end if
  mumps_par%nz_loc = nz_glob
  if (gmres) then
     if (.not. use_murge_element) then
        CALL MPI_Allreduce(m_loc,mumps_par%N,1,MPI_INTEGER,MPI_SUM,MPI_COMM_N,ierr)
        CALL MPI_Allreduce(mumps_par%NZ_loc,mumps_par%nz,1,MPI_INTEGER,MPI_SUM,MPI_COMM_N,ierr)
     end  if
  else
     CALL MPI_Allreduce(m_loc,mumps_par%N,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     CALL MPI_Allreduce(mumps_par%NZ_loc,mumps_par%nz,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if
  !------------------------------------------------------- colunm scaling of global distributed matrix
  CALL CPU_TIME(t_scale_0)


  if (.not. use_murge_element) then
     IF (ALLOCATED(column_scaling))  call tr_deallocate(column_scaling,"column_scaling")
     IF (ALLOCATED(column_local))    call tr_deallocate(column_local,"column_local")
     call tr_allocate(column_scaling,1,mumps_par%N,"column_scaling")
     call tr_allocate(column_local,1,mumps_par%N,"column_local")

     column_local = 1.d-20;   column_scaling = 1.d-20
     DO k=1,nz_glob
        j = jcn_glob(k)
        column_local(j) = MAX(column_local(j),ABS(A_glob(k)))
     ENDDO

     CALL MPI_AllReduce(column_local,column_scaling,mumps_par%N,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
     DO k=1,nz_glob
        j = jcn_glob(k)
        A_glob(k) = A_glob(k) / column_scaling(j)
     ENDDO

  end if

  CALL CPU_TIME(t_scale_1)

  IF (my_id .EQ. 0)  WRITE(*,'(A,f8.3)') ' PASTIX, scale     : ',t_scale_1-t_scale_0


  CALL CPU_TIME(t_comm_0)
  IF (.NOT. use_murge_element) THEN
     !------------------------------------------------------ collect the distributed matrix onto all procs
     IF (ALLOCATED(counts))        call tr_deallocate(counts,"counts")
     IF (ALLOCATED(displacements)) call tr_deallocate(displacements,"displacements")
     
     call tr_allocate(counts,1,n_cpu,"counts")
     call tr_allocate(displacements,1,n_cpu,"displacements")
     
     CALL MPI_Allgather(mumps_par%nz_loc,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
     
     displacements(1) = 0
     DO i=2,n_cpu
        displacements(i) = displacements(i-1) + counts(i-1)
     ENDDO
     
     IF (ASSOCIATED(mumps_par%IRN)) call tr_deallocatep(mumps_par%IRN,"mumps_par%IRN")
     IF (ASSOCIATED(mumps_par%JCN)) call tr_deallocatep(mumps_par%JCN,"mumps_par%JCN")
     IF (ASSOCIATED(mumps_par%A) )  call tr_deallocatep(mumps_par%A,"mumps_par%A")

     
     call tr_allocatep(mumps_par%IRN,1,mumps_par%nz,"mumps_par%IRN")
     call tr_allocatep(mumps_par%JCN,1,mumps_par%nz,"mumps_par%JCN")
     call tr_allocatep(mumps_par%A,1,mumps_par%nz,"mumps_par%A")


     CALL MPI_AllgatherV(IRN_glob,mumps_par%nz_loc,MPI_INTEGER,mumps_par%IRN, &
          counts,displacements,MPI_INTEGER,MPI_COMM_WORLD,ierr)
     
     CALL MPI_AllgatherV(JCN_glob,mumps_par%nz_loc,MPI_INTEGER,mumps_par%JCN, &
          counts,displacements,MPI_INTEGER,MPI_COMM_WORLD,ierr)
     
     CALL MPI_AllgatherV(A_glob,mumps_par%nz_loc,MPI_DOUBLE_PRECISION,mumps_par%A, &
          counts,displacements,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     

  END IF
  !IF (ASSOCIATED(mumps_par%rhs)) DEALLOCATE(mumps_par%rhs)     
  !ALLOCATE(mumps_par%rhs(mumps_par%n))
  !CALL MPI_AllReduce(RHS_glob,mumps_par%RHS,mumps_par%N,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  IF (ALLOCATED(sparskit_work)) DEALLOCATE(sparskit_work)
  ALLOCATE(sparskit_work(mumps_par%N + 1))
  IF ((.NOT. use_murge_element) .AND. (.NOT. murge_initialised)) THEN

     CALL coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)

     CALL CPU_TIME(t_comm_1)

     IF (my_id .EQ. 0)  WRITE(*,'(A,f8.3)') ' PASTIX, comm      : ',t_comm_1-t_comm_0

     CALL MURGE_Initialize(1, ierr)
     murge_id = 0;
     CALL MURGE_GetSolver(murge_solver, ierr)
     IF (murge_solver == MURGE_SOLVER_PASTIX) THEN
        CALL MURGE_SetDefaultOptions(murge_id, 0, ierr)
        CALL MURGE_SetOptionINT(murge_id, IPARM_VERBOSE,             API_VERBOSE_YES,  ierr)
        CALL MURGE_SetOptionINT(murge_id, IPARM_MATRIX_VERIFICATION, API_YES,         ierr)
        CALL MURGE_SetOptionINT(murge_id, IPARM_ITERMAX,             murge_iter,      ierr) ! refinement : max number of iterations
        !    CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_DOF     , n_tor * n_var,  ierr)   ! degrees of freedom per node (not correct)

!$omp parallel default(none) shared(pastix_nthrd, murge_nthrd)    
        murge_nthrd = omp_get_num_threads()
!$omp end parallel

        CALL MURGE_SetOptionINT(murge_id, IPARM_THREAD_NBR,          murge_nthrd,     ierr)
        CALL MURGE_SetOptionINT(murge_id, IPARM_LEVEL_OF_FILL,       murge_iluk,      ierr)
        CALL MURGE_SetOptionINT(murge_id, IPARM_INCOMPLETE,          murge_ricar,     ierr)
        CALL MURGE_SetOptionINT(murge_id, IPARM_AMALGAMATION_LEVEL,  murge_amalg,     ierr)
        CALL MURGE_SetOptionINT(murge_id, IPARM_MATRIX_VERIFICATION, API_YES, ierr)

        CALL MURGE_SetOptionREAL(murge_id, DPARM_EPSILON_MAGN_CTRL,    murge_pivot,   ierr)

     ENDIF


     CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_SYM,         murge_sym,   ierr)
     CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_BASEVAL,     1,   ierr)

     CALL MURGE_SetOptionREAL(murge_id, MURGE_RPARAM_EPSILON_ERROR, murge_epsilon, ierr)

     pastix_initialised = .TRUE.

  ENDIF

  IF ((.NOT. use_murge_element) .AND. (.NOT. pastix_analysed)) THEN

     CALL CPU_TIME(t_analysis_0)

     WRITE(*,*) '***********************************'
     WRITE(*,*) '* analyse Murge                   *'
     WRITE(*,*) '***********************************'
     ! this processor enters the A(myfirstrow:mylastrow, *) 
     ! part of the matrix non-zero pattern
     !CALL MURGE_GraphGlobalIJV(murge_id, n_glob, nz_glob, IRN_glob, JCN_glob, -1, ierr);
     !CALL MURGE_MatrixGlobalIJV(murge_id, n_glob, nz_glob, IRN_glob, JCN_glob, A_glob, -1, MURGE_ASSEMBLY_OVW, murge_sym, ierr);

     call system_clock(count=t0)
     CALL MURGE_GraphGlobalCSC(murge_id, mumps_par%n, mumps_par%jcn, mumps_par%irn, -1, ierr)
     call system_clock(count=t1)
     nb_periods = t1-t0
     if (t1<t0) nb_periods = nb_periods + nb_periodes_max   
     write(*,FMT_TIMING) my_id, ' system_clock elapsed time in MURGE_GraphGlobalCSC ',REAL(nb_periods)/nb_periodes_sec

     call system_clock(count=t0)
     CALL MURGE_MatrixGlobalCSC(murge_id, mumps_par%n, mumps_par%jcn, mumps_par%irn, mumps_par%A, -1, MURGE_ASSEMBLY_OVW, murge_sym, ierr)
     call system_clock(count=t1)
     nb_periods = t1-t0
     if (t1<t0) nb_periods = nb_periods + nb_periodes_max   
     write(*,FMT_TIMING) my_id, ' system_clock elapsed time in MURGE_MatrixGlobalCSC ',REAL(nb_periods)/nb_periodes_sec

     pastix_analysed = .TRUE.

     IF (my_id .EQ. 0)  WRITE(*,'(A,f8.3)') ' PASTIX, analysis  : ',t_analysis_1-t_analysis_0

  ENDIF

  CALL CPU_TIME(t_fact_0)


  WRITE(*,*) '***********************************'
  WRITE(*,*) '* call MURGE_SetGlobalRhs         *'
  WRITE(*,*) '***********************************'

  if (gmres) then
     !if (my_id_n .eq. 0) then
        
     if (allocated(deltas)) call tr_deallocate(deltas,"deltas")
     call tr_allocate(deltas,1,ndof_glob,"deltas")
     deltas = 0.d0

  end if
  call system_clock(count=t0)
  if (gmres) then
     CALL MURGE_SetGlobalRhs(murge_id, mumps_par%rhs, 0, MURGE_ASSEMBLY_OVW , ierr)
!      call tr_deallocate(rhs_tmp,"rhs_tmp")
  else
     CALL MURGE_SetGlobalRhs(murge_id, rhs_glob, 0, MURGE_ASSEMBLY_OVW , ierr)
  end if
  call system_clock(count=t1)
  nb_periods = t1-t0
  if (t1<t0) nb_periods = nb_periods + nb_periodes_max   
  write(*,FMT_TIMING) my_id, ' system_clock elapsed time in MURGE_SetGlobalRhs ',REAL(nb_periods)/nb_periodes_sec
  
  IF (ASSOCIATED(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs")     
  call tr_allocatep(mumps_par%rhs,1,murge_global_n*murge_ndof,"mumps_par%rhs")

  call system_clock(count=t0) 
  mumps_par%rhs = 0.0
  CALL MURGE_GetGlobalSolution(murge_id, mumps_par%rhs, 0, ierr)  
  call system_clock(count=t1)
  nb_periods = t1-t0
  if (t1<t0) nb_periods = nb_periods + nb_periodes_max   
  write(*,FMT_TIMING) my_id, ' system_clock elapsed time in MURGE_GetGlobalSolution ',REAL(nb_periods)/nb_periodes_sec

  write (*,*) "mumps_par%n", mumps_par%n, ndof_glob
  write (*,*) "size(deltas)", size(deltas)
  if (.not. gmres) then
     DO k=1,mumps_par%n
        deltas(k) =  mumps_par%rhs(k)  / column_scaling(k)
        !  write(*,*) k,deltas(k)
     ENDDO
  else
     if (my_id_n .eq. 0) then
        
        if (allocated(deltas)) call tr_deallocate(deltas,"deltas")
        call tr_allocate(deltas,1,ndof_glob,"deltas")
        deltas = 0.d0
        
        call tr_allocate(rhs_tmp,1,ndof_glob,"rhs_tmp")

        rhs_tmp = 0.d0
        
        if (my_id .eq. 0 ) then
        
           rhs_tmp(1:ndof_glob:n_tor) = mumps_par%rhs(1:murge_global_n*murge_ndof)
           
        else
        
           rhs_tmp(2*murge_harmonic-2:ndof_glob:n_tor) = mumps_par%rhs(1:murge_global_n*murge_ndof:2)
           rhs_tmp(2*murge_harmonic-1:ndof_glob:n_tor) = mumps_par%rhs(2:murge_global_n*murge_ndof:2)
           
        endif
        
        call MPI_AllReduce(RHS_tmp,deltas,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_MASTER,ierr)

        call tr_deallocate(rhs_tmp,"rhs_tmp")

     endif
  end if
  call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs")
#else
  print *, "Binary built without murge"
  call abort()
#endif
  RETURN
END SUBROUTINE solve_murge_all
