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
!*   method          - Solve method ('gmres' or 'direct')                      *
!*   my_id_n         - Identifier of the node in solver communicator.          *
!*   mpi_comm_master - masters MPI communicator.                               *
!*                                                                             *
!* Authors:                                                                    *
!*   Xavier Lacoste - xavier.lacoste@inria.fr                                  *
!*                                                                             *
!*******************************************************************************
SUBROUTINE solve_murge_all(n_cpu,my_id,index_min,index_max, i_tor,  method, & 
     my_id_n, mpi_comm_master)
  !---------------------------------------------------------------------
  ! subroutine solves the complete system of equation using pastix with
  ! distributed matrix on the main group mpi_comm_world
  !---------------------------------------------------------------------
  USE parameters
  USE mumps_module
  USE pastix_module
  USE murge_module
  USE global_distributed_matrix
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  
  ! Subroutine parameters:
  INTEGER                  :: n_cpu
  INTEGER                  :: my_id
  INTEGER                  :: index_min, index_max       ! global index_min, index_max for this cpu
  INTEGER                  :: i_tor(n_cpu)
  CHARACTER*8              :: method
  INTEGER                  :: my_id_n
  INTEGER                  :: mpi_comm_master

  ! local variables:
  REAL*8,ALLOCATABLE       :: column_local(:)
  REAL*8                   :: t_analysis_0, t_analysis_1, t_fact_0, t_fact_1, t_comm_0, t_comm_1
  REAL*8                   :: t_scale_0, t_scale_1
  INTEGER                  :: i, k, j, ierr, m_loc
  INTEGER,ALLOCATABLE      :: counts(:), displacements(:)
  REAL*8, ALLOCATABLE      :: rhs_tmp(:)

  WRITE(*,*) my_id,'*********************************'
  WRITE(*,*) my_id,'*  solve global matrix (PastiX) *'
  WRITE(*,*) my_id,'*********************************'

  if (use_murge_element) then
     m_loc = local_n * n_tor * n_var
  else
     m_loc = (index_max - index_min + 1) * n_tor * n_var
  end if
  mumps_par%nz_loc = nz_glob

  CALL MPI_Allreduce(m_loc,mumps_par%N,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_Allreduce(mumps_par%NZ_loc,mumps_par%nz,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

  !------------------------------------------------------- colunm scaling of global distributed matrix
  CALL CPU_TIME(t_scale_0)


  if (.not. use_murge_element) then
     IF (ALLOCATED(column_scaling))  DEALLOCATE(column_scaling)
     IF (ALLOCATED(column_local))    DEALLOCATE(column_local)
     ALLOCATE(column_scaling(mumps_par%N),column_local(mumps_par%N))

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
     IF (ALLOCATED(counts))        DEALLOCATE(counts)
     IF (ALLOCATED(displacements)) DEALLOCATE(displacements)
     
     ALLOCATE(counts(n_cpu),displacements(n_cpu))
     
     CALL MPI_Allgather(mumps_par%nz_loc,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
     
     displacements(1) = 0
     DO i=2,n_cpu
        displacements(i) = displacements(i-1) + counts(i-1)
     ENDDO
     
     IF (ASSOCIATED(mumps_par%IRN)) DEALLOCATE(mumps_par%IRN)
     IF (ASSOCIATED(mumps_par%JCN)) DEALLOCATE(mumps_par%JCN)
     IF (ASSOCIATED(mumps_par%A) )  DEALLOCATE(mumps_par%A)

     
     ALLOCATE(mumps_par%IRN(mumps_par%nz),mumps_par%JCN(mumps_par%nz),mumps_par%A(mumps_par%nz))


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
     id = 0;
     IF (solver == MURGE_SOLVER_PASTIX) THEN
        CALL MURGE_SetDefaultOptions(id, 0, ierr)
        CALL MURGE_SetOptionINT(id, IPARM_VERBOSE,             API_VERBOSE_YES,  ierr)
        CALL MURGE_SetOptionINT(id, IPARM_MATRIX_VERIFICATION, API_YES,         ierr)
        CALL MURGE_SetOptionINT(id, IPARM_ITERMAX,             murge_iter,      ierr) ! refinement : max number of iterations
        !    CALL MURGE_SetOptionINT(id, MURGE_IPARAM_DOF     , n_tor * n_var,  ierr)   ! degrees of freedom per node (not correct)
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

  ENDIF

  IF ((.NOT. use_murge_element) .AND. (.NOT. pastix_analysed)) THEN

     CALL CPU_TIME(t_analysis_0)

     WRITE(*,*) '***********************************'
     WRITE(*,*) '* analyse Murge                   *'
     WRITE(*,*) '***********************************'
     ! this processor enters the A(myfirstrow:mylastrow, *) 
     ! part of the matrix non-zero pattern
     !CALL MURGE_GraphGlobalIJV(id, n_glob, nz_glob, IRN_glob, JCN_glob, -1, ierr);
     !CALL MURGE_MatrixGlobalIJV(id, n_glob, nz_glob, IRN_glob, JCN_glob, A_glob, -1, MURGE_ASSEMBLY_OVW, murge_sym, ierr);

     write (*,*) 'colptr[n]-1 ', mumps_par%jcn(mumps_par%n+1)-1
     CALL MURGE_GraphGlobalCSC(id, mumps_par%n, mumps_par%jcn, mumps_par%irn, -1, ierr)
     write (*,*) "ierr", ierr
     CALL MURGE_MatrixGlobalCSC(id, mumps_par%n, mumps_par%jcn, mumps_par%irn, mumps_par%A, -1, MURGE_ASSEMBLY_OVW, murge_sym, ierr)
     CALL CPU_TIME(t_analysis_1)

     pastix_analysed = .TRUE.

     IF (my_id .EQ. 0)  WRITE(*,'(A,f8.3)') ' PASTIX, analysis  : ',t_analysis_1-t_analysis_0

  ENDIF

  CALL CPU_TIME(t_fact_0)


  WRITE(*,*) '***********************************'
  WRITE(*,*) '* call MURGE_SetGlobalRhs         *'
  WRITE(*,*) '***********************************'

  CALL MURGE_SetGlobalRhs(id, rhs_glob, -1,MURGE_ASSEMBLY_OVW , ierr)
  if (method .ne. "direct") then
     if (my_id_n .eq. 0) then
        
        if (allocated(deltas)) deallocate(deltas)
        allocate(deltas(ndof_glob))
        deltas = 0.d0
        
        allocate(rhs_tmp(ndof_glob))
        
        rhs_tmp = 0.d0
        
        if (my_id .eq. 0 ) then
           
           rhs_tmp(1:ndof_glob:n_tor) = rhs_glob(1:mumps_par%n)
           
        else
           
           rhs_tmp(2*i_tor(my_id+1)-2:ndof_glob:n_tor) = rhs_glob(1:mumps_par%n:2)
           rhs_tmp(2*i_tor(my_id+1)-1:ndof_glob:n_tor) = rhs_glob(2:mumps_par%n:2)
           
        endif
        
        deallocate(rhs_tmp)
     endif
  end if
  WRITE(*,*) '***********************************'
  WRITE(*,*) '* call MURGE_GetGlobalSolution    *'
  WRITE(*,*) '***********************************'
  IF (ASSOCIATED(mumps_par%rhs)) DEALLOCATE(mumps_par%rhs)     
  ALLOCATE(mumps_par%rhs(mumps_par%n))

  CALL MURGE_GetGlobalSolution(id, mumps_par%rhs, -1, ierr)
  
  CALL CPU_TIME(t_fact_1)

  IF (my_id .EQ. 0) WRITE(*,'(A,f8.3)')  ' PASTIX, fact/solv : ',t_fact_1-t_fact_0
  write (*,*) "mumps_par%n", mumps_par%n, ndof_glob
  write (*,*) "size(deltas)", size(deltas)
  if (method == "direct") then
     DO k=1,mumps_par%n
        deltas(k) =  mumps_par%rhs(k)  / column_scaling(k)
        !  write(*,*) k,deltas(k)
     ENDDO
  else
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
  end if
  deallocate(mumps_par%rhs)
  RETURN
END SUBROUTINE solve_murge_all
