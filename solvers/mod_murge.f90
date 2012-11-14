!> Variables and and routines related to the MURGE solver interface.
MODULE murge_module
  USE tr_module
  IMPLICIT NONE

  INCLUDE "murge.inc"

  !include "hips.inc"
  ! Indicate which solver is used
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_solver

  ! Solver identification number
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_id
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_id_prod
  INTEGER                                     :: murge_harmonic
  INTEGER                                     :: murge_ntor
  INTEGER                                     :: murge_first_tor
  INTEGER                                     :: murge_last_tor
  ! Local number of element
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_local_n
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_local_n_prod
  ! Global number of element
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_global_n
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_global_n_prod
  ! Number of dof by node
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_ndof
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_ndof_prod

  ! Local element list
  INTEGER(KIND=MURGE_INTS_KIND), ALLOCATABLE  :: murge_loc2glob(:)
  INTEGER(KIND=MURGE_INTS_KIND), ALLOCATABLE  :: murge_loc2glob_prod(:)
  ! Global element list
  INTEGER(KIND=MURGE_INTS_KIND), ALLOCATABLE  :: murge_glob2loc(:)
  INTEGER(KIND=MURGE_INTS_KIND), ALLOCATABLE  :: murge_glob2loc_prod(:)

  ! Indicate if we want to use murge or classical interface
  LOGICAL                                     :: use_murge

  ! Indicate if we want to use murge element building
  LOGICAL                                     :: use_murge_element
  LOGICAL                                     :: use_hips
  LOGICAL                                     :: murge_need_rebuild_sequence(2)
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_sequence_id(2)
  ! Indicate if murge has been initialized
  LOGICAL                                     :: murge_initialised

  ! Murge right-hand-side member
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_rhs
  PARAMETER (murge_rhs=0)

  ! Indicate if the matrix is symmetric
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_sym
  PARAMETER (murge_sym=MURGE_BOOLEAN_FALSE)

  ! Number of threads
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_nthrd =1

  ! Number of iteration in refinement
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_iter
  PARAMETER (murge_iter=10)

  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_ricar
  PARAMETER (murge_ricar=0)

  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_iluk
  PARAMETER (murge_iluk=3)

  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_amalg
  PARAMETER (murge_amalg=5)

  REAL*8                                      :: murge_epsilon
  PARAMETER (murge_epsilon=1.d-12)

  REAL*8                                      :: murge_pivot
  PARAMETER (murge_pivot=1.d-64)
  INTEGER                                     :: murge_comm
  INTEGER(KIND=MURGE_INTS_KIND), POINTER      :: MURGE_ROWS(:), MURGE_COLS(:)
  REAL(KIND=MURGE_COEF_KIND),    POINTER      :: MURGE_VALS(:)
  INTEGER,           ALLOCATABLE, TARGET      :: murge_assembly_first_entry(:)
  INTEGER                                     :: murge_assembly_step
  INTEGER                                     :: murge_elem_block_size
CONTAINS
  !>
  !! Subroutine: murge_add_one_entry
  !!
  !! Add one entry to the product and/or harminic matrix.
  !!
  !! @param index_node  row node index
  !! @param k           row var index
  !! @param in          row tor index
  !! @param index_node2 col node index
  !! @param k2          col var index
  !! @param in2         col tor index
  !! @param zbig        value
  !! @param solve_only  Do not add to harmonic matrix if .true.
  !! @param gmres       Do not add to product matrix if .false.
  !!
  SUBROUTINE murge_add_one_entry( index_node, k, in, index_node2, k2,          &
       &                          in2, zbig, murge_ntor, solve_only, gmres )
    USE parameters

    INTEGER :: index_node,  k,  in
    INTEGER :: index_node2, k2, in2
    REAL*8  :: zbig
    INTEGER :: murge_ntor
    INTEGER(kind=MURGE_INTS_KIND) :: ierr
    LOGICAL :: solve_only, gmres
    INTEGER(KIND=MURGE_INTS_KIND) :: row_idx, col_idx
#ifdef USE_MURGE
    IF (.NOT. solve_only) THEN
       row_idx = murge_ndof * (index_node - 1) + (k -1)*murge_ntor + 1
       col_idx = murge_ndof * (index_node2- 1) + (k2-1)*murge_ntor + 1
       IF (in /= 1) row_idx = row_idx + MOD(in, 2)
       IF (in2 /= 1) col_idx = col_idx + MOD(in2, 2)
       CALL MURGE_ASSEMBLYSETVALUE( murge_id, row_idx, col_idx, zbig,          &
            &                       ierr )
       IF (ierr /= MURGE_SUCCESS) THEN
          WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ",                     &
               "I", index_node, murge_ndof, k, in,                             &
               "J", index_node2, k2, in2
          STOP
       END IF
    END IF
    IF (gmres) THEN
       row_idx = n_tor*n_var * (index_node - 1) + (k -1)*n_tor + in
       col_idx = n_tor*n_var * (index_node2- 1) + (k2-1)*n_tor + in2
       CALL MURGE_ASSEMBLYSETVALUE( murge_id_prod, row_idx, col_idx, zbig,     &
            &                       ierr )

       IF (ierr /= MURGE_SUCCESS) THEN
          WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE(prod) ",               &
               "I", index_node, n_tor, n_var, k, in,                           &
               "J", index_node2, k2, in2
          STOP
       END IF
    END IF
#else
       print *, "Binary built without murge"
       call abort()
#endif
  END SUBROUTINE murge_add_one_entry


  !> Subroutine: murge_initialization
  !!
  !! Initialize murge instances.
  !!
  !! @param gmres       boolean indicating if we are running in GMRES mode.
  !! @param my_id       Rank of the process in MPI_COMM_WORLD
  !! @param mpi_comm_n  By harmonic communicator.
  !! @param i_tor       Harmonic ID.
  SUBROUTINE murge_initialization(gmres, my_id, mpi_comm_n, i_tor)
    USE parameters, ONLY : n_tor, n_var
    USE mpi_mod

    LOGICAL, INTENT(IN) :: gmres
    INTEGER, INTENT(IN) :: MPI_COMM_N
    INTEGER, INTENT(IN) :: my_id
    INTEGER, INTENT(IN) :: i_tor(:)

    INTEGER, EXTERNAL :: omp_get_num_threads
    INTEGER(KIND=MURGE_INTS_KIND) :: ierr
    INTEGER(KIND=MURGE_INTS_KIND) :: ndof

#ifdef USE_MURGE

    murge_need_rebuild_sequence(1) = .TRUE.
    murge_need_rebuild_sequence(2) = .TRUE.
    murge_harmonic  = 1
    murge_ntor      = n_tor
    murge_first_tor = 1
    murge_last_tor  = n_tor
    if (gmres) then
       murge_harmonic  = i_tor(my_id+1)
       murge_ntor      = 2
       murge_first_tor = 2*(murge_harmonic-1)
       murge_last_tor  = 2*(murge_harmonic-1)+1
       if (murge_harmonic .eq. 1) then
          murge_ntor      = 1
          murge_first_tor = 1
          murge_last_tor  = 1
       end if
    end if
    CALL tr_debug_write("murge_initialised begin")
    IF (.NOT. murge_initialised) THEN
       murge_id = 0
       !$omp PARALLEL shared(murge_nthrd)
       murge_nthrd = omp_get_num_threads()
       !$omp end PARALLEL
       CALL MURGE_GetSolver(murge_solver, ierr)

       IF (gmres) THEN
          murge_id_prod = 1
          CALL MURGE_Initialize(2, ierr)
          IF (murge_solver == MURGE_SOLVER_PASTIX) THEN
             CALL MURGE_SetDefaultOptions(murge_id_prod, 0, ierr)
             CALL MURGE_SetOptionINT(murge_id_prod, IPARM_VERBOSE,             &
                  &                  API_VERBOSE_YES, ierr)
             CALL MURGE_SetOptionINT(murge_id_prod,                            &
                  &                  IPARM_MATRIX_VERIFICATION,                &
                  &                  API_YES,         ierr)
             CALL MURGE_SetOptionINT(murge_id_prod, IPARM_THREAD_NBR,          &
                  &                  murge_nthrd,     ierr)
             CALL MURGE_SetOptionINT(murge_id_prod, MURGE_IPARAM_SYM,          &
                  &                  murge_sym,       ierr)
             CALL MURGE_SetOptionINT(murge_id_prod, MURGE_IPARAM_BASEVAL,      &
                  &                  1,               ierr)
             murge_ndof_prod = n_tor*n_var
!#define MURGE_PROD_NODE
#  ifdef MURGE_PROD_NODE
             ndof = murge_ndof_prod
#  else
             ndof = 1
#  endif
             CALL MURGE_SetOptionINT(murge_id_prod, MURGE_IPARAM_DOF,          &
                  &                  ndof,     ierr)
             CALL MURGE_SetCommunicator(murge_id_prod, MPI_COMM_WORLD, ierr)
             CALL MURGE_SetOptionREAL(murge_id_prod, DPARM_EPSILON_MAGN_CTRL,  &
                  &                   murge_pivot,     ierr)

             CALL MURGE_SetDefaultOptions(murge_id,      0, ierr)
!             CALL MURGE_SetOptionINT(murge_id, IPARM_MURGE_REFINMENT,         &
!                  &                  API_NO,     ierr)
             CALL MURGE_SetCommunicator(murge_id, MPI_COMM_N, ierr)
          END IF
       ELSE
          CALL MURGE_Initialize(1, ierr)
          IF (murge_solver == MURGE_SOLVER_PASTIX) THEN
             CALL MURGE_SetDefaultOptions(murge_id,      0, ierr)
             CALL MURGE_SetCommunicator(murge_id, MPI_COMM_WORLD, ierr)
          END IF
       END IF

       IF (murge_solver == MURGE_SOLVER_PASTIX) THEN

          CALL MURGE_SetOptionINT(murge_id, IPARM_VERBOSE,                     &
               &                  API_VERBOSE_YES, ierr)
          CALL MURGE_SetOptionINT(murge_id, IPARM_MATRIX_VERIFICATION,         &
               &                  API_YES,         ierr)
          ! refinement : max number of iterations
          CALL MURGE_SetOptionINT(murge_id, IPARM_ITERMAX,                     &
               &                  murge_iter,      ierr)
          CALL MURGE_SetOptionINT(murge_id, IPARM_MURGE_REFINEMENT,            &
               &                  API_NO,      ierr)
          ! degrees of freedom per node (not correct)
          IF (gmres) THEN
             IF ( i_tor(my_id+1) == 1 ) THEN
                murge_ndof = n_var
             ELSE
                murge_ndof = 2*n_var
             END IF
             CALL MURGE_SetOptionINT(murge_id, IPARM_DOF_COST,                 &
                  &                  2*n_var,      ierr)
          ELSE
             murge_ndof = n_tor*n_var
          END IF
          CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_DOF,                  &
               &                  murge_ndof,      ierr)
          ! TODO : omp_num_thread
          CALL MURGE_SetOptionINT(murge_id, IPARM_THREAD_NBR,                  &
               &                  murge_nthrd,     ierr)
          CALL MURGE_SetOptionINT(murge_id, IPARM_LEVEL_OF_FILL,               &
               &                  murge_iluk,      ierr)
          CALL MURGE_SetOptionINT(murge_id, IPARM_INCOMPLETE,                  &
               &                  murge_ricar,     ierr)
          CALL MURGE_SetOptionINT(murge_id, IPARM_AMALGAMATION_LEVEL,          &
               &                  murge_amalg,     ierr)
          CALL MURGE_SetOptionINT(murge_id, IPARM_MATRIX_VERIFICATION,         &
               &                  API_YES,         ierr)

          CALL MURGE_SetOptionREAL(murge_id, DPARM_EPSILON_MAGN_CTRL,          &
               &                   murge_pivot,     ierr)

       ENDIF


       CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_SYM,                     &
            &                  murge_sym,   ierr)
       CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_BASEVAL,                 &
            &                  1,   ierr)

       CALL MURGE_SetOptionREAL(murge_id, MURGE_RPARAM_EPSILON_ERROR,          &
            &                   murge_epsilon, ierr)

       murge_initialised = .TRUE.
    END IF
    CALL tr_debug_write("murge_initialised end")
#else
       PRINT *, "Binary built without murge"
       CALL abort()
#endif
  END SUBROUTINE murge_initialization

  SUBROUTINE murge_termination(gmres)
    LOGICAL :: gmres

    integer ierr

#ifdef USE_MURGE
    IF ( use_murge_element ) THEN
       IF (ALLOCATED(murge_glob2loc))                                          &
            call tr_deallocate(murge_glob2loc,"murge_glob2loc",CAT_DMATRIX)
       IF (ALLOCATED(murge_loc2glob))                                          &
            call tr_deallocate(murge_loc2glob,"murge_loc2glob",CAT_DMATRIX)
       IF (ALLOCATED(murge_loc2glob_prod))                                     &
            call tr_deallocate(murge_loc2glob_prod,"murge_loc2glob_prod",      &
            &                  CAT_DMATRIX)
       IF (ALLOCATED(murge_assembly_first_entry))                              &
            call tr_deallocate(murge_assembly_first_entry,                     &
            &                  "murge_assembly_first_entry", CAT_DMATRIX)
    END IF
    CALL MURGE_Clean(murge_id, ierr)
    if (gmres) then
       CALL MURGE_Clean(murge_id_prod, ierr)
    end if
#endif

  END SUBROUTINE murge_termination

  SUBROUTINE murge_setGraph(gmres, n, local_elms, n_local_elms, element_list, &
       &                    node_list, n_aa, my_id, my_id_trans, n_cpu_trans, &
       &                    MPI_COMM_N, MPI_COMM_TRANS)

    USE parameters,                ONLY : n_order, n_vertex_max
    USE data_structure,            ONLY : type_element, type_element_list,     &
         &                                type_node_list
    USE parameters,                ONLY : n_tor, n_var
    USE global_distributed_matrix, ONLY : ndof_glob
    USE clock_module,              ONLY : clcktype, clck_time, clck_ldiff,     &
         &                                FMT_TIMING
    USE mpi_mod,                   ONLY : MPI_REAL8, MPI_MAX, MPI_COMM_WORLD,  &
         &                                MPI_MIN, MPI_INTEGER, MPI_SUM

    logical,                  INTENT(IN)    :: gmres
    integer,                  INTENT(IN)    :: n
    integer, ALLOCATABLE,     INTENT(INOUT) :: local_elms(:)
    integer,                  INTENT(INOUT) :: n_local_elms
    type (type_element_list), INTENT(IN)    :: element_list
    type (type_node_list),    INTENT(INOUT) :: node_list
    integer,                  INTENT(IN)    :: n_aa
    integer,                  INTENT(IN)    :: my_id
    integer,                  INTENT(IN)    :: my_id_trans
    integer,                  INTENT(IN)    :: n_cpu_trans
    integer,                  INTENT(IN)    :: MPI_COMM_N, MPI_COMM_TRANS

    integer(KIND=MURGE_INTS_KIND) :: ierr
    integer(KIND=MURGE_INTL_KIND) :: nnz
    integer :: start
    integer :: i_elem
    integer :: i, j, k, inode1, knode, i_order, k_order
    integer(KIND=MURGE_INTS_KIND) :: index_node1, index_node2
    type (type_element)      :: element
    type(clcktype)           :: t0, t1
    REAL*8                   :: max_time, min_time, tsecond
    logical                  :: is_local
    integer                  :: sum_n_local_elms, max_n_local_elms, min_n_local_elms
    integer                  :: index_total, inode 

#ifdef USE_MURGE
    call tr_debug_write("murge_setgraph begin")



    murge_global_n      = n
    murge_global_n_prod = n
#ifndef MURGE_PROD_NODE
    murge_global_n_prod = murge_global_n_prod*n_tor*n_var
#endif
    nnz = n_local_elms*(n_order+1)*(n_order+1)*n_vertex_max*n_vertex_max

    Call MURGE_GRAPHBEGIN(murge_id, murge_global_n, nnz, ierr)
    IF (ierr /= MURGE_SUCCESS) THEN
       write (*,*) "ERROR in MURGE_GRAPHBEGIN"
       STOP
    END IF
    call clck_time(t0)

    call tr_debug_write("murge_setgraph loop start")
    DO i_elem = 1, n_local_elms

       element = element_list%element(local_elms(i_elem))
       DO i = 1, n_vertex_max

          inode1         = element%vertex(i)

          DO i_order = 1, n_order+1

             index_node1 = node_list%node(inode1)%index(i_order)

             ! Build nodes Matrices
             DO k= 1, n_vertex_max

                knode         = element%vertex(k)

                DO k_order = 1, n_order+1

                   index_node2 = node_list%node(knode)%index(k_order)


                   CALL MURGE_GRAPHEDGE(murge_id,                              &
                        index_node1,                                           &
                        index_node2,                                           &
                        ierr)
                   IF (ierr /= MURGE_SUCCESS) THEN
                      write (*,*) "N", n, n_AA,&
                           "I", index_node1, &
                           "J", index_node2
                      STOP
                   END IF

                END DO
             END DO
          END DO
       END DO
    END DO

    call tr_debug_write("murge_setgraph loop end")
    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    CALL MPI_Reduce(tsecond, max_time, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce(tsecond, min_time, 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    if (my_id .eq. 0) then
       write(*,FMT_TIMING) my_id, '# Elapsed time entering graph :',min_time
       write(*,FMT_TIMING) my_id, '# Elapsed time entering graph :',max_time
    end if


    call clck_time(t0)
    CALL MURGE_GRAPHEND(murge_id, ierr)
    IF (ierr /= MURGE_SUCCESS) THEN
       write (*,*) "ERROR in MURGE_GRAPHEND"
       STOP
    END IF
    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    CALL MPI_Reduce(tsecond, max_time, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce(tsecond, min_time, 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    if (my_id .eq. 0) then
       write(*,FMT_TIMING) my_id, '# Elapsed time in MURGE_GRAPHEND :',min_time
       write(*,FMT_TIMING) my_id, '# Elapsed time in MURGE_GRAPHEND :',max_time
    end if

    call clck_time(t0)
    CALL MURGE_GETLOCALNODENBR(murge_id, murge_local_n, ierr)
    IF (ierr /= MURGE_SUCCESS) THEN
       write (*,*) "ERROR in MURGE_GETLOCALNODENBR"
       STOP
    END IF
    murge_local_n_prod =                                                       &
         (murge_local_n-MOD(murge_local_n, n_cpu_trans))/n_cpu_trans
    if (my_id_trans .lt. MOD(murge_local_n, n_cpu_trans)) then
       murge_local_n_prod = murge_local_n_prod + 1
    end if
#ifndef MURGE_PROD_NODE
    murge_local_n_prod = murge_local_n_prod * murge_ndof_prod
    print *, my_id, "murge_local_n_prod", murge_local_n_prod
#endif
    murge_need_rebuild_sequence(1) = .true.
    CALL MURGE_PRODUCTSETLOCALNODENBR(murge_id_prod, murge_local_n_prod, ierr)
    IF (ierr /= MURGE_SUCCESS) THEN
       write (*,*) "ERROR in MURGE_PRODUCTSETLOCALNODENBR"
       STOP
    END IF
    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    CALL MPI_Reduce(tsecond, max_time, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce(tsecond, min_time, 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    if (my_id .eq. 0) then
       write(*,FMT_TIMING) my_id, '# Elapsed time in MURGE_GETLOCALNODENBR :',min_time
       write(*,FMT_TIMING) my_id, '# Elapsed time in MURGE_GETLOCALNODENBR :',max_time
    end if

    call tr_allocate(murge_loc2glob, 1, murge_local_n, "murge_loc2glob",       &
         &           CAT_DMATRIX)
    call tr_allocate( murge_loc2glob_prod, 1, murge_local_n_prod,              &
         &            "murge_loc2glob_prod",CAT_DMATRIX)
    write (*,*) "Local number of nodes", murge_local_n, "global",  n

    call clck_time(t0)
    CALL MURGE_GETLOCALNODELIST(murge_id, murge_loc2glob, ierr)
    if (allocated(murge_glob2loc))                                             &
         call tr_deallocate(murge_glob2loc,"murge_glob2loc",CAT_DMATRIX)

    IF (ierr /= MURGE_SUCCESS) THEN
       write (*,*) "ERROR in MURGE_GETLOCALNODELIST"
       STOP
    END IF
    start = murge_local_n_prod*my_id_trans
#ifndef MURGE_PROD_NODE
    start = start/ murge_ndof_prod
#endif
    start = start + MIN(my_id_trans, MOD(murge_local_n,n_cpu_trans))

#ifdef MURGE_PROD_NODE
    do i = 1, murge_local_n_prod
       murge_loc2glob_prod(i) = murge_loc2glob(start + i)
    end do
#else
    do i = 1, murge_local_n_prod/murge_ndof_prod
       do j = 1, murge_ndof_prod
          murge_loc2glob_prod((i-1)*murge_ndof_prod+j) =                       &
               (murge_loc2glob(start + i)-1)*murge_ndof_prod+j
       end do
    end do
#endif
    CALL MURGE_PRODUCTSETLOCALNODELIST(murge_id_prod, murge_loc2glob_prod, ierr)
    if (allocated(murge_glob2loc_prod))                                        &
         call tr_deallocate(murge_glob2loc_prod,"murge_glob2loc_prod",CAT_DMATRIX)

    IF (ierr /= MURGE_SUCCESS) THEN
       write (*,*) "ERROR in MURGE_PRODUCTSETLOCALNODELIST"
       STOP
    END IF
    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    CALL MPI_Reduce(tsecond, max_time, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce(tsecond, min_time, 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    if (my_id .eq. 0) then
       write(*,FMT_TIMING) my_id, '# Elapsed time in MURGE_GETLOCALNODELIST :',min_time
       write(*,FMT_TIMING) my_id, '# Elapsed time in MURGE_GETLOCALNODELIST :',max_time
    end if

    call clck_time(t0)
    ! Build local_elms from loc2glob
    n_local_elms = 0

    DO i_elem = 1, element_list%n_elements
       element = element_list%element(i_elem)
       LOOP_VERTEX : DO i=1,n_vertex_max
          inode1 = element%vertex(i)
          DO i_order = 1, n_order+1
             index_node1 = node_list%node(inode1)%index(i_order)
             call vertex_is_local(index_node1, is_local)
             IF (is_local) THEN	
                n_local_elms = n_local_elms + 1
                EXIT LOOP_VERTEX
             END IF
          END DO
       END DO LOOP_VERTEX
    END DO

    IF (ALLOCATED(local_elms)) call tr_deallocate(local_elms,"local_elms",CAT_FEM)
    ! Build local_elms from loc2glob
    call tr_allocate(local_elms,1,n_local_elms,"local_elms",CAT_FEM)

    n_local_elms = 0
    DO i_elem = 1, element_list%n_elements

       element = element_list%element(i_elem)
       L_I: DO i=1,n_vertex_max

          inode1	    = element%vertex(i)

          DO i_order = 1, n_order+1

             index_node1 = node_list%node(inode1)%index(i_order)
      
             call vertex_is_local(index_node1, is_local)
             IF (is_local) THEN	
                n_local_elms = n_local_elms + 1
                local_elms(n_local_elms) = i_elem
                exit L_I
             END IF
          END DO
       END DO L_I
    END DO
    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    CALL MPI_Reduce(tsecond, max_time, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce(tsecond, min_time, 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    if (my_id .eq. 0) then
       write(*,FMT_TIMING) my_id, '# Elapsed time local element list :',min_time
       write(*,FMT_TIMING) my_id, '# Elapsed time local element list :',max_time
    end if

    CALL MPI_Reduce(n_local_elms, sum_n_local_elms, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_N, ierr)
    CALL MPI_Reduce(n_local_elms, max_n_local_elms, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_N, ierr)
    CALL MPI_Reduce(n_local_elms, min_n_local_elms, 1, MPI_INTEGER, MPI_MIN, 0, MPI_COMM_N, ierr)

    if (my_id .eq. 0) then
       write(*,"(A70,I12)") ' maximum number of elements computed on one cpu ', max_n_local_elms
       write(*,"(A70,I12)") ' minimum number of elements computed on one cpu ', min_n_local_elms
       write(*,"(A70,I12)") ' number of elements computed over all cpus ', sum_n_local_elms
       write(*,"(A70,I12)") ' total number of elements ', element_list%n_elements
    end if

    index_total = -1
    do inode=1, node_list%n_nodes
       index_total = max(index_total,maxval(node_list%node(inode)%index))
    enddo

    ndof_glob  = index_total * n_tor * n_var
    
    node_list%n_dof = ndof_glob

    IF (MOD(n_local_elms, ((n_tor+1)/2)) .EQ. 0) THEN
       murge_elem_block_size = n_local_elms/((n_tor+1)/2)
    ELSE
       murge_elem_block_size = ((n_local_elms - MOD(n_local_elms,       &
            &                   ((n_tor+1)/2)))/((n_tor+1)/2))+1
    END IF

    IF (gmres) THEN
       murge_assembly_step = murge_elem_block_size*(n_tor+1)/2
       CALL MPI_AllReduce(murge_elem_block_size, murge_assembly_step, 1,   &
            &             MPI_INTEGER, MPI_SUM, MPI_COMM_TRANS, ierr)
    ELSE
       murge_assembly_step = murge_elem_block_size
    END IF


    call tr_debug_write("murge_setgraph end")
#else
    print *, "Binary built without murge"
    call abort()
#endif
  END SUBROUTINE murge_setGraph
end module murge_module
