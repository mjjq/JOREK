#ifdef MURGE_USE_SEQUENCE_HARM
#  ifndef MURGE_USE_SEQUENCE
#    define MURGE_USE_SEQUENCE
#  endif
#endif
!> Variables and and routines related to the MURGE solver interface.
MODULE murge_module
  USE tr_module
  USE ISO_C_BINDING, ONLY : C_INT, C_PTR
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
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_sequence_id(2), murge_sequence_id_harm(2)
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
  INTEGER(KIND=MURGE_INTS_KIND), POINTER      :: MURGE_ROWS_HARM(:), MURGE_COLS_HARM(:)
  REAL(KIND=MURGE_COEF_KIND),    POINTER      :: MURGE_VALS(:), MURGE_VALS_HARM(:)
  INTEGER,                   ALLOCATABLE      :: murge_assembly_first_entry(:)
  INTEGER,                   ALLOCATABLE      :: murge_assembly_first_entry_harm(:)
  INTEGER                                     :: murge_assembly_step
  INTEGER                                     :: murge_elem_block_size
CONTAINS

#ifndef MURGE_INTERFACE_MAJOR_VERSION
#  define NEED_REDEFINE_ASSEMBLYBEGIN
#else
! MURGE_INTERFACE_MAJOR_VERSION == 0 does not exists
#  if MURGE_INTERFACE_MAJOR_VERSION == 1
#    if MURGE_INTERFACE_MINOR_VERSION < 1
#      define NEED_REDEFINE_ASSEMBLYBEGIN
#    endif
#  endif
#endif

#ifdef NEED_REDEFINE_ASSEMBLYBEGIN
  ! Backward compatibility
   SUBROUTINE MURGE_ASSEMBLYBEGIN(ID, N, COEFNBR, OP, OP2, MODE, SYM, IERROR)
     INTEGER(KIND=4),      INTENT(IN)  :: ID, OP, OP2, MODE, SYM, N
     INTEGER(KIND=4),      INTENT(IN)  :: COEFNBR
     INTEGER(KIND=4),      INTENT(OUT) :: IERROR
     MURGE_ASSEMBLYBEGIN(ID, COEFNBR, OP, OP2, MODE, SYM, IERROR)
   END SUBROUTINE MURGE_ASSEMBLYBEGIN
#endif
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
  !! @param cnt         Entry counter for precond murge problem.
  !! @param cnt_prod    Entry counter for product matrix.
  !! @param only_count  Indicate if we do a count or a real assembly.
  !!
  SUBROUTINE murge_add_one_entry( index_node,  k,  in,                 &
       &                          index_node2, k2, in2,                &
       &                          zbig, solve_only, gmres,             &
       &                          cnt, cnt_prod, only_count)
    USE parameters

    INTEGER, INTENT(IN)    :: index_node,  k,  in
    INTEGER, INTENT(IN)    :: index_node2, k2, in2
    REAL*8,  INTENT(IN)    :: zbig
    LOGICAL, INTENT(IN)    :: solve_only, gmres
    INTEGER, INTENT(INOUT) :: cnt, cnt_prod
    LOGICAL, INTENT(IN)    :: only_count

    INTEGER(kind=MURGE_INTS_KIND) :: ierr
    INTEGER(KIND=MURGE_INTS_KIND) :: row_idx, col_idx
    LOGICAL :: is_local
#ifdef USE_MURGE
    IF (.NOT. solve_only) THEN
       IF ( (in  .ge. murge_first_tor  .and. in  .le. murge_last_tor) .and. &
            (in2 .ge. murge_first_tor  .and. in2 .le. murge_last_tor) ) THEN
          cnt = cnt + 1
          IF (.not. only_count) THEN
             row_idx = murge_ndof * (index_node - 1) + (k -1)*murge_ntor + 1
             col_idx = murge_ndof * (index_node2- 1) + (k2-1)*murge_ntor + 1
             IF (in  /= 1) row_idx = row_idx + MOD(in,  2)
             IF (in2 /= 1) col_idx = col_idx + MOD(in2, 2)
             CALL MURGE_ASSEMBLYSETVALUE( murge_id, row_idx, col_idx, zbig,    &
                  &                       ierr )
             IF (ierr /= MURGE_SUCCESS) THEN
                WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", ierr, cnt,    &
                     "I", index_node, murge_ndof, k, in,                       &
                     "J", index_node2, k2, in2
                STOP
             END IF
          END IF
       END IF
    END IF
    IF (gmres) THEN
       row_idx = n_tor*n_var * (index_node - 1) + (k -1)*n_tor + in
       col_idx = n_tor*n_var * (index_node2- 1) + (k2-1)*n_tor + in2
       call vertex_is_local_prod(index_node2, is_local)
       IF (is_local) THEN
          cnt_prod = cnt_prod + 1
          IF (.not. only_count) THEN
             CALL MURGE_ASSEMBLYSETVALUE( murge_id_prod, row_idx, col_idx,  &
                  &                       zbig, ierr )

             IF (ierr /= MURGE_SUCCESS) THEN
                WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE(prod) ",      &
                     "I", index_node, n_tor, n_var, k, in,                  &
                     "J", index_node2, k2, in2, cnt_prod, ierr
                STOP
             END IF
          END IF
       END IF
    END IF
#else
    PRINT *, "Binary built without murge"
    CALL abort()
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
             ndof = murge_ndof_prod
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
  !        CALL MURGE_SetOptionINT(murge_id, IPARM_MURGE_REFINEMENT,            &
  !             &                  API_NO,      ierr)
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

  !> Subroutine murge_termination
  !! 
  !! Clean Murge solvers environnement.
  !!
  !! @param gmres Indicate if the GMREs solver id used.
  !!
  SUBROUTINE murge_termination(gmres)
    LOGICAL, INTENT(IN) :: gmres

    INTEGER ierr

#ifdef USE_MURGE
    IF ( use_murge_element ) THEN
       IF (ALLOCATED(murge_glob2loc))                                          &
            CALL tr_deallocate(murge_glob2loc,"murge_glob2loc",CAT_DMATRIX)
       IF (ALLOCATED(murge_loc2glob))                                          &
            CALL tr_deallocate(murge_loc2glob,"murge_loc2glob",CAT_DMATRIX)
       IF (ALLOCATED(murge_loc2glob_prod))                                     &
            CALL tr_deallocate(murge_loc2glob_prod,"murge_loc2glob_prod",      &
            &                  CAT_DMATRIX)
       IF (ALLOCATED(murge_assembly_first_entry))                              &
            CALL tr_deallocate(murge_assembly_first_entry,                     &
            &                  "murge_assembly_first_entry", CAT_DMATRIX)
    END IF
    CALL MURGE_Clean(murge_id, ierr)
    IF (gmres) THEN
       CALL MURGE_Clean(murge_id_prod, ierr)
    END IF
#endif

  END SUBROUTINE murge_termination

  !> Subroutine: murge_setGraph
  !!
  !! Gives the graph to Murge solver.
  !! Get the solver distribtion.
  !! Build new local element list.
  !! Construct the product matrix distribution.
  !!
  !! @param gmres          Indicate if the solver used is the GMRES.
  !! @param n              Global number of columns in the matrix.
  !! @param local_elms     List of local elements.
  !! @param n_local_elms   Number of local elements.
  !! @param element_list   List of all the elements.
  !! @param node_list
  !! @param n_aa
  !! @param my_id          Global MPI rank.
  !! @param my_id_trans    Rank in the trans-harmonic communicator.
  !! @param n_cpu_trans    Number of process in the trans-harmonic communicator.
  !! @param MPI_COMM_N     Harmonic communicator.
  !! @param MPI_COMM_TRANS Trans-harmonic communicator.
  !!
  SUBROUTINE murge_setGraph(gmres, n, local_elms, n_local_elms, element_list, &
       &                    node_list, n_aa, my_id, my_id_trans, n_cpu_trans, &
       &                    MPI_COMM_N, MPI_COMM_TRANS)

    USE parameters,                ONLY : n_order, n_vertex_max
    USE data_structure,            ONLY : type_element, type_element_list,     &
         &                                type_node_list, MURGE_UserData_t
    USE parameters,                ONLY : n_tor, n_var
    USE global_distributed_matrix, ONLY : ndof_glob
    USE clock_module,              ONLY : clcktype, clck_time, clck_ldiff,     &
         &                                FMT_TIMING
    USE mpi_mod,                   ONLY : MPI_REAL8, MPI_MAX, MPI_COMM_WORLD,  &
         &                                MPI_MIN, MPI_INTEGER, MPI_SUM

    LOGICAL,                  INTENT(IN)    :: gmres
    INTEGER,                  INTENT(IN)    :: n
    INTEGER, ALLOCATABLE,     INTENT(INOUT) :: local_elms(:)
    INTEGER,                  INTENT(INOUT) :: n_local_elms
    TYPE (type_element_list), INTENT(IN)    :: element_list
    TYPE (type_node_list),    INTENT(INOUT) :: node_list
    INTEGER,                  INTENT(IN)    :: n_aa
    INTEGER,                  INTENT(IN)    :: my_id
    INTEGER,                  INTENT(IN)    :: my_id_trans
    INTEGER,                  INTENT(IN)    :: n_cpu_trans
    INTEGER,                  INTENT(IN)    :: MPI_COMM_N
    INTEGER,                  INTENT(IN)    :: MPI_COMM_TRANS

    INTEGER(KIND=MURGE_INTS_KIND) :: murge_ierr
    INTEGER(KIND=MURGE_INTL_KIND) :: nnz
    INTEGER :: start, ierr, comm_rank, comm_size
    INTEGER :: i_elem, iter
    INTEGER :: i, j, k, inode1, knode, i_order, k_order
    INTEGER(KIND=MURGE_INTS_KIND) :: index_node1, index_node2
    TYPE (type_element)      :: element
    TYPE(clcktype)           :: t0, t1
    REAL*8                   :: max_time, min_time, tsecond
    LOGICAL                  :: is_local
    INTEGER                  :: sum_n_local_elms
    INTEGER                  :: max_n_local_elms, min_n_local_elms
    INTEGER                  :: index_total, inode
    INTEGER,                       ALLOCATABLE  :: tmp_glob2loc_prod(:)
#ifdef MURGE_PROD_NO_COMM
    INTEGER(KIND=MURGE_INTS_KIND), ALLOCATABLE  :: tmp_loc2glob_prod(:)
    INTEGER(KIND=MURGE_INTS_KIND)               :: tmp_local_n_prod, tmp
    INTEGER :: ELM_INDEX
#endif
#ifdef MURGE_USE_GETLOCALELEMENTLIST
    type(MURGE_UserData_t) :: d
    INTEGER(C_INT) :: c_ierr
    INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MURGE_DUPLICATE_ELEMENTS = 0
    INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MURGE_DISTRIBUTE_ELEMENTS = 1
#  ifdef MURGE_USE_DUPLICATE_ELEMENT
    INTEGER(C_INT) :: mode = MURGE_DUPLICATE_ELEMENTS
#  else
    INTEGER(C_INT) :: mode = MURGE_DISTRIBUTE_ELEMENTS
#  endif
#endif
#ifdef USE_MURGE
    INTERFACE
       INTEGER(C_INT) FUNCTION MURGE_GetLocalElementNbr(id, N, globalElementNbr,  &
            &                                           localElementNbr, mode, d) &
            BIND(C, name="MURGE_GetLocalElementNbr")
         USE ISO_C_BINDING, ONLY : C_INT, C_PTR
         USE data_structure, ONLY : MURGE_UserData_t

         INTEGER(C_INT),  VALUE :: id, N, globalElementNbr
         INTEGER(C_INT)         :: localElementNbr
         INTEGER(C_INT),  VALUE :: mode
         TYPE(MURGE_UserData_t) :: d
       END FUNCTION MURGE_GetLocalElementNbr
       
       INTEGER(C_INT) FUNCTION MURGE_GetLocalElementList(id, element_list) &
            BIND(C, name="MURGE_GetLocalElementList")
         USE ISO_C_BINDING, ONLY : C_INT
         INTEGER(C_INT), VALUE :: id
         INTEGER(C_INT)        :: element_list(*)
       END FUNCTION MURGE_GetLocalElementList
    END INTERFACE

    CALL tr_debug_write("murge_setgraph begin")


    murge_global_n      = n
    murge_global_n_prod = n
    nnz = n_local_elms*(n_order+1)*(n_order+1)*n_vertex_max*n_vertex_max
#ifdef MURGE_USE_GETLOCALELEMENTLIST
    d%nVertexMax = n_vertex_max*(n_order+1)
    c_ierr = MURGE_GetLocalElementNbr(murge_id, murge_global_n, &
         &                            element_list%n_elements,  &
         &                            n_local_elms, mode, d)

    if (n_local_elms .GT. element_list%n_elements) then
       print *, n_local_elms, "is greater than", element_list%n_elements
       call abort()
    end if
    IF (ALLOCATED(local_elms)) &
         CALL tr_deallocate(local_elms,"local_elms",CAT_FEM)
    ! Build local_elms from loc2glob
    CALL tr_allocate(local_elms,1,n_local_elms,"local_elms",CAT_FEM)

    c_ierr = MURGE_GetLocalElementList(murge_id, local_elms)

    CALL MURGE_GETLOCALNODENBR(murge_id, murge_local_n, murge_ierr)
    IF (murge_ierr /= MURGE_SUCCESS) THEN
       WRITE (*,*) "ERROR in MURGE_GETLOCALNODENBR"
       STOP
    END IF
    CALL tr_allocate( murge_loc2glob, 1, murge_local_n,           &
         &            "murge_loc2glob",CAT_DMATRIX)
    CALL MURGE_GETLOCALNODELIST(murge_id, murge_loc2glob, murge_ierr)
    IF (ALLOCATED(murge_glob2loc))                                             &
         CALL tr_deallocate(murge_glob2loc,"murge_glob2loc",CAT_DMATRIX)

    IF (murge_ierr /= MURGE_SUCCESS) THEN
       WRITE (*,*) "ERROR in MURGE_GETLOCALNODELIST"
       STOP
    END IF
#else

    ! Give the graph to MURGE
    CALL MURGE_GRAPHBEGIN(murge_id, murge_global_n, nnz, murge_ierr)
    IF (murge_ierr /= MURGE_SUCCESS) THEN
       WRITE (*,*) "ERROR in MURGE_GRAPHBEGIN"
       STOP
    END IF
    CALL clck_time(t0)

    CALL tr_debug_write("murge_setgraph loop start")
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
                        murge_ierr)
                   IF (murge_ierr /= MURGE_SUCCESS) THEN
                      WRITE (*,*) "N", n, n_AA,&
                           "I", index_node1, &
                           "J", index_node2
                      STOP
                   END IF

                END DO
             END DO
          END DO
       END DO
    END DO

    CALL tr_debug_write("murge_setgraph loop end")
    CALL clck_time(t1)
    CALL clck_ldiff(t0,t1,tsecond)
    CALL MPI_Reduce( tsecond, max_time, 1, MPI_REAL8, MPI_MAX, 0, &
         &           MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce( tsecond, min_time, 1, MPI_REAL8, MPI_MIN, 0, &
         &           MPI_COMM_WORLD, ierr)
    IF (my_id .EQ. 0) THEN
       WRITE(*,FMT_TIMING) my_id, '# Elapsed time entering graph :',min_time
       WRITE(*,FMT_TIMING) my_id, '# Elapsed time entering graph :',max_time
    END IF


    CALL clck_time(t0)
    CALL MURGE_GRAPHEND(murge_id, murge_ierr)
    IF (murge_ierr /= MURGE_SUCCESS) THEN
       WRITE (*,*) "ERROR in MURGE_GRAPHEND"
       STOP
    END IF
    CALL clck_time(t1)
    CALL clck_ldiff(t0,t1,tsecond)
    CALL MPI_Reduce( tsecond, max_time, 1, MPI_REAL8, MPI_MAX, 0, &
         &           MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce( tsecond, min_time, 1, MPI_REAL8, MPI_MIN, 0, &
         &           MPI_COMM_WORLD, ierr)
    IF (my_id .EQ. 0) THEN
       WRITE(*,FMT_TIMING) my_id, '# Elapsed time in MURGE_GRAPHEND :',min_time
       WRITE(*,FMT_TIMING) my_id, '# Elapsed time in MURGE_GRAPHEND :',max_time
    END IF

    ! Get the solver distribution.
    CALL clck_time(t0)
    CALL MURGE_GETLOCALNODENBR(murge_id, murge_local_n, murge_ierr)
    IF (murge_ierr /= MURGE_SUCCESS) THEN
       WRITE (*,*) "ERROR in MURGE_GETLOCALNODENBR"
       STOP
    END IF

    CALL clck_time(t1)
    CALL clck_ldiff(t0,t1,tsecond)
    CALL MPI_Reduce( tsecond, max_time, 1, MPI_REAL8, MPI_MAX, 0, &
         &           MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce( tsecond, min_time, 1, MPI_REAL8, MPI_MIN, 0, &
         &           MPI_COMM_WORLD, ierr)
    IF (my_id .EQ. 0) THEN
       WRITE(*,FMT_TIMING) &
            my_id, '# Elapsed time in MURGE_GETLOCALNODENBR :',min_time
       WRITE(*,FMT_TIMING) &
            my_id, '# Elapsed time in MURGE_GETLOCALNODENBR :',max_time
    END IF

    CALL clck_time(t0)
    CALL tr_allocate( murge_loc2glob, 1, murge_local_n,           &
         &            "murge_loc2glob",CAT_DMATRIX)

    CALL MURGE_GETLOCALNODELIST(murge_id, murge_loc2glob, murge_ierr)
    IF (ALLOCATED(murge_glob2loc))                                             &
         CALL tr_deallocate(murge_glob2loc,"murge_glob2loc",CAT_DMATRIX)

    IF (murge_ierr /= MURGE_SUCCESS) THEN
       WRITE (*,*) "ERROR in MURGE_GETLOCALNODELIST"
       STOP
    END IF

    CALL clck_time(t1)
    CALL clck_ldiff(t0,t1,tsecond)
    CALL MPI_Reduce( tsecond, max_time, 1, MPI_REAL8, MPI_MAX, 0, &
         &           MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce( tsecond, min_time, 1, MPI_REAL8, MPI_MIN, 0, &
         &           MPI_COMM_WORLD, ierr)
    IF (my_id .EQ. 0) THEN
       WRITE(*,FMT_TIMING) &
            my_id, '# Elapsed time in MURGE_GETLOCALNODELIST :',min_time
       WRITE(*,FMT_TIMING) &
            my_id, '# Elapsed time in MURGE_GETLOCALNODELIST :',max_time
    END IF


    CALL clck_time(t0)

#ifdef MURGE_USE_DUPLICATE_ELEMENT
    ! Build local_elms from loc2glob
    n_local_elms = 0

    DO i_elem = 1, element_list%n_elements
       element = element_list%element(i_elem)
       LOOP_VERTEX: DO i=1,n_vertex_max
          inode1 = element%vertex(i)
          DO i_order = 1, n_order+1
             index_node1 = node_list%node(inode1)%index(i_order)
             CALL vertex_is_local(index_node1, is_local)
             IF (is_local) THEN
                n_local_elms = n_local_elms + 1
                EXIT LOOP_VERTEX
             END IF
          END DO
       END DO LOOP_VERTEX
    END DO

    IF (ALLOCATED(local_elms)) &
         CALL tr_deallocate(local_elms,"local_elms",CAT_FEM)
    ! Build local_elms from loc2glob
    CALL tr_allocate(local_elms,1,n_local_elms,"local_elms",CAT_FEM)

    n_local_elms = 0
    DO i_elem = 1, element_list%n_elements

       element = element_list%element(i_elem)
       L_I: DO i=1,n_vertex_max

          inode1 = element%vertex(i)

          DO i_order = 1, n_order+1

             index_node1 = node_list%node(inode1)%index(i_order)

             CALL vertex_is_local(index_node1, is_local)
             IF (is_local) THEN
                n_local_elms = n_local_elms + 1
                local_elms(n_local_elms) = i_elem
                EXIT L_I
             END IF
          END DO
       END DO L_I
    END DO
#endif 

    CALL clck_time(t1)
    CALL clck_ldiff(t0,t1,tsecond)
    CALL MPI_Reduce( tsecond, max_time, 1, MPI_REAL8, MPI_MAX, 0, &
         &           MPI_COMM_WORLD, ierr)
    CALL MPI_Reduce( tsecond, min_time, 1, MPI_REAL8, MPI_MIN, 0, &
         &           MPI_COMM_WORLD, ierr)
    IF (my_id .EQ. 0) THEN
       WRITE(*,FMT_TIMING) my_id, '# Elapsed time local element list :',min_time
       WRITE(*,FMT_TIMING) my_id, '# Elapsed time local element list :',max_time
    END IF
#endif

    CALL MPI_Reduce( n_local_elms, sum_n_local_elms, 1, MPI_INTEGER, MPI_SUM, &
         &           0, MPI_COMM_TRANS, ierr)
    CALL MPI_Reduce( n_local_elms, max_n_local_elms, 1, MPI_INTEGER, MPI_MAX, &
         &           0, MPI_COMM_TRANS, ierr)
    CALL MPI_Reduce( n_local_elms, min_n_local_elms, 1, MPI_INTEGER, MPI_MIN, &
         &           0, MPI_COMM_TRANS, ierr)

    IF (my_id .EQ. 0) THEN
       WRITE(*,"(A70,I12)")                                     &
            ' maximum number of elements computed on one cpu ', &
            max_n_local_elms/((n_tor+1)/2) + 1
       WRITE(*,"(A70,I12)")                                     &
            ' minimum number of elements computed on one cpu ', &
            min_n_local_elms/((n_tor+1)/2)
       WRITE(*,"(A70,I12)")                                     &
            ' number of elements computed over all cpus ', sum_n_local_elms
       WRITE(*,"(A70,I12)")                                     &
            ' total number of elements ', element_list%n_elements
    END IF

    !murge_elem_block_size = 1*murge_nthrd
    IF (MOD(n_local_elms, ((n_tor+1)/2)) .EQ. 0) THEN
       murge_elem_block_size = n_local_elms/((n_tor+1)/2)
    ELSE
       murge_elem_block_size = (( n_local_elms - &
            MOD( n_local_elms, (n_tor+1)/2 ) )/((n_tor+1)/2))+1
    END IF

    IF (gmres) THEN
       CALL MPI_AllReduce(murge_elem_block_size, murge_assembly_step, &
            1, MPI_INTEGER, MPI_SUM, MPI_COMM_TRANS, ierr)
    ELSE
       murge_assembly_step = murge_elem_block_size
    END IF

#  ifdef MURGE_PROD_NO_COMM
    IF (gmres) THEN
       ! Counting maximum number of entries in tmp_loc2glob_prod
       tmp_local_n_prod = 0
       DO i_elem = 1, n_local_elms, murge_assembly_step
          DO ELM_INDEX = 1, murge_elem_block_size
             IF ( i_elem + ELM_INDEX-1 +             &
                  murge_elem_block_size*my_id_trans > n_local_elms) EXIT
             element = element_list%element(local_elms(i_elem + &
                  ELM_INDEX-1 + murge_elem_block_size*my_id_trans))
             DO i=1,n_vertex_max
                inode1 = element%vertex(i)
                DO i_order = 1, n_order+1
                   index_node1 = node_list%node(inode1)%INDEX(i_order)
                   CALL vertex_is_local(index_node1, is_local)
                   IF (is_local) THEN
                      tmp_local_n_prod = tmp_local_n_prod + 1
                   END IF
                END DO
             END DO
          END DO
       END DO

       CALL tr_allocate( tmp_loc2glob_prod, 1, tmp_local_n_prod,              &
            &            "tmp_loc2glob_prod",CAT_DMATRIX)

       ! Building tmp_loc2glob_prod
       murge_local_n_prod = 0
       DO i_elem = 1, n_local_elms, murge_assembly_step
          DO ELM_INDEX = 1, murge_elem_block_size
             IF ( i_elem + ELM_INDEX-1 +             &
                  murge_elem_block_size*my_id_trans > n_local_elms) EXIT
             element = element_list%element(local_elms(i_elem + &
                  ELM_INDEX-1 + murge_elem_block_size*my_id_trans))
             DO i=1,n_vertex_max
                inode1 = element%vertex(i)
                LOOP_ORDER: DO i_order = 1, n_order+1
                   index_node1 = node_list%node(inode1)%INDEX(i_order)
                   CALL vertex_is_local(index_node1, is_local)
                   IF (is_local) THEN
                      DO j = 1, murge_local_n_prod
                         IF (index_node1 .EQ. tmp_loc2glob_prod(j)) CYCLE LOOP_ORDER
                         IF (index_node1 .LT. tmp_loc2glob_prod(j)) THEN
                            tmp = tmp_loc2glob_prod(j)
                            tmp_loc2glob_prod(j) = index_node1
                            index_node1 = tmp
                         END IF
                      END DO
                      murge_local_n_prod = murge_local_n_prod + 1
                      tmp_loc2glob_prod(murge_local_n_prod) = index_node1

                   END IF
                END DO LOOP_ORDER
             END DO
          END DO
       END DO

       CALL tr_allocate( murge_loc2glob_prod, 1, murge_local_n_prod,           &
            &            "murge_loc2glob_prod",CAT_DMATRIX)

       DO i = 1, murge_local_n_prod
          murge_loc2glob_prod(i) = tmp_loc2glob_prod(i)
       END DO

       CALL tr_deallocate( tmp_loc2glob_prod, "tmp_loc2glob_prod", CAT_DMATRIX)
    ELSE
#  endif

       murge_local_n_prod =                                                    &
            (murge_local_n-MOD(murge_local_n, n_cpu_trans))/n_cpu_trans
       if (my_id_trans .lt. MOD(murge_local_n, n_cpu_trans)) then
          murge_local_n_prod = murge_local_n_prod + 1
       end if

       CALL tr_allocate( murge_loc2glob_prod, 1, murge_local_n_prod,           &
            &            "murge_loc2glob_prod",CAT_DMATRIX)
       start = murge_local_n_prod*my_id_trans

       start = start + MIN(my_id_trans, MOD(murge_local_n,n_cpu_trans))

       DO i = 1, murge_local_n_prod
          murge_loc2glob_prod(i) = murge_loc2glob(start + i)
       END DO
#  ifdef MURGE_PROD_NO_COMM
    ENDIF
#  endif
    murge_need_rebuild_sequence(1) = .TRUE.
!    CALL MURGE_PRODUCTSETLOCALNODENBR(murge_id_prod, murge_local_n_prod,      &
!         murge_ierr)
    IF (murge_ierr /= MURGE_SUCCESS) THEN
       WRITE (*,*) "ERROR in MURGE_PRODUCTSETLOCALNODENBR"
       STOP
    END IF

#  ifdef MURGE_PROD_NO_COMM
    CALL MURGE_PRODUCTSETGLOBALNODENBR(murge_id_prod, murge_global_n_prod, &
         murge_ierr)
    IF (murge_ierr /= MURGE_SUCCESS) THEN
       WRITE (*,*) "ERROR in MURGE_PRODUCTSETGLOBALNODENBR"
       STOP
    END IF
#  endif


!    CALL MURGE_PRODUCTSETLOCALNODELIST(murge_id_prod, murge_loc2glob_prod, &
!         murge_ierr)
    IF (murge_ierr /= MURGE_SUCCESS) THEN
       WRITE (*,*) "ERROR in MURGE_PRODUCTSETLOCALNODELIST"
       STOP
    END IF

    IF (ALLOCATED(murge_glob2loc_prod))                                        &
         CALL tr_deallocate( murge_glob2loc_prod,"murge_glob2loc_prod", &
         &                   CAT_DMATRIX)

    CALL tr_allocate(murge_glob2loc_prod, 1, murge_global_n_prod, &
         "murge_glob2loc_prod", CAT_DMATRIX)
    CALL tr_allocate(tmp_glob2loc_prod, 1, murge_global_n_prod, &
         "tmp_glob2loc_prod", CAT_DMATRIX)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, comm_rank, ierr)
    murge_glob2loc_prod = comm_size + 1
    DO iter = 1, murge_local_n_prod
       murge_glob2loc_prod(murge_loc2glob_prod(iter)) = comm_rank
    END DO
    call MPI_Allreduce(murge_glob2loc_prod, tmp_glob2loc_prod,      &
         murge_global_n_prod, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, &
         ierr)
    index_total = 0
    DO iter = 1, murge_global_n_prod
       if (tmp_glob2loc_prod(iter) .ne. comm_rank) then
          murge_glob2loc_prod(iter) = -1 - tmp_glob2loc_prod(iter)
       else
          index_total = index_total + 1
       end if
    END DO

    DO iter = 1, murge_local_n_prod
       if (tmp_glob2loc_prod(murge_loc2glob_prod(iter)) == comm_rank) &
            murge_glob2loc_prod(murge_loc2glob_prod(iter)) = iter
    END DO
    CALL tr_deallocate( tmp_glob2loc_prod,"tmp_glob2loc_prod", &
         CAT_DMATRIX)

    index_total = -1
    DO inode=1, node_list%n_nodes
       index_total = MAX(index_total,MAXVAL(node_list%node(inode)%index))
    ENDDO

    ndof_glob  = index_total * n_tor * n_var

    node_list%n_dof = ndof_glob


    CALL tr_debug_write("murge_setgraph end")
#else
    PRINT *, "Binary built without murge"
    CALL abort()
#endif
  END SUBROUTINE murge_setGraph

  INTEGER(C_INT) FUNCTION getVertices(e, idx) BIND(C, name="getVertices")
    USE parameters,                ONLY : n_order, n_vertex_max
    USE data_structure,            ONLY : type_element, type_element_list
    USE data_structure,            ONLY : type_node_list
    USE nodes_elements,            ONLY : element_list, node_list
    INTEGER(C_INT), VALUE :: e
    INTEGER(C_INT)        :: idx(*)
    INTEGER :: i, j, i_order, inode1
    TYPE(type_element) :: element
    j = 1
    element = element_list%element(e)

    DO i = 1, n_vertex_max
       inode1         = element%vertex(i)
       DO i_order = 1, n_order+1
          idx(j) = node_list%node(inode1)%index(i_order)
          j = j + 1
       END DO
    END DO
    getVertices = 0
    RETURN 
  END FUNCTION getVertices
  
END MODULE murge_module
