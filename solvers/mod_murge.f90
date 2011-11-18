!> Variables and and routines related to the MURGE solver interface.
MODULE murge_module
  use tr_module 
  IMPLICIT NONE

  INCLUDE "murge.inc"
  !include "hips.inc"
  ! Indicate which solver is used
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_solver  

  ! Solver identification number
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_id      
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_id_prod
  INTEGER                                     :: murge_harmonic
  ! Local number of element
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_local_n 
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_local_n_prod 
  ! Global number of element
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_global_n 
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_global_n_prod 
  ! Number of dof by node
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_ndof

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
  SUBROUTINE murge_add_one_entry( index_node, k, in, index_node2, k2,&
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
       CALL MURGE_ASSEMBLYSETVALUE( murge_id, row_idx, col_idx, zbig,&
            &                       ierr )
       IF (ierr /= MURGE_SUCCESS) THEN
          WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
               "I", index_node, murge_ndof, k, in, &
               "J", index_node2, k2, in2
          STOP
       END IF
    END IF
    IF (gmres) THEN
       row_idx = n_tor*n_var * (index_node - 1) + (k -1)*n_tor + in
       col_idx = n_tor*n_var * (index_node2- 1) + (k2-1)*n_tor + in2 
       CALL MURGE_ASSEMBLYSETVALUE( murge_id_prod, row_idx, col_idx, zbig,&
            &                       ierr )

       IF (ierr /= MURGE_SUCCESS) THEN
          WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE(prod) ", &
               "I", index_node, n_tor, n_var, k, in, &
               "J", index_node2, k2, in2
          STOP
       END IF
    END IF
#else
       print *, "Binary built without murge"
       call abort()
#endif
  END SUBROUTINE murge_add_one_entry

  SUBROUTINE murge_initialization(gmres, my_id, mpi_comm_n, i_tor)
    USE parameters, ONLY : n_tor, n_var
    INCLUDE 'mpif.h'

    LOGICAL :: gmres
    INTEGER :: MPI_COMM_N
    INTEGER :: my_id
    INTEGER :: i_tor(:)

    INTEGER, EXTERNAL :: omp_get_num_threads
    INTEGER(KIND=MURGE_INTS_KIND) :: ierr
    INTEGER(KIND=MURGE_INTS_KIND) :: ndof
    
#ifdef USE_MURGE


    call tr_debug_write("murge_initialised begin")
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
             CALL MURGE_SetOptionINT(murge_id_prod, IPARM_VERBOSE,        &
                  &                  API_VERBOSE_YES, ierr)
             CALL MURGE_SetOptionINT(murge_id_prod,                       &
                  &                  IPARM_MATRIX_VERIFICATION,           &
                  &                  API_YES,         ierr)
             CALL MURGE_SetOptionINT(murge_id_prod, IPARM_THREAD_NBR,     &
                  &                  murge_nthrd,     ierr)
             CALL MURGE_SetOptionINT(murge_id_prod, MURGE_IPARAM_SYM,     &
                  &                  murge_sym,       ierr)
             CALL MURGE_SetOptionINT(murge_id_prod, MURGE_IPARAM_BASEVAL, &
                  &                  1,               ierr)
             ndof = n_tor*n_var
             CALL MURGE_SetOptionINT(murge_id_prod, MURGE_IPARAM_DOF,     &
                  &                  ndof,     ierr)  
             CALL MURGE_SetCommunicator(murge_id_prod, MPI_COMM_WORLD, ierr)

             CALL MURGE_SetDefaultOptions(murge_id,      0, ierr)
             CALL MURGE_SetOptionINT(murge_id, IPARM_MURGE_REFINMENT,     &
                  &                  API_NO,     ierr)  
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

          CALL MURGE_SetOptionINT(murge_id, IPARM_VERBOSE,             &
               &                  API_VERBOSE_YES, ierr)
          CALL MURGE_SetOptionINT(murge_id, IPARM_MATRIX_VERIFICATION, &
               &                  API_YES,         ierr)
          ! refinement : max number of iterations
          CALL MURGE_SetOptionINT(murge_id, IPARM_ITERMAX,             &
               &                  murge_iter,      ierr) 
          ! degrees of freedom per node (not correct)
          IF (gmres) THEN
             IF ( i_tor(my_id+1) == 1 ) THEN
                murge_ndof = n_var
             ELSE
                murge_ndof = 2*n_var
             END IF
             CALL MURGE_SetOptionINT(murge_id, IPARM_DOF_COST,         &
                  &                  2*n_var,      ierr) 
          ELSE
             murge_ndof = n_tor*n_var
          END IF
          CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_DOF,          &
               &                  murge_ndof,      ierr) 
          ! TODO : omp_num_thread
          CALL MURGE_SetOptionINT(murge_id, IPARM_THREAD_NBR,          &
               &                  murge_nthrd,     ierr)
          CALL MURGE_SetOptionINT(murge_id, IPARM_LEVEL_OF_FILL,       &
               &                  murge_iluk,      ierr)
          CALL MURGE_SetOptionINT(murge_id, IPARM_INCOMPLETE,          &
               &                  murge_ricar,     ierr)
          CALL MURGE_SetOptionINT(murge_id, IPARM_AMALGAMATION_LEVEL,  &
               &                  murge_amalg,     ierr)
          CALL MURGE_SetOptionINT(murge_id, IPARM_MATRIX_VERIFICATION, &
               &                  API_YES,         ierr)

          CALL MURGE_SetOptionREAL(murge_id, DPARM_EPSILON_MAGN_CTRL,  murge_pivot,     ierr)

       ENDIF


       CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_SYM,         &
            &                  murge_sym,   ierr)
       CALL MURGE_SetOptionINT(murge_id, MURGE_IPARAM_BASEVAL,     &
            &                  1,   ierr)

       CALL MURGE_SetOptionREAL(murge_id, MURGE_RPARAM_EPSILON_ERROR, &
            &                   murge_epsilon, ierr)

       murge_initialised = .TRUE.
    END IF
    call tr_debug_write("murge_initialised end")
#else
       print *, "Binary built without murge"
       call abort()
#endif
  END SUBROUTINE murge_initialization

  SUBROUTINE murge_setGraph(gmres, n, local_elms, n_local_elms, element_list, &
       &                    node_list, n_aa, my_id)

    use parameters, only : n_order, n_vertex_max
    use data_structure, only : type_element, type_element_list, type_node_list
    logical :: gmres
    integer :: n
    integer :: n_local_elms
    type (type_element_list) :: element_list
    type (type_node_list)    :: node_list
    integer :: n_aa
    integer :: my_id

    integer :: local_elms(:)

    integer(KIND=MURGE_INTS_KIND) :: ierr
    integer(KIND=MURGE_INTL_KIND) :: nnz
    integer :: t0, t1
    integer :: i_elem
    integer :: i, k, inode1, knode, i_order, k_order
    integer(KIND=MURGE_INTS_KIND) :: index_node1, index_node2
    type (type_element)      :: element
    integer :: nb_periods, nb_periodes_max, nb_periodes_sec
    character(len=20), parameter :: FMT_TIMING = "(I2,A70,F7.2)"
    
#ifdef USE_MURGE
    call tr_debug_write("murge_setgraph begin")
    call system_clock(count_rate=nb_periodes_sec, count_max=nb_periodes_max)
    
    murge_global_n      = n
    murge_global_n_prod = n
    nnz = n_local_elms*(n_order+1)*(n_order+1)*n_vertex_max*n_vertex_max 

    Call MURGE_GRAPHBEGIN(murge_id, murge_global_n, nnz, ierr)
    IF (ierr /= MURGE_SUCCESS) THEN
       write (*,*) "ERROR in MURGE_GRAPHBEGIN"
       STOP
    END IF
    IF (gmres) THEN
       CALL MURGE_GRAPHBEGIN(murge_id_prod, murge_global_n_prod, nnz, ierr)
       IF (ierr /= MURGE_SUCCESS) THEN
          write (*,*) "ERROR in MURGE_GRAPHBEGIN"
          STOP
       END IF
    END IF
    call system_clock(count=t0)

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


                   CALL MURGE_GRAPHEDGE(murge_id,  &
                        index_node1,         &
                        index_node2,         &
                        ierr)
                   IF (ierr /= MURGE_SUCCESS) THEN
                      write (*,*) "N", n, n_AA,&
                           "I", index_node1, &
                           "J", index_node2
                      STOP
                   END IF

                   CALL MURGE_GRAPHEDGE(murge_id_prod,  &
                        index_node1,         &
                        index_node2,         &
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


    CALL MURGE_GRAPHEND(murge_id_prod, ierr)
    IF (ierr /= MURGE_SUCCESS) THEN
       write (*,*) "ERROR in MURGE_GRAPHEND"
       STOP
    END IF

    call system_clock(count=t0)
    CALL MURGE_GETLOCALNODENBR(murge_id, murge_local_n, ierr)
    IF (ierr /= MURGE_SUCCESS) THEN
       write (*,*) "ERROR in MURGE_GETLOCALNODENBR"
       STOP
    END IF
    CALL MURGE_GETLOCALNODENBR(murge_id_prod, murge_local_n_prod, ierr)
    IF (ierr /= MURGE_SUCCESS) THEN
       write (*,*) "ERROR in MURGE_GETLOCALNODENBR"
       STOP
    END IF
    call system_clock(count=t1)
    nb_periods = t1-t0
    if (t1<t0) nb_periods = nb_periods + nb_periodes_max
    write(*,FMT_TIMING) my_id, ' system_clock elapsed time in MURGE_GETLOCALNODENBR ',REAL(nb_periods)/nb_periodes_sec


    call tr_allocate(murge_loc2glob,1,murge_local_n,"murge_loc2glob")
    call tr_allocate(murge_loc2glob_prod,1,murge_local_n_prod,"murge_loc2glob_prod")
    write (*,*) "Local number of nodes", murge_local_n, "global",  n
    call system_clock(count=t0)
    CALL MURGE_GETLOCALNODELIST(murge_id, murge_loc2glob, ierr)
    if (allocated(murge_glob2loc)) call tr_deallocate(murge_glob2loc,"murge_glob2loc")

    IF (ierr /= MURGE_SUCCESS) THEN
       write (*,*) "ERROR in MURGE_GETLOCALNODELIST"
       STOP
    END IF

    CALL MURGE_GETLOCALNODELIST(murge_id_prod, murge_loc2glob_prod, ierr)
    if (allocated(murge_glob2loc)) call tr_deallocate(murge_glob2loc,"murge_glob2loc")

    IF (ierr /= MURGE_SUCCESS) THEN
       write (*,*) "ERROR in MURGE_GETLOCALNODELIST"
       STOP
    END IF
    call system_clock(count=t1)
    nb_periods = t1-t0
    if (t1<t0) nb_periods = nb_periods + nb_periodes_max
    write(*,FMT_TIMING) my_id, ' system_clock elapsed time in MURGE_GETLOCALNODELIST ',REAL(nb_periods)/nb_periodes_sec

    !!murge_local_n_prod = murge_local_n/((n_tor+1)/2)
    !!if (my_id == 0) then
    !!   murge_local_n_prod = murge_local_n_prod + MOD(murge_local_n, (n_tor+1)/2)
    !!end if
    !!allocate (murge_loc2glob_prod(murge_local_n_prod))
    !!iter2=1
    !!do iter = my_id_trans+1, murge_local_n, (n_tor+1)/2
    !!   murge_loc2glob_prod(iter2) = murge_loc2glob(iter)
    !!   iter2 = iter2 + 1
    !!end do
    !!call MURGE_SETLOCALNODELIST(murge_id_prod, murge_local_n_prod, murge_loc2glob_prod, ierr)
    !!IF (ierr /= MURGE_SUCCESS) THEN
    !!   write (*,*) "ERROR in MURGE_SETLOCALNODELIST"
    !!   STOP
    !!END IF

    call tr_debug_write("murge_setgraph end")
#else
       print *, "Binary built without murge"
       call abort()
#endif
  END SUBROUTINE murge_setGraph
end module murge_module
