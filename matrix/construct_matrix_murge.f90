!*******************************************************************************
!* Subroutine: construct_matrix_murge                                          *
!*******************************************************************************
!*                                                                             *
!* Construct the matrix for the resolution using Murge interface (PaStiX).     *
!*                                                                             *
!* Parameters:                                                                 *
!*   my_id            - Identifier of the node in MPI_COMM_WORLD               *
!*   node_list        - List of nodes                                          *
!*   element_list     - List of all elements                                   *
!*   local_elms       - List of local elements                                 *
!*   n_local_elms     - Number of local elements                               *
!*   xpoint2          - ???                                                    *
!*   xcase2           - ???                                                    *
!*   psi_axis         - ???                                                    *
!*   psi_bnd          - ???                                                    *
!*   Z_xpoint         - ???                                                    *
!*   gmres            - Solve method (.true. for gmres, .false for 'direct')   *
!*   i_tor            - Tor number                                             *
!*   n_cpu            - Number of cpus                                         *
!*   mpi_comm_n       - Solver MPI communicator                                *
!*   MPI_COMM_TRANS   - Transversal communicator                               *
!*   my_id_trans      - ID in transversal communicator                         *
!*   n_cpu_trans      - Size of transversal communicator                       *
!*                                                                             *
!* Authors:                                                                    *
!*   Xavier Lacoste - xavier.lacoste@inria.fr                                  *
!*                                                                             *
!*******************************************************************************
MODULE THREAD_DATA
  USE tr_module
  USE data_structure,  ONLY : type_element_list,                               &
       &                      type_node_list
  USE mod_elt_matrix_fft
  USE mod_elt_matrix

  TYPE THREAD_DATA_TYPE
     ! -- local variables   --
     INTEGER                            :: thread_num
     LOGICAL                            :: ok
!#define PLENTY_TIMERS
#ifdef PLENTY_TIMERS
     INTEGER                            :: nb_periods_ass
     INTEGER                            :: nb_periods_set
     INTEGER                            :: nb_periods_elem_mat
     INTEGER                            :: nb_periods_comm
     INTEGER                            :: nb_periods_max
#endif
     INTEGER                            :: ndof_glob
     ! Different memory area on each thread
     REAL*8,                    POINTER :: ELM(:,:)
     REAL*8,                    POINTER :: RHS(:)

     ! --  global variables --
     ! read only, no mutex
     INTEGER,                   POINTER :: thread_nbr
     INTEGER,                   POINTER :: my_id
     INTEGER,                   POINTER :: my_id_trans
     INTEGER,                   POINTER :: MPI_COMM_TRANS
     TYPE (type_node_list),     POINTER :: node_list
     TYPE (type_element_list),  POINTER :: element_list
     LOGICAL,                   POINTER :: gmres, solve_only
     LOGICAL,                   POINTER :: xpoint2
     REAL*8,                    POINTER :: minRad
     REAL*8,                    POINTER :: R_axis
     REAL*8,                    POINTER :: Z_axis
     REAL*8,                    POINTER :: psi_axis
     REAL*8,                    POINTER :: psi_bnd
     REAL*8,                    POINTER :: Z_xpoint(:)
     REAL*8,                    POINTER :: psi_xpoint(:)
     INTEGER,                   POINTER :: xcase2
     INTEGER,                   POINTER :: local_elms(:)
     INTEGER,                   POINTER :: n_local_elms
     INTEGER,                   POINTER :: step
     INTEGER,                   POINTER :: elem_block_size
     INTEGER,                   POINTER :: harm_size, my_harm_size
     ! each thread access different part, no mutex
     REAL*8,                    POINTER :: PROD_MATRICES(:,:,:)
     INTEGER,                   POINTER :: PROD_COLROW(:,:,:)
     REAL*8,                    POINTER :: SEND_MATRICES(:,:,:,:)
     REAL*8,                    POINTER :: RECV_MATRICES(:,:,:,:)
     INTEGER,                   POINTER :: RECV_COLROW(:, :, :,:)
     INTEGER,                   POINTER :: matrix_nbr(:)
     INTEGER,                   POINTER :: matrix_nbr_rcv(:,:)
     REAL*8,                    POINTER :: rhs_loc_thread(:,:)
     ! accessed on thread 1
     REAL*8,                    POINTER :: rhs_loc(:)
  END TYPE THREAD_DATA_TYPE
CONTAINS

  INTEGER FUNCTION LOOP(DATA)
    USE data_structure,  ONLY : type_node, type_element
    USE parameters,      ONLY : n_vertex_max , n_var, n_order, n_tor
    USE murge_module,    ONLY : MURGE_SUCCESS, MURGE_ASSEMBLYBEGIN,            &
         &                      MURGE_ASSEMBLYSETNODEVALUES,                   &
         &                      MURGE_ASSEMBLYEND, murge_id,                   &
         &                      murge_id_prod, MURGE_COEF_KIND,                &
         &                      murge_ndof_prod
    USE mpi_mod

    IMPLICIT NONE

    ! Function parameters
    TYPE(THREAD_DATA_TYPE),  INTENT(INOUT) :: DATA

    ! internal variables
    INTEGER                        :: i_father, inode_father
    INTEGER                        :: ielm
    INTEGER                        :: inode
    INTEGER                        :: iv, iv2
    INTEGER                        :: index_node1
    INTEGER                        :: index_node2
    INTEGER                        :: k_order
    INTEGER                        :: i_order
    INTEGER                        :: index_rhs
    INTEGER                        :: index_ij
    INTEGER                        :: index_kl
    INTEGER                        :: index_large_i
    INTEGER                        :: index_large_k
    INTEGER                        :: i, j, k, l
    INTEGER                        :: index_send_mtx
    INTEGER                        :: new_row_mat_elem
    INTEGER                        :: new_col_mat_elem
    INTEGER                        :: row_harm
    INTEGER                        :: col_harm
    INTEGER                        :: cnt, cnt2, index_mtx
    INTEGER                        :: knode, inode1, inode2
    REAL(KIND=MURGE_COEF_KIND)     :: coefmtx(n_tor*n_var*n_tor*n_var)
    REAL(KIND=MURGE_COEF_KIND)     :: coefmtx_prod(n_tor*n_var*n_tor*n_var)
    INTEGER                        :: ierr
    TYPE (type_element)            :: element
    TYPE (type_element)            :: element_father
    TYPE (type_node)               :: nodes(n_vertex_max)
    LOGICAL                        :: is_local
    TYPE (type_node)               :: nodes_father(n_vertex_max)
    REAL*8,  POINTER               :: ELM(:,:)
    REAL*8,  POINTER               :: RHS(:)
    INTEGER                        :: next_matrix
    INTEGER                        :: ret
    INTEGER                        :: elem_size
    INTEGER                        :: ife
    INTEGER                        :: prod_size, mat_size
    INTEGER                        :: ELM_INDEX
    INTEGER                        :: node, iter, total, t0, t1, tt0, tt1
    INTEGER                        :: thread, iter_dof_row, iter_dof_col
    INTEGER                        :: row, col
    REAL(KIND=MURGE_COEF_KIND)     :: val
    INTEGER,               POINTER :: matrix_nbr_rcv(:,:), pt_matrix_nbr

    elem_size = n_tor*n_vertex_max*(n_order+1)*n_var
    cnt = 0
    cnt2 = 0
    pt_matrix_nbr => data%matrix_nbr(data%thread_num)
    data%ok = .TRUE.
    matrix_nbr_rcv => data%matrix_nbr_rcv

    ELM => data%ELM
    RHS => data%RHS

    !TODO: anticiper l'allocation ou l'allouer une fois pour toute
    ELEM : DO ife =1, data%n_local_elms, data%step
       pt_matrix_nbr = 0
       !!  elem_block_size elements are computed each loop iteration,
       !! each thread is in charge of a part of theses elements
#ifdef PLENTY_TIMERS
       CALL SYSTEM_CLOCK(count=t0)
#endif
       BLOCK_ELEM : DO ELM_INDEX = data%thread_num, data%elem_block_size,      &
            &                      data%thread_nbr

          !! If we have not treated all the local elements
          !! (elem_block_size may not divide n_local_elms)
          IF ( ife + ELM_INDEX-1 +                                             &
               data%elem_block_size*data%my_id_trans <=                        &
               data%n_local_elms) THEN
             ELM = 0.0
             RHS = 0.0

             ielm = data%local_elms( ife + ELM_INDEX - 1+                      &
                  data%elem_block_size*data%my_id_trans)

             element = data%element_list%element(ielm)

             i_father= data%element_list%element(ielm)%father

             IF ( i_father.NE.0) THEN
                element_father = data%element_list%element(i_father)
             ENDIF

             DO iv = 1, n_vertex_max

                IF ( i_father.NE.0) THEN
                   inode_father=element_father%vertex(iv)
                   nodes_father(iv) = data%node_list%node(inode_father)
                ENDIF
                inode     = element%vertex(iv)
                nodes(iv) = data%node_list%node(inode)

             ENDDO

             !! Compute the element matrix
#ifdef PLENTY_TIMERS
             CALL SYSTEM_CLOCK(count=tt0)
#endif
             IF (n_tor .GT. 3) THEN
                ! use fft for toroidal integration
                CALL element_matrix_fft(element,nodes, data%xpoint2,            &
                     &                  data%xcase2, data%minRad,               &
                     &                  data%R_axis, data%Z_axis, data%psi_axis,&
                     &                  data%psi_bnd, Data%z_xpoint,            &
                     &                  ELM, RHS, data%thread_num)
             ELSE
                ! use direct integration
                CALL element_matrix(element,nodes, data%xpoint2,               &
                     &              data%xcase2, data%minRad,                  &
                     &              data%R_axis, data%Z_axis, data%psi_axis,   &
                     &              data%psi_bnd, Data%z_xpoint,               &
                     &              ELM, RHS, data%thread_num)
                DO iv = 1, n_vertex_max ! boundary integrals

                   iv2  = MOD(iv, n_vertex_max) + 1

                   inode1 = element%vertex(iv)
                   inode2 = element%vertex(iv2)

                ENDDO
             ENDIF
#ifdef PLENTY_TIMERS
             CALL SYSTEM_CLOCK(count=tt1)
             data%nb_periods_elem_mat = data%nb_periods_elem_mat + tt1-tt0
             IF (tt1<tt0)                                                      &
                  data%nb_periods_elem_mat = data%nb_periods_elem_mat +        &
                  &                          data%nb_periods_max
#endif
             IF (element%n_sons .EQ. 0) THEN
                VERTEX_COL : DO i=1,n_vertex_max

                   !inode1         =node_out(i)! element%vertex(i)
                   inode1 =  element%vertex(i)
                   ORDER_COL : DO i_order = 1, n_order+1

                      ! index_node1 is the column
                      index_node1 = data%node_list%node(inode1)%index(i_order)

                      index_large_i = n_tor * n_var * (index_node1 - 1)

                      CALL vertex_is_local(index_node1, is_local)
                      IF (is_local) THEN

                         ! Set RHS member
                         DO j = 1, n_var * n_tor

                            ! index in the ELM matrix
                            index_ij = n_tor * n_var * (n_order+1) * (i-1) +   &
                                 & n_tor * n_var * (i_order-1) + j
                            ! index in global matrix
                            index_rhs = index_large_i+j
                            data%rhs_loc_thread(index_rhs,data%thread_num) =   &
                                 data%rhs_loc_thread(index_rhs,                &
                                 &                   data%thread_num) +        &
                                 &                   RHS(index_ij)
                         END DO

                         ! Build nodes Matrices
                         VERTEX_ROW : DO k=1,n_vertex_max

                            knode         = element%vertex(k)

                            ORDER_ROW : DO k_order = 1, n_order+1

                               index_node2 =                                   &
                                    data%node_list%node(knode)%index(k_order)

                               index_large_k = n_tor * n_var * (index_node2 - 1)

                               !coefmtx = 0
                               pt_matrix_nbr = pt_matrix_nbr + 1
                               DOF_COL : DO j = 1, n_var * n_tor
                                  ! Row index in the ELM matrix
                                  index_ij = n_tor * n_var *                   &
                                       (n_order+1) * (i-1) +                   &
                                       n_tor * n_var * (i_order-1) + j

                                  DOF_ROW : DO l = 1, n_var * n_tor
                                     ! BUILD node Matrix
                                     index_kl  = n_tor * n_var * (n_order      &
                                          &+1) * (k-1) + n_tor * n_var *       &
                                          & (k_order-1) + l

                                     IF ( data%gmres .AND.                     &
                                          .NOT. data%solve_only) THEN
                                        col_harm =                             &
                                             INT((MOD(j-1, n_tor)+1)/2) + 1
                                        row_harm =                             &
                                             INT((MOD(l-1, n_tor)+1)/2) + 1
                                        IF (col_harm == row_harm) THEN

                                           IF (col_harm == 1) THEN
                                              new_col_mat_elem =               &
                                                   INT((j-1)/n_tor)
                                              new_row_mat_elem =               &
                                                   INT((l-1)/n_tor)
                                              index_send_mtx =                 &
                                                   n_var*new_col_mat_elem +    &
                                                   new_row_mat_elem+1

                                           ELSE
                                              ! (num_var-1)*2 + &
                                              !   num_tor_local (0 or 1)
                                              new_col_mat_elem =               &
                                                   INT((j-1)/n_tor)*2 +        &
                                                   MOD(MOD(j-1, n_tor)+1,2)
                                              new_row_mat_elem =               &
                                                   INT((l-1)/n_tor)*2 +        &
                                                   &MOD(MOD(l-1, n_tor)+1,2)
                                              index_send_mtx =                 &
                                                   n_var*2*new_col_mat_elem +  &
                                                   new_row_mat_elem+1

                                           END IF
                                           data%SEND_MATRICES(index_send_mtx,  &
                                                &             pt_matrix_nbr,   &
                                                &             data%thread_num, &
                                                &             col_harm) =      &
                                                &  ELM(index_kl,index_ij)
                                        END IF
                                     END IF
                                     new_col_mat_elem = j-1
                                     new_row_mat_elem = l-1
                                     index_mtx = new_row_mat_elem+1            &
                                          & +new_col_mat_elem                  &
                                          & *(n_tor*n_var)
                                     data%PROD_MATRICES(index_mtx,             &
                                          &             pt_matrix_nbr,         &
                                          &             data%thread_num) =     &
                                          & ELM(index_kl,index_ij)
                                     !coefmtx(index_mtx) =&
                                     !     & ELM(index_kl&
                                     !     &,index_ij)

                                  END DO DOF_ROW
                               END DO DOF_COL
                               data%PROD_COLROW(1, pt_matrix_nbr,              &
                                    &           data%thread_num) = index_node1
                               data%PROD_COLROW(2, pt_matrix_nbr,              &
                                    &           data%thread_num) = index_node2
                            END DO ORDER_ROW
                         END DO VERTEX_ROW
                      END IF
                   END DO ORDER_COL
                END DO VERTEX_COL
             END IF

          END IF
       END DO BLOCK_ELEM
#ifdef PLENTY_TIMERS
       CALL SYSTEM_CLOCK(count=t1)
       data%nb_periods_ass = data%nb_periods_ass + t1-t0
       IF (t1<t0) &
            data%nb_periods_ass = &data%nb_periods_ass + data%nb_periods_max
       !$omp critical
       PRINT *, "TEMPS ASS", data%my_id, data%thread_num, t1-t0
       !$omp end critical
#endif

!$OMP barrier
#ifdef PLENTY_TIMERS
       CALL SYSTEM_CLOCK(count=t0)
#endif
       IF (.NOT. data%gmres) THEN
          ! We work on the full problem with direct method
          IF (data%thread_num == 1) THEN
             DO thread = 1, data%thread_nbr
                DO j = 1, data%matrix_nbr(thread)

                   index_node1 = data%PROD_COLROW(1, j, thread)
                   index_node2 = data%PROD_COLROW(2, j, thread)
                   DO iter = 1, (n_tor*n_var)**2
                      coefmtx(iter)= data%PROD_MATRICES(iter,j,thread)
                   END DO
#ifdef USE_MURGE
                   CALL MURGE_ASSEMBLYSETNODEVALUES(murge_id,                  &
                        &                           index_node2,               &
                        &                           index_node1,               &
                        &                           coefmtx, ierr)
#else
                   PRINT *, "Binary built without murge"
                   data%ok = .FALSE.
                   RETURN
#endif

                   IF (ierr /= MURGE_SUCCESS) THEN
                      WRITE (*,*) data%my_id, ":::",                           &
                           "I", index_node2,                                   &
                           "J", index_node1, cnt
                      data%ok = .FALSE.
                      RETURN
                   END IF
                END DO
             END DO
          END IF
       ELSE
          ! We have to construct a problem for the product
          IF (data%thread_num == 1) THEN
             DO thread = 1, data%thread_nbr
                DO j = 1, data%matrix_nbr(thread)
                   index_node1 = data%PROD_COLROW(1, j, thread)
                   index_node2 = data%PROD_COLROW(2, j, thread)
#ifdef USE_MURGE
!#define MURGE_PROD_NODE
#ifdef MURGE_PROD_NODE
                   DO iter = 1, (n_tor*n_var)**2
                      coefmtx(iter)= data%PROD_MATRICES(iter,j, thread)
                   END DO
                   CALL MURGE_ASSEMBLYSETNODEVALUES(murge_id_prod,             &
                        &                           index_node2,               &
                        &                           index_node1,               &
                        &                           coefmtx, ierr)
                   IF (ierr /= MURGE_SUCCESS) THEN
                      WRITE (*,*) data%my_id, ":::",                           &
                           "I", index_node2,                                   &
                           "J", index_node1, cnt
                      data%ok = .FALSE.
                      RETURN
                   END IF
#else
                   DO iter_dof_row = 1, murge_ndof_prod
                      DO iter_dof_col = 1, murge_ndof_prod
                         row = (index_node2-1)*murge_ndof_prod+iter_dof_row
                         col = (index_node1-1)*murge_ndof_prod+iter_dof_col
                         val = data%PROD_MATRICES(                             &
                              &  (iter_dof_col-1)*murge_ndof_prod+iter_dof_row,&
                              &  j, thread)

                         CALL MURGE_ASSEMBLYSETVALUE(murge_id_prod,            &
                              &                      row, col, val, ierr)
                         IF (ierr /= MURGE_SUCCESS) THEN
                            WRITE (*,*) data%my_id, ":::::",                   &
                                 "I", row,                                     &
                                 "J", col, ierr
                            data%ok = .FALSE.
                            RETURN
                         END IF

                      END DO
                   END DO
#endif
#else
                   PRINT *, "Binary built without murge"
                   data%ok = .FALSE.
                   RETURN
#endif
                END DO
             END DO
          END IF
          ! Just so that the first thread doesn't do all the job.
          IF ((data%thread_nbr == 1 .OR. data%thread_num == 2 ) .AND.          &
               & .NOT. DATA%solve_only) THEN

             matrix_nbr_rcv = 0
#ifdef PLENTY_TIMERS
             CALL SYSTEM_CLOCK(count=tt0)
#endif
             CALL MPI_Allgather(data%matrix_nbr, data%thread_nbr, MPI_INTEGER, &
                  &             matrix_nbr_rcv,  data%thread_nbr, MPI_INTEGER, &
                  &             data%MPI_COMM_TRANS, ierr)
             total = total + SUM(matrix_nbr_rcv)


             mat_size = &
                  ((data%elem_block_size/data%thread_nbr+1)*data%thread_nbr)*  &
                  (n_vertex_max*(n_order+1)*data%harm_size)**2
             CALL MPI_Alltoall(data%SEND_MATRICES, mat_size,                   &
                  &            MPI_DOUBLE_PRECISION,                           &
                  &            data%RECV_MATRICES, mat_size,                   &
                  &            MPI_DOUBLE_PRECISION,                           &
                  &            data%MPI_COMM_TRANS, ierr)

             prod_size = &
                  2*((data%elem_block_size/data%thread_nbr+1)*data%thread_nbr)*&
                  (n_vertex_max*(n_order+1))**2
             CALL MPI_Allgather(data%PROD_COLROW, prod_size, MPI_INTEGER,      &
                  &            data%RECV_COLROW, prod_size, MPI_INTEGER,       &
                  &            data%MPI_COMM_TRANS, ierr)


#ifdef PLENTY_TIMERS
             CALL SYSTEM_CLOCK(count=tt1)
             data%nb_periods_comm = data%nb_periods_comm + tt1-tt0
             IF (tt1<tt0)                                                      &
                  data%nb_periods_comm =                                       &
                  &   data%nb_periods_comm + data%nb_periods_max
#endif
             DO node = 1, (n_tor+1)/2
               DO thread = 1, data%thread_nbr
                 DO i = 1, matrix_nbr_rcv(thread, node)
                      index_node1 = data%RECV_COLROW(1, i, thread, node)
                      index_node2 = data%RECV_COLROW(2, i, thread, node)
                      DO iter = 1, data%my_harm_size**2
                         coefmtx(iter) = data%RECV_MATRICES(iter, i,           &
                              &                             thread, node)
                      END DO
#ifdef USE_MURGE
                      CALL MURGE_ASSEMBLYSETNODEVALUES(murge_id,               &
                           &                          index_node2,             &
                           &                          index_node1,             &
                           &                          coefmtx, ierr)

                      IF (ierr /= MURGE_SUCCESS) THEN
                         WRITE (*,*) data%my_id, "::::",                       &
                              "I", index_node2,                                &
                              "J", index_node1, ierr
                         data%ok = .FALSE.
                         RETURN
                      END IF
#else
                      PRINT *, "Binary built without murge"
                      data%ok = .FALSE.
                      RETURN
#endif
                      cnt2 = cnt2 + 1
                   END DO
                END DO
             END DO
          END IF
       END IF
#ifdef PLENTY_TIMERS
       CALL SYSTEM_CLOCK(count=t1)
       data%nb_periods_set = data%nb_periods_set + t1-t0
       IF (t1<t0)                                                              &
            data%nb_periods_set = data%nb_periods_set + data%nb_periods_max
#endif
       !$OMP barrier
    END DO ELEM
    IF (data%thread_num .EQ. 1) THEN
       DO thread = 1, data%thread_nbr
          DO j = 1, data%ndof_glob
             data%rhs_loc(j) = data%rhs_loc(j) + data%rhs_loc_thread(j,thread)
          END DO
       END DO
    END IF
    LOOP = 0
  END FUNCTION LOOP

END MODULE THREAD_DATA

SUBROUTINE construct_matrix_murge(my_id,node_list,element_list,                &
     &                            bnd_node_list, local_elms,                   &
     &                            n_local_elms, xpoint2, xcase2,               &
     &                            minRad, R_axis, Z_axis, psi_axis,            &
     &                            psi_bnd, Z_xpoint, psi_xpoint,               &
     &                            gmres, i_tor, n_cpu,                         &
     &                            mpi_comm_n, MPI_COMM_TRANS,                  &
     &                            my_id_trans, n_cpu_trans, solve_only)
  !---------------------------------------------------------------
  ! collect the element matrices into one large sparse matrix
  ! in coordinate format
  !---------------------------------------------------------------

  USE tr_module
  USE data_structure, ONLY : type_node, type_element,                          &
       &                     type_element_list, type_bnd_node_list,            &
       &                     type_node_list, thread_struct
  USE parameters,     ONLY : n_vertex_max , n_var, n_order, n_tor
  USE murge_module,   ONLY : MURGE_SUCCESS, MURGE_ASSEMBLYBEGIN,               &
       &                     MURGE_ASSEMBLYSETNODEVALUES,                      &
       &                     MURGE_ASSEMBLYEND, murge_id,                      &
       &                     murge_id_prod, MURGE_COEF_KIND,                   &
       &                     MURGE_ASSEMBLY_ADD,                               &
       &                     MURGE_ASSEMBLY_FOOL,                              &
       &                     MURGE_ASSEMBLY_RESPECT,                           &
       &                     MURGE_SYM, MURGE_SCAL_COL,                        &
       &                     MURGE_NORM_MAX_COL

  USE global_distributed_matrix, ONLY : RHS_GLOB, column_scaling, ndof_glob
  USE mumps_module, ONLY : mumps_par
  USE thread_data,  ONLY : thread_data_type, LOOP
  USE phys_module, ONLY : index_now
  USE mpi_mod
  IMPLICIT NONE
#include "r3_info.h"
  
  ! Subroutine parameters:
  INTEGER, TARGET                :: my_id
  TYPE (type_node_list), TARGET  :: node_list
  TYPE (type_element_list), TARGET :: element_list
  TYPE (type_bnd_node_list), TARGET :: bnd_node_list
  INTEGER, TARGET                :: n_local_elms
  INTEGER, TARGET                :: local_elms(n_local_elms)
  LOGICAL, TARGET                :: xpoint2
  REAL*8, TARGET                 :: minRad
  REAL*8, TARGET                 :: R_axis
  REAL*8, TARGET                 :: Z_axis
  REAL*8, TARGET                 :: psi_axis
  REAL*8, TARGET                 :: psi_bnd
  REAL*8, TARGET                 :: Z_xpoint(:)
  REAL*8, TARGET                 :: psi_xpoint(:)
  LOGICAL, TARGET                :: gmres
  INTEGER, TARGET                :: xcase2
  INTEGER                        :: i_tor(n_cpu)
  INTEGER                        :: n_cpu
  INTEGER                        :: column_number
  INTEGER                        :: mpi_comm_n
  INTEGER, TARGET                :: MPI_COMM_TRANS, my_id_trans, n_cpu_trans
  LOGICAL, TARGET                :: Solve_only
  ! local variables
  TYPE (type_element)            :: element
  INTEGER                        :: coefnbr, coefnbr_max, coefnbr_min
  INTEGER                        :: coefnbr_prod, coefnbr_prod_max
  INTEGER                        :: coefnbr_prod_min
  REAL(KIND=MURGE_COEF_KIND)     :: coefmtx(n_tor*n_var*n_tor*n_var)
  INTEGER, TARGET                :: elem_block_size
  INTEGER, TARGET                :: thread_nbr
  INTEGER                        :: elem_size
  INTEGER                        :: ELM_INDEX
  REAL*8,  POINTER               :: rhs_loc(:), rhs_loc_thread(:,:)
  REAL*8,  POINTER               :: PROD_MATRICES(:,:,:)
  REAL*8,  POINTER               :: SEND_MATRICES(:,:,:,:)
  REAL*8,  POINTER               :: RECV_MATRICES(:,:,:,:)
  INTEGER, POINTER               :: RECV_COLROW(  :,:,:,:)
  INTEGER, POINTER               :: PROD_COLROW(  :,:,:)
  INTEGER, POINTER               :: matrix_nbr(:)
  INTEGER, POINTER               :: matrix_nbr_rcv(:,:)
  REAL*8,  POINTER               :: ELM(:,:)
  REAL*8,  POINTER               :: RHS(:)
  LOGICAL                        :: ok, ok_recv

  TYPE (type_element)       :: element_father
  TYPE (type_node)          :: nodes_father(n_vertex_max)
  TYPE (type_node)          :: nodes(n_vertex_max)

  INTEGER :: index_ij, index_kl
  INTEGER :: index_large_k, k_order, knode
  INTEGER :: index_large_i
  INTEGER :: index_rhs
  INTEGER :: i, j, k, l
  INTEGER :: iv, iv2
  INTEGER :: ife, inode1, inode2, inode
  INTEGER :: i_father, inode_father
  INTEGER :: iter, ios
  INTEGER :: index_node1, index_node2
  INTEGER :: i_order, ielm, ierr
  LOGICAL :: is_local
  INTEGER :: cnt, cnt2
  INTEGER :: t0,t1,nb_periodes_max,nb_periodes_sec
  INTEGER :: nb_periods, max_periods, min_periods
  INTEGER :: nb_periods_ass_max, nb_periods_ass_min
  INTEGER :: nb_periods_set_max, nb_periods_set_min
  INTEGER :: nb_periods_elem_mat_max, nb_periods_elem_mat_min
  INTEGER :: nb_periods_comm_max, nb_periods_comm_min
  INTEGER, TARGET :: my_harm_num, my_harm_size, harm_size, node
  INTEGER, TARGET :: step
  INTEGER :: index_send_mtx, index_mtx, new_row_mat_elem, new_col_mat_elem
  INTEGER :: col_harm, row_harm
  CHARACTER(LEN=20), PARAMETER :: FMT_TIMING = "(A70,F7.2)"
  TYPE(THREAD_DATA_TYPE), ALLOCATABLE :: datas(:)
  INTEGER*8,              ALLOCATABLE :: threads(:)
  INTEGER :: ret, retval
  INTEGER, EXTERNAL :: omp_get_num_threads, omp_get_thread_num

  ! elapsed time
  CALL SYSTEM_CLOCK(count_rate=nb_periodes_sec,count_max=nb_periodes_max)
  ! timing
  CALL r3_info_begin (r3_info_index_0, 'solve_matrix_n')

!$omp PARALLEL shared(thread_nbr)
!$OMP master
  thread_nbr = omp_get_num_threads()
!$omp end master
!$omp end PARALLEL

  IF (my_id .EQ. 0) PRINT *, "THREAD_NBR", thread_nbr

  elem_block_size = 1*thread_nbr
  IF (MOD(elem_block_size, ((n_tor+1)/2)) .EQ. 0) THEN
     elem_block_size = n_local_elms/((n_tor+1)/2)
  ELSE
     elem_block_size = ((n_local_elms - MOD(elem_block_size,                   &
          &                                 ((n_tor+1)/2)))/((n_tor+1)/2))+1
  END IF
  elem_size = n_tor*n_vertex_max*(n_order+1)*n_var
  ! We allocate too much for harm_0
  harm_size = 2*n_var

  my_harm_num = i_tor(my_id+1)
  IF (my_harm_num == 1) THEN
     my_harm_size = n_var
  ELSE
     my_harm_size = 2*n_var
  END IF

  ! Allocate communication arrays
  IF (gmres .AND. .NOT. solve_only) THEN
     !****************************************************************
     !
     ! In GMRES mode, each processor will compute an elementar matrix
     !  and then all computed matrices will be exchanged through an
     !   MPI_allgather().
     !
     ! The same process is performed for right-hand-side member.
     !****************************************************************
     CALL tr_ALLOCATEp(RECV_MATRICES,1,harm_size**2,                           &
          1, (n_vertex_max*(n_order+1))**2*(elem_block_size/thread_nbr+1),     &
          1, thread_nbr, &
          1, (n_tor+1)/2, "RECV_MATRICES",CAT_DMATRIX)
     CALL tr_ALLOCATEp(SEND_MATRICES, 1, harm_size**2,                         &
          1, (n_vertex_max*(n_order+1))**2*(elem_block_size/thread_nbr+1),     &
          1, thread_nbr,                                                       &
          1, (n_tor+1)/2, "SEND_MATRICES",CAT_DMATRIX)
     CALL tr_ALLOCATEp(RECV_COLROW, 1, 2,                                      &
          1, (n_vertex_max*(n_order+1))**2*(elem_block_size/thread_nbr+1),     &
          1, thread_nbr, &
          1, (n_tor+1)/2, "RECV_COLROW",CAT_DMATRIX)
  END IF
  CALL tr_ALLOCATEp(PROD_MATRICES,1,(n_tor*n_var)**2,                          &
       1, (n_vertex_max*(n_order+1))**2*(elem_block_size/thread_nbr+1),        &
       1, thread_nbr, "PROD_MATRICES",CAT_DMATRIX)
  CALL tr_ALLOCATEp(PROD_COLROW, 1, 2,                                         &
       1, (n_vertex_max*(n_order+1))**2*(elem_block_size/thread_nbr+1),        &
       1, thread_nbr, "PROD_COLROW",CAT_DMATRIX)
  CALL tr_ALLOCATEp(matrix_nbr,1,thread_nbr, "matrix_nbr",CAT_DMATRIX)
  CALL tr_ALLOCATEp(matrix_nbr_rcv,1,thread_nbr,1,(n_tor+1)/2,"matrix_nbr_rcv",&
       &            CAT_DMATRIX)

  IF (my_id .EQ. 0) THEN
     WRITE(*,*) '****************************************'
     WRITE(*,*) '*  construct matrix MURGE              *'
     WRITE(*,*) '****************************************'
  END IF

  IF (ALLOCATED(rhs_glob)) CALL tr_deallocate(rhs_glob,"rhs_glob",CAT_DMATRIX)
  IF (.NOT. gmres) THEN
     column_number = ndof_glob
  ELSE
     IF (my_harm_num == 1) THEN
        column_number = ndof_glob/n_tor
     ELSE
        column_number = ndof_glob*2/n_tor
     END IF
  END IF

  CALL tr_allocate(rhs_glob,1,ndof_glob,"rhs_glob",CAT_DMATRIX)

  CALL tr_allocatep(rhs_loc,1,ndof_glob,"rhs_loc",CAT_DMATRIX)
  CALL tr_allocatep(rhs_loc_thread,1,ndof_glob,1,thread_nbr,"rhs_loc_thread",  &
       &            CAT_DMATRIX)

  RHS_glob = 0.d0
  RHS_loc  = 0.d0
  Rhs_loc_thread = 0.d0

  IF (gmres) THEN
     step = elem_block_size*(n_tor+1)/2
     CALL MPI_AllReduce(elem_block_size, step, 1, MPI_INTEGER, MPI_SUM,        &
          &             MPI_COMM_TRANS, ierr)
  ELSE
     step = elem_block_size
  END IF
  !
  ! Count coefnbr
  !
#ifdef PLENTY_TIMERS
  CALL SYSTEM_CLOCK(count=t0)
#endif
  coefnbr = 0
  DO ife =1, n_local_elms

     ielm = local_elms(ife)

     element = element_list%element(ielm)
     IF (element%n_sons .EQ. 0) THEN
        DO i=1,n_vertex_max

           inode1         = element%vertex(i)

           DO i_order = 1, n_order+1

              index_node1 = node_list%node(inode1)%index(i_order)

              CALL vertex_is_local(index_node1, is_local)
              IF (is_local) THEN
                 coefnbr = coefnbr + 1
              ENDIF

           ENDDO
        ENDDO
     END IF

  ENDDO


  coefnbr_prod = 0
  DO ife =1, n_local_elms, step
     DO ELM_INDEX = 1, elem_block_size
        IF ( ife + ELM_INDEX-1 +                                               &
             elem_block_size*my_id_trans > n_local_elms) EXIT

        ielm = local_elms(ife+ELM_INDEX-1+elem_block_size*my_id_trans)

        element = element_list%element(ielm)

        IF (element%n_sons .EQ. 0) THEN
           DO i=1,n_vertex_max

              inode1         = element%vertex(i)

              DO i_order = 1, n_order+1

                 index_node1 = node_list%node(inode1)%index(i_order)

                 CALL vertex_is_local(index_node1, is_local)
                 IF (is_local) THEN
                    coefnbr_prod = coefnbr_prod + 1
                 ENDIF
              END DO
           END DO
        END IF
     END DO

  ENDDO
  !*******************************************************************
  ! In gmres mode, each sub-communicator works on two tors, except the
  !  firs one which works on the first tor.
  !
  ! In direct mode, there is only one communicator, all the processors
  !  share the whole matrix.
  !*******************************************************************
  IF (.NOT. gmres) THEN
     coefnbr = coefnbr * (n_vertex_max)*(n_order+1)*((n_var*n_tor)**2)
  ELSE
     coefnbr_prod = &
          coefnbr_prod * (n_vertex_max)*(n_order+1)*((n_var*n_tor)**2)
     IF (my_harm_num == 1) THEN
        coefnbr = coefnbr * (n_vertex_max)*(n_order+1)*((n_var)**2)
     ELSE
        coefnbr = coefnbr * (n_vertex_max)*(n_order+1)*((n_var *2)**2)
     END IF
  END IF

#ifdef PLENTY_TIMERS
  CALL SYSTEM_CLOCK(count=t1)
  nb_periods = t1-t0
  IF (t1<t0) nb_periods = nb_periods + nb_periodes_max
  CALL MPI_Reduce(nb_periods, max_periods, 1, MPI_INTEGER, MPI_MAX, 0,         &
       &          MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(nb_periods, min_periods, 1, MPI_INTEGER, MPI_MIN, 0,         &
       &          MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) THEN
     WRITE(*,FMT_TIMING) ' maximum elapsed time counting assembly entries ',   &
          REAL(max_periods)/nb_periodes_sec
     WRITE(*,FMT_TIMING) ' minimum elapsed time counting assembly entries ',   &
          REAL(min_periods)/nb_periodes_sec
  END IF
#endif

  CALL MPI_Reduce(coefnbr, coefnbr_max, 1, MPI_INTEGER, MPI_MAX, 0,            &
       &          MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(coefnbr, coefnbr_min, 1, MPI_INTEGER, MPI_MIN, 0,            &
       &          MPI_COMM_WORLD, ierr)
  IF (.NOT. gmres) THEN
     IF (my_id .EQ. 0) THEN
        WRITE (*,"(A30,I10,A30)") ":: Murge Assembly phase :: ", coefnbr_max,  &
             " maximum entries"
        WRITE (*,"(A30,I10,A30)") ":: Murge Assembly phase :: ", coefnbr_min,  &
             " minimum entries"
     END IF
  ELSE
     CALL MPI_Reduce(coefnbr_prod, coefnbr_prod_max, 1, MPI_INTEGER, MPI_MAX,  &
          &          0, MPI_COMM_WORLD, ierr)
     CALL MPI_Reduce(coefnbr_prod, coefnbr_prod_min, 1, MPI_INTEGER, MPI_MIN,  &
          &          0, MPI_COMM_WORLD, ierr)
     IF (my_id .EQ. 0) THEN
        WRITE (*,"(A30,I10,A30)") ":: Murge Assembly phase :: ", coefnbr_max,  &
             " maximum entries"
        WRITE (*,"(A30,I10,A30)") ":: Murge Assembly phase :: ", coefnbr_min,  &
             " minimum entries"
        WRITE (*,"(A30,I10,A30)") ":: Murge Assembly phase :: ",               &
             coefnbr_prod_max, " maximum entries for product"
        WRITE (*,"(A30,I10,A30)") ":: Murge Assembly phase :: ",               &
             coefnbr_prod_min, " minimum entries for product"
     END IF
  END IF

  cnt = 0
  cnt2  = 0
  CALL SYSTEM_CLOCK(count=t0)
#ifdef USE_MURGE
  IF (.NOT. gmres .OR. .NOT. solve_only) THEN
     CALL MURGE_MATRIXRESET(murge_id, ierr)
     CALL MURGE_ASSEMBLYBEGIN(murge_id, coefnbr, MURGE_ASSEMBLY_ADD,           &
          &                   MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_RESPECT,      &
          &                   murge_sym, ierr)
  END IF
  IF (gmres) THEN
     CALL MURGE_MATRIXRESET(murge_id_prod, ierr)
     CALL MURGE_ASSEMBLYBEGIN(murge_id_prod, coefnbr_prod, MURGE_ASSEMBLY_ADD, &
          &                   MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_FOOL,         &
          &                   murge_sym, ierr)
  END IF
#else
  PRINT *, "Binary built without murge"
  CALL abort()
#endif
  !
  ! GMRES :
  !   Each processor Pi of each harmonic communicator will process
  !    elem_block_size elements.
  !
  !   For each element, for each n_vertex_max/order+1, we check if the
  !    column block is local. If it's local we had each node (n_tor
  !    *n_var**2) of this column to the product matrix and we had each
  !     digonal block, corresponding to one harmonic to the data to
  !      send.
  !
  !   Then the data is sent to the Pi processor of the correct harmonic
  !    and the local harmonic is received from other Pi programms.
  !
  !   After that the received data is treated, we had each node
  !    (my_harm_size**2) to the matrix to solve.
  !
  ! Direct :
  !   Here there is only one communicator, the matrix to solve is
  !    constructed like the product matrix. No communication has to be
  !     performed.

  ALLOCATE(datas(thread_nbr))
  CALL tr_register_mem(sizeof(datas),"datas",CAT_DMATRIX)
  ALLOCATE(threads(thread_nbr))
  CALL tr_register_mem(sizeof(threads),"threads",CAT_DMATRIX)
  DO iter = 1, thread_nbr
     !print *, "iter", iter
     datas(iter)%thread_num          = iter
     datas(iter)%thread_nbr          => thread_nbr
     datas(iter)%my_id               => my_id
     datas(iter)%my_id_trans         => my_id_trans
     datas(iter)%MPI_COMM_TRANS      => MPI_COMM_TRANS
     datas(iter)%node_list           => node_list
     datas(iter)%element_list        => element_list
     datas(iter)%ELM                 => thread_struct(iter)%ELM
     datas(iter)%RHS                 => thread_struct(iter)%RHS
     datas(iter)%rhs_loc             => rhs_loc
     datas(iter)%rhs_loc_thread      => rhs_loc_thread
     datas(iter)%SEND_MATRICES       => SEND_MATRICES
     datas(iter)%PROD_MATRICES       => PROD_MATRICES
     datas(iter)%RECV_MATRICES       => RECV_MATRICES
     datas(iter)%PROD_COLROW         => PROD_COLROW
     datas(iter)%RECV_COLROW         => RECV_COLROW
     datas(iter)%matrix_nbr          => matrix_nbr
     datas(iter)%matrix_nbr_rcv      => matrix_nbr_rcv
     datas(iter)%gmres               => gmres
     datas(iter)%solve_only          => solve_only
     datas(iter)%xpoint2             => xpoint2
     datas(iter)%xcase2              => xcase2
     datas(iter)%minRad              => minRad
     datas(iter)%R_axis              => R_axis
     datas(iter)%Z_axis              => Z_axis
     datas(iter)%psi_axis            => psi_axis
     datas(iter)%psi_bnd             => psi_bnd
     datas(iter)%Z_xpoint            => Z_xpoint
     datas(iter)%psi_xpoint          => psi_xpoint
     datas(iter)%local_elms          => local_elms
     datas(iter)%n_local_elms        => n_local_elms
     datas(iter)%step                => step
     datas(iter)%elem_block_size     => elem_block_size
     datas(iter)%harm_size           => harm_size
     datas(iter)%my_harm_size        => my_harm_size
#ifdef PLENTY_TIMERS
     datas(iter)%nb_periods_max      = nb_periodes_max
     datas(iter)%nb_periods_ass      = 0
     datas(iter)%nb_periods_set      = 0
     datas(iter)%nb_periods_elem_mat = 0
     datas(iter)%nb_periods_comm     = 0
#endif
     datas(iter)%ndof_glob           = ndof_glob
  END DO
!$OMP parallel private(iter,ret) default(shared)
  iter = 1+omp_get_thread_num()
  ret = LOOP(datas(iter))
!$OMP barrier
!$OMP end parallel
  ok = .TRUE.
#ifdef PLENTY_TIMERS
  nb_periods_ass_max = 0
  nb_periods_ass_min = datas(1)%nb_periods_ass
  nb_periods_set_max = 0
  nb_periods_set_min = datas(1)%nb_periods_set
  nb_periods_elem_mat_max = 0
  nb_periods_elem_mat_min = datas(1)%nb_periods_elem_mat
  nb_periods_comm_max = 0
  nb_periods_comm_min = datas(1)%nb_periods_elem_mat
#endif
  DO iter = 1, thread_nbr
     !call tr_deallocatep(datas(iter)%RHS,"datas(iter)%RHS")
     !call tr_deallocatep(datas(iter)%ELM,"datas(iter)%ELM")
     !NULLIFY(datas(iter)%RHS)
     !NULLIFY(datas(iter)%ELM)
     ok = datas(iter)%ok .AND. ok
#ifdef PLENTY_TIMERS
     nb_periods_ass_max = MAX(nb_periods_ass_max, datas(iter)%nb_periods_ass)
     nb_periods_ass_min = MIN(nb_periods_ass_min, datas(iter)%nb_periods_ass)
     nb_periods_set_max = MAX(nb_periods_set_max, datas(iter)%nb_periods_set)
     nb_periods_set_min = MIN(nb_periods_set_min, datas(iter)%nb_periods_set)
     nb_periods_elem_mat_max = MAX(nb_periods_elem_mat_max,                    &
          &                        datas(iter)%nb_periods_elem_mat)
     nb_periods_elem_mat_min = MIN(nb_periods_elem_mat_min,                    &
          &                        datas(iter)%nb_periods_elem_mat)
     nb_periods_comm_max = MAX(nb_periods_comm_max, datas(iter)%nb_periods_comm)
     nb_periods_comm_min = MIN(nb_periods_comm_min, datas(iter)%nb_periods_comm)
#endif
  END DO


  CALL MPI_Allreduce(ok, ok_recv, 1, MPI_LOGICAL, MPI_LAND,                    &
       &             MPI_COMM_WORLD, ierr)
  IF (.NOT. ok_recv) THEN
     PRINT *, "ERROR in MURGE call(s)"
     CALL abort()
  END IF
  CALL tr_unregister_mem(sizeof(datas),"datas",CAT_DMATRIX)
  DEALLOCATE(datas)
  CALL tr_unregister_mem(sizeof(threads),"threads",CAT_DMATRIX)
  DEALLOCATE(threads)

  IF (gmres .AND. .NOT. solve_only) THEN
     CALL tr_deallocatep(send_matrices,"send_matrices",CAT_DMATRIX)
     CALL tr_deallocatep(recv_matrices,"recv_matrices",CAT_DMATRIX)
     CALL tr_deallocatep(recv_colrow,"recv_colrow",CAT_DMATRIX)
  END IF
  CALL tr_deallocatep(prod_matrices,"prod_matrices",CAT_DMATRIX)
  CALL tr_deallocatep(prod_colrow,"prod_colrow",CAT_DMATRIX)
  CALL tr_deallocatep(matrix_nbr, "matrix_nbr",CAT_DMATRIX)
  CALL tr_deallocatep(matrix_nbr_rcv, "matrix_nbr_rcv",CAT_DMATRIX)
  CALL SYSTEM_CLOCK(count=t1)
#ifdef PLENTY_TIMERS
  CALL MPI_Reduce(nb_periods_ass_max, max_periods, 1, MPI_INTEGER, MPI_MAX, 0, &
       &          MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(nb_periods_ass_min, min_periods, 1, MPI_INTEGER, MPI_MIN, 0, &
       &          MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) THEN
     WRITE(*,FMT_TIMING) ' maximum elapsed time assemblying matrices ',        &
          REAL(max_periods)/nb_periodes_sec
     WRITE(*,FMT_TIMING) ' minimum elapsed time assemblying matrices ',        &
          REAL(min_periods)/nb_periodes_sec
  END IF

  CALL MPI_Reduce(nb_periods_set_max, max_periods, 1, MPI_INTEGER, MPI_MAX, 0, &
       &          MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(nb_periods_set_min, min_periods, 1, MPI_INTEGER, MPI_MIN, 0, &
       &          MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) THEN
     WRITE(*,FMT_TIMING) ' maximum elapsed time setting matrices into murge ', &
          REAL(max_periods)/nb_periodes_sec
     WRITE(*,FMT_TIMING) ' minimum elapsed time setting matrices into murge ', &
          REAL(min_periods)/nb_periodes_sec
  END IF

  CALL MPI_Reduce(nb_periods_elem_mat_max, max_periods, 1, MPI_INTEGER,        &
       &          MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(nb_periods_elem_mat_min, min_periods, 1, MPI_INTEGER,        &
       &          MPI_MIN, 0, MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) THEN
     WRITE(*,FMT_TIMING) ' maximum elapsed time building element matrices ',   &
          REAL(max_periods)/nb_periodes_sec
     WRITE(*,FMT_TIMING) ' minimum elapsed time building element matrices ',   &
          REAL(min_periods)/nb_periodes_sec
  END IF

  CALL MPI_Reduce(nb_periods_comm_max, max_periods, 1, MPI_INTEGER, MPI_MAX,   &
       &          0, MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(nb_periods_comm_min, min_periods, 1, MPI_INTEGER, MPI_MIN,   &
       &          0, MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) THEN
     WRITE(*,FMT_TIMING) ' maximum elapsed time communicating ',               &
          REAL(max_periods)/nb_periodes_sec
     WRITE(*,FMT_TIMING) ' minimum elapsed time communicating ',               &
          REAL(min_periods)/nb_periodes_sec
  END IF
#endif
  nb_periods = t1-t0
  IF (t1<t0) nb_periods = nb_periods + nb_periodes_max
  CALL MPI_Reduce(nb_periods, max_periods, 1, MPI_INTEGER, MPI_MAX, 0,         &
       &          MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(nb_periods, min_periods, 1, MPI_INTEGER, MPI_MIN, 0,         &
       &          MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) THEN
     WRITE(*,FMT_TIMING) ' maximum elapsed time in assembly loop ',            &
          REAL(max_periods)/nb_periodes_sec
     WRITE(*,FMT_TIMING) ' minimum elapsed time in assembly loop ',            &
          REAL(min_periods)/nb_periodes_sec
  END IF


  CALL SYSTEM_CLOCK(count=t0)

  IF (.NOT. gmres .OR. .NOT. solve_only) THEN
#ifdef USE_MURGE
     CALL MURGE_ASSEMBLYEND(murge_id, ierr)
#else
     PRINT *, "Binary built without murge"
     CALL abort()
#endif
     IF (ierr /= MURGE_SUCCESS) THEN
        IF (gmres) THEN
           WRITE (*,*) my_id, "::: error in assemblyend",                      &
                cnt2*my_harm_size**2, cnt2
        ELSE
           WRITE (*,*) my_id, "::: error in assemblyend..",                    &
                cnt*(n_tor*n_var)**2, cnt
        END IF
        CALL abort()
     END IF
  END IF

  IF (gmres) THEN
#ifdef USE_MURGE
     CALL MURGE_ASSEMBLYEND(murge_id_prod, ierr)
#else
       PRINT *, "Binary built without murge"
       CALL abort()
#endif
     IF (ierr /= MURGE_SUCCESS) THEN
        WRITE (*,*) my_id, "::: error in assemblyend..",                       &
             cnt*(n_tor*n_var)**2, cnt
        CALL abort()
     END IF
  END IF

  CALL SYSTEM_CLOCK(count=t1)
  nb_periods = t1-t0
  IF (t1<t0) nb_periods = nb_periods + nb_periodes_max
  CALL MPI_Reduce(nb_periods, max_periods, 1, MPI_INTEGER, MPI_MAX, 0,         &
       &          MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(nb_periods, min_periods, 1, MPI_INTEGER, MPI_MIN, 0,         &
       &          MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) THEN
     WRITE(*,FMT_TIMING) ' maximum elapsed time in MURGE_ASSEMBLYEND(s) ',     &
          REAL(max_periods)/nb_periodes_sec
     WRITE(*,FMT_TIMING) ' minimum elapsed time in MURGE_ASSEMBLYEND(s) ',     &
          REAL(min_periods)/nb_periodes_sec
  END IF

  !----------------------- boundary conditions

  CALL SYSTEM_CLOCK(count=t0)

!! For debugging purpose
!  CALL MPI_Reduce(RHS_loc, RHS_glob, ndof_glob, MPI_DOUBLE_PRECISION, MPI_SUM, &
!       &          0, MPI_COMM_WORLD, ierr)
!
!  call tr_locvnorms("cm_Rhs",RHS_glob,ndof_glob)
!  if (my_id .eq. 0) then
!     write(fname,'(A,I6.6)')"rhs",index_now
!     call tr_vdump(fname,RHS_glob,ndof_glob)
!  end if

  CALL boundary_conditions(my_id, node_list, element_list, bnd_node_list,      &
       &                   local_elms, n_local_elms, 0,                        &
       &                   0,  rhs_loc,    xpoint2, xcase2, psi_axis,          &
       &                   psi_bnd, Z_xpoint, psi_xpoint, gmres, solve_only)
  CALL MPI_Reduce(RHS_loc, RHS_glob, ndof_glob, MPI_DOUBLE_PRECISION, MPI_SUM, &
       &          0, MPI_COMM_WORLD, ierr)

!  if (my_id .eq. 0) then
!     write(fname,'(A,I6.6)')"rhsbc",index_now
!     call tr_vdump(fname,RHS_glob,ndof_glob)
!  end if

  call tr_locvnorms("cm_BCRhs",RHS_glob,ndof_glob)
  call tr_debug_writei("ndof_glob",ndof_glob)

  CALL SYSTEM_CLOCK(count=t1)
  nb_periods = t1-t0
  IF (t1<t0) nb_periods = nb_periods + nb_periodes_max
  CALL MPI_Reduce(nb_periods, max_periods, 1, MPI_INTEGER, MPI_MAX, 0,         &
       &          MPI_COMM_WORLD, ierr)
  CALL MPI_Reduce(nb_periods, min_periods, 1, MPI_INTEGER, MPI_MIN, 0,         &
       &          MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) THEN
     WRITE(*,FMT_TIMING) ' maximum elapsed time in boundary_conditions ',      &
          REAL(max_periods)/nb_periodes_sec
     WRITE(*,FMT_TIMING) ' minimum elapsed time in boundary_conditions ',      &
          REAL(min_periods)/nb_periodes_sec
  END IF

  IF (gmres .AND. .NOT. solve_only) THEN
     IF (ASSOCIATED(mumps_par%rhs))                                            &
          CALL tr_deallocatep(mumps_par%rhs,"mumps_par%rhs",CAT_DMATRIX)
     CALL tr_allocatep(mumps_par%rhs,1,column_number,"mumps_par%rhs",          &
          &            CAT_DMATRIX)
     mumps_par%n = column_number
     mumps_par%rhs = 0.d0

     IF (my_harm_num .EQ. 1 ) THEN

        mumps_par%rhs(1:column_number) = rhs_glob(1:ndof_glob:n_tor)

     ELSE
        mumps_par%rhs(1:column_number:2) =                                     &
             rhs_glob(2*(my_harm_num-1):ndof_glob:n_tor)
        mumps_par%rhs(2:column_number:2) =                                     &
             rhs_glob(2*(my_harm_num-1)+1:ndof_glob:n_tor)
     ENDIF


  END IF

  IF (.NOT. gmres .OR. .NOT. solve_only) THEN
     IF(ALLOCATED(column_scaling))                                             &
          CALL tr_deallocate(column_scaling,"column_scaling",CAT_DMATRIX)
     CALL tr_allocate(column_scaling,1,column_number,"column_scaling",         &
          &           CAT_DMATRIX)
     column_scaling = 1.d-20
#ifdef USE_MURGE
     CALL MURGE_GetGlobalNorm(murge_id, column_scaling, -1,                    &
          &                   MURGE_NORM_MAX_COL, ierr)
     CALL MURGE_ApplyGlobalScaling(murge_id, column_scaling,                   &
          &                        -1, MURGE_SCAL_COL, ierr)
#else
     PRINT *, "Binary built without murge"
     CALL abort()
#endif
  END IF
  IF (my_id .EQ. 0)                                                            &
       WRITE(*,*) '******** end construct matrix murge **********'

  !OPEN(unit=10, file='RHS_glob.txt', iostat=ios)
  !WRITE (10,*) RHS_glob
  !CLOSE(10)

  CALL tr_deallocatep(RHS_loc,"RHS_loc",CAT_DMATRIX)
  CALL tr_deallocatep(RHS_loc_thread,"RHS_loc_thread",CAT_DMATRIX)

  RETURN
END SUBROUTINE construct_matrix_murge
