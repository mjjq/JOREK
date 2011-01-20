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
  USE data_structure,  only : type_element_list&
       &, type_node_list


  TYPE THREAD_DATA_TYPE
     ! -- local variables   --
     INTEGER                            :: thread_num
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
     logical,                   POINTER :: gmres, solve_only
     LOGICAL,                   POINTER :: xpoint2
     REAL*8,                    POINTER :: psi_axis
     REAL*8,                    POINTER :: psi_bnd
     REAL*8,                    POINTER :: Z_xpoint
     INTEGER,                   POINTER :: local_elms(:)
     INTEGER,                   POINTER :: n_local_elms
     INTEGER,                   POINTER :: step
     INTEGER,                   POINTER :: elem_block_size
     INTEGER,                   POINTER :: harm_size, my_harm_size
     ! each thread access different part, no mutex
     REAL*8,                    POINTER :: PROD_MATRICES(:,:)
     INTEGER,                   POINTER :: PROD_COLROW(:,:)
     REAL*8,                    POINTER :: SEND_MATRICES(:,:,:)
     REAL*8,                    POINTER :: RECV_MATRICES(:,:,:)
     INTEGER,                   POINTER :: RECV_COLROW(:, :, :)
     ! only used by master, no mutex
     INTEGER,                   POINTER :: matrix_nbr_rcv(:)
     ! read/write, mutex protected
     REAL*8,                    POINTER :: rhs_loc(:)
     INTEGER,                   POINTER :: matrix_nbr
     ! barrier / mutex
     INTEGER*8,                 POINTER :: barrier
     INTEGER*8,                 POINTER :: mutex
  END TYPE THREAD_DATA_TYPE
contains

  INTEGER FUNCTION LOOP(data)
    !USE thread_data,     only : thread_data_type, synchro_x_threads
    USE data_structure,  only : type_node, type_element
    USE parameters,      only : n_vertex_max , n_var, n_order, n_tor
    USE murge_module,    only : MURGE_SUCCESS, MURGE_ASSEMBLYBEGIN,&
         & MURGE_ASSEMBLYSETNODEVALUES, MURGE_ASSEMBLYEND, murge_id,&
         & murge_id_prod, MURGE_COEF_KIND

    IMPLICIT NONE
    INCLUDE 'mpif.h'

    ! Function parameters
    TYPE(THREAD_DATA_TYPE),  INTENT(INOUT) :: data

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
    logical                        :: is_local
    TYPE (type_node)               :: nodes_father(n_vertex_max)
    REAL*8,  ALLOCATABLE           :: ELM(:,:)
    REAL*8,  ALLOCATABLE           :: RHS(:)
    INTEGER                        :: next_matrix
    INTEGER                        :: ret
    !INTEGER, external              :: fortran_pthread_mutex_lock, fortran_pthread_mutex_unlock
    INTEGER                        :: elem_size
    INTEGER                        :: ife
    INTEGER                        :: ELM_INDEX
    INTEGER                        :: node, iter, total, my_murde_id
    elem_size = n_tor*n_vertex_max*(n_order+1)*n_var
    cnt = 0
    cnt2 = 0
    !TODO: anticiper l'allocation ou l'allouer une fois pour toute
    ALLOCATE(ELM(elem_size, elem_size))
    ALLOCATE(RHS(elem_size))
    DO ife =1, data%n_local_elms, data%step 

!$OMP barrier       
!$OMP critical(matrix_nbr)       
       data%matrix_nbr = 0
!$OMP end critical(matrix_nbr)       
       !! do 
       DO ELM_INDEX = data%thread_num, data%elem_block_size, data%thread_nbr
          !print * , "ELM_INDEX", ELM_INDEX
          !print *, "ife+ELM_INDEX-1+elm_index*data%my_id_trans", ife,ELM_INDEX-1,data%elem_block_size,data%my_id_trans, "/", data%n_local_elms   
          IF (ife+ELM_INDEX-1+data%elem_block_size*data%my_id_trans <= data%n_local_elms) THEN
             ELM = 0.0
             RHS = 0.0

             ielm = data%local_elms(ife+ELM_INDEX-1+data%elem_block_size*data%my_id_trans)  

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

             IF (n_tor .GT. 3) THEN
                CALL element_matrix_fft(element,nodes, data%xpoint2, data%psi_axis,  &
                     &                  data%psi_bnd, Data%z_xpoint, ELM(:,:), RHS(:))      ! use fft for toroidal integration
             ELSE
                CALL element_matrix(element,nodes, data%xpoint2, data%psi_axis,&
                     &              data%psi_bnd, Data%z_xpoint, ELM(:,:), RHS(:))      ! use direct integration      
                DO iv = 1, n_vertex_max                                                                     ! boundary integrals

                   iv2  = MOD(iv, n_vertex_max) + 1

                   inode1 = element%vertex(iv)
                   inode2 = element%vertex(iv2)

                   !      if (     ((data%node_list%node(inode1)%boundary .eq. 1) .or.(data%node_list%node(inode1)%boundary .eq. 3)) &
                   !         .and. ((data%node_list%node(inode2)%boundary .eq. 1) .or.(data%node_list%node(inode2)%boundary .eq. 3)) ) then

                   ! nodes(1)  = data%node_list%node(inode1)
                   ! nodes(2)  = data%node_list%node(inode2)
                   ! vertex    = (/ iv, iv2 /)
                   !        direction = (/  1, 2   /)

                   ! write(*,*) iv,iv2,'boundary_matrix_open : ',inode1,inode2

                   !        call boundary_matrix_open(vertex, direction, element,nodes, data%xpoint2, data%psi_axis, data%psi_bnd, Data%z_xpoint, ELM, RHS)    ! for open field lines

                   !      endif

                   !IF (     ((data%node_list%node(inode1)%boundary .EQ. 2) .OR.(data%node_list%node(inode1)%boundary .EQ. 3)) &
                   !     .AND. ((data%node_list%node(inode2)%boundary .EQ. 2) .OR.(data%node_list%node(inode2)%boundary .EQ. 3)) ) THEN
                   !
                   !   nodes(1)  = data%node_list%node(inode1)
                   !   nodes(2)  = data%node_list%node(inode2)
                   !   vertex    = (/ iv, iv2 /)
                   !   direction = (/ 1, 3    /)

                   ! write(*,*) iv,iv2,'boundary_matrix : ',inode1,inode2

                   !        call boundary_matrix(vertex, direction, element,nodes, data%xpoint2, data%psi_axis, data%psi_bnd, Data%z_xpoint, ELM, RHS)    ! for closed field lines

                   !ENDIF

                ENDDO
             ENDIF

             !CALL ch_nod_rhs_elm(ielm,element,nodes,element_father,nodes_father,ELM,RHS,node_out) 
             IF (element%n_sons .EQ. 0) THEN
                DO i=1,n_vertex_max

                   !inode1         =node_out(i)! element%vertex(i)
                   inode1 =  element%vertex(i)
                   DO i_order = 1, n_order+1

                      ! index_node1 is the column
                      index_node1 = data%node_list%node(inode1)%index(i_order)

                      index_large_i = n_tor * n_var * (index_node1 - 1)

                      CALL vertex_is_local(index_node1, is_local)
                      IF (is_local) THEN

                         ! Set RHS member
                         DO j = 1, n_var * n_tor

                            ! index in the ELM matrix
                            index_ij = n_tor * n_var * (n_order+1) * (i-1) +&
                                 & n_tor * n_var * (i_order-1) + j   
                            ! index in global matrix
                            index_rhs = index_large_i+j
!$OMP critical(rhs_member)                            
                            data%rhs_loc(index_rhs) = data%rhs_loc(index_rhs) + RHS(index_ij)
!$OMP end critical (rhs_member)                           
                         END DO
                         ! Build nodes Matrices
                         DO k=1,n_vertex_max

                            knode         = element%vertex(k)

                            DO k_order = 1, n_order+1

                               index_node2 = data%node_list%node(knode)%index(k_order)

                               index_large_k = n_tor * n_var * (index_node2 - 1)

                               !coefmtx = 0
!$OMP critical(matrix_nbr)                            
                               data%matrix_nbr = data%matrix_nbr+1
                               next_matrix = data%matrix_nbr
!$OMP end critical (matrix_nbr)       
                               DO j = 1, n_var * n_tor
                                  ! Row index in the ELM matrix
                                  index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j   

                                  DO l = 1, n_var * n_tor
                                     ! BUILD node Matrix
                                     index_kl  = n_tor * n_var * (n_order&
                                          &+1) * (k-1) + n_tor * n_var *&
                                          & (k_order-1) + l   

                                     IF (data%gmres .AND. .NOT. data%solve_only) THEN
                                        col_harm = INT((MOD(j-1, n_tor)+1)/2) + 1
                                        row_harm = INT((MOD(l-1, n_tor)+1)/2) + 1
                                        IF (col_harm == row_harm) THEN

                                           IF (col_harm == 1) THEN
                                              new_col_mat_elem = INT((j-1)/n_tor)
                                              new_row_mat_elem = INT((l-1)&
                                                   &/n_tor)
                                              index_send_mtx = n_var*new_col_mat_elem +&
                                                   & new_row_mat_elem+1

                                           ELSE
                                              ! (num_var-1)*2 + num_tor_local (0 or 1)
                                              new_col_mat_elem = INT((j-1)/n_tor)*2 + MOD(MOD(j-1, n_tor)+1,2)
                                              new_row_mat_elem = INT((l-1)/n_tor)*2 + MOD(MOD(l-1, n_tor)+1,2)
                                              index_send_mtx = n_var*2   &
                                                   &*new_col_mat_elem +  &
                                                   & new_row_mat_elem+1

                                           END IF
                                           data%SEND_MATRICES(index_send_mtx,&
                                                &(next_matrix),&
                                                & col_harm) = ELM(index_kl& 
                                                &,index_ij)
                                        END IF
                                     END IF
                                     new_col_mat_elem = j-1
                                     new_row_mat_elem = l-1
                                     index_mtx = new_row_mat_elem+1&
                                          &+new_col_mat_elem&
                                          &*(n_tor*n_var)
                                     data%PROD_MATRICES(index_mtx, next_matrix) = &
                                          & ELM(index_kl&
                                          &,index_ij)
                                     !coefmtx(index_mtx) =&
                                     !     & ELM(index_kl&
                                     !     &,index_ij)

                                  END DO
                               END DO
                               data%PROD_COLROW(1, next_matrix) = index_node1
                               data%PROD_COLROW(2, next_matrix) = index_node2
                               !print *, data%my_id, data%thread_num, "Zone protegee"
                               !call fortran_pthread_mutex_lock(data%mutex, ret)
                               !!print *, data%my_id, data%thread_num, "Zone protegee in"
                               !IF (data%gmres) THEN
                               !   CALL MURGE_ASSEMBLYSETNODEVALUES(murge_id_prod, index_node2, index_node1, &
                               !        coefmtx, ierr)
                               !ELSE
                               !   CALL MURGE_ASSEMBLYSETNODEVALUES(murge_id, index_node2, index_node1, &
                               !        coefmtx, ierr)
                               !END IF
                               cnt = cnt + 1
                               !print *, data%my_id, data%thread_num, "Zone protegee out"
                               !call fortran_pthread_mutex_unlock(data%mutex, ret)
                               !print *, data%my_id, data%thread_num, "Zone protegee out final"
                               !IF (ierr /= MURGE_SUCCESS) THEN
                               !   WRITE (*,*) data%my_id, ":::", &
                               !        "I", index_node2, &
                               !        "J", index_node1, cnt
                               !   STOP
                               !END IF


                            END DO
                         END DO
                      END IF
                   END DO
                END DO
             END IF
             !print *, data%my_id, ife, "next_matrix", next_matrix
          END IF
       END DO
!$OMP barrier       

       IF (.not. data%gmres) THEN
          ! We work on the full problem with direct method
          my_murde_id = murge_id
          
          DO j = 1, data%matrix_nbr 
             index_node1 = data%PROD_COLROW(1, j) 
             index_node2 = data%PROD_COLROW(2, j) 
             DO iter = 1, (n_tor*n_var)**2
                coefmtx(iter)= data%PROD_MATRICES(iter,j)
             END DO
             CALL MURGE_ASSEMBLYSETNODEVALUES(my_murde_id, index_node2, index_node1, &
                  coefmtx, ierr)  
             IF (ierr /= MURGE_SUCCESS) THEN
                WRITE (*,*) data%my_id, ":::", &
                     "I", index_node2, &
                     "J", index_node1, cnt
                STOP
             END IF
          END DO

       ELSE
          ! We have to construct a problem for the product
          my_murde_id = murge_id_prod

          IF (data%thread_num == 1) THEN
             DO j = 1, data%matrix_nbr 
                index_node1 = data%PROD_COLROW(1, j) 
                index_node2 = data%PROD_COLROW(2, j) 
                DO iter = 1, (n_tor*n_var)**2
                   coefmtx(iter)= data%PROD_MATRICES(iter,j)
                END DO
                CALL MURGE_ASSEMBLYSETNODEVALUES(my_murde_id, index_node2, index_node1, &
                     coefmtx, ierr)  
                IF (ierr /= MURGE_SUCCESS) THEN
                   WRITE (*,*) data%my_id, ":::", &
                        "I", index_node2, &
                        "J", index_node1, cnt
                   STOP
                END IF
             END DO
          END IF
          ! Just so that the first thread doesn't do all the job.
          IF ((data%thread_nbr == 1 .or. data%thread_num == 2 ) .and. & 
               & .NOT. DATA%solve_only) THEN
             data%matrix_nbr_rcv = 0
             CALL MPI_Allgather(data%matrix_nbr, 1, MPI_INTEGER, &
                  &            data%matrix_nbr_rcv, 1,   MPI_INTEGER, &
                  &            data%MPI_COMM_TRANS, ierr)
             total = total + sum(data%matrix_nbr_rcv)
             !print *, data%my_id, data%thread_num, "data%matrix_nbr", data%matrix_nbr, data%matrix_nbr_rcv, total, cnt
             CALL MPI_Alltoall(data%SEND_MATRICES, data%elem_block_size*(n_vertex_max&
                  &            *(n_order+1)*data%harm_size)**2, MPI_DOUBLE_PRECISION, &
                  &            data%RECV_MATRICES, data%elem_block_size*(n_vertex_max&
                  &            *(n_order+1)*data%harm_size)**2, MPI_DOUBLE_PRECISION, &
                  &            data%MPI_COMM_TRANS, ierr)
             
             CALL MPI_Allgather(data%PROD_COLROW, data%elem_block_size*2*(n_vertex_max&
                  &            *(n_order+1))**2, MPI_INTEGER, &
                  &            data%RECV_COLROW, data%elem_block_size*2*(n_vertex_max&
                  &            *(n_order+1))**2, MPI_INTEGER, &
                  &            data%MPI_COMM_TRANS, ierr)
          
             DO node = 1, (n_tor+1)/2
                DO i = 1, data%matrix_nbr_rcv(node)
                   index_node1 = data%RECV_COLROW(1, i, node)
                   index_node2 = data%RECV_COLROW(2, i, node)
                   DO iter = 1, data%my_harm_size**2
                      coefmtx(iter) = data%RECV_MATRICES(iter, i, node)
                   END DO
                   
                   CALL MURGE_ASSEMBLYSETNODEVALUES(murge_id,         &
                        & index_node2,       &
                        & index_node1,     &
                        & coefmtx, ierr)
                   IF (ierr /= MURGE_SUCCESS) THEN
                      WRITE (*,*) data%my_id, "::::", &
                           "I", index_node2, &
                           "J", index_node1,&
                           & cnt2, cnt, ife, i, node, data%matrix_nbr_rcv(node),&
                           & data%n_local_elms,data%element_list%n_elements
                      STOP
                   END IF
                   cnt2 = cnt2 + 1
                   
                END DO
             END DO
          END IF
       END IF
    END DO
    DEALLOCATE(RHS,ELM)
    LOOP = 0
    !print *, "cnt, cnt2", cnt, cnt2
  END FUNCTION LOOP

END MODULE THREAD_DATA

SUBROUTINE construct_matrix_murge(my_id,node_list,element_list, local_elms, &
     n_local_elms, xpoint2,psi_axis,psi_bnd,Z_xpoint, gmres, i_tor, n_cpu, &
     mpi_comm_n, MPI_COMM_TRANS, my_id_trans, n_cpu_trans, solve_only)
  !---------------------------------------------------------------
  ! collect the element matrices into one large sparse matrix
  ! in coordinate format
  !---------------------------------------------------------------

  USE data_structure, only : type_node, type_element, type_element_list&
       &, type_node_list
  USE parameters,     only : n_vertex_max , n_var, n_order, n_tor
  USE murge_module,   only : MURGE_SUCCESS, MURGE_ASSEMBLYBEGIN,   &
       & MURGE_ASSEMBLYSETNODEVALUES, MURGE_ASSEMBLYEND, murge_id, &
       & murge_id_prod, MURGE_COEF_KIND, MURGE_ASSEMBLY_ADD,       &
       & MURGE_ASSEMBLY_FOOL, MURGE_SYM, MURGE_SCAL_COL,           &
       & MURGE_NORM_MAX_COL
  !USE data_structure
  USE global_distributed_matrix, only : RHS_GLOB, column_scaling, ndof_glob
  !USE phys_module
  !USE murge_module
  USE mumps_module, only : mumps_par
  USE thread_data,  only : thread_data_type, LOOP
  IMPLICIT NONE
  INCLUDE 'mpif.h'
#include "r3_info.h"

  ! Subroutine parameters:
  INTEGER, TARGET                :: my_id
  TYPE (type_node_list), target  :: node_list
  TYPE (type_element_list), target :: element_list
  INTEGER, target                :: n_local_elms
  INTEGER, target                :: local_elms(n_local_elms)
  LOGICAL, target                :: xpoint2
  REAL*8, target                 :: psi_axis
  REAL*8, target                 :: psi_bnd
  REAL*8, target                 :: Z_xpoint
  LOGICAL, target                :: gmres
  INTEGER                        :: i_tor(n_cpu)
  INTEGER                        :: n_cpu
  INTEGER                        :: column_number
  INTEGER                        :: mpi_comm_n
  INTEGER, target                :: MPI_COMM_TRANS, my_id_trans, n_cpu_trans
  LOGICAL, target                :: Solve_only
  ! local variables
  TYPE (type_element)            :: element
  INTEGER                        :: coefnbr, coefnbr_prod
  REAL(KIND=MURGE_COEF_KIND)     :: coefmtx(n_tor*n_var*n_tor*n_var)
  INTEGER, TARGET                :: elem_block_size
  INTEGER, TARGET                :: thread_nbr
  INTEGER                        :: elem_size
  INTEGER                        :: ELM_INDEX
  REAL*8,  POINTER               :: rhs_loc(:)
  REAL*8,  POINTER               :: PROD_MATRICES(:,:)
  REAL*8,  POINTER               :: SEND_MATRICES(:,:,:)
  REAL*8,  POINTER               :: RECV_MATRICES(:,:,:)
  INTEGER, POINTER               :: RECV_COLROW(  :,:,:)
  INTEGER, POINTER               :: PROD_COLROW(  :,:)
  REAL*8,  POINTER               :: ELM(:,:)
  REAL*8,  POINTER               :: RHS(:)
  INTEGER, TARGET                :: matrix_nbr
  INTEGER, POINTER               :: matrix_nbr_rcv(:)

  type (type_element)       :: element_father
  type (type_node)          :: nodes_father(n_vertex_max)
  type (type_node)          :: nodes(n_vertex_max)

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
  INTEGER :: nb_periods
  INTEGER, TARGET :: my_harm_num, my_harm_size, harm_size, node
  INTEGER, target :: step
  INTEGER :: index_send_mtx, index_mtx, new_row_mat_elem, new_col_mat_elem
  INTEGER :: col_harm, row_harm
  CHARACTER(LEN=20), PARAMETER :: FMT_TIMING = "(I2,A70,F7.2)"
  INTEGER*8, target   :: barrier
  TYPE(THREAD_DATA_TYPE), ALLOCATABLE :: datas(:)
  INTEGER*8,              ALLOCATABLE :: threads(:)
  INTEGER :: ret, retval
  INTEGER *8, target :: mutex
  INTEGER, external :: omp_get_num_threads, omp_get_thread_num

  CALL SYSTEM_CLOCK(count_rate=nb_periodes_sec,count_max=nb_periodes_max) ! elapsed time
  CALL r3_info_begin (r3_info_index_0, 'solve_matrix_n')                  ! timing
!$omp PARALLEL shared(thread_nbr)
!$OMP master
  thread_nbr = omp_get_num_threads()
  print *, "THREAD_NBR", thread_nbr
!$omp end master
!$omp end PARALLEL
  elem_block_size = thread_nbr
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
     ALLOCATE(RECV_MATRICES(harm_size**2,     (n_vertex_max*(n_order+1))**2*elem_block_size, (n_tor+1)/2))
     ALLOCATE(SEND_MATRICES(harm_size**2,     (n_vertex_max*(n_order+1))**2*elem_block_size, (n_tor+1)/2)) 

     ALLOCATE(RECV_COLROW(2, (n_vertex_max*(n_order+1))**2* elem_block_size, (n_tor+1)/2))
     ALLOCATE(matrix_nbr_rcv((n_tor+1)/2))
  END IF
  ALLOCATE(PROD_MATRICES((n_tor*n_var)**2, (n_vertex_max*(n_order+1))**2*elem_block_size)) 
  ALLOCATE(PROD_COLROW(2, (n_vertex_max*(n_order+1))**2* elem_block_size))

  WRITE(*,*) '****************************************'
  WRITE(*,*) '*  construct matrix MURGE              *'
  WRITE(*,*) '****************************************'
  WRITE(*,*) ' n_elements (local)       : ',my_id,n_local_elms

  IF (ALLOCATED(rhs_glob))        DEALLOCATE(rhs_glob)
  IF (.NOT. gmres) THEN
     column_number = ndof_glob
  ELSE
     IF (my_harm_num == 1) THEN
        column_number = ndof_glob/n_tor
     ELSE
        column_number = ndof_glob*2/n_tor
     END IF
  END IF

  ALLOCATE(rhs_glob(ndof_glob))

  ALLOCATE(rhs_loc(ndof_glob))

  RHS_glob = 0.d0
  RHS_loc  = 0.d0

  IF (gmres) THEN
     step = elem_block_size*(n_tor+1)/2
  ELSE
     step = elem_block_size
  END IF
  !
  ! Count coefnbr 
  !
  CALL SYSTEM_CLOCK(count=t0)
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
        IF (ife+ELM_INDEX-1+elem_block_size*my_id_trans > n_local_elms) EXIT
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
     coefnbr_prod = coefnbr_prod * (n_vertex_max)*(n_order+1)*((n_var*n_tor)&
          &**2)
     IF (my_harm_num == 1) THEN 
        coefnbr = coefnbr * (n_vertex_max)*(n_order+1)*((n_var)**2)
     ELSE
        coefnbr = coefnbr * (n_vertex_max)*(n_order+1)*((n_var *2)**2)
     END IF
  END IF

  CALL SYSTEM_CLOCK(count=t1)
  nb_periods = t1-t0
  IF (t1<t0) nb_periods = nb_periods + nb_periodes_max   
  WRITE(*,FMT_TIMING) my_id, ' system_clock elapsed time counting assembly entries ',REAL(nb_periods)/nb_periodes_sec

  IF (.NOT. gmres) THEN
     WRITE (*,*) ":: Murge Assembly phase :: ", coefnbr, " entries on processor ", my_id 
  ELSE
     WRITE (*,*) ":: Murge Assembly phase :: ", coefnbr, "/",&
          & coefnbr_prod, " entries on processor ", my_id 
  END IF
  cnt = 0
  cnt2  = 0
  CALL SYSTEM_CLOCK(count=t0)
  IF (.NOT. solve_only) THEN
     CALL MURGE_MATRIXRESET(murge_id, ierr)
     CALL MURGE_ASSEMBLYBEGIN(murge_id, coefnbr, MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_ADD, &
          MURGE_ASSEMBLY_FOOL, murge_sym, ierr)
  END IF
  CALL MURGE_MATRIXRESET(murge_id_prod, ierr)
  CALL MURGE_ASSEMBLYBEGIN(murge_id_prod, coefnbr_prod, MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_ADD, &
       MURGE_ASSEMBLY_FOOL, murge_sym, ierr)

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

  Allocate(datas(thread_nbr))
  Allocate(threads(thread_nbr))
  Do iter = 1, thread_nbr
     !print *, "iter", iter
     datas(iter)%thread_num        = iter
     datas(iter)%thread_nbr        => thread_nbr
     datas(iter)%my_id             => my_id
     datas(iter)%my_id_trans       => my_id_trans
     datas(iter)%MPI_COMM_TRANS    => MPI_COMM_TRANS     
     datas(iter)%node_list         => node_list          
     datas(iter)%element_list      => element_list       
     ALLOCATE(datas(iter)%ELM(elem_size, elem_size))
     ALLOCATE(datas(iter)%RHS(elem_size))
     datas(iter)%rhs_loc           => rhs_loc          
     datas(iter)%SEND_MATRICES     => SEND_MATRICES    
     datas(iter)%PROD_MATRICES     => PROD_MATRICES    
     datas(iter)%RECV_MATRICES     => RECV_MATRICES    
     datas(iter)%PROD_COLROW       => PROD_COLROW     
     datas(iter)%RECV_COLROW       => RECV_COLROW      
     datas(iter)%matrix_nbr        => matrix_nbr         
     datas(iter)%gmres             => gmres            
     datas(iter)%solve_only        => solve_only       
     datas(iter)%xpoint2           => xpoint2            
     datas(iter)%psi_axis          => psi_axis           
     datas(iter)%psi_bnd           => psi_bnd            
     datas(iter)%Z_xpoint          => Z_xpoint           
     datas(iter)%local_elms        => local_elms      
     datas(iter)%n_local_elms      => n_local_elms       
     datas(iter)%step              => step               
     datas(iter)%elem_block_size   => elem_block_size    
     datas(iter)%barrier           => barrier            
     datas(iter)%mutex             => mutex              
     datas(iter)%matrix_nbr_rcv    => matrix_nbr_rcv
     datas(iter)%harm_size         => harm_size        
     datas(iter)%my_harm_size      => my_harm_size
  end Do
!$OMP parallel private(iter,ret) default(shared)
  iter = 1+omp_get_thread_num()
  ret = LOOP(datas(iter))
!$OMP barrier
!$OMP end parallel

  Do iter = 1, thread_nbr
     DEALLOCATE(datas(iter)%RHS)
     DEALLOCATE(datas(iter)%ELM)
     NULLIFY(datas(iter)%RHS)
     NULLIFY(datas(iter)%ELM)
  end Do
  deallocate(datas)
  deallocate(threads)

  IF (gmres .AND. .NOT. solve_only) THEN
     DEALLOCATE(SEND_MATRICES, RECV_MATRICES,&
          & RECV_COLROW, matrix_nbr_rcv)
  END IF
  DEALLOCATE(PROD_MATRICES, PROD_COLROW)
  CALL SYSTEM_CLOCK(count=t1)
  nb_periods = t1-t0
  IF (t1<t0) nb_periods = nb_periods + nb_periodes_max   
  WRITE(*,FMT_TIMING) my_id, ' system_clock elapsed time in assembly loop ',REAL(nb_periods)/nb_periodes_sec
  CALL SYSTEM_CLOCK(count=t0)

  IF (.NOT. gmres .OR. .NOT. solve_only) THEN
     CALL MURGE_ASSEMBLYEND(murge_id, ierr)
     IF (ierr /= MURGE_SUCCESS) THEN
        IF (gmres) THEN
           WRITE (*,*) my_id, "::: error in assemblyend", &
                cnt2*my_harm_size**2, cnt2
        ELSE
           WRITE (*,*) my_id, "::: error in assemblyend..", &
                cnt*(n_tor*n_var)**2, cnt
        END IF
        STOP
     END IF
  END IF

  IF (gmres) THEN
     CALL MURGE_ASSEMBLYEND(murge_id_prod, ierr)
     IF (ierr /= MURGE_SUCCESS) THEN
        WRITE (*,*) my_id, "::: error in assemblyend..", &
             cnt*(n_tor*n_var)**2, cnt
        STOP
     END IF
  END IF

  CALL SYSTEM_CLOCK(count=t1)
  nb_periods = t1-t0
  IF (t1<t0) nb_periods = nb_periods + nb_periodes_max   
  WRITE(*,FMT_TIMING) my_id, ' system_clock elapsed time in MURGE_ASSEMBLYEND(s)',REAL(nb_periods)/nb_periodes_sec

  !----------------------- boundary conditions

  CALL SYSTEM_CLOCK(count=t0)
  WRITE (*,*) MY_ID, " : Reduce..."
  CALL MPI_Reduce(RHS_loc, RHS_glob, ndof_glob, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  CALL boundary_conditions(my_id, node_list, element_list, local_elms, n_local_elms, 0, &
       0,         xpoint2, psi_axis, psi_bnd, Z_xpoint, gmres, solve_only)
  CALL SYSTEM_CLOCK(count=t1)
  nb_periods = t1-t0
  IF (t1<t0) nb_periods = nb_periods + nb_periodes_max   
  WRITE(*,FMT_TIMING) my_id, ' system_clock elapsed time in boundary_conditions ',REAL(nb_periods)/nb_periodes_sec

  IF (gmres .AND. .NOT. solve_only) THEN
     IF (ASSOCIATED(mumps_par%rhs)) DEALLOCATE(mumps_par%rhs)
     ALLOCATE(mumps_par%rhs(column_number))
     mumps_par%n = column_number
     mumps_par%rhs = 0.d0

     IF (my_harm_num .EQ. 1 ) THEN

        mumps_par%rhs(1:column_number) = rhs_glob(1:ndof_glob:n_tor)

     ELSE
        mumps_par%rhs(1:column_number:2) = rhs_glob(2*(my_harm_num-1):ndof_glob:n_tor)
        mumps_par%rhs(2:column_number:2) = rhs_glob(2*(my_harm_num-1)+1:ndof_glob:n_tor)
     ENDIF


  END IF

  IF (.NOT. gmres) THEN
     IF(ALLOCATED(column_scaling))   DEALLOCATE(column_scaling)
     ALLOCATE(column_scaling(column_number))
     CALL MURGE_GetGlobalNorm(murge_id, column_scaling, -1, MURGE_NORM_MAX_COL, ierr)
     CALL MURGE_ApplyGlobalScaling(murge_id, column_scaling, -1, MURGE_SCAL_COL, ierr)
  END IF

  WRITE(*,*) '******** end construct matrix murge **********'
  OPEN(unit=10, file='RHS_glob.txt', iostat=ios)
  WRITE (10,*) RHS_glob
  CLOSE(10)
  DEALLOCATE(RHS_loc)

  RETURN
END SUBROUTINE construct_matrix_murge


