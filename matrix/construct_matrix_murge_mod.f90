! If PROD_MATRICES_STORAGE is defined: 
!   The elementery matrices concerning the product matrix
!   are stored during the assembly loop and then entered
!   later (as for harmonics)
! else
!   The elementary matrices are entered during the assembly
!   loop using a MURGE_AssemblySetValues() protected in a 
!   critical section.
!
#define PROD_MATRICES_STORAGE

! If PLENTY_TIMERS is defined some more timing will be displayed.

! MURGE_USE_SEQUENCE_HARM requires MURGE_USE_SEQUENCE
! Uses sequence in Murge that will remember the order
! the values were entered (not tested since a while,
! may be removed)
#ifdef MURGE_USE_SEQUENCE_HARM
#  ifndef MURGE_USE_SEQUENCE
#    define MURGE_USE_SEQUENCE
#  endif
#endif

! Some tool functions used in construct matrix 
MODULE CONSTRUCT_MATRIX_MURGE_TOOLS
CONTAINS
  ! Set a value to NaN
  REAL*8 FUNCTION SET_NAN()
    real*8 :: arg = -1.0
    real*8 :: x
    set_nan = dsqrt(arg)
    return
  END FUNCTION SET_NAN

  ! Check if a value is NaN
  LOGICAL FUNCTION IS_NAN(x)
    real*8 :: x
    is_nan = (x .ne. x) 
    return 
  END FUNCTION IS_NAN

  ! This function is used for reduction on RHS.
  ! NaN will be present in OUTVEC(i) only if both OUTVEC(i) and INVEC(i) are NaN
  ! If INVEC(i) is zero and OUTVEC(i) is not NaN, we keep OUTVEC(i)
  ! Otherwise we keep OUTVEC(i)
  SUBROUTINE KEEP_NON_NAN(INVEC, INOUTVEC, LEN, TYPE)
    INTEGER:: LEN, TYPE, i
    REAL*8 :: INVEC(LEN), INOUTVEC(LEN)
    DO i = 1, LEN
       if (.not. IS_NAN(INVEC(i))) then
          if (IS_NAN(INOUTVEC(i)) .or. INVEC(i) .ne. 0.d0) then
             INOUTVEC(i) = INVEC(i)
          end if
       end if
    END DO
    return
  END SUBROUTINE KEEP_NON_NAN
END MODULE CONSTRUCT_MATRIX_MURGE_TOOLS

! Data that are used in the assembly loop and function LOOP that
! implements the assembly loop
MODULE THREAD_DATA
  USE tr_module
  USE data_structure,  ONLY : type_element_list, type_node_list
#ifdef MURGE_USE_SEQUENCE
  USE murge_module, ONLY : MURGE_INTS_KIND, MURGE_COEF_KIND
#endif
  USE mod_elt_matrix_fft
  USE mod_elt_matrix
  use phys_module, only: n_tor_fft_thresh

  TYPE THREAD_DATA_TYPE
     ! -- local variables   --
     INTEGER                            :: thread_num
#ifdef MURGE_USE_SEQUENCE
     INTEGER                            :: cnt_entries
     INTEGER                            :: first_entry
#endif
#ifdef MURGE_USE_SEQUENCE_HARM
     INTEGER                            :: cnt_entries_harm
     INTEGER                            :: first_entry_harm
#endif
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
#ifdef MURGE_USE_SEQUENCE
     INTEGER,                   POINTER :: mode
     INTEGER(MURGE_INTS_KIND),  POINTER :: ROWS(:), COLS(:)
     REAL(MURGE_COEF_KIND),     POINTER :: VALS(:)
#endif
#ifdef MURGE_USE_SEQUENCE_HARM
     INTEGER(MURGE_INTS_KIND),  POINTER :: ROWS_HARM(:), COLS_HARM(:)
     REAL(MURGE_COEF_KIND),     POINTER :: VALS_HARM(:)
#endif
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
     REAL*8,                    POINTER :: R_xpoint(:)
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
  
  ! Assembly loop
  INTEGER FUNCTION LOOP(DATA)
    USE data_structure,  ONLY : type_node, type_element
    use mod_parameters,      ONLY : n_vertex_max , n_var, n_order, n_tor, jorek_model
    USE murge_module,    ONLY : MURGE_SUCCESS,                                 &
         &                      MURGE_ASSEMBLYSETNODEVALUES,                   &
         &                      MURGE_ASSEMBLYEND, murge_id,                   &
         &                      murge_id_prod, MURGE_COEF_KIND,                &
         &                      murge_ndof_prod
    USE phys_module,     ONLY : refinement
    USE murge_module,    ONLY : MURGE_ASSEMBLYBEGIN => MURGE_ASSEMBLYBEGIN_WRAPPER, vertex_is_local
    USE mpi_mod
    use construct_matrix_mod, only : elementary_matrix_build
    use data_structure,   only : thread_struct
    use mod_ch_nod_rhs_elm, only : ch_nod_rhs_elm

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
#ifdef MURGE_USE_SEQUENCE
    INTEGER                        :: index
#endif
    REAL(KIND=MURGE_COEF_KIND)     :: val
    INTEGER,               POINTER :: matrix_nbr_rcv(:,:), pt_matrix_nbr
    integer                           :: node_out(n_vertex_max)

    integer                        :: murge_id_glob
    IF (.NOT. data%gmres) THEN
       murge_id_glob = murge_id_prod
    ELSE
       murge_id_glob = murge_id
    END IF
    elem_size = n_tor*n_vertex_max*(n_order+1)*n_var
    cnt = 0
    cnt2 = 0
    pt_matrix_nbr => data%matrix_nbr(data%thread_num)
    data%ok = .TRUE.
    matrix_nbr_rcv => data%matrix_nbr_rcv

    ELM => data%ELM
    RHS => data%RHS

#ifdef MURGE_USE_SEQUENCE
    ! Reset entries counter
    DATA%cnt_entries = 0
#endif

#ifdef MURGE_USE_SEQUENCE_HARM
    DATA%cnt_entries_harm = 0
#endif
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

             ! --- Get element
             ielm    = data%local_elms( ife + ELM_INDEX - 1+                      &
                  data%elem_block_size*data%my_id_trans)

             element = data%element_list%element(ielm)
    
             ! --- Define nodes (this depends on whether our element has been refined)
             if (refinement) then
      
                i_father = data%element_list%element(ielm)%father

                if (i_father .ne. 0) then
                   element_father = data%element_list%element(i_father)
                   do iv = 1, n_vertex_max
                      inode_father     = element_father%vertex(iv)
                      nodes_father(iv) = data%node_list%node(inode_father)
                   enddo
                endif
                
             else
    	 
                do iv = 1, n_vertex_max
                   inode     = element%vertex(iv)
                   nodes(iv) = data%node_list%node(inode)
                enddo

             endif


             !! Compute the element matrix
#ifdef PLENTY_TIMERS
             CALL SYSTEM_CLOCK(count=tt0)
#endif
#ifdef MURGE_USE_SEQUENCE
             IF (DATA%mode .EQ. 3) THEN
#endif
                call elementary_matrix_build(element, nodes, data%xpoint2, data%xcase2,       &
                     &                       data%minRad, data%R_axis, data%Z_axis,           &
                     &                       data%psi_axis, data%psi_bnd, data%R_xpoint,      &
                     &                       data%Z_xpoint,                                   &
                     &                       data%thread_num, ife,              &
                     &                       data%n_local_elms, data%node_list)

#ifdef MURGE_USE_SEQUENCE
             END IF
#endif
#ifdef PLENTY_TIMERS
             CALL SYSTEM_CLOCK(count=tt1)
             data%nb_periods_elem_mat = data%nb_periods_elem_mat + tt1-tt0
             IF (tt1<tt0)                                                      &
                  data%nb_periods_elem_mat = data%nb_periods_elem_mat +        &
                  &                          data%nb_periods_max
#endif
             ! --- Define element nodes (depends if it's refined)
             if (refinement) then   
                call ch_nod_rhs_elm(ielm,element,nodes,element_father,nodes_father, &
                    thread_struct(data%thread_num)%ELM, thread_struct(data%thread_num)%RHS, node_out)
                ! note that ELM points to thread_struct(data%thread_num)%ELM
             else
                do i=1, n_vertex_max
                   node_out(i) = element%vertex(i)   
                enddo
             endif

             ! --- We only look at non-refined elements
             if ((.not. refinement) .or. (refinement .and. (element%n_sons .eq. 0))) then
                VERTEX_COL : DO i=1,n_vertex_max

                   inode1 =  node_out(i)
                   ORDER_COL : DO i_order = 1, n_order+1

                      ! index_node1 is the column
                      index_node1 = data%node_list%node(inode1)%index(i_order)

                      index_large_i = n_tor * n_var * (index_node1 - 1)

                      CALL vertex_is_local(index_node1, is_local)
                      IF (is_local) THEN
#ifdef MURGE_USE_SEQUENCE
                         IF (DATA%mode .eq. 3) THEN
#endif
                            ! Set RHS member
                            DO j = 1, n_var * n_tor

                               ! index in the ELM matrix
                               index_ij = n_tor * n_var * (n_order+1) *        &
                                    (i-1) + n_tor * n_var * (i_order-1) + j
                               ! index in global matrix
                               index_rhs = index_large_i+j
                               !$omp critical
                               if (index_rhs .eq. 51622) then
                                  write(*,"(I2,1X,A16,I6,A3,I3,A3,E20.12,1X,E20.12,1X,I4)") data%my_id, "rhs_loc_thread( ", index_rhs, " , ", data%thread_num, " ) ", data%rhs_loc(index_rhs), RHS(index_ij),ielm
                               end if
                               data%rhs_loc(index_rhs) = data%rhs_loc(index_rhs) + RHS(index_ij)
                               !$omp end critical                                  
                               !!!   if (index_rhs .eq. 51622) then
                               !!!      !$omp critical
                               !!!      write(*,"(I2,1X,A16,I6,A3,I3,A3,E20.12,1X,E20.12,1X,I4)") data%my_id, "rhs_loc_thread( ", index_rhs, " , ", data%thread_num, " ) ", data%rhs_loc_thread(index_rhs, data%thread_num), RHS(index_ij),ielm
                               !!!      !$omp end critical                                  
                               !!!   end if
                               !!!   
                               !!!   data%rhs_loc_thread(index_rhs,                  &
                               !!!        &              data%thread_num) =          &
                               !!!        data%rhs_loc_thread(index_rhs,             &
                               !!!        &                   data%thread_num) +     &
                               !!!        &                   RHS(index_ij)
                            END DO
#ifdef MURGE_USE_SEQUENCE
                         END IF
#endif

                         ! Build nodes Matrices
                         VERTEX_ROW : DO k=1,n_vertex_max

                            knode         = node_out(k)

                            ORDER_ROW : DO k_order = 1, n_order+1

                               index_node2 =                                   &
                                    data%node_list%node(knode)%index(k_order)

                               index_large_k = n_tor * n_var * (index_node2 - 1)

                               !coefmtx = 0
                               pt_matrix_nbr = pt_matrix_nbr + 1
#ifdef MURGE_USE_SEQUENCE
                               index = DATA%first_entry + DATA%cnt_entries
                               SELECT CASE (DATA%mode)
                               CASE (1)
                               CASE (2)
                                  DATA%ROWS(index) = index_node2
                                  DATA%COLS(index) = index_node1
                               CASE (3)
#endif
                                  DOF_COL : DO j = 1, n_var * n_tor
                                     ! Row index in the ELM matrix
                                     index_ij = n_tor * n_var *                &
                                          (n_order+1) * (i-1) +                &
                                          n_tor * n_var * (i_order-1) + j

                                     DOF_ROW : DO l = 1, n_var * n_tor
                                        ! BUILD node Matrix
                                        index_kl  = n_tor * n_var * (n_order   &
                                             &+1) * (k-1) + n_tor * n_var *    &
                                             & (k_order-1) + l

                                        IF ( data%gmres .AND.                  &
                                             .NOT. data%solve_only) THEN
                                           col_harm =                          &
                                                INT((MOD(j-1, n_tor)+1)/2) + 1
                                           row_harm =                          &
                                                INT((MOD(l-1, n_tor)+1)/2) + 1
                                           IF (col_harm == row_harm) THEN

                                              IF (col_harm == 1) THEN
                                                 new_col_mat_elem =            &
                                                      INT((j-1)/n_tor)
                                                 new_row_mat_elem =            &
                                                      INT((l-1)/n_tor)
                                                 index_send_mtx =              &
                                                      n_var*new_col_mat_elem + &
                                                      new_row_mat_elem+1

                                              ELSE
                                                 ! (num_var-1)*2 + &
                                                 !   num_tor_local (0 or 1)
                                                 new_col_mat_elem =            &
                                                      INT((j-1)/n_tor)*2 +     &
                                                      MOD(MOD(j-1, n_tor)+1,2)
                                                 new_row_mat_elem =            &
                                                      INT((l-1)/n_tor)*2 +     &
                                                      &MOD(MOD(l-1, n_tor)+1,2)
                                                 index_send_mtx =              &
                                                      n_var*2*new_col_mat_elem &
                                                      + new_row_mat_elem+1

                                              END IF

                                              data%SEND_MATRICES(              &
                                                   &          index_send_mtx,  &
                                                   &          pt_matrix_nbr,   &
                                                   &          data%thread_num, &
                                                   &          col_harm) =      &
                                                   &  ELM(index_kl,index_ij)
                                           END IF
                                        END IF
                                        new_col_mat_elem = j-1
                                        new_row_mat_elem = l-1
                                        index_mtx = new_row_mat_elem+1         &
                                             & +new_col_mat_elem               &
                                             & *(n_tor*n_var)
#ifdef MURGE_USE_SEQUENCE
                                        data%VALS((index-1)*(n_tor*n_var)**2+  &
                                             &     index_mtx) =                &
                                             &  ELM(index_kl,index_ij)
#else
#ifdef PROD_MATRICES_STORAGE
                                        data%PROD_MATRICES(index_mtx,          &
                                             &             pt_matrix_nbr,      &
                                             &             data%thread_num) =  &
                                             & ELM(index_kl,index_ij)
#else
                                        !$omp critical
                                        CALL MURGE_ASSEMBLYSETVALUE(murge_id_prod,             &
                                             &                      (index_node2-1)*n_tor*n_var + l,               &
                                             &                      (index_node1-1)*n_tor*n_var + j,               &
                                             &                       ELM(index_kl,index_ij), ierr)
                                        !$omp end critical
#endif
#endif
                                     END DO DOF_ROW
                                  END DO DOF_COL
#ifdef MURGE_USE_SEQUENCE
                               CASE DEFAULT
                                  !$omp critical
                                  WRITE (0,*), __FILE__,__LINE__,              &
                                       "Unknown mode", DATA%mode
                                  !$omp end critical
                                  CALL ABORT()
                               END SELECT
                               data%cnt_entries = data%cnt_entries+1
#endif
#ifdef MURGE_USE_SEQUENCE
                               IF (data%mode .eq. &
#  ifdef MURGE_USE_SEQUENCE_HARM
                                    2 &
#  else
                                    3 &
#  endif
                                    ) THEN
#endif
                                  data%PROD_COLROW(1, pt_matrix_nbr,           &
                                       &        data%thread_num) = index_node1
                                  data%PROD_COLROW(2, pt_matrix_nbr,           &

                                       &        data%thread_num) = index_node2
                                  !! HERE We could insert the entry without needs of this data%PROD_MATRICES
                                  !! within a omp critical
                                  !!CALL MURGE_ASSEMBLYSETNODEVALUES(murge_id_glob,             &
                                  !!     &                           index_node2,               &
                                  !!     &                           index_node1,               &
                                  !!     &                           coefmtx, ierr)

                                  !!!!!omp critical
                                  !!!!if (index_node2 .eq. 1 .and. index_node1 .eq. 1) then
                                  !!!!   print *, "PROD, 1, 1",  data%PROD_MATRICES(1,          &
                                  !!!!        &             pt_matrix_nbr,      &
                                  !!!!        &             data%thread_num)
                                  !!!!end if
                                  !!!!!omp end critical
#ifdef MURGE_USE_SEQUENCE
                               END IF
#endif
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
            data%nb_periods_ass = data%nb_periods_ass + data%nb_periods_max
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
                       !$omp critical
                       WRITE (*,*) data%my_id, ":::",                           &
                            "I", index_node2,                                   &
                            "J", index_node1, cnt
                       !$omp end critical
                       data%ok = .FALSE.
                       RETURN
                    END IF
                 END DO
              END DO
           END IF
       ELSE
#ifdef PROD_MATRICES_STORAGE
          ! We have to construct a problem for the product
          IF (data%thread_num == 1) THEN
             DO thread = 1, data%thread_nbr
                DO j = 1, data%matrix_nbr(thread)
                   index_node1 = data%PROD_COLROW(1, j, thread)
                   index_node2 = data%PROD_COLROW(2, j, thread)
#ifdef USE_MURGE
#  ifndef MURGE_USE_SEQUENCE
                   DO iter = 1, (n_tor*n_var)**2
                      coefmtx(iter)= data%PROD_MATRICES(iter,j, thread)
                   END DO
                   cnt = cnt + 1
                   if ( index_node2 .eq. 1+(85126-1)/(n_tor*n_var) .and. &
                        index_node1 .eq. 1+(85117-1)/(n_tor*n_var)) then
                      print *, data%my_id, "..PROD, 18, 17", coefmtx((mod(85117-1,n_tor*n_var))*(n_tor*n_var)+mod(85126-1,n_tor*n_var)+1)
                   end if
                   CALL MURGE_ASSEMBLYSETNODEVALUES(murge_id_prod,             &
                        &                           index_node2,               &
                        &                           index_node1,               &
                        &                           coefmtx, ierr)
                   !if (index_node1 .eq.1) then
                   !   do iter = 1, n_tor*n_var
                   !      print *, data%my_id, "coefmtx_prod(",iter, ")", coefmtx(iter)
                   !   end do
                   !end if
                   IF (ierr /= MURGE_SUCCESS) THEN
                      !$omp critical
                      WRITE (*,*) data%my_id, ":::",                           &
                           "I", index_node2,                                   &
                           "J", index_node1, cnt
                      !$omp end critical
                      data%ok = .FALSE.
                      RETURN
                   END IF
#  endif  
#else     
                   PRINT *, "Binary built without murge"
                   data%ok = .FALSE.
                   RETURN
#endif    
                END DO
             END DO
          END IF
#endif

          ! Just so that the first thread doesn't do all the job.
          IF ((data%thread_nbr == 1 .OR. data%thread_num == 2 ) .AND.      &
               & .NOT. DATA%solve_only) THEN
#ifdef MURGE_USE_SEQUENCE_HARM
             IF (data%mode .eq. 1 ) THEN
                ! for harmonics all is done by thread 1
                ! Needs to be improved...
                index = 0
                DO thread = 1, data%thread_nbr
                   index = index + data%matrix_nbr(thread)
                END DO

                CALL MPI_Allreduce(index,      &
                     &             data%cnt_entries_harm, &
                     &             1, MPI_INTEGER,        &
                     &             MPI_SUM,               &
                     &             data%mpi_comm_trans, ierr)
             END IF
#endif
#ifdef MURGE_USE_SEQUENCE
             IF (data%mode &
#  ifdef MURGE_USE_SEQUENCE_HARM
                  .ge. 2 &
#  else
                  .eq. 3 &
#  endif
                  ) THEN
#endif

                matrix_nbr_rcv = 0
#ifdef PLENTY_TIMERS
                CALL SYSTEM_CLOCK(count=tt0)
#endif
                CALL MPI_Allgather(data%matrix_nbr, data%thread_nbr, &
                     &             MPI_INTEGER,                      &
                     &             matrix_nbr_rcv,  data%thread_nbr, &
                     &             MPI_INTEGER,                      &
                     &             data%MPI_COMM_TRANS, ierr)
                total = total + SUM(matrix_nbr_rcv)

#ifdef MURGE_USE_SEQUENCE_HARM
                IF (data%mode .eq.3) THEN
#endif
                   mat_size = &
                        ((data%elem_block_size/data%thread_nbr+1)*     &
                        data%thread_nbr)*                              &
                        (n_vertex_max*(n_order+1)*data%harm_size)**2
                   CALL MPI_Alltoall(data%SEND_MATRICES, mat_size,             &
                        &            MPI_DOUBLE_PRECISION,                     &
                        &            data%RECV_MATRICES, mat_size,             &
                        &            MPI_DOUBLE_PRECISION,                     &
                        &            data%MPI_COMM_TRANS, ierr)
#ifdef MURGE_USE_SEQUENCE_HARM
                END IF
#endif
#ifdef MURGE_USE_SEQUENCE_HARM
                IF (data%mode .eq. 2) THEN
#endif
                   prod_size = &
                        2*((data%elem_block_size/data%thread_nbr+1)*     &
                        data%thread_nbr)*                                &
                        (n_vertex_max*(n_order+1))**2
                   CALL MPI_Allgather(data%PROD_COLROW, prod_size, MPI_INTEGER,&
                        &            data%RECV_COLROW, prod_size, MPI_INTEGER, &
                        &            data%MPI_COMM_TRANS, ierr)
#ifdef MURGE_USE_SEQUENCE_HARM
                END IF
#endif

#ifdef PLENTY_TIMERS
                CALL SYSTEM_CLOCK(count=tt1)
                data%nb_periods_comm = data%nb_periods_comm + tt1-tt0
                IF (tt1<tt0)                                                   &
                     data%nb_periods_comm =                                    &
                     &   data%nb_periods_comm + data%nb_periods_max
#endif
                DO node = 1, (n_tor+1)/2
                   DO thread = 1, data%thread_nbr
                      DO i = 1, matrix_nbr_rcv(thread, node)
                         index_node1 = data%RECV_COLROW(1, i, thread, node)
                         index_node2 = data%RECV_COLROW(2, i, thread, node)
#ifdef USE_MURGE
#  ifdef MURGE_USE_SEQUENCE_HARM
                         index = DATA%first_entry_harm + DATA%cnt_entries_harm
                         SELECT CASE (DATA%mode)
                         CASE (1)
                         CASE (2)
                            DATA%COLS_HARM(index) = index_node1
                            DATA%ROWS_HARM(index) = index_node2
                         CASE (3)
                            DO iter = 1, data%my_harm_size**2
                               DATA%VALS_HARM((index-1)*data%my_harm_size**2 +    &
                                    iter) = data%RECV_MATRICES(iter, i,        &
                                    &                          thread, node)
                            END DO
                         CASE DEFAULT
                            !$omp critical
                            WRITE (0,*), __FILE__,__LINE__,                 &
                                 "Unknown mode", DATA%mode
                            !$omp end critical
                            CALL ABORT()
                         END SELECT
                         data%cnt_entries_harm = data%cnt_entries_harm+1
#  else
                         DO iter = 1, data%my_harm_size**2
                            coefmtx(iter) = data%RECV_MATRICES(iter, i,        &
                                 &                             thread, node)
                         END DO
                         CALL MURGE_ASSEMBLYSETNODEVALUES(murge_id,            &
                              &                          index_node2,          &
                              &                          index_node1,          &
                              &                          coefmtx, ierr)

                         IF (ierr /= MURGE_SUCCESS) THEN
                            WRITE (*,*) data%my_id, "::::",                    &
                                 "I", index_node2,                             &
                                 "J", index_node1, ierr
                            data%ok = .FALSE.
                            RETURN
                         END IF
#  endif
#else
                         PRINT *, "Binary built without murge"
                         data%ok = .FALSE.
                         RETURN
#endif
                         cnt2 = cnt2 + 1
                      END DO
                   END DO
                END DO
#ifdef MURGE_USE_SEQUENCE
             END IF
#endif
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
#ifdef MURGE_USE_SEQUENCE
    IF (data%mode .EQ. 3) THEN
#endif
       !!!   IF (data%thread_num .EQ. 1) THEN
       !!!      DO thread = 1, data%thread_nbr
       !!!         DO j = 1, data%ndof_glob
       !!!            data%rhs_loc(j) = data%rhs_loc(j) + &
       !!!                 data%rhs_loc_thread(j,thread)
       !!!            if (j .eq. 51622) THEN
       !!!               !$omp critical
       !!!               write(*,"(I2,1X,A16,I6,A3,E20.12,1X,E20.12,1X,I4)") data%my_id, "rhs_loc( ", j,  " ) ", data%rhs_loc(j), data%rhs_loc_thread(j, thread)
       !!!               !$omp end critical  
       !!!            end if
       !!!         END DO
       !!!      END DO
       !!!   END IF
#ifdef MURGE_USE_SEQUENCE
    END IF
#endif
    LOOP = 0
  END FUNCTION LOOP

END MODULE THREAD_DATA

MODULE construct_matrix_murge_mod
contains
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
SUBROUTINE construct_matrix_murge(my_id,node_list,element_list,                &
     &                            bnd_node_list, local_elms,                   &
     &                            n_local_elms, xpoint2, xcase2,               &
     &                            minRad, R_axis, Z_axis, psi_axis,            &
     &                            psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,     &
     &                            gmres, i_tor, n_cpu,                         &
     &                            mpi_comm_n, MPI_COMM_TRANS,                  &
     &                            my_id_trans, n_cpu_trans, solve_only)
  !---------------------------------------------------------------
  ! collect the element matrices into one large sparse matrix
  ! in coordinate format
  !---------------------------------------------------------------

  USE tr_module
  USE Construct_matrix_murge_tools, only : keep_non_nan, is_nan, set_nan
  USE data_structure, ONLY : type_node, type_element,                          &
       &                     type_element_list, type_bnd_node_list,            &
       &                     type_node_list, thread_struct
  use mod_parameters,     ONLY : n_vertex_max , n_var, n_order, n_tor
  USE murge_module,   ONLY : MURGE_SUCCESS,                                    &
       &                     MURGE_ASSEMBLYSETNODEVALUES,                      &
       &                     MURGE_ASSEMBLYEND, murge_id,                      &
       &                     murge_id_prod, MURGE_COEF_KIND,                   &
       &                     MURGE_ASSEMBLY_ADD,                               &
       &                     MURGE_ASSEMBLY_FOOL,                              &
       &                     MURGE_ASSEMBLY_RESPECT,                           &
       &                     MURGE_SYM, MURGE_SCAL_COL,                        &
       &                     MURGE_NORM_MAX_COL, murge_need_rebuild_sequence,  &
       &                     MURGE_ROWS, MURGE_COLS, MURGE_VALS,               &
       &                     MURGE_ROWS_HARM, MURGE_COLS_HARM,                 &
       &                     MURGE_VALS_HARM, MURGE_BOOLEAN_TRUE,              &
       &                     murge_sequence_id, murge_assembly_first_entry,    &
       &                     murge_sequence_id_harm,                           &
       &                     murge_assembly_first_entry_harm
  USE murge_module,    ONLY : MURGE_ASSEMBLYBEGIN => MURGE_ASSEMBLYBEGIN_WRAPPER
  USE murge_module, ONLY : MURGE_MatrixReset
  USE murge_module, ONLY : MURGE_GetGlobalNorm
  USE murge_module, ONLY : MURGE_ApplyGlobalScaling

  USE global_distributed_matrix, ONLY : RHS_GLOB, column_scaling, ndof_glob
  USE mumps_module, ONLY : mumps_par
  USE thread_data,  ONLY : thread_data_type, LOOP
  USE phys_module,  ONLY : index_now, freeboundary, freeboundary_equil, resistive_wall
  USE mpi_mod
  USE murge_module, ONLY : murge_assembly_step, murge_elem_block_size, &
       murge_global_n, murge_global_n_prod, vertex_is_local
  use mod_boundary_conditions, only : boundary_conditions
  !$ use omp_lib
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
  REAL*8, TARGET                 :: R_xpoint(:)
  REAL*8, TARGET                 :: Z_xpoint(:)
  REAL*8, TARGET                 :: psi_xpoint(:)
  LOGICAL, TARGET                :: gmres
  INTEGER, TARGET                :: xcase2
  INTEGER                        :: n_cpu
  INTEGER                        :: i_tor(n_cpu)
  INTEGER                        :: column_number
  INTEGER                        :: mpi_comm_n
  INTEGER, TARGET                :: MPI_COMM_TRANS, my_id_trans, n_cpu_trans
  LOGICAL, TARGET                :: Solve_only
  ! local variables
  TYPE (type_element)            :: element
  INTEGER                        :: n_local_elms_min
  INTEGER                        :: n_local_elms_max
  INTEGER                        :: n_local_elms_sum
  INTEGER                        :: coefnbr, coefnbr_max, coefnbr_min
  INTEGER                        :: coefnbr_prod, coefnbr_prod_max
  INTEGER                        :: coefnbr_prod_min
  REAL(KIND=MURGE_COEF_KIND)     :: coefmtx(n_tor*n_var*n_tor*n_var)
  INTEGER, TARGET                :: elem_block_size
  INTEGER, TARGET                :: thread_nbr
  INTEGER                        :: elem_size
  INTEGER                        :: ELM_INDEX
  REAL*8,  POINTER               :: rhs_loc(:), rhs_loc_thread(:,:), rhs_loc_bc(:), rhs_loc_bc_recv(:)
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
  CHARACTER(LEN=128) :: fname
  INTEGER :: MPI_KEEP_NON_NAN
#ifdef MURGE_USE_SEQUENCE
  INTEGER, TARGET   :: mode
  INTEGER :: seq_coefnbr, cpu
#endif
#ifdef MURGE_USE_SEQUENCE_HARM
  INTEGER :: seq_coefnbr_harm
#endif
  call tr_print_memsize("BeginConstM")

  ! elapsed time
  CALL SYSTEM_CLOCK(count_rate=nb_periodes_sec,count_max=nb_periodes_max)
  ! Timing call
  CALL r3_info_begin (r3_info_index_0, 'construct_matrix_murge')

  IF (my_id .EQ. 0) THEN
     WRITE(*,*) '****************************************'
     WRITE(*,*) '*  construct matrix MURGE              *'
     WRITE(*,*) '****************************************'
     !write(*,*) ' solve_only : ', solve_only
  END IF

  ! --- Memory tracking
  call tr_print_memsize("DebConstM")
  
!$omp PARALLEL shared(thread_nbr)
!$OMP master
  thread_nbr = omp_get_num_threads()
!$omp end master
!$omp end PARALLEL

#ifdef MURGE_USE_SEQUENCE
  if (.not. gmres) then
     write(*,*) 'FATAL : MURGE_USE_SEQUENCE only works with GMRES=.true.'
     call MPI_FINALIZE(IERR)
     stop
  end if
#endif
  IF (my_id .EQ. 0) PRINT *, "THREAD_NBR", thread_nbr

  elem_block_size = murge_elem_block_size

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
#ifndef MURGE_USE_SEQUENCE
#ifdef PROD_MATRICES_STORAGE
  CALL tr_ALLOCATEp(PROD_MATRICES,1,(n_tor*n_var)**2,                          &
       1, (n_vertex_max*(n_order+1))**2*(elem_block_size/thread_nbr+1),        &
       1, thread_nbr, "PROD_MATRICES",CAT_DMATRIX)
#endif
#endif
  CALL tr_ALLOCATEp(PROD_COLROW, 1, 2,                                         &
       1, (n_vertex_max*(n_order+1))**2*(elem_block_size/thread_nbr+1),        &
       1, thread_nbr, "PROD_COLROW",CAT_DMATRIX)
  CALL tr_ALLOCATEp(matrix_nbr,1,thread_nbr, "matrix_nbr",CAT_DMATRIX)
  CALL tr_ALLOCATEp(matrix_nbr_rcv,1,thread_nbr,1,(n_tor+1)/2,"matrix_nbr_rcv",&
       &            CAT_DMATRIX)


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

  step = murge_assembly_step

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
     CALL MPI_Reduce(n_local_elms, n_local_elms_max, 1, MPI_INTEGER, MPI_MAX,  &
          &          0, MPI_COMM_WORLD, ierr)
     CALL MPI_Reduce(n_local_elms, n_local_elms_min, 1, MPI_INTEGER, MPI_MIN,  &
          &          0, MPI_COMM_WORLD, ierr)
     CALL MPI_Reduce(n_local_elms, n_local_elms_sum, 1, MPI_INTEGER, MPI_SUM,  &
          &          0, MPI_COMM_WORLD, ierr)
     CALL MPI_Reduce(coefnbr_prod, coefnbr_prod_max, 1, MPI_INTEGER, MPI_MAX,  &
          &          0, MPI_COMM_WORLD, ierr)
     CALL MPI_Reduce(coefnbr_prod, coefnbr_prod_min, 1, MPI_INTEGER, MPI_MIN,  &
          &          0, MPI_COMM_WORLD, ierr)
     IF (my_id .EQ. 0) THEN
        WRITE (*,"(A30,I10,A30)") ":: Murge Assembly phase :: ",               &
             n_local_elms_max/n_cpu_trans, " maximum local elements"
        WRITE (*,"(A30,I10,A30)") ":: Murge Assembly phase :: ",               &
             n_local_elms_min/n_cpu_trans, " minimum local elements"
        WRITE (*,"(A30,I10,A30)") ":: Murge Assembly phase :: ",               &
             n_local_elms_sum/n_cpu_trans, " total elements computed"
        WRITE (*,"(A30,I10,A30)") ":: Murge Assembly phase :: ",               &
             element_list%n_elements, " elements"
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
     IF (ierr .ne. MURGE_SUCCESS) then
        PRINT *, "ERROR in MURGE_MATRIXRESET :", ierr
        call abort()
     END IF
#  ifndef MURGE_USE_SEQUENCE_HARM
#    ifndef MURGE_USE_DUPLICATE_ELEMENT
     CALL MURGE_ASSEMBLYBEGIN(murge_id, murge_global_n, coefnbr,               &
          &                   MURGE_ASSEMBLY_ADD,                              &
          &                   MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_FOOL,         &
          &                   murge_sym, ierr)
#    else
     CALL MURGE_ASSEMBLYBEGIN(murge_id, murge_global_n, coefnbr,               &
          &                   MURGE_ASSEMBLY_ADD,                              &
          &                   MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_RESPECT,      &
          &                   murge_sym, ierr)
#    endif
     IF (ierr .ne. MURGE_SUCCESS) then
        PRINT *, "ERROR in MURGE_ASSEMBLYBEGIN :", ierr
        call abort()
     END IF
#  endif
  END IF
  IF (gmres) THEN
     CALL MURGE_MATRIXRESET(murge_id_prod, ierr)
     IF (ierr .ne. MURGE_SUCCESS) then
        PRINT *, "ERROR in MURGE_MATRIXRESET :", ierr
        call abort()
     END IF
     print *, "murge_global_n, murge_global_n_prod", murge_global_n, murge_global_n_prod
#  ifndef MURGE_USE_SEQUENCE
#    ifdef MURGE_PROD_NO_COMM
     CALL MURGE_ASSEMBLYBEGIN(murge_id_prod, murge_global_n_prod, coefnbr_prod,      &
          &                   MURGE_ASSEMBLY_ADD,                              &
          &                   MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_RESPECT,      &
          &                   murge_sym, ierr)
#    else
     CALL MURGE_ASSEMBLYBEGIN(murge_id_prod, murge_global_n_prod, coefnbr_prod,      &
          &                   MURGE_ASSEMBLY_ADD,                              &
          &                   MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_FOOL,         &
          &                   murge_sym, ierr)
#    endif
     IF (ierr .ne. MURGE_SUCCESS) then
        PRINT *, "ERROR in MURGE_ASSEMBLYBEGIN :", ierr
        call abort()
     END IF

#  endif
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
#ifdef MURGE_USE_SEQUENCE
  if (.not. ALLOCATED(murge_assembly_first_entry) ) then
     ALLOCATE(murge_assembly_first_entry(thread_nbr))
     CALL tr_register_mem(sizeof(murge_assembly_first_entry),                  &
          "murge_assembly_first_entry",CAT_DMATRIX)
     murge_assembly_first_entry = -1
  end if
#endif
#ifdef MURGE_USE_SEQUENCE_HARM
  if (.not. ALLOCATED(murge_assembly_first_entry_harm) ) then
     ALLOCATE(murge_assembly_first_entry_harm(thread_nbr))
     CALL tr_register_mem(sizeof(murge_assembly_first_entry_harm),             &
          "murge_assembly_first_entry_harm",CAT_DMATRIX)
     murge_assembly_first_entry = -1
  end if
#endif
  DO iter = 1, thread_nbr
     !print *, "iter", iter
     datas(iter)%thread_num          = iter
#ifdef MURGE_USE_SEQUENCE
     datas(iter)%cnt_entries         = 0
     datas(iter)%first_entry         = murge_assembly_first_entry(iter)
     datas(iter)%mode                => mode
     datas(iter)%VALS                => MURGE_VALS
     datas(iter)%ROWS                => MURGE_ROWS
     datas(iter)%COLS                => MURGE_COLS
#endif
#ifdef MURGE_USE_SEQUENCE_HARM
     datas(iter)%VALS_HARM           => MURGE_VALS_HARM
     datas(iter)%ROWS_HARM           => MURGE_ROWS_HARM
     datas(iter)%COLS_HARM           => MURGE_COLS_HARM
     datas(iter)%cnt_entries_harm    = 0
     datas(iter)%first_entry_harm    = murge_assembly_first_entry_harm(iter)
#endif
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
#ifndef MURGE_USE_SEQUENCE
#ifdef PROD_MATRICES_STORAGE
     datas(iter)%PROD_MATRICES       => PROD_MATRICES
#endif
#endif
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
     datas(iter)%R_xpoint            => R_xpoint
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
#ifdef MURGE_USE_SEQUENCE
  IF (murge_need_rebuild_sequence(1)) THEN
     murge_need_rebuild_sequence(1) = .FALSE.
     IF (ASSOCIATED(MURGE_VALS)) THEN
        DEALLOCATE(MURGE_VALS)
        NULLIFY(MURGE_VALS)
     END IF
     mode = 1
#endif
     !$OMP parallel private(iter,ret) default(shared)
     iter = 1+omp_get_thread_num()
     ret = LOOP(datas(iter))
     !$OMP barrier
     !$OMP end parallel

#ifdef MURGE_USE_SEQUENCE
     datas(1)%first_entry = 1
     DO iter = 2, thread_nbr
        datas(iter)%first_entry = datas(iter-1)%first_entry +                  &
             datas(iter-1)%cnt_entries
     END DO
     seq_coefnbr = datas(thread_nbr)%first_entry +                             &
          &        datas(thread_nbr)%cnt_entries - 1

     ALLOCATE(MURGE_ROWS(seq_coefnbr))
     ALLOCATE(MURGE_COLS(seq_coefnbr))
     DO iter = 1, thread_nbr
        datas(iter)%ROWS                => MURGE_ROWS
        datas(iter)%COLS                => MURGE_COLS
        murge_assembly_first_entry(iter) = datas(iter)%first_entry
     END DO
#  ifdef MURGE_USE_SEQUENCE_HARM
     datas(1)%first_entry_harm = 1
     DO iter = 2, thread_nbr
        datas(iter)%first_entry_harm = datas(iter-1)%first_entry_harm +   &
             datas(iter-1)%cnt_entries_harm
     END DO
     seq_coefnbr_harm = datas(thread_nbr)%first_entry_harm +                   &
          &        datas(thread_nbr)%cnt_entries_harm - 1

     ALLOCATE(MURGE_ROWS_HARM(seq_coefnbr_harm))
     ALLOCATE(MURGE_COLS_HARM(seq_coefnbr_harm))
     !print *, seq_coefnbr_harm
     DO iter = 1, thread_nbr
        datas(iter)%ROWS_HARM                => MURGE_ROWS_HARM
        datas(iter)%COLS_HARM                => MURGE_COLS_HARM
        murge_assembly_first_entry_harm(iter) = datas(iter)%first_entry_harm
     END DO
#  endif
     mode = 2
     !$OMP parallel private(iter,ret) default(shared)
     iter = 1+omp_get_thread_num()
     ret = LOOP(datas(iter))
     !$OMP barrier
     !$OMP end parallel

#  ifdef MURGE_PROD_NO_COMM
     CALL MURGE_AssemblySetSequence(murge_id_prod, seq_coefnbr,                &
          &                         MURGE_ROWS, MURGE_COLS,                    &
          &                         MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_ADD,    &
          &                         MURGE_ASSEMBLY_RESPECT, MURGE_BOOLEAN_TRUE,&
          &                         murge_sequence_id(1), ierr)
#  else
     CALL MURGE_AssemblySetSequence(murge_id_prod, seq_coefnbr,                &
          &                         MURGE_ROWS, MURGE_COLS,                    &
          &                         MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_ADD,    &
          &                         MURGE_ASSEMBLY_FOOL, MURGE_BOOLEAN_TRUE,   &
          &                         murge_sequence_id(1), ierr)
#  endif
#  ifdef MURGE_USE_SEQUENCE_HARM
#    ifndef MURGE_USE_DUPLICATE_ELEMENT
     CALL MURGE_AssemblySetSequence(murge_id, seq_coefnbr_harm,                &
          &                         MURGE_ROWS_HARM, MURGE_COLS_HARM,          &
          &                         MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_ADD,    &
          &                         MURGE_ASSEMBLY_FOOL, MURGE_BOOLEAN_TRUE,&
          &                         murge_sequence_id_harm(1), ierr)
#    else
     CALL MURGE_AssemblySetSequence(murge_id, seq_coefnbr_harm,                &
          &                         MURGE_ROWS_HARM, MURGE_COLS_HARM,          &
          &                         MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_ADD,    &
          &                         MURGE_ASSEMBLY_RESPECT, MURGE_BOOLEAN_TRUE,&
          &                         murge_sequence_id_harm(1), ierr)
#    endif
#  endif
     DEALLOCATE(MURGE_ROWS, MURGE_COLS)
     NULLIFY(MURGE_ROWS, MURGE_COLS)
     ALLOCATE(MURGE_VALS(seq_coefnbr*(n_tor*n_var)**2))
#  ifdef MURGE_USE_SEQUENCE_HARM
     DEALLOCATE(MURGE_ROWS_HARM, MURGE_COLS_HARM)
     NULLIFY(MURGE_ROWS_HARM, MURGE_COLS_HARM)
     ALLOCATE(MURGE_VALS_HARM(seq_coefnbr_harm*(my_harm_size**2)))
#  endif
  END IF
  mode = 3
#  ifdef MURGE_USE_SEQUENCE_HARM
  MURGE_VALS_HARM = 0.
#  endif
  MURGE_VALS = 0.
  DO iter = 1, thread_nbr
     datas(iter)%VALS                => MURGE_VALS
#  ifdef MURGE_USE_SEQUENCE_HARM
     datas(iter)%VALS_HARM           => MURGE_VALS_HARM
#  endif
  END DO
  !$OMP parallel private(iter,ret) default(shared)
  iter = 1+omp_get_thread_num()
  ret = LOOP(datas(iter))
  !$OMP barrier
  !$OMP end parallel
  CALL MURGE_AssemblyUseSequence(murge_id_prod, murge_sequence_id(1),          &
       &                         MURGE_VALS, ierr)
#  ifdef MURGE_USE_SEQUENCE_HARM
  CALL MURGE_AssemblyUseSequence(murge_id, murge_sequence_id_harm(1),     &
       &                         MURGE_VALS_HARM, ierr)
#  endif
#endif
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
#ifndef MURGE_USE_SEQUENCE
#ifdef PROD_MATRICES_STORAGE
  CALL tr_deallocatep(prod_matrices,"prod_matrices",CAT_DMATRIX)
#endif
#endif
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

#ifndef MURGE_USE_SEQUENCE_HARM
  IF (.NOT. gmres .OR. .NOT. solve_only) THEN
#  ifdef USE_MURGE
     CALL MURGE_ASSEMBLYEND(murge_id, ierr)
#  else
     PRINT *, "Binary built without murge"
     CALL abort()
#  endif
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
#endif

#ifndef MURGE_USE_SEQUENCE
  IF (gmres) THEN
#  ifdef USE_MURGE
     CALL MURGE_ASSEMBLYEND(murge_id_prod, ierr)
#  else
     PRINT *, "Binary built without murge"
     CALL abort()
#  endif
     IF (ierr /= MURGE_SUCCESS) THEN
        WRITE (*,*) my_id, "::: error in assemblyend..",                       &
             cnt*(n_tor*n_var)**2, cnt
        CALL abort()
     END IF
  END IF
#endif
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

  ! --- Add vacuum response (boundary integral) for free boundary computations
  if (freeboundary) then
     print *, "NOT IMPLEMENTED CORRECTLY"
     CALL ABORT()
     !call global_matrix_structure_vacuum(node_list, bnd_node_list, index_min, index_max) !###TODO### move somewhere else
     !call vacuum_boundary_integral(my_id, bnd_node_list, node_list, bnd_elm_list,                   &
     ! freeboundary_equil, resistive_wall, index_min, index_max, rhs_loc, tstep, index_now)
  end if

  CALL SYSTEM_CLOCK(count=t0)

#ifdef NORMTRACE
! For debugging purpose
  CALL MPI_Reduce(RHS_loc, RHS_glob, ndof_glob, MPI_DOUBLE_PRECISION, MPI_SUM,&
       &          0, MPI_COMM_WORLD, ierr)

  call tr_locvnorms("cm_Rhs",RHS_glob,ndof_glob)
  if (my_id .eq. 0) then
     !print *, my_id, "rhs_glob(1)", rhs_glob(1), index_now
     write(fname,'(A,I6.6)')"rhs",index_now
     call tr_vdump(fname,RHS_glob,ndof_glob)
  end if
  !print *, my_id, "rhs_loc(1)", rhs_loc(1), index_now
  !print *, my_id, "rhs_gloab(1)", rhs_glob(1), index_now
#endif
  CALL tr_deallocatep(RHS_loc_thread,"RHS_loc_thread",CAT_DMATRIX)
#ifndef MURGE_USE_DUPLICATE_ELEMENT
  ! One element can appear on two MPI process, thus the boundary conditions can be set several time, we have to take care of that.
  CALL tr_allocatep(rhs_loc_bc,1,ndof_glob,"rhs_loc_bc",CAT_DMATRIX)
  CALL tr_allocatep(rhs_loc_bc_recv,1,ndof_glob,"rhs_loc_bc_recv",CAT_DMATRIX)
  ! right-hand-side is filled with NaN before receiving boundary conditions.
  ! This way we can detect which values are set in boundary conditions.
  do i = 1, ndof_glob
     rhs_loc_bc(i)      = set_nan()
  end do
  do i = 1, ndof_glob
     rhs_loc_bc_recv(i)      = set_nan()
  end do
  CALL boundary_conditions(my_id, node_list, element_list, bnd_node_list, local_elms, n_local_elms,&
    0, 0, rhs_loc_bc, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint,         &
    psi_xpoint, gmres, solve_only)
  ! Reduce only non NaN values from boundary conditions.
  ! KEEP_NON_NAN will keep either NaN if only NaN were entered, zeros if only NaN and zeros, or last non-zero value otherwise
  CALL MPI_Op_create(keep_non_nan, .true., MPI_KEEP_NON_NAN, ierr)
  CALL MPI_Allreduce(RHS_loc_bc, RHS_loc_bc_recv, ndof_glob, MPI_DOUBLE_PRECISION, MPI_KEEP_NON_NAN, MPI_COMM_WORLD,  &
       ierr)
  CALL tr_deallocatep(RHS_loc_bc,"rhs_loc_bc",CAT_DMATRIX)
  CALL MPI_Op_free(MPI_KEEP_NON_NAN, ierr)
  ! Copy non NaN values from boundary conditions into rhs_loc
  if (my_id .eq. 0) then
     DO i = 1, ndof_glob
        if (.not. IS_NAN(rhs_loc_bc_recv(i))) then
           rhs_loc(i) = rhs_loc_bc_recv(i)
        end if
     END DO
  else
     DO i = 1, ndof_glob
        if (.not. IS_NAN(rhs_loc_bc_recv(i))) then
           rhs_loc(i) = 0.d0
        end if
     END DO
  end if
  CALL tr_deallocatep(RHS_loc_bc_recv,"rhs_loc_bc_recv",CAT_DMATRIX)
#else
  CALL boundary_conditions(my_id, node_list, element_list, bnd_node_list, local_elms, n_local_elms,&
       0, 0, rhs_loc,    xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint,         &
       psi_xpoint, gmres, solve_only)
#endif

  CALL MPI_Reduce(RHS_loc, RHS_glob, ndof_glob, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,  &
    ierr)
  call tr_deallocatep(RHS_loc,"RHS_loc",CAT_DMATRIX)

  ! --- For debugging purpose
  if (my_id .eq. 0) then
     write(fname,'(A,I6.6)')"rhsbc",index_now
     call tr_vdump(fname,RHS_glob,ndof_glob)
  end if

  ! --- Memory tracking
  call tr_locvnorms("cm_BCRhs",RHS_glob,ndof_glob)
  call tr_debug_writei("ndof_glob",ndof_glob)


  ! -- Timer
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
  ! ICI c'est OK
  !CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)

  ! --- Construct per harmonic RHS
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

  ! --- Apply column scaling
  IF (.NOT. gmres .OR. .NOT. solve_only) THEN
     IF(ALLOCATED(column_scaling))                                             &
          CALL tr_deallocate(column_scaling,"column_scaling",CAT_DMATRIX)
     CALL tr_allocate(column_scaling,1,column_number,"column_scaling",         &
          &           CAT_DMATRIX)
     column_scaling = 1.d-20
#ifdef USE_MURGE
     CALL MURGE_GetGlobalNorm(murge_id, column_scaling, -1,                    &
          &                   MURGE_NORM_MAX_COL, ierr)
     if (my_id .eq. 0) then
        write(fname,'(A,I6.6)')"column_scaling",index_now
        call tr_vdump(fname,column_scaling,column_number)
     end if
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

  ! --- Timing
  call r3_info_end(r3_info_index_0)
  call tr_print_memsize("EndConstM")

  RETURN
END SUBROUTINE construct_matrix_murge
END MODULE construct_matrix_murge_mod
