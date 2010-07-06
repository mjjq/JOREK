!*******************************************************************************
!* Subroutine: construct_matrix_murge                                          *
!*******************************************************************************
!*                                                                             *
!* Construct the matrix for the resolution using Murge interface (PaStiX).     *
!*                                                                             *
!* Parameters:                                                                 *
!*   my_id        - Identifier of the node in MPI_COMM_WORLD                   *
!*   node_list    - List of nodes                                              *
!*   element_list - List of all elements                                       *
!*   local_elms   - List of local elements                                     *
!*   n_local_elms - Number of local elements                                   *
!*   xpoint2      - ???                                                        *
!*   psi_axis     - ???                                                        *
!*   psi_bnd      - ???                                                        *
!*   Z_xpoint     - ???                                                        *
!*   gmres        - Solve method (.true. for gmres, .false for 'direct')       *
!*   i_tor        - Tor number                                                 *
!*   n_cpu        - Number of cpus                                             *
!*   mpi_comm_n   - Solver MPI communicator                                    *
!*                                                                             *
!* Authors:                                                                    *
!*   Xavier Lacoste - xavier.lacoste@inria.fr                                  *
!*                                                                             *
!*******************************************************************************
SUBROUTINE construct_matrix_murge(my_id,node_list,element_list, local_elms, &
     n_local_elms, xpoint2,psi_axis,psi_bnd,Z_xpoint, gmres, i_tor, n_cpu, &
     mpi_comm_n)
  !---------------------------------------------------------------
  ! collect the element matrices into one large sparse matrix
  ! in coordinate format
  !---------------------------------------------------------------
  USE data_structure
  USE global_distributed_matrix
  USE phys_module
  USE murge_module
  USE mumps_module

  IMPLICIT NONE
  INCLUDE 'mpif.h'
  
  ! Subroutine parameters:
  INTEGER                        :: my_id
  TYPE (type_node_list)          :: node_list
  TYPE (type_element_list)       :: element_list
  INTEGER                        :: local_elms(*)
  INTEGER                        :: n_local_elms
  LOGICAL                        :: xpoint2
  REAL*8                         :: psi_axis
  REAL*8                         :: psi_bnd
  REAL*8                         :: Z_xpoint
  LOGICAL                        :: gmres
  INTEGER                        :: i_tor(n_cpu)
  INTEGER                        :: n_cpu
  INTEGER                        :: column_number
  INTEGER                        :: mpi_comm_n
  INTEGER                        :: new_col_mat_elem, new_row_mat_elem

  ! local variables
  TYPE (type_element)      :: element
  TYPE (type_node)         :: nodes(n_vertex_max)
  REAL*8, ALLOCATABLE :: rhs_loc(:)
  INTEGER :: coefnbr
  REAL(KIND=MURGE_COEF_KIND) :: coefmtx(n_tor*n_var*n_tor*n_var)
  REAL*8  :: ELM(n_tor*n_vertex_max*(n_order+1)*n_var,n_tor*n_vertex_max*(n_order+1)*n_var)

  REAL*8  :: RHS(n_tor*n_vertex_max*(n_order+1)*n_var)
  INTEGER :: i, ife, iv, inode, inode1, knode, j, k, l, index_ij, index_kl
  INTEGER :: index_large_i, index_large_k, index_node1, index_node2, i_order, k_order, ielm, ierr

  LOGICAL :: is_local
  INTEGER :: index_rhs
  INTEGER :: index_mtx
  INTEGER :: cnt
  WRITE(*,*) '****************************************'
  WRITE(*,*) '*  construct matrix MURGE              *'
  WRITE(*,*) '****************************************'
  WRITE(*,*) ' n_elements (local)       : ',my_id,n_local_elms

  IF (ALLOCATED(rhs_glob))        DEALLOCATE(rhs_glob)
  IF (.NOT. gmres) THEN
     column_number = ndof_glob
  ELSE
     IF (murge_id == 0) THEN
        column_number = ndof_glob/n_tor
     ELSE
        column_number = ndof_glob*2/n_tor
     END IF    
  END IF

  ALLOCATE(rhs_glob(ndof_glob))
  IF (.NOT. ALLOCATED(rhs_loc))   ALLOCATE(rhs_loc(ndof_glob))
  RHS_glob = 0.d0
  RHS_loc  = 0.d0

  !
  ! Count coefnbr 
  !
  coefnbr = 0
  DO ife =1, n_local_elms
     
     ielm = local_elms(ife)
     
     element = element_list%element(ielm)
     
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
     
  ENDDO
  IF (.NOT. gmres) THEN
     coefnbr = coefnbr * (n_vertex_max)*(n_order+1)* (n_var * n_tor)* (n_var * n_tor)
  ELSE
     IF (murge_id == 0) THEN 
        coefnbr = coefnbr * (n_vertex_max)*(n_order+1)* (n_var )* (n_var )
     ELSE
        coefnbr = coefnbr * (n_vertex_max)*(n_order+1)* (n_var )* (n_var ) * 4
     END IF
  END IF


  WRITE (*,*) ":: Murge Assembly phase :: ", coefnbr, " entries on processor ", my_id 
  cnt = 0
  CALL MURGE_ASSEMBLYBEGIN(murge_id, coefnbr, MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_ADD, &
       MURGE_ASSEMBLY_FOOL, murge_sym, ierr)

  DO ife =1, n_local_elms

     ielm = local_elms(ife)

     element = element_list%element(ielm)

     DO iv = 1, n_vertex_max

        inode     = element%vertex(iv)
        nodes(iv) = node_list%node(inode)

     ENDDO

     IF (n_tor .GT. 3) THEN
        CALL element_matrix_fft(element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)      ! use fft for toroidal integration
     ELSE
        CALL element_matrix(element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)           ! use direct integration
     ENDIF

     DO i=1,n_vertex_max

        inode1         = element%vertex(i)

        DO i_order = 1, n_order+1
           
           ! index_node1 is the column
           index_node1 = node_list%node(inode1)%index(i_order)

           index_large_i = n_tor * n_var * (index_node1 - 1)

           CALL vertex_is_local(index_node1, is_local)
           IF (is_local) THEN
              
              ! Set RHS member
              DO j = 1, n_var * n_tor
                    
                 index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j   ! index in the ELM matrix
                    index_rhs = index_large_i+j

                 IF (index_rhs > ndof_glob .OR. index_rhs < 1) THEN
                    WRITE (*,*) "index_rhs", index_rhs
                    WRITE (*,*) "column_number", column_number
                    WRITE (*,*) "index_large_i", index_large_i
                    WRITE (*,*) "j", j
                    STOP
                 END IF
                 rhs_loc(index_rhs) = rhs_loc(index_rhs) + RHS(index_ij)
              END DO
              ! Build nodes Matrices
              DO k=1,n_vertex_max
                 
                 knode         = element%vertex(k)
                 
                 DO k_order = 1, n_order+1

                    index_node2 = node_list%node(knode)%index(k_order)
                    
                    index_large_k = n_tor * n_var * (index_node2 - 1)

                    coefmtx = 0
                    DO j = 1, n_var * n_tor
                       ! Row index in the ELM matrix
                       index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j   
                       
                       DO l = 1, n_var * n_tor
                          ! BUILD node Matrix
                          index_kl  = n_tor * n_var * (n_order+1) * (k-1) + n_tor * n_var * (k_order-1) + l   
                          IF (.NOT. gmres) THEN
                             ! Column index in the ELM matrix
                             index_mtx = l+(j-1)*(n_var * n_tor)
                             coefmtx(index_mtx) = ELM(index_kl,index_ij)
                          ELSE
                             IF (INT((MOD(j-1, n_tor)+1)/2) == murge_id .AND. &
                                  INT((MOD(l-1, n_tor)+1)/2) == murge_id) THEN
                                IF (murge_id == 0) THEN
                                   new_col_mat_elem = INT((j-1)/n_tor)
                                   new_row_mat_elem = INT((l-1)/n_tor)
                                   index_mtx = new_row_mat_elem+1+new_col_mat_elem*(n_var) 
                                ELSE
                                   new_col_mat_elem = INT((j-1)/n_tor)*2 + MOD(MOD(j-1, n_tor)+1,2)
                                   new_row_mat_elem = INT((l-1)/n_tor)*2 + MOD(MOD(l-1, n_tor)+1,2)
                                   index_mtx = new_row_mat_elem+1+new_col_mat_elem*(n_var*2) 
                                END IF
                                coefmtx(index_mtx) = ELM(index_kl,index_ij)
                             END IF
                          END IF
                       ENDDO
                          
                    ENDDO

                    IF (.NOT. gmres) THEN
                       cnt = cnt + n_var*n_var * n_tor * n_tor
                    ELSE
                       IF (murge_id == 0) THEN
                          cnt = cnt + n_var*n_var
                       ELSE
                          cnt = cnt + n_var*n_var * 4
                       END IF
                    END IF

                    CALL MURGE_ASSEMBLYSETNODEVALUES(murge_id, index_node2, index_node1, &
                         coefmtx, ierr)
                    IF (ierr /= MURGE_SUCCESS) THEN
                       WRITE (*,*) my_id, ":::", &
                            "I", index_node2, &
                            "J", index_node1, cnt
                       STOP
                    END IF

                 ENDDO
              ENDDO
           ENDIF

        ENDDO
     ENDDO

  ENDDO
  WRITE (*,*) my_id, " : MURGE_ASSEMBLYEND in"
  CALL MURGE_ASSEMBLYEND(murge_id, ierr)
  WRITE (*,*) my_id, " : MURGE_ASSEMBLYEND out"

  !----------------------- boundary conditions
  
  CALL boundary_conditions_murge(my_id,node_list,element_list,local_elms,n_local_elms, &
       psi_bnd)
  
  IF (.NOT. gmres) THEN
     WRITE (*,*) MY_ID, " : Reduce..."
     CALL MPI_Allreduce(RHS_loc,RHS_glob,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_N,ierr)

  ELSE
     IF (murge_id == 0) THEN
        WRITE (*,*) my_id, " : Reduce...."
        CALL MPI_Allreduce(RHS_loc,RHS_glob,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_N, ierr)
     ELSE
        WRITE (*,*) my_id, " : Reduce.."
        CALL MPI_AllReduce(RHS_loc,RHS_glob,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_N,ierr)
     END IF

     if (associated(mumps_par%rhs)) deallocate(mumps_par%rhs)
     allocate(mumps_par%rhs(column_number))
     mumps_par%n = column_number
     mumps_par%rhs = 0.d0

     if (murge_id .eq. 0 ) then
        
        mumps_par%rhs(1:column_number) = rhs_glob(1:ndof_glob:n_tor)
        
     else
           
        mumps_par%rhs(1:column_number:2) = rhs_glob(2*murge_id:ndof_glob:n_tor)
        mumps_par%rhs(2:column_number:2) = rhs_glob(2*murge_id+1:ndof_glob:n_tor)
        
     endif
        
     
  END IF
  
  IF (.NOT. gmres) THEN
     IF(ALLOCATED(column_scaling))   DEALLOCATE(column_scaling)
     ALLOCATE(column_scaling(column_number))
     CALL MURGE_GetGlobalNorm(murge_id, column_scaling, -1, MURGE_NORM_MAX_COL, ierr)
     CALL MURGE_ApplyGlobalScaling(murge_id, column_scaling, -1, MURGE_SCAL_COL, ierr)
  END IF
  WRITE(*,*) '******** end construct matrix murge **********'
  
  DEALLOCATE(RHS_loc)

  RETURN
END SUBROUTINE construct_matrix_murge
