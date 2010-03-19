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
  IF (.not. gmres) THEN
     column_number = ndof_glob
  ELSE
     IF (i_tor(my_id+1) == 1) THEN
        column_number = ndof_glob/n_tor
     ELSE
        column_number = ndof_glob*2/n_tor
     END IF    
  END IF

  ALLOCATE(rhs_glob(column_number))
  IF (.NOT. ALLOCATED(rhs_loc))   ALLOCATE(rhs_loc(column_number))
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
  IF (.not. gmres) THEN
     coefnbr = coefnbr * (n_vertex_max)*(n_order+1)* (n_var * n_tor)* (n_var * n_tor)
  ELSE
     IF (i_tor(my_id+1) == 1) THEN 
        coefnbr = coefnbr * (n_vertex_max)*(n_order+1)* (n_var )* (n_var )
     ELSE
        coefnbr = coefnbr * (n_vertex_max)*(n_order+1)* (n_var )* (n_var ) * 4
     END IF
  END IF


  WRITE (*,*) ":: Murge Assembly phase :: ", coefnbr, " entries on processor ", my_id 
  cnt = 0
  CALL MURGE_ASSEMBLYBEGIN(id, coefnbr, MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_ADD, &
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
                 IF (.not. gmres) THEN
                    index_rhs = index_large_i+j
                 ELSE
                    IF (i_tor(my_id+1) == 1) THEN
                       index_rhs = index_large_i/n_tor+INT((j-1)/n_tor)+1
                    ELSE
                       index_rhs = index_large_i*2/n_tor+INT((j-1)/n_tor)*2+MOD(MOD(j-1, n_tor)+1,2)+1 
                    ENDIF
                 END IF

                 if (index_rhs > column_number .or. index_rhs < 1) then
                    write (*,*) "index_rhs", index_rhs
                    write (*,*) "column_number", column_number
                    write (*,*) "index_large_i", index_large_i
                    write (*,*) "j", j
                    stop
                 end if
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
                          IF (.not. gmres) THEN
                             ! Column index in the ELM matrix
                             index_mtx = l+(j-1)*(n_var * n_tor)
                             coefmtx(index_mtx) = ELM(index_kl,index_ij)
                          ELSE
                             IF (INT((MOD(j-1, n_tor)+1)/2)+1 == i_tor(my_id+1) .and. &
                                  INT((MOD(l-1, n_tor)+1)/2)+1 == i_tor(my_id+1)) THEN
                                IF (i_tor(my_id+1) == 1) THEN
                                   index_mtx = INT((l-1)/n_tor)+1+(INT((j-1)/n_tor))*(n_var) 
                                ELSE
                                   index_mtx = INT((l-1)/n_tor)+MOD(MOD(l-1,n_tor),2)+1+ &
                                        (INT((j-1)/n_tor)*2+ MOD(MOD(j-1,n_tor),2))*(n_var*2)
                                END IF
                                coefmtx(index_mtx) = ELM(index_ij,index_kl)
                             END IF
                          END IF
                       ENDDO
                          
                    ENDDO

                    IF (.not. gmres) THEN
                       cnt = cnt + n_var*n_var * n_tor * n_tor
                    ELSE
                       IF (i_tor(my_id+1) == 1) THEN
                          cnt = cnt + n_var*n_var
                       ELSE
                          cnt = cnt + n_var*n_var * 4
                       END IF
                    END IF

                    CALL MURGE_ASSEMBLYSETNODEVALUES(id, index_node2, index_node1, &
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
  write (*,*) my_id, " : MURGE_ASSEMBLYEND in"
  CALL MURGE_ASSEMBLYEND(id, ierr)
  write (*,*) my_id, " : MURGE_ASSEMBLYEND out"

  !----------------------- boundary conditions
  
  IF (.not. gmres) THEN
     CALL boundary_conditions_murge(my_id,node_list,element_list,local_elms,n_local_elms, &
          psi_bnd)
  END IF

  IF (.not. gmres) THEN
     write (*,*) MY_ID, " : Reduce..."
     CALL MPI_Allreduce(RHS_loc,RHS_glob,column_number,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_N,ierr)
  ELSE
     IF (i_tor(my_id+1) == 1) THEN
        write (*,*) my_id, " : Reduce...."
        CALL MPI_Allreduce(RHS_loc,RHS_glob,column_number,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_N, ierr)
     ELSE
        write (*,*) my_id, " : Reduce.."
        CALL MPI_AllReduce(RHS_loc,RHS_glob,column_number,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_N,ierr)
     END IF    
  END IF
  
  IF(ALLOCATED(column_scaling))   DEALLOCATE(column_scaling)
  ALLOCATE(column_scaling(column_number))

  CALL MURGE_GetGlobalNorm(id, column_scaling, -1, MURGE_NORM_MAX_COL, ierr)
  CALL MURGE_ApplyGlobalScaling(id, column_scaling, -1, MURGE_SCAL_COL, ierr)
  write(*,*) '******** end construct matrix murge **********'
  
  DEALLOCATE(RHS_loc)

  RETURN
END SUBROUTINE construct_matrix_murge
