SUBROUTINE construct_matrix_murge(my_id,node_list,element_list, local_elms, &
     n_local_elms, xpoint2,psi_axis,psi_bnd,Z_xpoint, method, i_tor, n_cpu, mpi_comm_n)
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

  ! Parameters
  INTEGER                        :: my_id
  TYPE (type_node_list)          :: node_list
  TYPE (type_element_list)       :: element_list
  INTEGER                        :: local_elms(*)
  INTEGER                        :: n_local_elms
  REAL*8                         :: psi_axis
  REAL*8                         :: psi_bnd
  REAL*8                         :: Z_xpoint
  CHARACTER*8                    :: method
  INTEGER                        :: i_tor(n_cpu)
  INTEGER                        :: n_cpu
  INTEGER                        :: mpi_comm_n

  ! local variables
  TYPE (type_element)      :: element
  TYPE (type_node)         :: nodes(n_vertex_max)
  REAL*8, ALLOCATABLE :: rhs_loc(:)
  INTEGER :: coefnbr
  REAL(KIND=MURGE_COEF_KIND) :: coefmtx(n_tor*n_var*n_tor*n_var)
  REAL(KIND=MURGE_COEF_KIND) :: coef
  REAL*8  :: ELM(n_tor*n_vertex_max*(n_order+1)*n_var,n_tor*n_vertex_max*(n_order+1)*n_var)

  REAL*8  :: RHS(n_tor*n_vertex_max*(n_order+1)*n_var)
  REAL*8  :: ELM2(n_tor*n_vertex_max*(n_order+1)*n_var,n_tor*n_vertex_max*(n_order+1)*n_var)
  REAL*8  :: RHS2(n_tor*n_vertex_max*(n_order+1)*n_var)
  REAL*8  :: zbig, T0, Vpar0, bigR, dT0_ds, dVpar0_ds, dBigR_ds
  REAL*8  :: R_s, R_t, Z_s, Z_t, ps0_s, ps0_t, ps0_x, ps0_y, direction, xjac
  REAL*8  :: Vpar0_pol_R, Vpar0_pol_Z, Vpol_R, Vpol_Z, znormal_R, znormal_Z
  REAL*8  :: Vpar0_perp, Vpol_perp, Btot, cs_fraction, ratio
  REAL*8  :: grad_s, grad_psi, u0_s, u0_t, u0_x, u0_y
  INTEGER :: i, in, ife, iv, inode, inode1, inode2, knode, j, k, l, index_ij, index_kl
  INTEGER :: index_i, index_large_i, index_large_k, index_node, index_node1, index_node2, i_order, k_order, ic, ielm, ierr
  INTEGER :: ijA_position,ijA_position2, nz_AA2, n_AA2, ilarge2, kv, kT, ku, ilarge_vv, ilarge_vT, ilarge_vus
  LOGICAL :: xpoint2
  LOGICAL :: is_local
  INTEGER :: index_rhs
  INTEGER :: index_mtx
  INTEGER :: cnt
  WRITE(*,*) '****************************************'
  WRITE(*,*) '*  construct matrix MURGE              *'
  WRITE(*,*) '****************************************'
  WRITE(*,*) ' n_elements (local)       : ',my_id,n_local_elms


  IF (ALLOCATED(rhs_glob))        DEALLOCATE(rhs_glob)
  IF (method .EQ. "direct") THEN
     ALLOCATE(rhs_glob(ndof_glob))
     IF (.NOT. ALLOCATED(rhs_loc))   ALLOCATE(rhs_loc(ndof_glob))
  ELSE
     IF (i_tor(my_id+1) == 1) THEN
        ALLOCATE(rhs_glob(ndof_glob/n_tor))
        IF (.NOT. ALLOCATED(rhs_loc))   ALLOCATE(rhs_loc(ndof_glob/n_tor))
     ELSE
        ALLOCATE(rhs_glob(ndof_glob*2/n_tor))
        IF (.NOT. ALLOCATED(rhs_loc))   ALLOCATE(rhs_loc(ndof_glob*2/n_tor))
     END IF    
  END IF
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
  IF (method .EQ. "direct") THEN
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
              IF (method == "direct") THEN
                 DO j = 1, n_var * n_tor
                    
                    index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j   ! index in the ELM matrix
                    index_rhs = index_large_i+j
                    if (index_rhs > ndof_glob .or. index_rhs < 1) then
                       write (*,*) "index_rhs", index_rhs
                       write (*,*) "ndof_glob", ndof_glob*2/n_tor
                       write (*,*) "index_large_i", index_large_i
                       write (*,*) "j", j
                             stop
                          end if
                    rhs_loc(index_rhs) = rhs_loc(index_rhs) + RHS(index_ij)
                 END DO
              ELSE
                 IF (i_tor(my_id+1) == 1) THEN
                    DO j = 1, n_var * n_tor
                    
                       IF (INT((MOD(j-1, n_tor)+1)/2)+1 == i_tor(my_id+1)) THEN
                          ! index in the RHS
                          index_ij  = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j 
                          index_rhs = index_large_i/n_tor+INT((j-1)/n_tor)+1
                          if (index_rhs > ndof_glob/n_tor .or. index_rhs < 1) then
                             write (*,*) "index_rhs", index_rhs
                             write (*,*) "ndof_glob/n_tor", ndof_glob/n_tor
                             write (*,*) "index_large_i", index_large_i
                             write (*,*) "j", j
                             stop
                          end if
                          rhs_loc(index_rhs) = rhs_loc(index_rhs) + RHS(index_ij)
                       END IF
                       
                    END DO
                 ELSE
                    DO j = 1, n_var * n_tor
                    
                       IF (INT((MOD(j-1, n_tor)+1)/2)+1 == i_tor(my_id+1)) THEN
                          ! index in the RHS
                          index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j 
                          index_rhs = index_large_i*2/n_tor+INT((j-1)/n_tor)*2+MOD(MOD(j-1, n_tor)+1,2)+1 
                          if (index_rhs > ndof_glob*2/n_tor .or. index_rhs < 1) then
                             write (*,*) "index_rhs", index_rhs
                             write (*,*) "ndof_glob*2/n_tor", ndof_glob*2/n_tor
                             write (*,*) "index_large_i", index_large_i
                             write (*,*) "j", j, INT((j-1)/n_tor), MOD(MOD(j-1, n_tor)+1,2)+1 
                             stop
                          end if
                          rhs_loc(index_rhs) = rhs_loc(index_rhs) + RHS(index_ij)
                       END IF
                       
                    END DO
                 END IF
              END IF
              ! Build nodes Matrices
              DO k=1,n_vertex_max
                 
                 knode         = element%vertex(k)
                 
                 DO k_order = 1, n_order+1

                    index_node2 = node_list%node(knode)%index(k_order)
                    
                    index_large_k = n_tor * n_var * (index_node2 - 1)

                    coefmtx = 0
                    ! BUILD node Matrix
                    IF (method == "direct") THEN
                       DO j = 1, n_var * n_tor
                          ! Row index in the ELM matrix
                          index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j   
                          
                          DO l = 1, n_var * n_tor
                             ! Column index in the ELM matrix
                             index_kl  = n_tor * n_var * (n_order+1) * (k-1) + n_tor * n_var * (k_order-1) + l   
                             index_mtx = j+(l-1)*(n_var * n_tor)
                             coefmtx(index_mtx) = ELM(index_ij,index_kl)
                          ENDDO
                          
                       ENDDO
                       cnt = cnt + n_var*n_var * n_tor * n_tor
                       CALL MURGE_ASSEMBLYSETNODEVALUES(id, index_node2, index_node1, &
                            coefmtx, ierr)
                       IF (ierr /= MURGE_SUCCESS) THEN
                          WRITE (*,*) my_id, ":::", &
                               "I", index_node1, &
                               "J", index_node2, cnt
                          STOP
                       END IF
                    ELSE
                       IF (i_tor(my_id+1) == 1) THEN
                          
                          DO j = 1, n_var * n_tor
                             ! Row index in the ELM matrix
                             index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j   
                          
                             DO l = 1, n_var * n_tor
                                ! Column index in the ELM matrix
                                
                                IF (INT((MOD(j-1, n_tor)+1)/2)+1 == i_tor(my_id+1) .and. &
                                     INT((MOD(l-1, n_tor)+1)/2)+1 == i_tor(my_id+1)) THEN
                                   
                                   index_kl  = n_tor * n_var * (n_order+1) * (k-1) + &
                                        n_tor * n_var * (k_order-1) + l     
                                   index_mtx = INT((j-1)/n_tor)+1+(INT((l-1)/n_tor))*(n_var) 
                                   coefmtx(index_mtx) = ELM(index_ij,index_kl)
                                END IF
                             ENDDO
                             
                          ENDDO
                          CALL MURGE_ASSEMBLYSETNODEVALUES(id, index_node2, index_node1, &
                               coefmtx, ierr)
                          cnt = cnt + n_var*n_var
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*) my_id, "::", &
                                  "I", index_node1, &
                                  "J", index_node2, cnt
                             STOP
                          END IF
                       ELSE
                          DO j = 1, n_var * n_tor
                             ! Row index in the ELM matrix
                             index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j   
                          
                             DO l = 1, n_var * n_tor
                                ! Column index in the ELM matrix
                                
                                IF (INT((MOD(j-1, n_tor)+1)/2)+1 == i_tor(my_id+1) .and. &
                                     INT((MOD(l-1, n_tor)+1)/2)+1 == i_tor(my_id+1)) THEN
                                   
                                   index_kl = n_tor * n_var * (n_order+1) * (k-1) + n_tor * n_var * (k_order-1) + l
                                   index_mtx = INT((j-1)/n_tor)+MOD(MOD(j-1,n_tor),2)+1+ &
                                        (INT((l-1)/n_tor)*2+ MOD(MOD(l-1,n_tor),2))*(n_var*2)
                                   coefmtx(index_mtx) = ELM(index_ij,index_kl)
                                END IF
                             ENDDO
                             
                          ENDDO
                          CALL MURGE_ASSEMBLYSETNODEVALUES(id, index_node2, index_node1, &
                               coefmtx, ierr)
                          cnt = cnt + n_var*n_var * 4
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*) my_id, ":", &
                                  "I", index_node1, &
                                  "J", index_node2, cnt
                             STOP
                          END IF
                       END IF
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
  
  IF (method == "direct") THEN
     CALL boundary_conditions_murge(my_id,node_list,element_list,local_elms,n_local_elms, &
          xpoint2,psi_axis,psi_bnd,Z_xpoint)
  END IF
  
  IF (method .EQ. "direct") THEN
     write (*,*) MY_ID, " : Reduce..."
     CALL MPI_Allreduce(RHS_loc,RHS_glob,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_N,ierr)
  ELSE
     IF (i_tor(my_id+1) == 1) THEN
        write (*,*) my_id, " : Reduce...."
        CALL MPI_Allreduce(RHS_loc,RHS_glob,ndof_glob/n_tor,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_N, ierr)
     ELSE
        write (*,*) my_id, " : Reduce.."
        CALL MPI_AllReduce(RHS_loc,RHS_glob,ndof_glob*2/n_tor,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_N,ierr)
     END IF    
  END IF


   write(*,*) '******** end construct matrix **********'

  DEALLOCATE(RHS_loc)

  RETURN
END SUBROUTINE construct_matrix_murge
