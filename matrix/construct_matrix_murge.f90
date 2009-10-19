SUBROUTINE construct_matrix_murge(my_id,node_list,element_list, local_elms, &
     n_local_elms, index_min,index_max, xpoint2,psi_axis,psi_bnd,Z_xpoint)
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
  INTEGER                        :: index_min
  INTEGER                        :: index_max
  REAL*8                         :: psi_axis
  REAL*8                         :: psi_bnd
  REAL*8                         :: Z_xpoint
  INTEGER                        :: cnt
  INTEGER                        :: cnt2

  ! local variables
  TYPE (type_element)      :: element
  TYPE (type_node)         :: nodes(n_vertex_max)
  REAL*8, ALLOCATABLE :: rhs_loc(:)
  INTEGER :: index_min_loc, index_max_loc, coefnbr
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
  INTEGER :: i_bnd, i, in, ife, iv, inode, inode1, inode2, knode, j, k, l, index_ij, index_kl
  INTEGER :: index_i, index_large_i, index_large_k, index_node, index_node1, index_node2, i_order, k_order, ic, ielm, ierr
  INTEGER :: ijA_position,ijA_position2, nz_AA2, n_AA2, ilarge2, kv, kT, ku, ilarge_vv, ilarge_vT, ilarge_vus
  LOGICAL :: xpoint2
  LOGICAL :: is_local

  !write(*,*) '****************************************'
  !write(*,*) '*  construct matrix                    *'
  !write(*,*) '****************************************'
  !write(*,*) ' n_elements (local)       : ',my_id,n_local_elms
  !write(*,*) ' index_min,index_max      : ',my_id,index_min,index_max

  i_bnd = 0

  DO i=1, n_local_elms

     ielm = local_elms(i)

     DO iv=1,n_vertex_max

        inode = element_list%element(ielm)%vertex(iv)

        IF (node_list%node(inode)%boundary .EQ. 1) i_bnd = i_bnd + 1
        IF (node_list%node(inode)%boundary .EQ. 2) i_bnd = i_bnd + 1
        IF (node_list%node(inode)%boundary .EQ. 3) i_bnd = i_bnd + 2

        index_min_loc = MIN(index_min_loc, MINVAL(node_list%node(iv)%index))
        index_max_loc = MAX(index_max_loc, MAXVAL(node_list%node(iv)%index))

     ENDDO

  ENDDO


  IF (ALLOCATED(rhs_glob))        DEALLOCATE(rhs_glob)
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
              coefnbr = coefnbr +  (n_vertex_max)*(n_order+1)* (n_var * n_tor)* (n_var * n_tor)
           ENDIF

        ENDDO
     ENDDO

  ENDDO
  cnt = 0;
  cnt2 = 0;
  write (*,*) ":: Murge Assembly phase :: ", coefnbr, " entries"

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

                 rhs_loc(index_large_i+j) = rhs_loc(index_large_i+j) + RHS(index_ij)
              END DO
              
              ! Build nodes Matrices
              DO k=1,n_vertex_max
                 
                 knode         = element%vertex(k)
                 
                 DO k_order = 1, n_order+1

                    index_node2 = node_list%node(knode)%index(k_order)
                    
                    index_large_k = n_tor * n_var * (index_node2 - 1)

                    coefmtx = 0
                    ! BUILD node Matrix
                    DO j = 1, n_var * n_tor
                       ! Row index in the ELM matrix
                       index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j   

                       DO l = 1, n_var * n_tor
                          ! Column index in the ELM matrix
                          index_kl = n_tor * n_var * (n_order+1) * (k-1) + n_tor * n_var * (k_order-1) + l   
                          coefmtx(j+(l-1)*(n_var * n_tor)) = ELM(index_ij,index_kl)

                       ENDDO

                    ENDDO
                    
                    
                    CALL MURGE_ASSEMBLYSETNODEVALUES(id, index_node2, index_node1, &
                         coefmtx, ierr)
                    IF (ierr /= MURGE_SUCCESS) THEN
                       write (*,*) &
                            "I", index_node1, &
                            "J", index_node2
                       STOP
                    END IF
                    cnt = cnt + n_var * n_tor * n_var * n_tor 

                 ENDDO
              ENDDO
              cnt2 = cnt2 + (n_vertex_max)*(n_order+1)* n_var * n_tor * n_var * n_tor 
           ENDIF

        ENDDO
     ENDDO

  ENDDO
  CALL MURGE_ASSEMBLYEND(id, ierr)

  !write(*,*) ' nz_aa : ',my_id, nz_aa

  !----------------------- boundary conditions

  call boundary_conditions_murge(my_id,node_list,element_list,local_elms,n_local_elms, &
       xpoint2,psi_axis,psi_bnd,Z_xpoint)

  CALL MPI_Reduce(RHS_loc,RHS_glob,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  !write(*,*) '******** end construct matrix **********'

  DEALLOCATE(RHS_loc)

  RETURN
END SUBROUTINE construct_matrix_murge
