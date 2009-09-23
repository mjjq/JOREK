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
  TYPE (type_node_list)    :: nodes
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

  IF (.NOT. ALLOCATED(A_glob))    ALLOCATE(A_glob(nz_glob))
  IF (.NOT. ALLOCATED(irn_glob))  ALLOCATE(irn_glob(nz_glob))
  IF (.NOT. ALLOCATED(jcn_glob))  ALLOCATE(jcn_glob(nz_glob))

  IF (ALLOCATED(rhs_glob))        DEALLOCATE(rhs_glob)

  ALLOCATE(rhs_glob(ndof_glob))

  IF (.NOT. ALLOCATED(rhs_loc))   ALLOCATE(rhs_loc(ndof_glob))

  irn_glob = 0
  jcn_glob = 0
  A_glob   = 0.d0
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

           CALL vertex_is_local(index_node1, loc2glob, local_n, is_local)
           IF (is_local) THEN 
              coefnbr = coefnbr +  (n_vertex_max)*(n_order+1)* (n_var * n_tor)* (n_var * n_tor)
           ENDIF

        ENDDO
     ENDDO

  ENDDO
  cnt = 0;
  cnt2 = 0;
  write (*,*) ":: Premiere phase d'assemblage ::", n_local_elms, coefnbr

  CALL MURGE_ASSEMBLYBEGIN(id, coefnbr, MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_ADD, &
       MURGE_ASSEMBLY_FOOL, murge_sym, ierr)

  DO ife =1, n_local_elms

     ielm = local_elms(ife)

     element = element_list%element(ielm)

     DO iv = 1, n_vertex_max

        inode          = element%vertex(iv)
        nodes%node(iv) = node_list%node(inode)

     ENDDO

     IF (n_tor .GT. 3) THEN
        CALL element_matrix_fft(element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)      ! use fft for toroidal integration
     ELSE
        CALL element_matrix(element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)           ! use direct integration
     ENDIF

     DO i=1,n_vertex_max

        inode1         = element%vertex(i)

        DO i_order = 1, n_order+1

           index_node1 = node_list%node(inode1)%index(i_order)

           index_large_i = n_tor * n_var * (index_node1 - 1)

           CALL vertex_is_local(index_node1, loc2glob, local_n, is_local)
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
                       
                       index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j   ! index in the ELM matrix

                       DO l = 1, n_var * n_tor

                          index_kl = n_tor * n_var * (n_order+1) * (k-1) + n_tor * n_var * (k_order-1) + l   ! index in the ELM matrix
                          coefmtx(j+l*(n_var * n_tor)) = ELM(index_ij,index_kl)

                       ENDDO

                    ENDDO
                    
                    CALL MURGE_ASSEMBLYSETNODEVALUES(id, index_node1, index_node2, &
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

  zbig = 1.d10
  ! Compute coefnbr for boundary values
  coefnbr = 0
  DO i=1, n_local_elms

     ielm = local_elms(i)
     
     DO iv=1, n_vertex_max
        
        inode = element_list%element(ielm)%vertex(iv)

        IF (node_list%node(inode)%boundary .NE. 0) THEN
           
           DO in=1, n_tor

              DO k=1, n_var

                 !------------------------------------ the open field lines (in case of x-point grid)
                 IF ((node_list%node(inode)%boundary .EQ. 1) .OR. (node_list%node(inode)%boundary .EQ. 3)) THEN

                    !            if ((k .eq.   1) .or. (k .eq. 92) .or. (k .eq. 3)  .or. &
                    !                (k .eq.  4)  .or. (k .eq. 95)  .or. (k .eq. 96) .or. (k .eq.97) ) then

                    IF ((k .EQ.   1) .OR. (k .EQ. 2) .OR. (k .EQ. 3) .OR. &
                         (k .EQ.  4)  .OR. (k .EQ. 5) .OR. (k .EQ. 6) .OR. (k .EQ.97) ) THEN

                       index_node = node_list%node(inode)%index(1)

                       CALL vertex_is_local(index_node, loc2glob, local_n, is_local)
                       IF (is_local) THEN
                          
                          coefnbr = coefnbr + 1

                       ENDIF

                       index_node = node_list%node(inode)%index(2)
                       
                       CALL vertex_is_local(index_node, loc2glob, local_n, is_local)
                       IF (is_local) THEN

                          coefnbr = coefnbr + 1

                       ENDIF

                    ENDIF ! Test on k
                    
                    ! Other value of k
                    IF (k .EQ. 7) THEN

                       index_node  = node_list%node(inode)%index(1)             ! position of value

                       xjac  =  R_s*Z_t - R_t*Z_s
                       ps0_x = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
                       ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac

                       u0_x = (   Z_t * u0_s - Z_s * u0_t ) / xjac
                       u0_y = ( - R_t * u0_s + R_s * u0_t ) / xjac

                       direction = + ps0_x / ABS(ps0_x)             ! temporary solution for lower x-point only

                       grad_psi = SQRT(ps0_x**2 + ps0_y**2)

                       Btot = SQRT(F0**2 + ps0_x**2 + ps0_y**2) / BigR

                       CALL vertex_is_local(index_node, loc2glob, local_n, is_local)
                       IF (is_local) THEN
                         
                          coefnbr = coefnbr + 1

                       ENDIF

                       index_node = node_list%node(inode)%index(2)

                       CALL vertex_is_local(index_node, loc2glob, local_n, is_local)
                       IF (is_local) THEN

                          coefnbr = coefnbr + 1
  
                       ENDIF

                    ENDIF ! Value of k

                 ENDIF !boundary 1 or 3

                 !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                 IF ((node_list%node(inode)%boundary .EQ. 2) .OR. (node_list%node(inode)%boundary .EQ. 3)) THEN

                    !            if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
                    !                (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 96) .or. (k .eq. 7) ) then
                    IF ((k .EQ. 1) .OR. (k .EQ. 2) .OR. (k .EQ. 3) .OR. &
                         (k .EQ. 4) .OR. (k .EQ. 5) .OR. (k .EQ. 6) .OR. (k .EQ. 7) ) THEN

                       index_node = node_list%node(inode)%index(1)

                       CALL vertex_is_local(index_node, loc2glob, local_n, is_local)
                       IF (is_local) THEN

     
                          coefnbr = coefnbr + 1

                       ENDIF

                       index_node = node_list%node(inode)%index(3)

                       CALL vertex_is_local(index_node, loc2glob, local_n, is_local)
                       IF (is_local) THEN
                         
                          coefnbr = coefnbr + 1

                       ENDIF

                    ENDIF ! value of k

                 ENDIF ! boundary 2 or 3

              ENDDO

           ENDDO
        ENDIF
     ENDDO
  ENDDO
  coefnbr = coefnbr *  (n_var * n_tor) *  (n_var * n_tor)

  write (*,*) ":: Seconde phase d'assemblage ::", n_local_elms, coefnbr
  CALL MURGE_ASSEMBLYBEGIN(id, coefnbr, MURGE_ASSEMBLY_ADD, MURGE_ASSEMBLY_ADD, &
       MURGE_ASSEMBLY_FOOL, murge_sym, ierr)
  DO i=1, n_local_elms

     ielm = local_elms(i)

     DO iv=1, n_vertex_max

        inode = element_list%element(ielm)%vertex(iv)

        IF (node_list%node(inode)%boundary .NE. 0) THEN
     
           
           DO in=1, n_tor

              DO k=1, n_var

                 !------------------------------------ the open field lines (in case of x-point grid)
                 IF ((node_list%node(inode)%boundary .EQ. 1) .OR. (node_list%node(inode)%boundary .EQ. 3)) THEN

                    !            if ((k .eq.   1) .or. (k .eq. 92) .or. (k .eq. 3)  .or. &
                    !                (k .eq.  4)  .or. (k .eq. 95)  .or. (k .eq. 96) .or. (k .eq.97) ) then

                    IF ((k .EQ.   1) .OR. (k .EQ. 2) .OR. (k .EQ. 3) .OR. &
                         (k .EQ.  4)  .OR. (k .EQ. 5) .OR. (k .EQ. 6) .OR. (k .EQ.97) ) THEN

                       index_node = node_list%node(inode)%index(1)

                       CALL vertex_is_local(index_node, loc2glob, local_n, is_local)
                       IF (is_local) THEN
                                
                          coefmtx = 0
                          coefmtx(((k-1)*n_tor + in)*(n_tor * n_var) + ((k-1)*n_tor + in)) = zbig;
                          CALL MURGE_ASSEMBLYSETNODEVALUES(id, index_node, index_node, &
                               coefmtx, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             write (*,*) 376, &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF


                       ENDIF

                       index_node = node_list%node(inode)%index(2)
                       
                       CALL vertex_is_local(index_node, loc2glob, local_n, is_local)
                       IF (is_local) THEN

                          coefmtx = 0
                          coefmtx(((k-1)*n_tor + in)*(n_tor * n_var) + ((k-1)*n_tor + in)) = zbig;
                          CALL MURGE_ASSEMBLYSETNODEVALUES(id, index_node, index_node, &
                               coefmtx, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             write (*,*) 396, &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF
                       ENDIF

                    ENDIF ! Test on k
                    
                    ! Other value of k
                    IF (k .EQ. 7) THEN

                       index_node  = node_list%node(inode)%index(1)             ! position of value
                       index_node2 = node_list%node(inode)%index(2)             ! position of first deriative

                       T0        = node_list%node(inode)%values(1,1,6)
                       Vpar0     = node_list%node(inode)%values(1,1,7)
                       BigR      = node_list%node(inode)%x(1,1)
                       dT0_ds    = node_list%node(inode)%values(1,2,6)
                       dVpar0_ds = node_list%node(inode)%values(1,2,7)
                       dBigR_ds  = node_list%node(inode)%x(2,1)

                       ps0_s     = node_list%node(inode)%values(1,2,1)
                       ps0_t     = node_list%node(inode)%values(1,3,1)

                       U0_s      = node_list%node(inode)%values(1,2,2)
                       U0_t      = node_list%node(inode)%values(1,3,2)

                       R_s       = node_list%node(inode)%x(2,1)
                       R_t       = node_list%node(inode)%x(3,1)
                       Z_s       = node_list%node(inode)%x(2,2)
                       Z_t       = node_list%node(inode)%x(3,2)

                       xjac  =  R_s*Z_t - R_t*Z_s
                       ps0_x = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
                       ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac

                       u0_x = (   Z_t * u0_s - Z_s * u0_t ) / xjac
                       u0_y = ( - R_t * u0_s + R_s * u0_t ) / xjac

                       direction = + ps0_x / ABS(ps0_x)             ! temporary solution for lower x-point only

                       grad_psi = SQRT(ps0_x**2 + ps0_y**2)

                       Btot = SQRT(F0**2 + ps0_x**2 + ps0_y**2) / BigR

                       IF (in .EQ. 1) THEN

                          !                write(*,'(A,3e14.6,A,e14.6)') ' Boundary : ',Vpar0, -BigR**2 * u0_s/ps0_s, direction*sqrt(GAMMA*T0)/Btot,&
                          !                                              ' error : ',Vpar0 - BigR**2 * u0_s/ps0_s - direction*sqrt(GAMMA*T0)/Btot

                       ENDIF

                       CALL vertex_is_local(index_node, loc2glob, local_n, is_local)
                       IF (is_local) THEN

                          index_large_i = n_tor * n_var * (index_node - 1)

                          ku = 2
                          kv = 7
                          kT = 6

                          ilarge_vv  = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kv-1)*n_tor + in
                          ilarge_vT  = ijA_position  - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in
                          ilarge_vus = ijA_position2 - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (ku-1)*n_tor + in
                          coefmtx = 0
                          coefmtx(((kv-1)*n_tor + in)*(n_tor * n_var) + ((kv-1)*n_tor + in)) = zbig

                          coefmtx(((kv-1)*n_tor + in)*(n_tor * n_var) + ((kT-1)*n_tor + in)) = &
                               - zbig / Btot * 0.5d0 * GAMMA / SQRT(GAMMA*T0) * direction

                          coefmtx(((kv-1)*n_tor + in)*(n_tor * n_var) + ((ku-1)*n_tor + in)) = &
                               - zbig * BigR**2 / ps0_s
                          CALL MURGE_ASSEMBLYSETNODEVALUES(id, index_node, index_node, &
                               coefmtx, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             write (*,*) 471, &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF
                          RHS_glob(n_tor*n_var * (index_node-1) + (kv-1)*n_tor + in) = &
                               Zbig * ( - Vpar0 + BigR**2 * U0_s /ps0_s + direction*SQRT(GAMMA*T0))

                       ENDIF

                       index_node = node_list%node(inode)%index(2)

                       CALL vertex_is_local(index_node, loc2glob, local_n, is_local)
                       IF (is_local) THEN


                          index_large_i = n_tor * n_var * (index_node - 1)

                          kv = 7
                          kT = 6

                          ilarge_vv = ijA_position - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kv-1)*n_tor + in
                          ilarge_vT = ijA_position - 1 + ((kv-1)*n_tor + in-1) * n_var*n_tor + (kT-1)*n_tor + in

                          !coefmtx = 0
                          !coefmtx(((kv-1)*n_tor + in)*(n_tor * n_var) + ((kv-1)*n_tor + in)) = zbig

                          !coefmtx(((kT-1)*n_tor + in)*(n_tor * n_var) + ((kv-1)*n_tor + in)) = &
                          ! - zbig / Btot * 0.5d0 * GAMMA / sqrt(GAMMA*T0) * direction
                          !CALL MURGE_ASSEMBLYSETNODEVALUES(id, index_node, index_node, &
                          !     coefmtx, ierr)
                          !               RHS_glob(n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in) = &
                          !                       Zbig*(-dVpar0_ds +  0.5d0 / Btot * GAMMA / sqrt(GAMMA*T0) * dT0_ds * direction)

                       ENDIF

                    ENDIF ! Value of k

                 ENDIF !boundary 1 or 3

                 !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                 IF ((node_list%node(inode)%boundary .EQ. 2) .OR. (node_list%node(inode)%boundary .EQ. 3)) THEN

                    !            if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
                    !                (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 96) .or. (k .eq. 7) ) then
                    IF ((k .EQ. 1) .OR. (k .EQ. 2) .OR. (k .EQ. 3) .OR. &
                         (k .EQ. 4) .OR. (k .EQ. 5) .OR. (k .EQ. 6) .OR. (k .EQ. 7) ) THEN

                       index_node = node_list%node(inode)%index(1)
                       CALL vertex_is_local(index_node, loc2glob, local_n, is_local)
                       IF (is_local) THEN
                          coefmtx = 0
                          coefmtx(((k-1)*n_tor + in)*(n_tor * n_var) + ((k-1)*n_tor + in)) = zbig
                          CALL MURGE_ASSEMBLYSETNODEVALUES(id, index_node, index_node, &
                               coefmtx, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             write (*,*) 527, &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF
                       ENDIF

                       index_node = node_list%node(inode)%index(3)
                       CALL vertex_is_local(index_node, loc2glob, local_n, is_local)
                       IF (is_local) THEN
                          coefmtx = 0
                          coefmtx(((k-1)*n_tor + in)*(n_tor * n_var) + ((k-1)*n_tor + in)) = zbig
                          CALL MURGE_ASSEMBLYSETNODEVALUES(id, index_node, index_node, &
                               coefmtx, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             write (*,*) 542, &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF
                       ENDIF

                    ENDIF ! value of k

                 ENDIF ! boundary 2 or 3

              ENDDO

           ENDDO
        ENDIF
     ENDDO
  ENDDO
  CALL MURGE_ASSEMBLYEND(id, ierr)
  CALL MPI_Reduce(RHS_loc,RHS_glob,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  !write(*,*) '******** end construct matrix **********'

  DEALLOCATE(RHS_loc)

  RETURN
END SUBROUTINE construct_matrix_murge
