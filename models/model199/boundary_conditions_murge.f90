SUBROUTINE boundary_conditions_murge(my_id, node_list, element_list, local_elms, &
     n_local_elms, xpoint2,psi_axis,psi_bnd,Z_xpoint)			    
  !---------------------------------------------------------------
  ! add the boundary condition to the global matrix
  !---------------------------------------------------------------
  USE data_structure
  USE global_distributed_matrix
  USE phys_module
  USE murge_module

  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE (type_node_list)    :: node_list
  TYPE (type_element_list) :: element_list
  TYPE (type_element)      :: element
  TYPE (type_node_list)    :: nodes

  INTEGER :: my_id, local_elms(*), n_local_elms
  REAL*8  :: zbig, psi_axis, psi_bnd, Z_xpoint, T0, Vpar0, bigR, dT0_ds, dVpar0_ds, dBigR_ds
  REAL*8  :: R_s, R_t, Z_s, Z_t, ps0_s, ps0_t, ps0_x, ps0_y, direction, xjac
  REAL*8  :: Vpar0_pol_R, Vpar0_pol_Z, Vpol_R, Vpol_Z, znormal_R, znormal_Z
  REAL*8  :: Vpar0_perp, Vpol_perp, Btot, cs_fraction, ratio
  REAL*8  :: grad_s, grad_psi, u0_s, u0_t, u0_x, u0_y
  INTEGER :: i_bnd, i, in, ife, iv, inode, inode1, inode2, knode, j, k, l, index_ij, index_kl
  INTEGER :: index_i, index_large_i, index_large_k, index_node, index_node1, index_node2, i_order, k_order, ic, ielm, ierr
  INTEGER :: ijA_position,ijA_position2, nz_AA2, n_AA2, ilarge2, kv, kT, ku, ilarge_vv, ilarge_vT, ilarge_vus
  LOGICAL :: xpoint2
  INTEGER :: coefnbr
  LOGICAL :: is_local
  zbig = 1.d12

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

                    IF ((k .EQ. 1) .OR. (k .EQ. 2) .OR. (k .EQ. 3) .OR. &
                         (k .EQ. 4) .OR. (k .EQ. 5) .OR. (k .EQ. 6) ) THEN

                       index_node = node_list%node(inode)%index(1)

                       CALL vertex_is_local(index_node, is_local)
                       IF (is_local) THEN

                          coefnbr = coefnbr + 1
                       ENDIF

                       index_node = node_list%node(inode)%index(2)

                       CALL vertex_is_local(index_node, is_local)
                       IF (is_local) THEN

                          coefnbr = coefnbr + 1

                       ENDIF

                    ENDIF

                 ENDIF

                 !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                 IF ((node_list%node(inode)%boundary .EQ. 2) .OR. (node_list%node(inode)%boundary .EQ. 3)) THEN

                    IF ((k .EQ. 1) .OR. (k .EQ. 2) .OR. (k .EQ. 3) .OR. &
                         (k .EQ. 4) .OR. (k .EQ. 5) .OR. (k .EQ. 6) ) THEN

                       index_node = node_list%node(inode)%index(1)
                       CALL vertex_is_local(index_node, is_local)
                       IF (is_local) THEN
                          coefnbr = coefnbr + 1
                       ENDIF

                       index_node = node_list%node(inode)%index(3)
                       CALL vertex_is_local(index_node, is_local)
                       IF (is_local) THEN
                          coefnbr = coefnbr + 1
                       ENDIF

                    ENDIF

                 ENDIF

              ENDDO

           ENDDO
        ENDIF
     ENDDO
  ENDDO
  write (*,*) ":: Murge Boundary Assembly phase :: ", coefnbr, " entries"

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

                    IF ((k .EQ. 1) .OR. (k .EQ. 2) .OR. (k .EQ. 3) .OR. &
                         (k .EQ. 4) .OR. (k .EQ. 5) .OR. (k .EQ. 6) ) THEN

                       index_node = node_list%node(inode)%index(1)

                       CALL vertex_is_local(index_node, is_local)
                       IF (is_local) THEN
                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node - 1) + (k-1)*n_tor + in, &
                               n_tor * n_var * (index_node - 1) + (k-1)*n_tor + in, &
                               zbig, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF
                       END IF

                       index_node = node_list%node(inode)%index(2)

                       CALL vertex_is_local(index_node, is_local)
                       IF (is_local) THEN
                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node - 1) + (k-1)*n_tor + in, &
                               n_tor * n_var * (index_node - 1) + (k-1)*n_tor + in, &
                               zbig, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF
                       END IF

                    ENDIF

                 ENDIF

                 !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                 IF ((node_list%node(inode)%boundary .EQ. 2) .OR. (node_list%node(inode)%boundary .EQ. 3)) THEN

                    IF ((k .EQ. 1) .OR. (k .EQ. 2) .OR. (k .EQ. 3) .OR. &
                         (k .EQ. 4) .OR. (k .EQ. 5) .OR. (k .EQ. 6) ) THEN

                       index_node = node_list%node(inode)%index(1)
                       CALL vertex_is_local(index_node, is_local)
                       IF (is_local) THEN
                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node - 1) + (k-1)*n_tor + in, &
                               n_tor * n_var * (index_node - 1) + (k-1)*n_tor + in, &
                               zbig, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF
                       END IF

                       index_node = node_list%node(inode)%index(3)
                       CALL vertex_is_local(index_node, is_local)
                       IF (is_local) THEN
                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node - 1) + (k-1)*n_tor + in, &
                               n_tor * n_var * (index_node - 1) + (k-1)*n_tor + in, &
                               zbig, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF
                       END IF
                    ENDIF

                 ENDIF

              ENDDO

           ENDDO
        ENDIF
     ENDDO
  ENDDO

  CALL MURGE_ASSEMBLYEND(id, ierr)

  RETURN
END SUBROUTINE boundary_conditions_murge
