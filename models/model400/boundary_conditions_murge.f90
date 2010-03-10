!*******************************************************************************
!* Subroutine: boundary_condition_murge                                        *
!*******************************************************************************
!*                                                                             *
!* Add boundary condition on the matrix given to murge (PaStiX) solver.        *
!*                                                                             *
!* Important: If this file is modified, boundary_conditions.f90 should also    *
!*            be modified.                                                     *
!*            boundary_conditions_murge.f90 in other models folders may also   *
!*            need modifications.                                              *
!*                                                                             *
!* Parameters:                                                                 *
!*   my_id        - Identifier of the node in MPI_COMM_WORLD                   *
!*   node_list    - List of nodes                                              *
!*   element_list - List of all elements                                       *
!*   local_elms   - List of local elements                                     *
!*   n_local_elms - Number of local elements                                   *
!*   psi_bnd      - UNUSED                                                     *
!*                                                                             *
!*******************************************************************************
SUBROUTINE boundary_conditions_murge(my_id,node_list,element_list,local_elms, &
     n_local_elms, psi_bnd)			    
  !---------------------------------------------------------------
  ! add the boundary condition to the global matrix
  !---------------------------------------------------------------
  USE data_structure
  USE global_distributed_matrix
  USE phys_module
  USE murge_module

  IMPLICIT NONE
  INCLUDE 'mpif.h'

  ! Subroutine parameters
  INTEGER                  :: my_id
  INTEGER                  :: local_elms(*)
  INTEGER                  :: n_local_elms
  TYPE (type_node_list)    :: node_list
  TYPE (type_element_list) :: element_list
  REAL*8                   :: psi_bnd

  ! Internal parameters

  REAL*8  :: zbig, Ti0, Te0, Vpar0, bigR, dTi0_ds, dTe0_ds, dVpar0_ds, dBigR_ds
  REAL*8  :: R_s, R_t, Z_s, Z_t, ps0_s, ps0_t, ps0_x, ps0_y, direction, xjac
  REAL*8  :: Btot
  REAL*8  :: grad_psi, u0_s, u0_t, u0_x, u0_y
  INTEGER :: i, in, iv, inode, k
  INTEGER :: index_node, index_node2, ielm, ierr
  INTEGER :: kv, kTi, kTe, ku
  INTEGER :: coefnbr
  LOGICAL :: is_local

  coefnbr = 0

  zbig = 1.d10

  DO i=1, n_local_elms

     ielm = local_elms(i)

     DO iv=1, n_vertex_max

        inode = element_list%element(ielm)%vertex(iv)

        IF (node_list%node(inode)%boundary .NE. 0) THEN

           DO in=1, n_tor

              DO k=1, n_var

                 !------------------------------------ the open field lines (in case of x-point grid)
                 IF ((node_list%node(inode)%boundary .EQ. 1) .OR. (node_list%node(inode)%boundary .EQ. 3)) THEN

                    IF ((k .EQ.   1) .OR. (k .EQ. 2) .OR. (k .EQ. 3)  .OR. &
                         (k .EQ.  4)  .OR. (k .EQ. 95)  .OR. (k .EQ. 96) .OR. (k .EQ.97) .OR. (k .EQ.98) ) THEN

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

                    IF (k .EQ. 7) THEN

                       index_node  = node_list%node(inode)%index(1)             ! position of value
                       index_node2 = node_list%node(inode)%index(2)             ! position of first deriative

                       CALL vertex_is_local(index_node, is_local)
                       IF (is_local) THEN

                          coefnbr = coefnbr + 4
                       ENDIF

                       CALL vertex_is_local(index_node2, is_local)
                       IF (is_local) THEN

                          coefnbr = coefnbr + 3
                       ENDIF

                    ENDIF

                 ENDIF

                 !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                 IF ((node_list%node(inode)%boundary .EQ. 2) .OR. (node_list%node(inode)%boundary .EQ. 3)) THEN

                    IF ((k .EQ. 1) .OR. (k .EQ. 2) .OR. (k .EQ. 3) .OR. &
                         (k .EQ. 4) .OR. (k .EQ. 95) .OR. (k .EQ. 96) .OR. (k .EQ. 7) .OR. (k .EQ. 98) ) THEN

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


  write (*,*) my_id, ":: Murge Boundary Assembly phase :: ", coefnbr, " entries"

  CALL MURGE_ASSEMBLYBEGIN(id, coefnbr, MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW, &
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

                    IF ((k .EQ.   1) .OR. (k .EQ. 2) .OR. (k .EQ. 3)  .OR. &
                         (k .EQ.  4)  .OR. (k .EQ. 95)  .OR. (k .EQ. 96) .OR. (k .EQ.97) .OR. (k .EQ.98) ) THEN

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

                    IF (k .EQ. 7) THEN

                       index_node  = node_list%node(inode)%index(1)             ! position of value
                       index_node2 = node_list%node(inode)%index(2)             ! position of first deriative

                       Ti0       = node_list%node(inode)%values(1,1,6)
                       Te0       = node_list%node(inode)%values(1,1,8)
                       Vpar0     = node_list%node(inode)%values(1,1,k)
                       BigR      = node_list%node(inode)%x(1,1)
                       dTi0_ds    = node_list%node(inode)%values(1,2,6)
                       dTe0_ds    = node_list%node(inode)%values(1,2,8)
                       dVpar0_ds = node_list%node(inode)%values(1,2,k)
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


                       CALL vertex_is_local(index_node, is_local)
                       IF (is_local) THEN
                        
                          ku  = 2
                          kv  = 7
                          kTi = 6
                          kTe = 8

                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node - 1) + (kv-1)*n_tor + in, &
                               n_tor * n_var * (index_node - 1) + (kv-1)*n_tor + in, &
                               zbig, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF

                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node - 1) + (kv -1)*n_tor + in, &
                               n_tor * n_var * (index_node - 1) + (kTi-1)*n_tor + in, &
                               - zbig / Btot * 0.5d0 * GAMMA / SQRT(GAMMA*(Ti0 + Te0)) * direction, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF

                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node - 1) + (kv -1)*n_tor + in, &
                               n_tor * n_var * (index_node - 1) + (kTe-1)*n_tor + in, &
                               - zbig / Btot * 0.5d0 * GAMMA / SQRT(GAMMA*(Ti0 + Te0)) * direction, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF

                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node  - 1) + (kv -1)*n_tor + in, &
                               n_tor * n_var * (index_node2 - 1) + (ku -1)*n_tor + in, &
                               - zbig * BigR**2 / ps0_s, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF


                          RHS_glob(n_tor*n_var * (index_node-1) + (kv-1)*n_tor + in) = &
                               Zbig * ( - Vpar0 + BigR**2 * U0_s /ps0_s + direction*SQRT(GAMMA*(Ti0 + Te0))/ Btot)

                       END IF


                       CALL vertex_is_local(index_node2, is_local)
                       IF (is_local) THEN

                          kv  = 7
                          kTi = 6
                          kTe = 8

                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node2 - 1) + (kv-1)*n_tor + in, &
                               n_tor * n_var * (index_node2 - 1) + (kv-1)*n_tor + in, &
                               zbig, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF

                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node2 - 1) + (kv-1)*n_tor + in, &
                               n_tor * n_var * (index_node2 - 1) + (kTi-1)*n_tor + in, &
                               - zbig / Btot * 0.5d0 * GAMMA / SQRT(GAMMA*(Ti0 + Te0)) * direction, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF

                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node2 - 1) + (kv-1)*n_tor + in, &
                               n_tor * n_var * (index_node2 - 1) + (kTe-1)*n_tor + in, &
                               - zbig / Btot * 0.5d0 * GAMMA / SQRT(GAMMA*(Ti0 + Te0)) * direction, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF

                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node2 - 1) + (kv-1)*n_tor + in, &
                               n_tor * n_var * (index_node  - 1) + (kTi-1)*n_tor + in, &
                               + zbig / Btot * 0.25d0 * GAMMA**(2.d0) / ((GAMMA*(Ti0 + Te0))**1.5d0) * dTi0_ds * direction, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF

                          CALL MURGE_ASSEMBLYSETVALUE(id, &
                               n_tor * n_var * (index_node2 - 1) + (kv-1)*n_tor + in, &
                               n_tor * n_var * (index_node  - 1) + (kTe-1)*n_tor + in, &
                               + zbig / Btot * 0.25d0 * GAMMA**(2.d0) / ((GAMMA*(Ti0 + Te0))**1.5d0) * dTi0_ds * direction, ierr)
                          IF (ierr /= MURGE_SUCCESS) THEN
                             WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                                  "I", index_node, &
                                  "J", index_node
                             STOP
                          END IF

                          RHS_glob(n_tor * n_var * (index_node-1) + (kv-1)*n_tor + in) = &
                               Zbig*(-dVpar0_ds +  0.5d0 / Btot * GAMMA / SQRT(GAMMA*(Ti0 + Te0)) * (dTi0_ds + dTe0_ds) * direction)

                       ENDIF

                    ENDIF

                 ENDIF

                 !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                 IF ((node_list%node(inode)%boundary .EQ. 2) .OR. (node_list%node(inode)%boundary .EQ. 3)) THEN

                    IF ((k .EQ. 1) .OR. (k .EQ. 2) .OR. (k .EQ. 3) .OR. &
                         (k .EQ. 4) .OR. (k .EQ. 95) .OR. (k .EQ. 96) .OR. (k .EQ. 7) .OR. (k .EQ. 98) ) THEN

                       index_node = node_list%node(inode)%index(1)

                       CALL MURGE_ASSEMBLYSETVALUE(id, &
                            n_tor * n_var * (index_node - 1) + (k -1)*n_tor + in, &
                            n_tor * n_var * (index_node - 1) + (k -1)*n_tor + in, &
                            zbig, ierr)
                       IF (ierr /= MURGE_SUCCESS) THEN
                          WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                               "I", index_node, &
                               "J", index_node
                          STOP
                       END IF

                       index_node = node_list%node(inode)%index(3)

                       CALL MURGE_ASSEMBLYSETVALUE(id, &
                            n_tor * n_var * (index_node - 1) + (k -1)*n_tor + in, &
                            n_tor * n_var * (index_node - 1) + (k -1)*n_tor + in, &
                            zbig, ierr)
                       IF (ierr /= MURGE_SUCCESS) THEN
                          WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                               "I", index_node, &
                               "J", index_node
                          STOP
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
