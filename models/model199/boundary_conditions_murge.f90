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
!* Authors:                                                                    *
!*   Xavier Lacoste - xavier.lacoste@inria.fr                                  *
!*                                                                             *
!*******************************************************************************
SUBROUTINE boundary_conditions_murge(my_id, node_list, element_list, & 
     local_elms, n_local_elms, psi_bnd)			    
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
  REAL*8  :: zbig
  INTEGER :: i, in, iv, inode, k
  INTEGER :: index_node, ielm, ierr
  INTEGER :: coefnbr
  LOGICAL :: is_local
  zbig = 1.d12

  ! Count the number of values that will be entered
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


  ! Add the boundary entries to murge Matrix.
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
