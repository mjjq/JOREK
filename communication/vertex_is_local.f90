!
! Routine: vertex_is_local
!
! check if vertex is in the SORTED in ascendent order list 
! loc2glob of size n and return it in islocal.
!
SUBROUTINE vertex_is_local(vertex, islocal)

  USE murge_module

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: vertex
  LOGICAL, INTENT(OUT) :: islocal

  INTEGER              :: iter
  INTEGER              :: imin
  INTEGER              :: imax
  INTEGER              :: imid

  imin = 1
  imax = local_n

  IF (ALLOCATED(glob2loc)) THEN
     IF (vertex > global_n) THEN
        write (*,*) "vertex_is_local: Vertex out of range"
        STOP
     END IF
        
     IF (glob2loc(vertex) > 0) THEN 
        islocal = .true.
     ELSE
        islocal = .false.
     END IF
  ELSE
     ALLOCATE(glob2loc(global_n))
     glob2loc = -1
     DO iter = 1, local_n
        glob2loc(loc2glob(iter)) = iter
     END DO
     IF (glob2loc(vertex) > 0) THEN 
        islocal = .true.
     ELSE
        islocal = .false.
     END IF
  END IF

END SUBROUTINE vertex_is_local
