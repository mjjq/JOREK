!
! Routine: vertex_is_local
!
! check if vertex is in the SORTED in ascendent order list 
! loc2glob of size n and return .true. or .false. in islocal.
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
  imax = murge_local_n

  IF (ALLOCATED(murge_glob2loc)) THEN

  ELSE
     ALLOCATE(murge_glob2loc(murge_global_n))
     murge_glob2loc = -1
     DO iter = 1, murge_local_n
        murge_glob2loc(murge_loc2glob(iter)) = iter
     END DO
  END IF
  IF (vertex > murge_global_n) THEN
     write (*,*) "vertex_is_local: Vertex out of range"
     STOP
  END IF
  
  IF (murge_glob2loc(vertex) > 0) THEN 
     islocal = .true.
  ELSE
     islocal = .false.
  END IF
END SUBROUTINE vertex_is_local
