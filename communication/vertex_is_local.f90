!
! Routine: vertex_is_local
!
! check if vertex is in the SORTED in ascendent order list 
! loc2glob of size n and return it in islocal.
!
SUBROUTINE vertex_is_local(vertex, loc2glob, n, islocal)

  INTEGER, INTENT(IN)  :: vertex
  INTEGER, INTENT(IN)  :: loc2glob(n)
  INTEGER, INTENT(IN)  :: n
  LOGICAL, INTENT(OUT) :: islocal
  
  INTEGER              :: imin
  INTEGER              :: imax
  INTEGER              :: imid

  imin = 1
  imax = n
  
  DO WHILE(imin < imax)
     imid = imin + (imax - imin +1)/2
     IF (loc2glob(imin) == vertex) THEN
        islocal = .true.
        RETURN
     END IF
     IF (loc2glob(imid) == vertex) THEN
        islocal = .true.
        RETURN
     END IF
     IF (loc2glob(imid) < vertex) THEN
        if (imin == imid) then
           islocal = .false.
           return
        end if
        imin = imid
     ELSE
        if (imax == imid) then
           islocal = .false.
           return
        end if
        imax = imid
     END IF
  END DO
  
  islocal = .false.

END SUBROUTINE vertex_is_local
