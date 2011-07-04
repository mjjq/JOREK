SUBROUTINE split_brodcast(type,MPI_COMM_N)
  !
  ! Split MPI_BCAST if MPI buffer beyond 2Go
  ! 
  USE mumps_module
  !  
  IMPLICIT NONE
  !
  INCLUDE 'mpif.h'
  !
  !IN
  CHARACTER*8           :: type
  INTEGER               :: MPI_COMM_N
  !LOCAL
  INTEGER               :: buff_max
  INTEGER               :: ierr
  INTEGER               :: i,crit,nz_split,nz_split_end,ie,is
  !
  IF ( trim(type) .EQ. 'intIRN' .OR. trim(type) .EQ. 'intJCN')  THEN
     buff_max = 500000000
  ELSE IF ( trim(type) .EQ. 'double' )  THEN
     buff_max = 250000000
  ELSE
     WRITE(*,*) 'ERROR : type in split_broadcast.f90'
  ENDIF
  !
  IF ( mumps_par%nz > buff_max ) THEN 
     !
     crit =  mumps_par%nz / buff_max      
     DO i=1,crit
        is = (i-1)*buff_max+1
        ie = i*buff_max
        IF (trim(type).EQ.'intIRN') THEN
           CALL MPI_BCAST(mumps_par%IRN(is:ie),buff_max,MPI_INTEGER,0,MPI_COMM_N,ierr)
        ELSE IF (trim(type).EQ.'intJCN') THEN
           CALL MPI_BCAST(mumps_par%JCN(is:ie),buff_max,MPI_INTEGER,0,MPI_COMM_N,ierr)   
        ELSE IF (trim(type).EQ.'double') THEN
           CALL MPI_BCAST(mumps_par%A(is:ie),buff_max,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
        ENDIF
     ENDDO
     !
     nz_split_end = mod(mumps_par%nz,buff_max)
     is = crit * buff_max + 1
     ie = mumps_par%nz
     IF (trim(type).EQ.'intIRN') THEN
        CALL MPI_BCAST(mumps_par%IRN(is:ie),nz_split_end,MPI_INTEGER,0,MPI_COMM_N,ierr)
     ELSE IF (trim(type).EQ.'intJCN') THEN
        CALL MPI_BCAST(mumps_par%JCN(is:ie),nz_split_end,MPI_INTEGER,0,MPI_COMM_N,ierr)   
     ELSE IF (trim(type).EQ.'double') THEN
        CALL MPI_BCAST(mumps_par%A(is:ie),nz_split_end,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
     ENDIF
  ELSE
     IF (trim(type).EQ.'intIRN') THEN
        CALL MPI_BCAST(mumps_par%IRN,mumps_par%nz,MPI_INTEGER,0,MPI_COMM_N,ierr)
     ELSE IF (trim(type).EQ.'intJCN') THEN   
        CALL MPI_BCAST(mumps_par%JCN,mumps_par%nz,MPI_INTEGER,0,MPI_COMM_N,ierr)
     ELSE IF (trim(type).EQ.'double') THEN
        CALL MPI_BCAST(mumps_par%A,mumps_par%nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
     ENDIF
     !
  ENDIF
  !
END SUBROUTINE split_brodcast
