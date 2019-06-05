subroutine initialise_mumps(MPI_COMM_in)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use mumps_module
use mpi_mod

implicit none

integer :: MPI_COMM_in
#ifdef USE_MUMPS
mumps_par%COMM = MPI_COMM_in                   ! Define a communicator for mumps

mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1
mumps_par%ICNTL(13) = -1

call  DMUMPS(mumps_par)                        ! Initialize an instance of mumps

mumps_par%ICNTL(2)  = -1                       ! Verbosity levels
mumps_par%ICNTL(3)  = -1
mumps_par%ICNTL(4)  = 6

mumps_par%ICNTL(14) = 20                       ! memory working space increase
!mumps_par%ICNTL(15) = 0                       ! memory balance only, 1: flops
#endif
return
end
