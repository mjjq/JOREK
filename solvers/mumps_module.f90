module mumps_module
  
  save
  
#ifdef USE_MUMPS
  include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure
#else
  include 'no_dmumps_struct.h'
#endif
  
  type (DMUMPS_STRUC) :: mumps_par
  logical             :: use_mumps,no_zeros_mumps, use_mumps_BLR
  REAL*8              :: mumps_BLR_eps
  integer             :: MPI_COMM_MUMPS_EQUIL, MPI_GROUP_MUMPS_EQUIL
  
end module mumps_module
