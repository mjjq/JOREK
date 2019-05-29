module mumps_module
  
  save
  
#ifdef USE_MUMPS
  include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure
#else
  include 'no_dmumps_struct.h'
#endif
  
  type (DMUMPS_STRUC) :: mumps_par
  logical             :: use_mumps
  integer             :: mumps_ordering        ! ordering option (7:automatic, 3:Scotch, 4:PORD, 5:METIS), default: 7
  logical             :: no_zeros_mumps
  logical             :: use_mumps_BLR         ! switch  for BLR compression, default: .false.
  REAL*8              :: mumps_BLR_eps         ! accuracy of BLR compression, default: 0
  integer             :: MPI_COMM_MUMPS_EQUIL, MPI_GROUP_MUMPS_EQUIL
  
end module mumps_module
