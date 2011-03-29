module mumps_module
  
  save
  
#ifdef USE_MUMPS
  include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure
#else
  include 'no_dmumps_struct.h'
#endif
  
  type (DMUMPS_STRUC) :: mumps_par
  logical             :: use_mumps,use_matrix_whitout_zeros_mumps
  integer             :: MPI_COMM_MUMPS_EQUIL, MPI_GROUP_MUMPS_EQUIL
  
end module mumps_module
