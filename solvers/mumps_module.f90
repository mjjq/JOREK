module mumps_module
#ifdef USE_MUMPS  
  save
  

  include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure
  
  type (DMUMPS_STRUC) :: mumps_par
  !integer             :: mumps_ordering        ! ordering option (7:automatic, 3:Scotch, 4:PORD, 5:METIS), default: 7
  !logical             :: no_zeros_mumps
  !integer             :: MPI_COMM_MUMPS_EQUIL, MPI_GROUP_MUMPS_EQUIL
  
#endif
end module mumps_module
