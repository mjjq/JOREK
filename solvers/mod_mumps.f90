module mumps_module
 include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure
 type (DMUMPS_STRUC) :: mumps_par
 logical             :: use_mumps,use_matrix_whitout_zeros_mumps
endmodule mumps_module
