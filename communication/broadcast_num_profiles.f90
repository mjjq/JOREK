subroutine broadcast_num_profiles(my_id)
!----------------------------------------------------------
! Broadcast the numerical input profiles from MPI thread 0 to the others
!----------------------------------------------------------
use tr_module
use phys_module

implicit none

include 'mpif.h'               ! MPI fortran include file

integer, intent(in) :: my_id

integer :: ierr

if ( num_rho ) then
  call MPI_BCAST(num_rho_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_rho_x,1,num_rho_len,"num_rho_x")
     call tr_allocate(num_rho_y0,1,num_rho_len,"num_rho_y0")
  end if
  call MPI_BCAST(num_rho_x,num_rho_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rho_y0,num_rho_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_T ) then
  call MPI_BCAST(num_T_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_T_x,1,num_T_len,"num_T_x")
     call tr_allocate(num_T_y0,1,num_T_len,"num_T_y0")
  end if
  call MPI_BCAST(num_T_x,num_T_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_T_y0,num_T_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_ffprime ) then
  call MPI_BCAST(num_ffprime_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_ffprime_x,1,num_ffprime_len,"num_ffprime_x")
     call tr_allocate(num_ffprime_y0,1,num_ffprime_len,"num_ffprime_y0")
  end if
  call MPI_BCAST(num_ffprime_x,num_ffprime_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_ffprime_y0,num_ffprime_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

return
end subroutine broadcast_num_profiles
