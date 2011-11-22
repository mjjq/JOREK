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
     call tr_allocate(num_rho_x,1,num_rho_len,"num_rho_x",CAT_UNKNOWN)
     call tr_allocate(num_rho_y0,1,num_rho_len,"num_rho_y0",CAT_UNKNOWN)
     call tr_allocate(num_rho_y1,1,num_rho_len,"num_rho_y1",CAT_UNKNOWN)
     call tr_allocate(num_rho_y2,1,num_rho_len,"num_rho_y2",CAT_UNKNOWN)
     call tr_allocate(num_rho_y3,1,num_rho_len,"num_rho_y3",CAT_UNKNOWN)
  end if
  call MPI_BCAST(num_rho_x,num_rho_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rho_y0,num_rho_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rho_y1,num_rho_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rho_y2,num_rho_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rho_y3,num_rho_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_T ) then
  call MPI_BCAST(num_T_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_T_x,1,num_T_len,"num_T_x",CAT_UNKNOWN)
     call tr_allocate(num_T_y0,1,num_T_len,"num_T_y0",CAT_UNKNOWN)
     call tr_allocate(num_T_y1,1,num_T_len,"num_T_y1",CAT_UNKNOWN)
     call tr_allocate(num_T_y2,1,num_T_len,"num_T_y2",CAT_UNKNOWN)
     call tr_allocate(num_T_y3,1,num_T_len,"num_T_y3",CAT_UNKNOWN)
  end if
  call MPI_BCAST(num_T_x,num_T_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_T_y0,num_T_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_T_y1,num_T_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_T_y2,num_T_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_T_y3,num_T_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_ffprime ) then
  call MPI_BCAST(num_ffprime_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_ffprime_x,1,num_ffprime_len,"num_ffprime_x",CAT_UNKNOWN)
     call tr_allocate(num_ffprime_y0,1,num_ffprime_len,"num_ffprime_y0",CAT_UNKNOWN)
     call tr_allocate(num_ffprime_y1,1,num_ffprime_len,"num_ffprime_y1",CAT_UNKNOWN)
     call tr_allocate(num_ffprime_y2,1,num_ffprime_len,"num_ffprime_y2",CAT_UNKNOWN)
  end if
  call MPI_BCAST(num_ffprime_x,num_ffprime_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_ffprime_y0,num_ffprime_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_ffprime_y1,num_ffprime_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_ffprime_y2,num_ffprime_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_d_perp ) then
  call MPI_BCAST(num_d_perp_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_d_perp_x,1,num_d_perp_len,"num_d_perp_x")
     call tr_allocate(num_d_perp_y,1,num_d_perp_len,"num_d_perp_y")
  end if
  call MPI_BCAST(num_d_perp_x,num_d_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_d_perp_y,num_d_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_zk_perp ) then
  call MPI_BCAST(num_zk_perp_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_zk_perp_x,1,num_zk_perp_len,"num_zk_perp_x")
     call tr_allocate(num_zk_perp_y,1,num_zk_perp_len,"num_zk_perp_y")
  end if
  call MPI_BCAST(num_zk_perp_x,num_zk_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_zk_perp_y,num_zk_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( jorek_model == 400 ) then
  if ( num_zk_e_perp ) then
    call MPI_BCAST(num_zk_e_perp_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( my_id /= 0 ) then
       call tr_allocate(num_zk_e_perp_x,1,num_zk_e_perp_len,"num_zk_e_perp_x")
       call tr_allocate(num_zk_e_perp_y,1,num_zk_e_perp_len,"num_zk_e_perp_y")
    end if
    call MPI_BCAST(num_zk_e_perp_x,num_zk_e_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(num_zk_e_perp_y,num_zk_e_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end if
  
  if ( num_zk_i_perp ) then
    call MPI_BCAST(num_zk_i_perp_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( my_id /= 0 ) then
       call tr_allocate(num_zk_i_perp_x,1,num_zk_i_perp_len,"num_zk_i_perp_x")
       call tr_allocate(num_zk_i_perp_y,1,num_zk_i_perp_len,"num_zk_i_perp_y")
    end if
    call MPI_BCAST(num_zk_i_perp_x,num_zk_i_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(num_zk_i_perp_y,num_zk_i_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end if
end if

return
end subroutine broadcast_num_profiles
