!> Initialize parameters and broadcast them to all MPI procs.
subroutine initialise_and_broadcast_parameters(my_id)
  
  use parameters,  only: n_tor, n_period
  use phys_module, only: mode, mode_type
  
  implicit none
  
  ! --- Routine parameters
  integer, intent(in) :: my_id
  
  ! --- Local parameters
  integer :: itor
  
  call initialise_parameters(my_id)
  
  ! --- Broadcast input parameters from MPI thread 0 to the others.
  call broadcast_phys(my_id)
  
  ! --- Broadcast numerical input profiles from MPI thread 0 to the others.
  call broadcast_num_profiles(my_id)
  
  ! --- Fill the arrays mode (toroidal mode number n) and mode_type (cos or sin).
  do itor=1, n_tor
    mode(itor)        = int(itor / 2) * n_period
    if ( (itor==1) .or. (mod(itor,2)==0) ) then
      mode_type(itor) = 'cos'
    else
      mode_type(itor) = 'sin'
    end if
  end do
  
end subroutine initialise_and_broadcast_parameters
