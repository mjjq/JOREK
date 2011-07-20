subroutine initialise_and_broadcast_parameters(my_id)
  
  implicit none
  
  ! --- Routine parameters
  integer, intent(in) :: my_id
  
  call initialise_parameters(my_id)
  
  ! --- Broadcast input parameters from MPI thread 0 to the others.
  call broadcast_phys(my_id)
  
  ! --- Broadcast numerical input profiles from MPI thread 0 to the others.
  call broadcast_num_profiles(my_id)
  
end subroutine initialise_and_broadcast_parameters
