!> Initialize parameters and broadcast them to all MPI procs.
subroutine initialise_and_broadcast_parameters(my_id, filename)
  
  use constants, only: mu_zero
  use mod_parameters,  only: n_tor, n_period
  use phys_module
  
  implicit none
  
  ! --- Routine parameters
  integer,                      intent(in) :: my_id
  character(len=*),             intent(in) :: filename

  call initialise_parameters(my_id, filename)
  
  ! --- Broadcast input parameters from MPI thread 0 to the others.
  call broadcast_phys(my_id)
  
  ! --- Broadcast numerical input profiles from MPI thread 0 to the others.
  call broadcast_num_profiles(my_id)

  ! --- Initialize the time-stepping parameters.
  call update_time_evol_params()
  
  ! --- Calculate normalization factors.
  sqrt_mu0_rho0      = sqrt( mu_zero * ( central_density * 1.d20 * central_mass * mass_proton ) )
  sqrt_mu0_over_rho0 = sqrt( mu_zero / ( central_density * 1.d20 * central_mass * mass_proton ) )
  
  ! --- Deprecated input parameters ---
  if ( use_murge ) then
    write(*,*) 'ERROR: use_murge=.true. is not supported any more. Remove this parameter from the namelist input file.'
    stop
  else if ( use_murge_element ) then
    write(*,*) 'ERROR: use_murge_element=.true. is not supported any more. Remove this parameter from the namelist input file.'
    stop
  end if
  ! -----------------------------------
  
end subroutine initialise_and_broadcast_parameters
