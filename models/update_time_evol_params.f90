!> Update the time evolution parameters time_evol_theta and time_evol_zeta.
subroutine update_time_evol_params()
  
  use phys_module, only: time_evol_scheme, time_evol_theta, time_evol_zeta
  
  implicit none
  
  if ( time_evol_scheme == 'implicit Euler' ) then
    time_evol_theta = 1.0d0
    time_evol_zeta  = 0.0d0
  else if ( time_evol_scheme == 'Crank-Nicholson' ) then
    time_evol_theta = 0.5d0
    time_evol_zeta  = 0.0d0
  else if ( time_evol_scheme == 'Gears' ) then
    time_evol_theta = 1.0d0
    time_evol_zeta  = 0.5d0
  else
    write(*,*) 'ERROR: Illegal value for time_evol_scheme: ', trim(time_evol_scheme)
    stop
  end if
  
end subroutine update_time_evol_params
