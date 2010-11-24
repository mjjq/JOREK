subroutine derive_num_profiles()
!-----------------------------------------------------------------------
! Determine the derivatives of the numerical input profiles
!-----------------------------------------------------------------------
  
  use phys_module
  
  implicit none
  
  if ( num_rho ) then
    call derivProf(num_rho_x, num_rho_y0, num_rho_len, num_rho_y1)
    call derivProf(num_rho_x, num_rho_y1, num_rho_len, num_rho_y2)
    call derivProf(num_rho_x, num_rho_y2, num_rho_len, num_rho_y3)
  end if
  
  if ( num_T ) then
    call derivProf(num_T_x, num_T_y0, num_T_len, num_T_y1)
    call derivProf(num_T_x, num_T_y1, num_T_len, num_T_y2)
    call derivProf(num_T_x, num_T_y2, num_T_len, num_T_y3)
  end if
  
  if ( num_ffprime ) then
    call derivProf(num_ffprime_x, num_ffprime_y0, num_ffprime_len, num_ffprime_y1)
    call derivProf(num_ffprime_x, num_ffprime_y1, num_ffprime_len, num_ffprime_y2)
  end if
  
end subroutine derive_num_profiles
