!> Read numerical input profiles
subroutine read_num_profiles()
  
  use phys_module
  use profiles
  
  implicit none
  
  num_rho = ( rho_file /= 'none' )
  if ( num_rho ) then
    call readProf(num_rho_x, num_rho_y0, num_rho_len, rho_file)
    if ( num_rho_len <= 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(rho_file)//'".'
      stop
    end if
    rho_1 = num_rho_y0(num_rho_len)
    num_rho_y0 = num_rho_y0 - rho_1
  end if
  
  num_T = ( T_file /= 'none' )
  if ( num_T ) then
    call readProf(num_T_x, num_T_y0, num_T_len, T_file)
    if ( num_T_len <= 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(T_file)//'".'
      stop
    end if
    T_1 = num_T_y0(num_T_len)
    num_T_y0 = num_T_y0 - T_1
  end if
  
  num_ffprime = ( ffprime_file /= 'none' )
  if ( num_ffprime ) then
    call readProf(num_ffprime_x, num_ffprime_y0, num_ffprime_len, ffprime_file)
    if ( num_ffprime_len <= 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(ffprime_file)//'".'
      stop
    end if
    FF_1 = 0.
  end if
  
  num_d_perp = ( d_perp_file /= 'none' )
  if ( num_d_perp ) then
    call readProf(num_d_perp_x, num_d_perp_y, num_d_perp_len, d_perp_file)
    if ( num_d_perp_len <= 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(d_perp_file)//'".'
      stop
    end if
  end if

  num_zk_perp = ( zk_perp_file /= 'none' )
  if ( num_zk_perp ) then
    call readProf(num_zk_perp_x, num_zk_perp_y, num_zk_perp_len, zk_perp_file)
    if ( num_zk_perp_len <= 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(zk_perp_file)//'".'
      stop
    end if
  end if

  if ( jorek_model == 400 ) then
    num_zk_e_perp = ( zk_e_perp_file /= 'none' )
    if ( num_zk_e_perp ) then
      call readProf(num_zk_e_perp_x, num_zk_e_perp_y, num_zk_e_perp_len, zk_e_perp_file)
      if ( num_zk_e_perp_len <= 2 ) then 
        write(*,*) '  ERROR: Could not read the numerical profile "'//trim(zk_e_perp_file)//'".'
        stop
      end if
    end if
    
    num_zk_i_perp = ( zk_i_perp_file /= 'none' )
    if ( num_zk_i_perp ) then
      call readProf(num_zk_i_perp_x, num_zk_i_perp_y, num_zk_i_perp_len, zk_i_perp_file)
      if ( num_zk_i_perp_len <= 2 ) then 
        write(*,*) '  ERROR: Could not read the numerical profile "'//trim(zk_i_perp_file)//'".'
        stop
      end if
    end if
  end if

return
end subroutine read_num_profiles
