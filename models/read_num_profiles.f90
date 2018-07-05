!> Read numerical input profiles
subroutine read_num_profiles(my_id)
  
  use phys_module
  use profiles
  
  implicit none
  
  integer, intent(in) :: my_id
  
  num_rho = ( rho_file /= 'none' )
  if ( num_rho .and. ( my_id == 0 ) ) then
    call readProf(num_rho_x, num_rho_y0, num_rho_len, rho_file)
    if ( num_rho_len < 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(rho_file)//'".'
      stop
    end if
    rho_1 = num_rho_y0(num_rho_len)
    num_rho_y0 = num_rho_y0 - rho_1
  end if
  
  num_rhon = ( rhon_file /= 'none' )
  if ( num_rhon .and. ( my_id == 0 ) ) then
    call readProf(num_rhon_x, num_rhon_y0, num_rhon_len, rhon_file)
    if ( num_rhon_len < 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(rhon_file)//'".'
      stop
    end if
    rhon_1 = num_rhon_y0(num_rhon_len)
    num_rhon_y0 = num_rhon_y0 - rhon_1
  end if
  
  num_T = ( T_file /= 'none' )
  if ( num_T .and. ( my_id == 0 ) ) then
    call readProf(num_T_x, num_T_y0, num_T_len, T_file)
    if ( num_T_len < 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(T_file)//'".'
      stop
    end if
    T_1 = num_T_y0(num_T_len)
    num_T_y0 = num_T_y0 - T_1
    T_0 = num_T_y0(1)
  end if
  
  num_Te = ( Te_file /= 'none' )
  if ( num_Te .and. ( my_id == 0 ) ) then
    call readProf(num_Te_x, num_Te_y0, num_Te_len, Te_file)
    if ( num_Te_len < 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(Te_file)//'".'
      stop
    end if
    Te_1 = num_Te_y0(num_Te_len)
    num_Te_y0 = num_Te_y0 - Te_1
  end if
  
  num_Ti = ( Ti_file /= 'none' )
  if ( num_Ti .and. ( my_id == 0 ) ) then
    call readProf(num_Ti_x, num_Ti_y0, num_Ti_len, Ti_file)
    if ( num_Ti_len < 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(Ti_file)//'".'
      stop
    end if
    Ti_1 = num_Ti_y0(num_Ti_len)
    num_Ti_y0 = num_Ti_y0 - Ti_1
  end if
  
!  num_rhon = ( rhon_file /= 'none' )
!  if ( num_rhon .and. ( my_id == 0 ) ) then
!    call readProf(num_rhon_x, num_rhon_y0, num_rhon_len, rhon_file)
!    if ( num_rhon_len < 2 ) then 
!      write(*,*) '  ERROR: Could not read the numerical profile !"'//trim(rhon_file)//'".'
!      stop
!    end if
!    rhon_1 = num_rhon_y0(num_rhon_len)
!    num_rhon_y0 = num_rhon_y0 - rhon_1
!  end if
  
  
  num_ffprime = ( ffprime_file /= 'none' )
  if ( num_ffprime .and. ( my_id == 0 ) ) then
    call readProf(num_ffprime_x, num_ffprime_y0, num_ffprime_len, ffprime_file)
    if ( num_ffprime_len < 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(ffprime_file)//'".'
      stop
    end if
    FF_0 = num_ffprime_y0(1)
    FF_1 = num_ffprime_y0(num_ffprime_len)
  end if
  
  num_d_perp = ( d_perp_file /= 'none' )
  if ( num_d_perp .and. ( my_id == 0 ) ) then
    call readProf(num_d_perp_x, num_d_perp_y, num_d_perp_len, d_perp_file)
    if ( num_d_perp_len < 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(d_perp_file)//'".'
      stop
    end if
  end if

  num_zk_perp = ( zk_perp_file /= 'none' )
  if ( num_zk_perp .and. ( my_id == 0 ) ) then
    call readProf(num_zk_perp_x, num_zk_perp_y, num_zk_perp_len, zk_perp_file)
    if ( num_zk_perp_len < 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile "'//trim(zk_perp_file)//'".'
      stop
    end if
  end if

  if ( jorek_model == 400 ) then
    num_zk_e_perp = ( zk_e_perp_file /= 'none' )
    if ( num_zk_e_perp .and. ( my_id == 0 ) ) then
      call readProf(num_zk_e_perp_x, num_zk_e_perp_y, num_zk_e_perp_len, zk_e_perp_file)
      if ( num_zk_e_perp_len < 2 ) then 
        write(*,*) '  ERROR: Could not read the numerical profile "'//trim(zk_e_perp_file)//'".'
        stop
      end if
    end if
    
    num_zk_i_perp = ( zk_i_perp_file /= 'none' )
    if ( num_zk_i_perp .and. ( my_id == 0 ) ) then
      call readProf(num_zk_i_perp_x, num_zk_i_perp_y, num_zk_i_perp_len, zk_i_perp_file)
      if ( num_zk_i_perp_len < 2 ) then 
        write(*,*) '  ERROR: Could not read the numerical profile "'//trim(zk_i_perp_file)//'".'
        stop
      end if
    end if
 end if

if (NEO) then
   num_neo_file= ( neo_file /= 'none')
   if ( num_neo_file .and. ( my_id == 0 ) ) then
      write(*,*) 'using ki and mui profiles from file "'//trim(neo_file)//'"'
      call readProfNeo(num_neo_psi, num_amu_value, num_aki_value, num_neo_len, neo_file)
      
      if ( num_neo_len <= 2 ) then 
         write(*,*) '  ERROR: Could not read the numerical profile.'
         stop
      end if
   end if
endif

num_rot = ( rot_file /= 'none' )
if ( (num_rot) .and. ( my_id == 0 ) ) then
  call readProf(num_rot_x, num_rot_y0, num_rot_len, rot_file)
  if ( num_rot_len < 2 ) then 
    write(*,*) '  ERROR: Could not read the numerical profile "'//trim(rot_file)//'".'
    stop
  end if
end if


return
end subroutine read_num_profiles
