!> Module for settings which can be modified with the 'set' command (used by jorek2_postproc)
module settings
  
  use convert_character
  
  implicit none
  
  !> Derived datatype for one setting
  type type_setting
    character(len=128)  :: name   !< Name of the setting, e.g., 'linepoints'
    character(len=1024) :: value  !< Value of the setting, e.g., '200'
    character(len=1024) :: descr  !< Description
  end type type_setting
  
  integer, parameter :: n_max_settings = 100    !< Maximum number of settings
  integer            :: n_settings     = 0      !< Current number of settings
  type(type_setting) :: setting(n_max_settings) !< List of all settings
  
  
  
  private
  public set_setting, get_setting, get_int_setting, get_float_setting, get_log_setting,            &
    print_settings
  
  
  
  save
  
  
  
  contains
   
  
  
  !> Sets the value of the setting name to a given value
  subroutine set_setting(name, value, error, descr)
    ! --- Routine parameters
    integer,                     intent(out)    :: error   !< Error flag
    character(len=*),            intent(in)     :: name    !< Name of setting to be changed
    character(len=*),            intent(in)     :: value   !< (New) value for the setting
    character(len=*), optional,  intent(in)     :: descr   !< (New) description for the setting
      
    ! --- Local variables
    integer     :: i
    
    error = 0
    
    do i = 1, n_settings
      if (setting(i)%name == name) then
        ! Entry already exists, change its value
        setting(i)%value = value
        if ( present(descr) ) setting(i)%descr = descr
        return
      end if
    end do

    if (n_settings == n_max_settings) then
      write(*,*) 'ERROR: maximum of possible settings (', n_max_settings, ') reached.'
      error = 1
      return
    end if
    
    ! --- Entry does not exist yet, create it.
    n_settings = n_settings + 1
    setting(n_settings)%name  = name
    setting(n_settings)%value = value
    setting(n_settings)%descr = ''
    if ( present(descr) ) setting(n_settings)%descr = descr

  end subroutine set_setting
  

  
  !> Returns the character value of a requested setting
  character(len=1024) function get_setting(name, error)
    
    ! --- Routine parameters
    character(len=*),   intent(in)     :: name    !< Name of requested setting
    integer,            intent(out)    :: error   !< Error flag
    
    integer     :: i
    
    error = 0
    
    do i = 1, n_settings
      if (setting(i)%name == name) then
        get_setting = setting(i)%value
        return
      end if
    end do
    
    error = 1 ! (if setting was not found)
  end function get_setting
  
  
  
  !> Returns the integer value of a setting
  integer function get_int_setting(name, error)
    
    ! --- Routine parameters
    character(len=*),   intent(in)     :: name    !< Name of requested setting
    integer,            intent(out)    :: error   !< Error flag
    
    character(len=1024) :: value
    
    value = get_setting(name, error)
    if ( error /= 0 ) return
    
    get_int_setting = to_int(value, error)
  end function get_int_setting
  
  
  
  !> Returns the float value of a setting
  real*8 function get_float_setting(name, error)
    
    ! --- Routine parameters
    character(len=*),   intent(in)     :: name     !< Name of requested setting
    integer,            intent(out)    :: error   !< Error flag
    
    character(len=1024) :: value
    
    value = get_setting(name, error)
    if ( error /= 0 ) return
    
    get_float_setting = to_float(value, error)
  end function get_float_setting
  
  
  
  !> Returns the logical value of a setting
  logical function get_log_setting(name, error)
    
    ! --- Routine parameters
    character(len=*),   intent(in)     :: name     !< Name of requested setting
    integer,            intent(out)    :: error   !< Error flag
    
    character(len=1024) :: value
    
    value = get_setting(name, error)
    if ( error /= 0 ) return
    
    get_log_setting = to_log(value, error)
  end function get_log_setting
  
  
  
  !> List all existing settings
  subroutine print_settings()

    integer  :: i
    
    do i = 1, n_settings
      write(*,'(1x,5a)') trim(setting(i)%name), '=', trim(setting(i)%value), '          !', trim(setting(i)%descr)
    end do
  
  end subroutine print_settings

end module settings
