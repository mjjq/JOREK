!> Module for converting character strings to integer or float values
!! (used by jorek2_postproc)
module convert_character

  use parameters, only: n_var, variable_names

  implicit none
  
  
  
  private
  public lower_case, to_int, to_float, get_variable_number
  
  
  
  contains
  
  
  
  !> Convert a character string to lower case
  function lower_case(s)
    
    ! --- Routine parameters
    character(len=*), intent(in) :: s          !< Characters string to be converted
    character(len=len(s))        :: lower_case !< Return value (lower case character string)
    
    ! --- Local variables
    integer :: i, ic, nlen
    
    lower_case = s
    
    nlen = len(lower_case)
    
    do i = 1, nlen
      ic = ichar(lower_case(i:i))
      if (ic >= 65 .and. ic <= 90) lower_case(i:i) = char(ic+32)
    end do
    
  end function lower_case
  
  
  
  !> Convert a character string to an integer
  integer function to_int(s, error)
    character(len=*), intent(in)  :: s     !< Character string to be converted
    integer,          intent(out) :: error !< Error flag
    
    read(s,*,iostat=error) to_int
    if ( error /= 0 ) write(*,*) 'ERROR: Parameter "', trim(s), '" is not an integer number.'
  end function to_int
  
  
  
  !> Convert a character string to a floating point number
  real*8 function to_float(s, error)
    character(len=*), intent(in)  :: s     !< Character string to be converted
    integer,          intent(out) :: error !< Error flag
    
    read(s,*,iostat=error) to_float
    if ( error /= 0 ) write(*,*) 'ERROR: Parameter "', trim(s), '" is not a floating point number.'
  end function to_float
  
  
  
  !> Get the variable number for a given character string that can contain the variable name or
  !! its index.
  integer function get_variable_number(s, error)
    Character(len=*), intent(in)     :: s
    integer,          intent(out)    :: error
    
    integer  :: i
    
    error = 0
    
    do i = 1, n_var
      if ( lower_case(trim(s)) == lower_case(trim(variable_names(i))) ) then
        get_variable_number = i
        return
      end if
    end do
    get_variable_number = to_int(s, error)
    
  end function get_variable_number
  
  
  
end module convert_character
