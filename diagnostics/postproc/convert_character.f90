!> Module for converting character strings to integer or float values (used by jorek2_postproc)
module convert_character

  use mod_parameters, only: n_var, variable_names

  implicit none
  
  
  
  private
  public lower_case, to_int, to_log, to_float, real2str
  
  
  
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
  
  
  
  !> Convert a character string to a logical
  logical function to_log(s, error)
    character(len=*), intent(in)  :: s     !< Character string to be converted
    integer,          intent(out) :: error !< Error flag
    
    read(s,*,iostat=error) to_log
    if ( error /= 0 ) write(*,*) 'ERROR: Parameter "', trim(s), '" is not a logical value.'
  end function to_log
  
  
  
  !> Convert a character string to a floating point number
  real*8 function to_float(s, error)
    character(len=*), intent(in)  :: s     !< Character string to be converted
    integer,          intent(out) :: error !< Error flag
    
    read(s,*,iostat=error) to_float
    if ( error /= 0 ) write(*,*) 'ERROR: Parameter "', trim(s), '" is not a floating point number.'
  end function to_float
  
  
  
  !> Convert a real number into a character string (for filename generation).
  character(len=15) function real2str(r,f) result(s)
    real*8, intent(in) :: r
    character(len=*), optional, intent(in) :: f ! format
    
    character(len=30) :: form
    
    form = '(f7.3)'
    if ( present(f) ) form=f
    
    write(s,trim(form)) r
    s = adjustl(s)
    
  end function real2str
  
  
  
end module convert_character
