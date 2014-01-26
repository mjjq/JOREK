!> This module allows to apply a Fourier analysis or a Fourier filter to the result array
!! of the routine eval_expr() from the module mod_expression.
module mod_four_filter
  
  
  
  
  
  use parameters
  
  
  
  
  
  implicit none
  
  
  
  
  
  public
  
  
  
  
  
  ! --- Constants
  character(len=15), parameter, private :: THIS_MOD_NAME = 'mod_four_filter'

  
  
  
  
  
  type t_four_filter
    !### store somehow which harmonics to keep
    !### maybe even smooth "window"?
  end type t_four_filter
  
  
  
  
  
  contains
  
  
  
  
  
  subroutine perform_four_transform(input, output, forward, poloidal, toroidal)
    !### r2r bidirectional
    !### pol, tor, both
  end subroutine perform_four_transform
  
  
  
  
  
  subroutine cleanup_four_filter(filter)
    !###
  end subroutine cleanup_four_filter
  
  
  
  
  subroutine create_four_filter(filter, some-filter-specification)
    call cleanup_four_filter(filter)
    !###
  end subroutine create_four_filter
  
  
  
  
  subroutine apply_four_filter(inoutput)
    !### apply filter, i.e., set parts of the coeffs to zero
  end subroutine apply_four_filter
  
  
  
  
  
  !### all-in-one-routine: trafo+create filter+apply filter+backtrafo
  
  
  
  
  
end module mod_four_filter
