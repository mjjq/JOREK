module mod_catalyst_adaptor
#ifdef USE_CATALYST

  use iso_c_binding

  implicit none

  public :: catalyst_adaptor_initialise
  public :: catalyst_adaptor_finalise
  public :: catalyst_adaptor

  private

  character(kind=c_char, len=4096), public :: catalyst_script

  interface

    subroutine catalyst_adaptor_initialise(a_catalyst_script) bind(C)
      use iso_c_binding
      character(kind=c_char), intent(in) :: a_catalyst_script 
    end subroutine catalyst_adaptor_initialise

    subroutine catalyst_adaptor_finalise() bind(C)
    end subroutine catalyst_adaptor_finalise

    ! empty function to get the dependency generator to work with catalyst_adaptor.cpp
    subroutine catalyst_adaptor() bind(C)
    end subroutine catalyst_adaptor

  end interface

#endif 
end module mod_catalyst_adaptor