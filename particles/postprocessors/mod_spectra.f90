!> mod_spectra contains base variables and procedures for
!> the spectra modules
module mod_spectra
implicit none

private
public :: spectrum_base

!> Variable and type definitions ---------------------
!> spectrum: abstract class containing the basic types
!> and procedures for generating radiation spectra
type,abstract :: spectrum_base
  integer :: n_points  !< number of wave lengths or colors
  integer :: n_spectra !< number of wave length or color intervals
  !> wave lengths or colors
  real*8,dimension(:,:),allocatable :: points
  contains
  procedure,pass(spectrum) :: allocate_spectrum_base
  procedure,pass(spectrum) :: deallocate_spectrum_base
end type spectrum_base

contains
!> Procedures spectrum base -------------------------
!> allocate spectrum base and set counters
subroutine allocate_spectrum_base(spectrum,n_points,n_spectra)
  implicit none
  !> input-outputs
  class(spectrum_base),intent(inout) :: spectrum
  !> inputs
  integer,intent(in) :: n_points,n_spectra
  if(.not.allocated(spectrum%points).and.((n_points.gt.0).and.&
  (n_spectra.gt.0))) &
  allocate(spectrum%points(n_points,n_spectra))
  spectrum%n_points  = n_points
  spectrum%n_spectra = n_spectra
end subroutine allocate_spectrum_base

!> deallocate spectrum base and set counters to 0
subroutine deallocate_spectrum_base(spectrum)
  implicit none
  !> inputs-outputs
  class(spectrum_base),intent(inout) :: spectrum
  !> cleanup
  if(allocated(spectrum%points)) deallocate(spectrum%points)
  spectrum%n_points  = -1
  spectrum%n_spectra = -1
end subroutine deallocate_spectrum_base

!>---------------------------------------------------

end module mod_spectra
