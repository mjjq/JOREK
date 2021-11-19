!> mod_spectra contains variables and procedures for
!> generating radiation spectra and spectral lines
module mod_spectra
implicit none

private
public :: spectrum

!> Variable and type definitions ---------------------
!> spectrum: abstract class containing the basic types
!> and procedures for generating radiation spectra
type,abstract :: spectrum_base
  integer :: n_points !< number of wave lengths or colors
  !> wave lengths or colors
  real*8,dimension(:),allocatable :: points
  contains
  procedure,pass(spectrum_base) :: allocate_spectrum => &
  allocate_spectrum_base
  procedure,pass(spectrum_base) :: deallocate_spectrum => &
  deallocate_spectrum_base
end type spectrum_base

!> generate a set of spectral points for integration
!> using a uniform distribution
type,extends(spectrum_base) :: spectrum_rng_uniform
  real*8 :: i_pdf !< 1 over probability density function
  contains
  procedure,pass(spectrum_rng_uniform) :: deallocate_spectrum => &
  deallocate_spectrum_rng_uniform
end type spectrum_rng_uniform

!> Interfaces ---------------------------------------
interface spectrum_base
  module procedure construct_spectrum_base
end interface

contains
!> Constructors -------------------------------------
function construct_spectrum_base(n_points) result(spectrum)
  implicit none
  integer,intent(in) :: n_points
  type(spectrum_base) :: spectrum
  call spectrum%allocate_spectrum(n_points)
end function construct_spectrum_base

!> Procedures spectrum base -------------------------
!> allocate spectrum base and set counters
subroutine allocate_spectrum_base(spectrum,n_points)
  implicit none
  !> input-outputs
  type(spectrum_base),intent(inout) :: spectrum
  !> inputs
  integer,intent(in) :: n_points
  if(.not.allocated(spectrum%points).and.n_points.gt.0) &
  allocate(spectrum%points(n_points))
  spectrum%n_points = n_points
end subroutine allocate_spectrum_base

!> deallocate spectrum base and set counters to 0
subroutine deallocate_spectrum_base(spectrum)
  implicit none
  !> inputs-outputs
  type(spectrum_base),intent(inout) :: spectrum
  !> cleanup
  if(allocated(spectrum%points)) deallocate(spectrum%points)
  spectrum%n_points = 0
end subroutine deallocate_spectrum_base

!> Procedures spectrum rng uniform ------------------

!> deallocate spectrum rng uniform and cleanup
subroutine deallocate_spectrum_rng_uniform(spectrum)
  implicit none
  !> inputs-outputs
  type(spectrum_rng_uniform),intent(inout) :: spectrum
  call deallocate_spectrum_base(spectrum)
  spectrum%i_pdf = -1.d0
end subroutine deallocate_spectrum_rng_uniform

!>---------------------------------------------------
end module mod_spectra
