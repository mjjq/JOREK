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
  integer :: n_points  !< number of wave lengths or colors
  integer :: n_spectra !< number of wave length or color intervals
  !> wave lengths or colors
  real*8,dimension(:,:),allocatable :: points
  contains
  procedure,pass(spectrum_base) :: allocate_spectrum   => allocate_spectrum_base
  procedure,pass(spectrum_base) :: deallocate_spectrum => deallocate_spectrum_base
end type spectrum_base

!> generate a set of spectral points for integration
!> using a uniform distribution
type,extends(spectrum_base) :: spectrum_rng_uniform
  real*8,dimension(:),allocatable :: min_wlen !< lower wavelenght of the interval
  real*8,dimension(:),allocatable :: i_pdf !< 1 over probability density function
  contains
  procedure,pass(spectrum_rng_uniform) :: allocate_spectrum      => allocate_spectrum_rng_uniform
  procedure,pass(spectrum_rng_uniform) :: set_spectrum_variables => set_uniform_spectrum
  procedure,pass(spectrum_rng_uniform) :: generate_spectrum      => generate_uniform_rng_spectrum
  procedure,pass(spectrum_rng_uniform) :: integrate_data         => integrate_rng_uniform
  procedure,pass(spectrum_rng_uniform) :: deallocate_spectrum    => deallocate_spectrum_rng_uniform
end type spectrum_rng_uniform

!> Interfaces ---------------------------------------
interface spectrum_base
  module procedure construct_spectrum_base
end interface

interface spectrum_rng_uniform
  module procedure construct_spectrum_rng_uniform
end interface

contains
!> Constructors -------------------------------------
function construct_spectrum_base(n_points,n_spectra) &
result(spectrum)
  implicit none
  integer,intent(in) :: n_points
  class(spectrum_base) :: spectrum
  call spectrum%allocate_spectrum(n_points,n_spectra)
end function construct_spectrum_base

!> construct a uniform random spectrum
!> inputs:
!>   n_points:  (integer) number of discrete random variables
!>   n_spectra: (integer) number of spectra
!>   min_wlen:  (real8) minimum wavelength (interval lowerbound)
!>   max_wlen:  (real8) maximum wavelength (interval upperbound)
!> outputs:
!>   spectrum: (spectrum_rng_uniform) uniform random spectrum
function construct_spectrum_rng_uniform(n_points,n_spectra,&
min_wlen,max_wlen) result(spectrum)
  implicit none
  integer,intent(in)                     :: n_points,n_spectra
  real*8,dimension(n_spectra),intent(in) :: min_wlen,max_wlen
  class(spectrum_rng_uniform)            :: spectrum 
  call spectrum%allocate_spectrum(n_points,n_spectra,&
  min_wlen,max_wlen)
end function construct_spectrum_rng_uniform

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
  spectrum&n_spectra = -1
end subroutine deallocate_spectrum_base

!> Procedures spectrum rng uniform ------------------

!> allocate the spectrum_rng_uniform datatype
!> inputs:
!>   spectrum:  (spectrum_rng_uniform) generates and integrates
!>              variables along a uniform spectral distribution
!>   n_points:  (integer) number of spectral points
!>   n_spectra: (integer) number of spectral intervals
!>   min_wlen:  (real8)(n_spectra) minimum wavelength
!>   max_wlen:  (real8)(n_spectra) maximum wavelength
!> outputs:
!>   spectrum: (spectrum_rng_uniform) generates and integrates
!>             variables along a uniform spectral distribution
subroutine allocate_spectrum_rng_uniform(spectrum,n_points,&
n_spectra,min_wlen,max_wlen)
  implicit none
  !> inputs
  integer,intent(in) :: n_points,n_spectra
  real*8,dimension(n_spectra) :: min_wlen,max_wlen
  !> inputs-outputs
  class(spectrum_rng_uniform) :: spectrum
  !> allocated all variables
  call allocate_spectrum_base(spectrum,n_poins,n_spectra)
  call spectrum%set_uniform_spectrum(n_spectra,min_wlen,max_wlen)
end subroutine allocate_spectrum_rng_uniform

!> set the spectrum properties
!> inputs:
!>   spectrum:  (spectrum_rng_uniform) generates and integrates
!>              variables along a uniform spectral distribution
!>   n_spectra: (integer) number of spectral intervals
!>   min_wlen:  (real8)(n_spectra) minimum wavelength
!>   max_wlen:  (real8)(n_spectra) maximum wavelength
!> outputs:
!>   spectrum: (spectrum_rng_uniform) generates and integrates
!>             variables along a uniform spectral distribution
subroutine set_uniform_spectrum(spectrum,n_spectra,min_wlen,max_wlen)
  implicit none
  !> inputs-outputs
  class(spectrum_rng_uniform),intent(inout) :: spectrum
  !> inputs
  integer,intent(in) :: n_spectra
  real*8,dimension(n_spectra),intent(in) :: min_wlen,max_wlen
  !> set values
  if(.not.allocated(spectrum%min_wlen)) allocate(spectrum%min_wlen(n_spectra))
  if(.not.allocated(spectrum%i_pdf))    allocate(spectrum%i_pdf)
  if(spectrum%n_spectra.ne.n_spectra) then
    deallocate(spectrum%min_wlen); deallocate(specrun%i_pdf);
    allocate(spectrum%min_wlen(n_spectrum))
    allocate(spectrum%i_pdf(n_spectrum))
    spectrum%n_spectra = n_spectra
    write(*,*) 'WARNING: n_spectra is changed: regenerate spectrum points!'
  endif
  spectrum%min_wlen = min_wlen
  spectrum%i_pdf = max_wlen - min_wlen
end subroutine set_uniform_spectrum

!> generate uniform random spectrum
!> inputs:
!>   spectrum: (spectrum_rng_uniform) generates and integrates 
!>             variables along a uniform spectral distribution
!>   n_rng:
!>   rng:      (type_rng)(n_rng) random number generators
!> outputs:
!>   spectrum: (spectrum_rng_uniform) generates and integrates 
!>             variables along a uniform spectral distribution
subroutine generate_uniform_rng_spectrum(spcetrum,n_rng,rng)
  use mod_rng
  implicit none
  !> inputs-outpus
  class(spectrum_rng_uniform),intent(inout) :: spectrum
  !> inputs
  type(type_rng)
  
end subroutine generate_uniform_rng_spectrum

!> integrate data computed using uniformly distributed spectra
!> inputs:
!>   spectrum:     (spectrum_rng_uniform) generates and integrates 
!>                 variables along a uniform spectral distribution
!>   uniform_data: (real8)(n_points,n_spectra) data obtained from
!>                 a uniform spectral distribution
!> outputs:
!>   spectrum:  (spectrum_rng_uniform) generates and integrates 
!>              variables along a uniform spectral distribution
!>   integrals: (n_spectra) integrals of uniform_data for each spectrum
subroutine integrate_rng_uniform(spectrum,uniform_data,integrals)
  implicit none
  !> inputs-outputs
  class(spectrum_uniform_rng),intent(inout) :: spectrum
  !> inputs
  real*8,dimension(spectrum%n_points,spectrum%n_spectra),intent(in) :: uniform_data
  !> outputs
  real*8,dimension(spectrum%n_spectra),intent(out) :: integrals
  !> variables
  integer :: ii
  !> integrate
  do ii=1,spectrum%n_spectra
    integrals(ii) = sum(uniform_data(:,ii))
  enddo
  integrals = integrals*spectrum%i_pdf/spectrum%n_points
end subroutine integrate_rng_uniform

!> deallocate spectrum rng uniform and cleanup
subroutine deallocate_spectrum_rng_uniform(spectrum)
  implicit none
  !> inputs-outputs
  class(spectrum_rng_uniform),intent(inout) :: spectrum
  call deallocate_spectrum_base(spectrum)
  if(allocated(spectrum%min_wlen)) deallocate(spectrum%min_wlen)
  if(allocated(spectrum%i_pdf))    deallocate(spectrum%i_pdf)
end subroutine deallocate_spectrum_rng_uniform

!>---------------------------------------------------
end module mod_spectra
