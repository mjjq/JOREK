!> mod_spectra contains variables and procedures for
!> generating radiation spectra and spectral lines
module mod_spectra
implicit none

private
public :: spectrum_base,spectrum_rng_uniform

!> Variable and type definitions ---------------------
!> spectrum: abstract class containing the basic types
!> and procedures for generating radiation spectra
type,abstract :: spectrum_base
  integer :: n_points  !< number of wave lengths or colors
  integer :: n_spectra !< number of wave length or color intervals
  !> wave lengths or colors
  real*8,dimension(:,:),allocatable :: points
end type spectrum_base

!> generate a set of spectral points for integration
!> using a uniform distribution
type,extends(spectrum_base) :: spectrum_rng_uniform
  real*8,dimension(:),allocatable :: min_wlen !< lower wavelenght of the interval
  real*8,dimension(:),allocatable :: i_pdf !< 1 over probability density function
  contains
  procedure,pass(spectrum) :: allocate_spectrum      => allocate_spectrum_rng_uniform
  procedure,pass(spectrum) :: set_spectrum_interval  => set_uniform_spectrum_interval
  procedure,pass(spectrum) :: generate_spectrum      => generate_uniform_rng_spectrum
  procedure,pass(spectrum) :: integrate_data         => integrate_rng_uniform
  procedure,pass(spectrum) :: deallocate_spectrum    => deallocate_spectrum_rng_uniform
end type spectrum_rng_uniform

!> Interfaces ---------------------------------------

interface spectrum_rng_uniform
  module procedure construct_spectrum_rng_uniform
end interface

contains
!> Constructors -------------------------------------
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
  real*8,dimension(n_spectra),optional,intent(in) :: min_wlen,max_wlen
  type(spectrum_rng_uniform),target      :: spectrum
  real*8,dimension(2*n_spectra)          :: real8_param 
  if(present(min_wlen).and.present(max_wlen)) then
    real8_param(1:n_spectra) = min_wlen
    real8_param(n_spectra+1:2*n_spectra) = max_wlen
    call spectrum%allocate_spectrum(n_points,n_spectra,real8_param)
  else
    call spectrum%allocate_spectrum(n_points,n_spectra)
  endif
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
  spectrum%n_spectra = -1
end subroutine deallocate_spectrum_base

!> Procedures spectrum rng uniform ------------------

!> allocate the spectrum_rng_uniform datatype
!> inputs:
!>   spectrum:  (spectrum_rng_uniform) generates and integrates
!>              variables along a uniform spectral distribution
!>   n_points:  (integer) number of spectral points
!>   n_spectra: (integer) number of spectral intervals
!>   real8_param: (real8)(2*n_spectra)(optional) minimum:
!>                1:n_spectral-> minimum wavelengths
!>                n_spectra+1:2*n_spectra->maximum wavelengths
!>   max_wlen:  (integer)(0)(optional) maximum wavelength
!> outputs:
!>   spectrum: (spectrum_rng_uniform) generates and integrates
!>             variables along a uniform spectral distribution
subroutine allocate_spectrum_rng_uniform(spectrum,n_points,&
n_spectra,real8_param,int_param)
  implicit none
  !> inputs
  integer,intent(in) :: n_points,n_spectra
  real*8,dimension(2*n_spectra),intent(in),optional  :: real8_param
  integer,dimension(0),intent(in),optional :: int_param
  !> inputs-outputs
  class(spectrum_rng_uniform),intent(inout) :: spectrum
  !> allocated all variables
  call allocate_spectrum_base(spectrum,n_points,n_spectra)
  if(present(real8_param)) &
  call spectrum%set_spectrum_interval(n_spectra,real8_param(1:n_spectra),&
  real8_param(n_spectra+1:2*n_spectra))
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
subroutine set_uniform_spectrum_interval(spectrum,n_spectra,min_wlen,max_wlen)
  implicit none
  !> inputs-outputs
  class(spectrum_rng_uniform),intent(inout) :: spectrum
  !> inputs
  integer,intent(in) :: n_spectra
  real*8,dimension(n_spectra),intent(in) :: min_wlen,max_wlen
  !> set values
  if(.not.allocated(spectrum%min_wlen)) allocate(spectrum%min_wlen(n_spectra))
  if(.not.allocated(spectrum%i_pdf))    allocate(spectrum%i_pdf(n_spectra))
  if(spectrum%n_spectra.ne.n_spectra) then
    deallocate(spectrum%min_wlen); deallocate(spectrum%i_pdf);
    allocate(spectrum%min_wlen(n_spectra))
    allocate(spectrum%i_pdf(n_spectra))
    spectrum%n_spectra = n_spectra
    write(*,*) 'WARNING: n_spectra is changed: regenerate spectrum points!'
  endif
  spectrum%min_wlen = min_wlen
  spectrum%i_pdf = max_wlen - min_wlen
end subroutine set_uniform_spectrum_interval

!> generate uniform random spectrum
!> inputs:
!>   spectrum: (spectrum_rng_uniform) generates and integrates 
!>             variables along a uniform spectral distribution
!>   rngs:     (type_rng)(n_omp_threads,n_spectra) random number generators
!> outputs:
!>   spectrum: (spectrum_rng_uniform) generates and integrates 
!>             variables along a uniform spectral distribution
subroutine generate_uniform_rng_spectrum(spectrum,rngs)
  use mod_rng
  !$ use omp_lib
  implicit none
  !> inputs-outpus
  class(spectrum_rng_uniform),intent(inout)   :: spectrum
  class(type_rng),dimension(:),allocatable,intent(inout) :: rngs
  !> variables
  integer :: ii,thread_id
  real*8,dimension(spectrum%n_spectra) :: rands
  real*8,dimension(spectrum%n_spectra,spectrum%n_points) :: local_points
  !> generate spectrum from uniform random number distribution
  thread_id = 0
  !$omp parallel default(private) shared(spectrum,rngs,local_points)
  !$ thread_id = omp_get_thread_num()
  !$omp do
  do ii=1,spectrum%n_points
    call rngs(thread_id+1)%next(rands)
    local_points(:,ii) = spectrum%min_wlen+spectrum%i_pdf*rands
  enddo
  !$omp end do
  !$omp end parallel
  spectrum%points = transpose(local_points)
  
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
  class(spectrum_rng_uniform),intent(inout) :: spectrum
  !> inputs
  real*8,dimension(spectrum%n_points,spectrum%n_spectra),intent(in) :: uniform_data
  !> outputs
  real*8,dimension(spectrum%n_spectra),intent(out) :: integrals
  !> variables
  integer :: ii
!$ integer :: jj
  !> integrate
#ifdef _OPENMP
  !$omp parallel do default(private) shared(spectrum,uniform_data) &
  !$omp collapse(2) reduction(+:integrals)
  do jj=1,spectrum%n_spectra
    do ii=1,spectrum%n_points
      integrals(jj) = integrals(jj) + uniform_data(ii,jj)
    enddo
  enddo
  !$omp end parallel do
#else
  do ii=1,spectrum%n_spectra
    integrals(ii) = sum(uniform_data(:,ii))
  enddo
#endif
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
