!> the mod_pinhole_lens contains variables and procedures defining
!> the pinhole lens model
module mod_pinhole_lens
use mod_lens, only: lens
implicit none

private
public :: pinhole_lens

!> Variables ----------------------------------------------------------
type,extends(lens) :: pinhole_lens
  contains
  procedure,pass(lens_inout) :: init_pinhole
  procedure,pass(lens_inout) :: lens_sampling => pinhole_sampling
  procedure,pass(lens_inout) :: lens_pdf      => pinhole_pdf
end type pinhole_lens

!> Interfaces ---------------------------------------------------------
contains

!> Procedures ---------------------------------------------------------
!> initialise the pinhole camera
!> inputs:
!>   lens_inout:     (pinhole_lens) pinhole lens to be initialised
!>   n_x:            (integer) size of the pinhole position array
!>   pinhole_center: (n_x) position of the pinhole cartesian coord.
!> outputs:
!>   lens_inout: (pinhole_lens) initialised pinhole lens
subroutine init_pinhole(lens_inout,n_x,pinhole_center)
  implicit none
  !> inputs-outputs
  class(pinhole_lens),intent(inout) :: lens_inout
  !> inputs
  integer,intent(in) :: n_x
  real*8,dimension(n_x),intent(in) :: pinhole_center
  call lens_inout%allocate_lens(n_x,pinhole_center)
end subroutine init_pinhole

!> return a set of samples on the pinhole lens (the pinhole position)
!> inputs:
!>   lens_inout: (pinhole_lens) the pinhole lens class
!>   n_samples:  (integer) number of samples to be generated
!> outputs:
!>   lens_inout: (pinhole_lens) the pinhole lens class
!>   x_pos:      (n_x,n_samples) arrays of the sampled lens points
subroutine pinhole_sampling(lens_inout,n_samples,x_pos)
  implicit none
  !> inputs-outputs
  class(pinhole_lens),intent(inout) :: lens_inout
  !> inputs
  integer,intent(in) :: n_samples
  !> outputs:
  real*8,dimension(lens_inout%n_x,n_samples),intent(out) :: x_pos
  !>
  integer :: ii
  !$omp parallel do default(shared) firstprivate(n_samples) private(ii)
  do ii=1,n_samples
    x_pos(:,ii) = lens_inout%center
  enddo
  !$omp end parallel do
end subroutine pinhole_sampling

!> return return the pdf for a set of sampled positions
!> Due to the fact that the pdf of a pinhole camera is a 
!> delta-Dirac distribution, 1.d0 is returned by default
!> inputs:
!>   lens_inout: (pinhole_lens) the pinhole lens class
!>   n_samples:  (integer) number of samples to be generated
!>   x_pos:      (n_x,n_samples) arrays of the sampled lens points
!> outputs:
!>   lens_inout: (pinhole_lens) the pinhole lens class
!>   pdf:        (n_samples) arrays contaning the sample pdfs
subroutine pinhole_pdf(lens_inout,n_samples,x_pos,pdf)
  implicit none
  !> inputs-outputs
  class(pinhole_lens),intent(inout) :: lens_inout
  !> inputs
  integer,intent(in) :: n_samples
  real*8,dimension(lens_inout%n_x,n_samples),intent(in)  :: x_pos
  !> outputs:
  real*8,dimension(n_samples),intent(out) :: pdf
  !>
  integer :: ii
  pdf = 1.d0
end subroutine pinhole_pdf 

!>---------------------------------------------------------------------
end module mod_pinhole_lens
