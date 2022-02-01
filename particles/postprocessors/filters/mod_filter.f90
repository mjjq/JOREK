!> mod_filter contains all variables and procedures
!> defining filters
module mod_filter
implicit none
private
public :: filter

!> Variable and data type definitions ----------------------
!> basic type for the filter class
type,abstract :: filter
  integer                          :: n_dimensions !< filter dimension
  integer,dimension(:),allocatable :: stencil_shape !< shape of the stencil
  contains
  procedure,pass(filter_inout)     :: allocate_filter
  procedure,pass(filter_inout)     :: deallocate_filter
  procedure(int_init_filter),pass(filter_inout),deferred  :: init_filter
  procedure(int_comp_fil_pos),pass(filter_inout),deferred :: compute_filter_from_position
  procedure(int_comp_fil_pos_vect),pass(filter_inout),deferred :: compute_filter_from_position_vectorial 
end type filter

!> Interfaces ----------------------------------------------
interface
  !> initialise the filter data types
  !> inputs:
  !>   filter_inout:     (filter) filter to be initialised
  !>   n_dimensions:     (integer) filter dimension
  !>   stencil_shape_in: (integer)(n_dimension)(optional) filter stencil shape
  !> outputs:
  !>   filter_inout:     (filter) initialised filter
  subroutine int_init_filter(filter_inout,n_dimensions,stencil_shape_in)
    IMPORT :: filter
    implicit none
    !> inputs-outputs
    class(filter),intent(inout) :: filter_inout
    !> inputs
    integer,intent(in) :: n_dimensions
    integer,dimension(n_dimensions),intent(in),optional :: stencil_shape_in
  end subroutine int_init_filter

  !> compute_filter_from_position compute the value of the filter
  !> given a normalised position
  !> inputs:
  !>   filter_inout:  (filter) filter type
  !>   pos:           (real8)(n_dimension) position where the filter weight is computed
  !> outputs:
  !>   filter_inout:  (filter) filter type
  !>   weight:        (real8) filter weight
  subroutine int_comp_fil_pos(filter_inout,pos,weight)
    IMPORT :: filter
    implicit none
    !> inputs-outputs
    class(filter),intent(inout) :: filter_inout
    !> inputs
    real*8,dimension(filter_inout%n_dimensions),intent(in) :: pos
    !> outputs:
    real*8,intent(out) :: weight
  end subroutine int_comp_fil_pos

  !> compute_filter_from_position compute the value of the filter
  !> given a normalised position: vectorial version.
  !> inputs:
  !>   filter_inout:  (filter) filter type
  !>   n_positions:   (integer) number of positions
  !>   pos:           (real8)(n_dimension,n_positions) positions where 
  !>                         the filter weights are computed
  !> outputs:
  !>   filter_inout:  (filter) filter type
  !>   weights:       (real8)(n_positions) filter weights
  subroutine int_comp_fil_pos_vect(filter_inout,n_positions,pos,weights)
    IMPORT :: filter
    implicit none
    !> inputs-outputs
    class(filter),intent(inout) :: filter_inout
    !> inputs
    integer,intent(in) :: n_positions
    real*8,dimension(filter_inout%n_dimensions,n_positions),intent(in) :: pos
    !> outputs:
    real*8,dimension(n_positions),intent(out) :: weights
  end subroutine int_comp_fil_pos_vect
end interface

contains
!> Procedures ----------------------------------------------
!> allocate filter variables
!> inputs:
!>   filter_inout: (filter) filter to be allocated
!>   n_dimensions: (integer) filter dimensionality (1D,2D,...)
!> outputs:
!>   filter_inout: (filter) allocate filter
subroutine allocate_filter(filter_inout,n_dimensions)
  implicit none
  class(filter),intent(inout) :: filter_inout
  integer,intent(in)          :: n_dimensions 
  if(allocated(filter_inout%stencil_shape)) deallocate(filter_inout%stencil_shape)
  allocate(filter_inout%stencil_shape(n_dimensions))
  filter_inout%stencil_shape = 0; filter_inout%n_dimensions = n_dimensions;
end subroutine allocate_filter

!> allocate filter variables
!> deallocate filter variables
!> inputs:
!>   filter_inout: (filter) filter to be deallocated
!> outputs:
!>   filter_inout: (filter) deallocate filter
subroutine deallocate_filter(filter_inout)
  implicit none
  class(filter),intent(inout) :: filter_inout
  if(allocated(filter_inout%stencil_shape)) deallocate(filter_inout%stencil_shape)
  filter_inout%n_dimensions = 0
end subroutine deallocate_filter

!>----------------------------------------------------------
end module mod_filter
