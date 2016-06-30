!> Module containing abstract type for multi-dimensional (Q)RNGs
!! subclass type_rng to implement other generators
module mod_rng
  implicit none

  type, abstract :: type_rng
    contains
      procedure (initialize), public, deferred :: initialize
      procedure (next),       public, deferred :: next
  end type

  interface
    subroutine initialize(rng, n_dims, seed, n_streams, i_stream, ifail)
      import :: type_rng
      implicit none
      class(type_rng), intent(inout) :: rng
      integer, intent(in)  :: n_dims !< Dimension of the generated output vector
      integer, intent(in)  :: seed !< Seed for the RNG if required
      integer, intent(in)  :: n_streams !< Number of output streams needed
      integer, intent(in)  :: i_stream !< Index of this output stream
      integer, intent(out) :: ifail !< Error code
    end subroutine initialize

    subroutine next(rng, out)
      import :: type_rng
      implicit none
      class(type_rng), intent(inout) :: rng
      real*8, dimension(:), intent(out) :: out
    end subroutine next
  end interface
end module mod_rng
