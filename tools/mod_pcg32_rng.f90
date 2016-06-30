!> Wrapper module for the pcg32 rng in the mod_rng type
module mod_pcg32_rng
use mod_pcg32
use mod_rng
implicit none

private

type, public, extends(type_rng) :: pcg32_rng
  private
    type(pcg_state_setseq_64), dimension(:), allocatable :: state
  contains
    procedure, public :: initialize => initialize_pcg32_rng
    procedure, public :: next => next_pcg32_rng
end type

contains
  subroutine initialize_pcg32_rng(rng, n_dims, seed, n_streams, i_stream, ifail)
    implicit none
    class(pcg32_rng), intent(inout) :: rng
    integer, intent(in)  :: n_dims !< Dimension of the generated output vector (not used)
    integer, intent(in)  :: seed !< Seed for the RNG if required
    integer, intent(in)  :: n_streams !< Number of output streams needed (not used)
    integer, intent(in)  :: i_stream !< Index of this output stream (sequence number)
    integer, intent(out) :: ifail !< Error code

    ifail = 0
    if (n_dims .le. 0) ifail = 1
    if (seed .eq. 0) ifail = 2
    if (ifail .ne. 0) return

    ! n_streams and i_stream are not used
    if (allocated(rng%state)) deallocate(rng%state)
    allocate(rng%state(1))

    call pcg32_srandom_r(rng%state(1), int(seed, 8), int(i_stream, 8))
  end subroutine initialize_pcg32_rng

  subroutine next_pcg32_rng(rng, out)
    implicit none
    class(pcg32_rng), intent(inout) :: rng
    real*8, dimension(:), intent(out) :: out
    call pcg32_random_doubles_r(rng%state(1), out)
  end subroutine next_pcg32_rng
end module mod_pcg32_rng
