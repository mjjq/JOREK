module mod_sobseq_rng
  use mod_sobseq
  use mod_rng
  implicit none
  private

  type, public, extends(type_rng) :: sobseq_rng
    private
      type(sobol_state), dimension(:), allocatable :: state
    contains
      procedure, public :: initialize => initialize_sobseq_rng
      procedure, public :: next => next_sobseq_rng
  end type

  ! Direction numbers taken from Joe and Kuo
  integer, parameter, dimension(1:8)   :: s = (/1,2,3,3,4,4,5,5/)
  integer, parameter, dimension(1:8)   :: a = (/0,1,1,2,1,4,2,4/)
  integer, parameter, dimension(5,1:8) :: m = reshape((/1,0,0,0,0, &
                                                1,3,0,0,0, &
                                                1,3,1,0,0, &
                                                1,1,1,0,0, &
                                                1,1,3,3,0, &
                                                1,3,5,13,0,&
                                                1,1,5,5,17,&
                                                1,1,5,5,5/), (/5,8/))
  ! TODO why is there no 1?
contains
  subroutine initialize_sobseq_rng(rng, n_dims, seed, n_streams, i_stream, ifail)
    implicit none
    class(sobseq_rng), intent(inout) :: rng
    integer, intent(in)  :: n_dims !< Dimension of the generated output vector (not used)
    integer, intent(in)  :: seed !< Seed for the RNG if required
    integer, intent(in)  :: n_streams !< Number of output streams needed (not used)
    integer, intent(in)  :: i_stream !< Index of this output stream (sequence number)
    integer, intent(out) :: ifail !< Error code

    integer :: i
    real*8 :: DUMMY_REAL

    ifail = 0
    if (n_dims .le. 0) ifail = 1
    if (seed .eq. 0) ifail = 2
    if (n_streams .le. 0) ifail = 3
    if (i_stream .le. 0) ifail = 4
    if (i_stream .gt. n_streams) ifail = 5
    if (ifail .ne. 0) return

    if (allocated(rng%state)) deallocate(rng%state)
    allocate(rng%state(n_dims))


    ! Seed with default values
    do i=1,n_dims
      call rng%state(i)%initialize(s(i), a(i), m(:,i), stride=ilog2_b(n_streams))
      DUMMY_REAL = rng%state(i)%skip_ahead(i_stream)
    enddo
  end subroutine initialize_sobseq_rng

  subroutine next_sobseq_rng(rng, out)
    implicit none
    class(sobseq_rng), intent(inout)  :: rng
    real*8, dimension(:), intent(out) :: out

    integer :: i

    do i=1,size(rng%state,1)
      out(i) = rng%state(i)%next()
    enddo
  end subroutine next_sobseq_rng



  !> Integer logarithm in base 2
  function ilog2_b(val) result(res)
    integer, intent(in) :: val
    integer             :: res
    integer             :: tmp

    res = -1
    ! Negative values not allowed
    if (val < 1) return

    tmp = val
    do while (tmp > 0)
      res = res + 1
      tmp = shiftr(tmp, 1)
    enddo
  end function ilog2_b
end module mod_sobseq_rng
