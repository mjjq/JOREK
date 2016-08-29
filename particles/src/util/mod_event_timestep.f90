!> calculate timesteps and event times to match fixed timestep pushers
module mod_event_timestep
implicit none
contains

subroutine fix_event_timestep(pusher_timesteps, event_start, event_step, constraints, ierr)
  use mod_constants, only: TICK
  real*8, intent(inout), dimension(:)                                      :: pusher_timesteps !< steps of each of the used pushers
  real*8, intent(inout), dimension(:)                                      :: event_start
  real*8, intent(inout), dimension(size(event_start))                      :: event_step
  logical, intent(in), dimension(size(event_start),size(pusher_timesteps)) :: constraints !< whether event i constrains pusher j
  integer, intent(out) :: ierr
  integer :: i, j, k
  integer :: n, p !< dimension for GLSE problem. See http://www.netlib.org/lapack/lug/node28.html (with m=n)
  integer :: num_pushers, num_events

  real*8, dimension(:), allocatable   :: x !< timesteps of pushers and events
  real*8, dimension(:), allocatable   :: c, d !< constraints and reference value
  real*8, dimension(:,:), allocatable :: A, B_real, B, B_tmp
  real*8, dimension(:), allocatable :: work
  integer :: lwork, info
  integer, dimension(:), allocatable :: ipiv

  real*8, parameter :: tolerance = 1d-8 ! for comparing numbers of order 1000, no problem

  ierr = 0
  num_pushers = size(pusher_timesteps)
  num_events  = size(event_start)

  ! Set the number of variables, the weighting matrix A and the target vector c
  n = num_pushers + 2*num_events
  allocate(x(n), c(n), A(n,n))
  ! Set optimization parameters
  x = [pusher_timesteps, event_start, event_step] ! concatenate for initial guess
  c(:) = 1.d0 ! reference value = A x0 = 1
  A(:,:) = 0.d0 ! A contains a normalization by the current timestep size
  do i=1,n
    if (i .gt. num_pushers .and. i .le. num_pushers + num_events .and. abs(x(i)) .le. TICK) then
      ! if x(i) is zero normalize by the event step, and set c(i) to zero
      A(i,i) = 1.d0/x(i+num_events)
      c(i) = 0.d0
    else
      A(i,i) = 1.d0/x(i)
    end if
  end do


  ! convert the constraints into B_real
  p = 2*count(constraints)
  k = 0
  allocate(B_real(p,n))
  B_real(:,:) = 0.d0
  do i=1,size(constraints,1) ! i numbers the event
    do j=1,size(constraints,2) ! j numbers the pusher
      if (constraints(i,j)) then
        k = k+2
        ! equation: m * t_j = T_step,i (where m is the number of steps to fit)
        B_real(k-1,j) = event_step(i)/pusher_timesteps(j)
        B_real(k-1,num_pushers+num_events+i) = -1.d0
        ! equation: l * t_j = T_start,i (where l is the number of steps to fit)
        B_real(k,j) = event_start(i)/pusher_timesteps(j)
        B_real(k,num_pushers+i) = -1.d0
      end if
    end do
  end do

  ! TODO: heuristic algorithm to match timesteps (rows in the matrix) if needed (if p > n)
  ! so we can find a matrix B with rank p <= n so it has a solution

  allocate(B_tmp(p, n), ipiv(min(p,n)))
  B_tmp = real(nint(B_real),8) ! round all values
  call dgetrf(p, n, B_tmp, p, ipiv, info)
  if (info .lt. 0) then
    write(*,*) "ERROR: dgetrf info: ", info
  end if
  ! B_tmp now contains the LU factorisation of nint(B_real)
  ! copy these rows to B to get the row-reduced echelon form (i.e. a full-rank constraint matrix)

  if (info .gt. 0) then
    p = info - 1
  else
    p = min(p,n)
  end if
  ! test for the last zero (workaround)
  if (abs(B_tmp(p,p) - 0.d0) .le. tolerance) p = p-1

  allocate(B(p,n))
  B = 0.d0
  do i=1,p
    B(i,i:n) = B_tmp(i,i:n)
  end do
  allocate(d(p))
  d(:) = 0.d0

  if (.false.) then ! TODO add debug logging flag
    write(*,"(A,100f9.3)") "c=", c
    write(*,"(A,100f9.3)") "d=", d
    write(*,"(A,100f9.3)") "A(i,i)=", [(A(i,i), i=1, size(A,1))]
    do i=1,size(B_real,1)
      write(*,"(A,i1,A,100f8.4)") "B1(",i,",:)=", real(nint(B_real(i,:)),8)
    end do
    do i=1,p
      write(*,"(A,i1,A,100f8.4)") "B2(",i,",:)=", B(i,:)
    end do
    write(*,*) "x0=", x
  end if

  ! check if the system is solvable (don't think this will occur at all)
  if (p .gt. n) then
    write(*,"(A,i2,A,i2,A)") "ERROR: Too many constraints (", p, ">", n, ") to find a proper timestep!"
    ierr = 3
    return
  else if (p .eq. 0) then ! no work to do
    return
  end if

  ! Let lapack solve the system (alters A, c, d,x,work)
  ! get the optimum size of the work array
  allocate(work(1))
  call dgglse(n,n,p,A,n,B,p,c,d,x,work,-1,info)
  lwork = nint(work(1))
  deallocate(work);allocate(work(lwork))
  call dgglse(n,n,p,A,n,B,p,c,d,x,work,lwork,info)
  if (info .ne. 0) then
    write(*,*) "ERROR: dgglse info: ", info
    ierr = info
  else
    ! Save values
    if (.false.) write(*,*) "x=", x
    do i=1,n
      if (i .le. num_pushers) then
        pusher_timesteps(i) = x(i)
      else if (i .le. num_pushers + num_events) then
        event_start(i-num_pushers) = x(i)
      else
        event_step(i-num_pushers-num_events) = x(i)
      end if
    end do

    ! if debug: verify whether all timesteps 'fit' in integer values into events
  end if
end subroutine fix_event_timestep
end module mod_event_timestep
