!> calculate timesteps and event times to match fixed timestep pushers
module mod_event_timestep
implicit none
private
public fix_event_timestep
contains

!> Calculate the optimum times for pushers, event start and event steps to fit constraints
!> while staying as close as possible to the original values.
!>
!> Each of the constraints represents two equations:
!> \[ k*t_i = T_{start,j} \] i.e. the pusher timestep \(t_i\) must fit an integer
!> \(k\) number of times between \(0\) and \(T_{start,j}\)
!> and \[ l*t_i = T_{step,j} \]
!> i.e. the pusher timestep \(t_i\) must fit an integer $l$ number of times
!> between \(0\) and \(T_{step,j}\) where pushers are numbered with \(i\)
!> and events with \(j\). \(k\) and \(l\) must be some integer value, determined by dividing the 
!> requested timesteps and rounding.
!>
!> We gather the timesteps into a vector \(x = \) [pusher_timesteps, event_start, event_step] of
!> size \(n\). The resulting system of constraints we can write as \(Bx=d\) where \(B\) is a full-rank
!> matrix of size \(p\) by \(n\) and \(d = 0\) (and has size \(p\)).
!> To calculate \(B\) the LAPACK routine DGETRF is used on the full constraints
!> matrix \(B_{tmp}\) (which has size \(k\) by \(n\)).
!> The upper triangular part returned is in the row-reduced echelon form.
!> From this, rows having a nonzero diagonal element are selected into \(B\).
!> 
!> The problem reduces to a constrained linear least-squares optimization problem
!> which can be solved by the LAPACK routine DGGLSE. This minimizes \(|c-Ax|_2\) subject to \(Bx=d\).
!> The weight matrix \(A\) is given by \(diag(1/x)\), i.e. 1/x put on the diagonal of \(A\).
!> This ensures the minimization of the relative change. There are two special cases to consider
!> here, the one where an initial value is 0 and where it is huge.
!>
!> For the first case, which should only occur for the event_start time, we take event_step
!> as normalization. In the second case the weight factor is close to zero, and the result for
!> this event_step is ignored.
subroutine fix_event_timestep(pusher_timesteps, event_start, event_step, constraints, ierr)
  real*8, parameter :: TICK = 1d-12
  real*8, intent(inout), dimension(:)                                      :: pusher_timesteps !< steps of each of the used pushers
  real*8, intent(inout), dimension(:)                                      :: event_start
  real*8, intent(inout), dimension(size(event_start))                      :: event_step
  logical, intent(in), dimension(size(event_start),size(pusher_timesteps)) :: constraints !< whether event i constrains pusher j
  integer, intent(out) :: ierr
  integer :: i, j, k
  integer :: n, p !< dimension for GLSE problem. See http://www.netlib.org/lapack/lug/node28.html (with m=n)
  integer :: num_pushers, num_event_start, num_event_step
  integer :: i_constrain_start, i_constrain_step

  real*8, dimension(:), allocatable   :: x !< timesteps of pushers and events
  real*8, dimension(:), allocatable   :: c, d !< constraints and reference value
  real*8, dimension(:,:), allocatable :: A, B_real, B, B_tmp
  real*8, dimension(:), allocatable :: work
  integer :: lwork, info
  integer, dimension(:), allocatable :: ipiv

  real*8, parameter :: tolerance = 1d-8 ! for comparing numbers of order 1000, no problem

  ierr = 0
  num_pushers     = size(pusher_timesteps)
  num_event_start = count(abs(event_start-0.d0) .gt. TICK) ! number of nonzeros here
  num_event_step  = count(event_step .lt. huge(event_step)*(1.d0-tolerance)) ! number of < huge numbers

  ! Matrix construction:
  ! All timesteps can be changed by default.
  ! Starting times can be changed if they are nonzero
  ! Step times can be changed if they are not huge(0.d0)

  ! Set the number of variables, the weighting matrix A and the target vector c
  ! count the number of non-huge event_steps
  n = num_pushers + num_event_start + num_event_step
  allocate(x(n), c(n), A(n,n))
  x = [pusher_timesteps, pack(event_start, mask=abs(event_start-0.d0) .gt. TICK), &
      pack(event_step, mask=event_step .lt. huge(event_step)*(1.d0-tolerance))]

  ! Weight matrix and reference
  c(:) = 1.d0 ! reference value = A x0 = 1
  A(:,:) = 0.d0 ! A contains a normalization by the current timestep size
  do i=1,n
    A(i,i) = 1.d0/x(i)
  end do

  ! convert the constraints into B_real
  k = 2*count(constraints) ! maximum number of constraints
  p = 0
  allocate(B_real(k,n))
  B_real(:,:) = 0.d0
  do i=1,size(constraints,1) ! i numbers the event
    do j=1,size(constraints,2) ! j numbers the pusher
      if (constraints(i,j)) then
        ! the index below is given by the number of constrained events before this
        ! plus the number of pushers, (optionally plus the number of start constraints, for step constraints)
        if (event_step(i) .lt. huge(event_step(i))*(1.d0-tolerance)) then
          p = p+1
          i_constrain_step  = count(event_step(1:i) .lt. huge(event_step(1:i))*(1.d0-tolerance))
          ! equation: m * t_j = T_step,i (where m is the number of steps to fit)
          B_real(p,j) = event_step(i)/pusher_timesteps(j)
          B_real(p,num_pushers+num_event_start+i_constrain_step) = -1.d0
        end if
        if (abs(event_start(i)-0.d0) .gt. TICK) then
          p = p+1
          i_constrain_start = count(abs(event_start(1:i)-0.d0) .gt. TICK)
          ! equation: l * t_j = T_start,i (where l is the number of steps to fit)
          B_real(p,j) = event_start(i)/pusher_timesteps(j)
          B_real(p,num_pushers+i_constrain_start) = -1.d0
        end if
      end if
    end do
  end do

  ! TODO: heuristic algorithm to match timesteps (rows in the matrix) if needed (if p > n)
  ! so we can find a matrix B with rank p <= n so it has a solution


  allocate(B_tmp(k, n), ipiv(min(k,n)))
  B_tmp = real(nint(B_real),8) ! round all values
  call dgetrf(k, n, B_tmp, k, ipiv, info)
  if (info .lt. 0) then
    write(*,*) "ERROR: dgetrf info: ", info
  end if
  ! B_tmp now contains the LU factorisation of nint(B_real)
  ! copy these rows to B to get the row-reduced echelon form (i.e. a full-rank constraint matrix)

  ! count the number of non-zero diagonal elements in k
  p = 0
  do i=1,min(k,n)
    if (abs(B_tmp(i,i)) .gt. tolerance) p = p+1
  end do

  allocate(B(p,n))
  B = 0.d0; p = 0
  do i=1,min(k,n)
    if (abs(B_tmp(i,i)) .gt. tolerance) then ! do not copy rows that have a zero on the diagonal
      p = p+1
      B(p,i:n) = B_tmp(i,i:n)
    end if
  end do
  allocate(d(p))
  d(:) = 0.d0

  if (.false.) then ! TODO add debug logging flag
    write(*,"(A,100g10.3)") "c=", c
    write(*,"(A,100g10.3)") "d=", d
    write(*,"(A,100g10.3)") "A(i,i)=", [(A(i,i), i=1, size(A,1))]
    do i=1,size(B_real,1)
      write(*,"(A,i1,A,100g10.3)") "B1(",i,",:)=", real(nint(B_real(i,:)),8)
    end do
    do i=1,p
      write(*,"(A,i1,A,100g10.3)") "B2(",i,",:)=", B(i,:)
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
  if (info .ne. 0) then
    write(*,*) "ERROR: dgglse setup info: ", info
    ierr = info
    return
  end if
  lwork = nint(work(1))
  deallocate(work);allocate(work(lwork))
  call dgglse(n,n,p,A,n,B,p,c,d,x,work,lwork,info)
  if (info .ne. 0) then
    write(*,*) "ERROR: dgglse info: ", info
    ierr = info
  else
    ! Save values
    if (.false.) write(*,*) "x=", x
    do i=1,num_pushers
      pusher_timesteps(i) = x(i)
    end do
    j=1
    do i=1,num_event_start
      if (abs(event_start(j)-0.d0) .gt. TICK) then
        event_start(j) = x(i+num_pushers)
        j = j+1
      end if
    end do
    j=1
    do i=1,num_event_step
      if (event_step(i) .lt. huge(event_step(i))*(1.d0-tolerance)) then ! if it is not huge
        event_step(j) = x(i+num_pushers+num_event_start)
        j = j+1
      end if
    end do
    ! if debug: verify whether all timesteps 'fit' in integer values into events
  end if
end subroutine fix_event_timestep
end module mod_event_timestep
