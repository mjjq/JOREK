!> This module allows to apply a Fourier analysis or a Fourier filter to the result array
!! of the routine eval_expr() from the module mod_expression.
!!
!! ### take into account n_period!!!
module mod_four_filter
  
  
  
  
  
  use iso_c_binding
  use parameters
  
  
  
  
  
  implicit none
  
  
  
  
  
#ifdef USE_FFTW
  include 'fftw3.f03'
#endif
  
  
  
  
  
  public
  private aux
  
  
  
  
  
  ! --- Constants
  character(len=15), parameter, private :: THIS_MOD_NAME     = 'mod_four_filter'
  integer,           parameter, private :: MAX_MODE_LIST_LEN = 300
  integer,           parameter, private :: PLUS_INF          = +99999
  
  
  
  
  
  !> Datatype specifying a certain mode number or range of mode numbers.
  type t_mode_specification
    integer :: m_start
    integer :: m_end
    integer :: n_start
    integer :: n_end
  end type t_mode_specification
  
  
  
  
  
  !> Datatype describing a Fourier filter, i.e., which harmonics to keep.
  type t_four_filter
    integer                    :: n_keep                         !< Length of keep_list
    type(t_mode_specification) :: keep_list(MAX_MODE_LIST_LEN)   !< Keep when applying filter.
  end type t_four_filter
  
  
  
  
  
  contains
  
  
  
  
  
  !> Perform a fast Fourier transform of the result array
  subroutine perform_four_trafo(result, forward)
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: result(:,:,:,:)
    logical,             intent(in)    :: forward
    
    ! --- Local variables
    integer   :: i, j
    integer*8 :: plan, n(4)
    
    !### use advanced fftw interface
    !### check result allocated
    !### allow to perform only 1d transforms
    
    ! --- Dimensionality of result array
    n(:) = (/ size(result,1), size(result,2), size(result,3), size(result,4) /)
    
    ! --- Normalization for back transform (Fourier to real space)
    !     (this is done in order to have really the coefficients of the sin and cos functions)
    if ( .not. forward ) then
      !### optimize this
      result(:,:,:,:)             = result(:,:,:,:) / 4
      result(1,:,:,:)             = result(1,:,:,:) * 2
      result(:,1,:,:)             = result(:,1,:,:) * 2
      result(n(1)/2+1:n(1),:,:,:) = -result(n(1)/2+1:n(1),:,:,:)
      result(:,n(2)/2+1:n(2),:,:) = -result(:,n(2)/2+1:n(2),:,:)
    end if
    
    ! --- Perform the Fast Fourier Transformation
    if ( forward ) then
      call dfftw_plan_r2r_2d(plan, n(1), n(2), result, result, FFTW_R2HC, FFTW_R2HC, FFTW_ESTIMATE)
    else
      call dfftw_plan_r2r_2d(plan, n(1), n(2), result, result, FFTW_HC2R, FFTW_HC2R, FFTW_ESTIMATE)
    end if
    do j = 1, size(result,4)
      do i = 1, size(result,3)
        call dfftw_execute_r2r(plan, result(:,:,i,j), result(:,:,i,j))
      end do
    end do
    call dfftw_destroy_plan(plan)
    
    ! --- Normalization for forward transform (real to Fourier space)
    !     (this is done to recover the original array after one forward and back transform)
    if ( forward ) then
      !### optimize this
      result(:,:,:,:)             = ( result(:,:,:,:) * 4 ) / ( n(2) * n(1) )
      result(1,:,:,:)             = result(1,:,:,:) / 2
      result(:,1,:,:)             = result(:,1,:,:) / 2
      result(n(1)/2+1:n(1),:,:,:) = -result(n(1)/2+1:n(1),:,:,:)
      result(:,n(2)/2+1:n(2),:,:) = -result(:,n(2)/2+1:n(2),:,:)
    end if
    
  end subroutine perform_four_trafo
  
  
  
  
  
  !> Initialize a Fourier filter data structure.
  !!
  !! The returned filter will remove all harmonics. Call filter_add to specify harmonics to keep.
  subroutine init_four_filter(filter)
    
    ! --- Routine parameters
    type (t_four_filter), intent(inout) :: filter
    
    filter%n_keep   = 0
    
  end subroutine init_four_filter
  
  
  
  
  !> Add a mode specification to the keep-list of a filter.
  subroutine filter_add(filter, ierr, m, n, m_start, m_end, n_start, n_end)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':eval_expr'
    
    ! --- Routine parameters
    type (t_four_filter), intent(inout) :: filter
    integer,              intent(inout) :: ierr
    integer, optional,    intent(in)    :: m, n, m_start, m_end, n_start, n_end
    
    integer :: m_start2, m_end2, n_start2, n_end2
    
    ierr = 0
    
    ! --- Poloidal mode number(s).
    if      ( (      present(m)) .and. (.not. present(m_start)) .and. (.not. present(m_end)) ) then
      m_start2 = m
      m_end2 = m
    else if ( (.not. present(m)) .and. (      present(m_start)) .and. (.not. present(m_end)) ) then
      m_start2 = m_start
      m_end2 = PLUS_INF
    else if ( (.not. present(m)) .and. (.not. present(m_start)) .and. (      present(m_end)) ) then
      m_start2 = 0
      m_end2 = m_end
    else if ( (.not. present(m)) .and. (      present(m_start)) .and. (      present(m_end)) ) then
      m_start2 = m_start
      m_end2 = m_end
    else if ( (.not. present(m)) .and. (.not. present(m_start)) .and. (.not. present(m_end)) ) then
      m_start2 = 0
      m_end2 = PLUS_INF
    else
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': Invalid combination of m, m_start, m_end.'
      ierr = 100
      return
    end if
    
    ! --- Toroidal mode number(s).
    if      ( (      present(n)) .and. (.not. present(n_start)) .and. (.not. present(n_end)) ) then
      n_start2 = n
      n_end2 = n
    else if ( (.not. present(n)) .and. (      present(n_start)) .and. (.not. present(n_end)) ) then
      n_start2 = n_start
      n_end2 = PLUS_INF
    else if ( (.not. present(n)) .and. (.not. present(n_start)) .and. (      present(n_end)) ) then
      n_start2 = 0
      n_end2 = n_end
    else if ( (.not. present(n)) .and. (      present(n_start)) .and. (      present(n_end)) ) then
      n_start2 = n_start
      n_end2 = n_end
    else if ( (.not. present(n)) .and. (.not. present(n_start)) .and. (.not. present(n_end)) ) then
      n_start2 = 0
      n_end2 = PLUS_INF
    else
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': Invalid combination of n, n_start, n_end.'
      ierr = 101
      return
    end if
    
    if ( filter%n_keep >= MAX_MODE_LIST_LEN ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': n_keep >= MAX_MODE_LIST_LEN.'
      ierr = 102
      return
    end if
    filter%n_keep = filter%n_keep + 1
    filter%keep_list(filter%n_keep)%m_start = m_start2
    filter%keep_list(filter%n_keep)%m_end = m_end2
    filter%keep_list(filter%n_keep)%n_start = n_start2
    filter%keep_list(filter%n_keep)%n_end = n_end2
    
  end subroutine filter_add
  
  
  
  
  
  !> Create a simple filter on the fly.
  function simple_filter(m, n, m_start, m_end, n_start, n_end) result(filter)
    type (t_four_filter) :: filter
    
    ! --- Routine parameters
    integer, optional,    intent(in)    :: m, n, m_start, m_end, n_start, n_end
    
    integer :: ierr
    
    call init_four_filter(filter)
    call filter_add(filter, ierr, m, n, m_start, m_end, n_start, n_end)
    
  end function simple_filter
  
  
  
  
  
  !> Print information on a Fourier filter for diagnostic purposes.
  subroutine print_filter(filter)
  
    ! --- Routine parameters
    type (t_four_filter), intent(in) :: filter
    
    ! --- Local variables
    integer :: i, m_start, m_end, n_start, n_end
    
    write(*,*)
    write(*,*) '=== FOURIER FILTER: KEEP THE FOLLOWING HARMONICS ===='
    do i = 1, filter%n_keep
      m_start = filter%keep_list(i)%m_start
      m_end = filter%keep_list(i)%m_end
      n_start = filter%keep_list(i)%n_start
      n_end = filter%keep_list(i)%n_end
      write(*,*) '(m,n) = ('//trim(aux(m_start,m_end))//', '//trim(aux(n_start,n_end))//')'
    end do
    write(*,*) '====================================================='
    write(*,*)
    
  end subroutine print_filter
  
  
  
  
  
  !> Auxilliary routine for print_filter: Returns mode number or range of mode numbers like 3...5.
  function aux(mi,ma) result(descr)
    character(len=32)   :: descr
    integer, intent(in) :: mi,ma
    
    character(len=16) :: s
    
    descr = ''
    if ( mi == ma ) then
      write(s,*) mi
      descr = adjustl(s)
      return
    end if
    
    write(s,*) mi
    descr = trim(descr) // trim(adjustl(s)) // '...'
    
    if ( ma == PLUS_INF ) then
      s = '+inf'
    else
      write(s,*) ma
    end if
    descr = trim(descr) // trim(adjustl(s))
    
  end function aux
  
  
  
  
  
  !> Apply a Fourier filter to a Fourier transformed result array.
  subroutine apply_four_filter(result, filter)
    
    ! --- Routine parameters
    real*8, allocatable,            intent(inout) :: result(:,:,:,:)
    type (t_four_filter),           intent(in)    :: filter
    
    ! --- Local variables
    real*8, allocatable :: filt_fact(:,:) !< Filter factors (currently 1 or 0)
    integer :: nn(4), num_pol_cos, num_pol_sin, num_tor_cos, num_tor_sin, m_max, n_max, i, m_start,&
      m_end, n_start, n_end, m, n
    
    nn(:) = (/ size(result,1), size(result,2), size(result,3), size(result,4) /)
    
    num_pol_cos = nn(2) / 2 + 1
    num_pol_sin = (nn(2)-1) / 2
    m_max       = num_pol_cos - 1
    
    num_tor_cos = nn(1) / 2 + 1
    num_tor_sin = (nn(1)-1) / 2
    n_max       = num_tor_cos - 1
    
    allocate( filt_fact(0:n_max,0:m_max) )
    filt_fact(:,:) = 0.d0 ! (remove all harmonics by default)
    
    do i = 1, filter%n_keep
      m_start = max(filter%keep_list(i)%m_start, 0    )
      m_end   = min(filter%keep_list(i)%m_end,   m_max)
      n_start = max(filter%keep_list(i)%n_start, 0    )
      n_end   = min(filter%keep_list(i)%n_end,   n_max)
      filt_fact(n_start:n_end,m_start:m_end) = 1.d0
    end do
    
    do m = 0, m_max
      do n = 0, n_max
        
        if ( filt_fact(n,m) == 1.d0 ) cycle ! nothing to be done in this case
        
        ! --- toroidal cos, poloidal cos component
        result(n+1,m+1,:,:) = result(n+1,m+1,:,:) * filt_fact(n,m)
        
        ! --- toroidal cos, poloidal sin component
        if ( (m /= 0) .and. (m <= num_pol_sin) ) then
          result(n+1,nn(2)+1-m,:,:) = result(n+1,nn(2)+1-m,:,:) * filt_fact(n,m)
        end if
        
        ! --- toroidal sin, poloidal cos component
        if ( (n /= 0) .and. (n <= num_tor_sin) ) then
          result(nn(1)+1-n,m+1,:,:) = result(nn(1)+1-n,m+1,:,:) * filt_fact(n,m)
        end if
        
        ! --- toroidal sin, poloidal sin component
        if ( (m /= 0) .and. (n /= 0 ) .and. (m <= num_pol_sin) .and. (n <= num_tor_sin) ) then
          result(nn(1)+1-n,nn(2)+1-m,:,:) = result(nn(1)+1-n,nn(2)+1-m,:,:) * filt_fact(n,m)
        end if
        
      end do
    end do
    
    deallocate(filt_fact)
    
  end subroutine apply_four_filter
  
  
  
  
  
  !> All-in-one-routine: Transform, filter, transform back.
  subroutine transform_and_filter(result, filter)
    
    ! --- Routine parameters
    real*8, allocatable,            intent(inout) :: result(:,:,:,:)
    type (t_four_filter),           intent(in)    :: filter
    
    call perform_four_trafo(result, .true.)
    call apply_four_filter(result, filter)
    call perform_four_trafo(result, .false.)
    
  end subroutine
  
  
  
  
  
end module mod_four_filter
