!> This module allows to apply a Fourier analysis or a Fourier filter to the result array
!! of the routine eval_expr() from the module mod_expression.
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
  integer,           parameter          :: TOROIDAL_TRAFO    = 0
  integer,           parameter          :: POLOIDAL_TRAFO    = 1
  integer,           parameter          :: POLTOR_TRAFO      = 2
  logical,           parameter          :: FORWARD_TRAFO     = .true.
  logical,           parameter          :: BACKWARD_TRAFO    = .false.
  
  
  
  
  
  !> Datatype specifying a certain mode number or range of mode numbers.
  type t_mode_specification
    integer :: m_start
    integer :: m_end
    integer :: n_start
    integer :: n_end
  end type t_mode_specification
  
  
  
  
  
  !> Datatype describing a Fourier filter, i.e., which harmonics to keep.
  type t_four_filter
    integer                    :: n_keep = 0                     !< Length of keep_list
    type(t_mode_specification) :: keep_list(MAX_MODE_LIST_LEN)   !< Keep when applying filter.
  end type t_four_filter
  
  
  
  
  
  contains
  
  
  
  
  
  !> Perform a fast Fourier transform of the result array
  subroutine perform_pol_trafo(re, cp, ierr)
    
    ! --- Routine parameters
    real*8,     allocatable, intent(inout) :: re(:,:,:,:)
    complex*16, allocatable, intent(inout) :: cp(:,:,:,:)
    integer,                 intent(inout) :: ierr
    
    call perform_four_trafo(re, cp, POLOIDAL_TRAFO, ierr)
    
  end subroutine perform_pol_trafo
  
  
  
  
  
  !> Perform a fast Fourier transform of the result array
  subroutine perform_tor_trafo(re, cp, ierr)
    
    ! --- Routine parameters
    real*8,     allocatable, intent(inout) :: re(:,:,:,:)
    complex*16, allocatable, intent(inout) :: cp(:,:,:,:)
    integer,                 intent(inout) :: ierr
    
    call perform_four_trafo(re, cp, TOROIDAL_TRAFO, ierr)
    
  end subroutine perform_tor_trafo
  
  
  
  
  
  !> Perform a fast Fourier transform of the result array
  subroutine perform_poltor_trafo(re, cp, ierr)
    
    ! --- Routine parameters
    real*8,     allocatable, intent(inout) :: re(:,:,:,:)
    complex*16, allocatable, intent(inout) :: cp(:,:,:,:)
    integer,                 intent(inout) :: ierr
    
    call perform_four_trafo(re, cp, POLTOR_TRAFO, ierr)
    
  end subroutine perform_poltor_trafo
  
  
  
  
  
  !> Perform a fast Fourier transform of the result array
  subroutine perform_four_trafo(re, cp, trafo_type, ierr)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':perform_four_trafo'
    
    ! --- Routine parameters
    real*8,     allocatable, intent(inout) :: re(:,:,:,:)
    complex*16, allocatable, intent(inout) :: cp(:,:,:,:)
    integer,                 intent(in)    :: trafo_type !< see constants TOROIDAL_TRAFO etc.
    integer,                 intent(inout) :: ierr
    
    ! --- Local variables
    integer   :: i, j, k
    integer*8 :: plan, n(4)
    logical   :: forward
    real*8, allocatable :: tmp(:)
    
    ierr = 0
    
#ifdef USE_FFTW
    
    if ( allocated(re) .and. (.not. allocated(cp)) ) then
      forward = .true.
    else if ( (.not. allocated(re)) .and. allocated(cp) ) then
      forward = .false.
    else
      ierr = 100
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': Exactly one of the arrays re and cp ' // &
        'must be allocated (re for forward transform; cp for backward transform).'
      return
    end if
    
    if ( trafo_type == POLTOR_TRAFO ) then
      
      if ( forward ) then
        
        n(:) = (/ size(re,1), size(re,2), size(re,3), size(re,4) /)
        allocate( cp(n(1)/2+1,n(2),n(3),n(4)) )
        
        call dfftw_plan_dft_r2c_2d(plan, n(1), n(2), re, cp, FFTW_ESTIMATE)
        
        do i = 1, n(4)
          do j = 1, n(3)
            call dfftw_execute_dft_r2c(plan, re(:,:,j,i), cp(:,:,j,i))
          end do
        end do
        
        cp(:,:,:,:) = cp(:,:,:,:) / ( n(1) * n(2) )
        
        call dfftw_destroy_plan(plan)
        
      else ! backward
        
        n(:) = (/ size(cp,1), size(cp,2), size(cp,3), size(cp,4) /)
        allocate( re(2*(n(1)-1),n(2),n(3),n(4)) )
        
        call dfftw_plan_dft_r2c_2d(plan, n(1), n(2), cp, re, FFTW_ESTIMATE)
        
        do i = 1, n(4)
          do j = 1, n(3)
            call dfftw_execute_dft_r2c(plan, cp(:,:,j,i), re(:,:,j,i))
          end do
        end do
        
        call dfftw_destroy_plan(plan)
        
      end if
      
    else if ( trafo_type == TOROIDAL_TRAFO ) then
      
      stop!###
      
    else if ( trafo_type == POLOIDAL_TRAFO ) then
      
      stop!###
      
    else
      
      ierr = 100
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': parameter trafo_type has illegal value'
      return
      
    end if
    
#else
#error "You need fftw library and to define -DUSE_FFTW"
#endif
  
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
  
  
  
  
  
  !> Apply a Fourier filter.
  subroutine apply_four_filter(result, filter, n_skip, ierr)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':apply_four_filter'
    
    ! --- Routine parameters
    real*8, allocatable,        intent(inout) :: result(:,:,:,:)
    type (t_four_filter),       intent(in)    :: filter
    integer,                    intent(in)    :: n_skip !< Don't filter first n_skip expressions
    integer,                    intent(inout) :: ierr
    
    ! --- Local variables
    real*8,     allocatable :: filt_fact(:,:) !< Filter factors (currently 1 or 0)
    complex*16, allocatable :: cp(:,:,:,:)
    integer :: i, j, m, n_in_period, n, m_max, n_max, n_in_period_max, nn(4), m_start, m_end,      &
      n_start, n_end
    
    ierr = 0
    
    if ( .not. allocated(result) ) then
      ierr = 100
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': array result not allocated'
      return
    end if
    
    nn(:) = (/ size(result,1), size(result,2), size(result,3), size(result,4) /)
    
    ! --- Create Fourier filter (filter factors)
    m_max = nn(2)/2
    n_in_period_max = nn(1)/2
    n_max = n_in_period_max * n_period
    
    if ( allocated(filt_fact) ) deallocate(filt_fact)
    allocate( filt_fact(0:n_max+1,-m_max-1:m_max+1) )
    filt_fact(:,:) = 0.d0 ! (remove all harmonics by default)
    
    do i = 1, filter%n_keep
      m_start = max(filter%keep_list(i)%m_start, -m_max)
      m_end   = min(filter%keep_list(i)%m_end,   +m_max)
      n_start = max(filter%keep_list(i)%n_start,  0    )
      n_end   = min(filter%keep_list(i)%n_end,    n_max)
      filt_fact(n_start:n_end,m_start:m_end) = 1.d0
    end do
    
    ! --- Transform data
    call perform_four_trafo(result, cp, POLTOR_TRAFO, ierr)
    
    ! --- Apply the Fourier filter (multiply by filter factors)
    do m = -m_max, m_max
      do n_in_period = 0, n_in_period_max
        
        n = n_in_period * n_period ! toroidal mode number
        
        write(37,'(18i7)') m, n, m_max, n_max, size(filt_fact,2), size(filt_fact,1)
        
        if ( filt_fact(n,m) == 1.d0 ) cycle ! nothing to be done in this case
        
        if ( m >= 0 ) then
          cp(n_in_period+1,m+1,:,n_skip+1:nn(4))       =                                           &
            cp(n_in_period+1,m+1,:,n_skip+1:nn(4))       * filt_fact(n,m)
        else ! m < 0
          cp(n_in_period+1,nn(2)+1+m,:,n_skip+1:nn(4)) =                                           &
            cp(n_in_period+1,nn(2)+1+m,:,n_skip+1:nn(4)) * filt_fact(n,m)
        end if
        
      end do
    end do
    
    ! --- Transform data back
    deallocate( result )
    call perform_four_trafo(result, cp, POLTOR_TRAFO, ierr)
    
    ! --- Clean up
    if ( allocated(filt_fact) ) deallocate(filt_fact)
    if ( allocated(cp) ) deallocate(cp)
    
  end subroutine apply_four_filter
  
  
  
  
  
end module mod_four_filter
