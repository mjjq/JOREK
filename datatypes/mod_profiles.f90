module profiles
!-----------------------------------------------------------------------
! Routines for reading/writing/manipulating 1D data profiles
!-----------------------------------------------------------------------
  
  
  
  implicit none
  
  
  
  private
  
  ! --- Public routines
  public constructProf, destructProf, resizeProf, readProf, writProf, derivProf
  
  
  
  contains
  
  
  
  ! Construct a profile.
  recursive subroutine constructProf(x, y, len)

    real, allocatable, intent(inout) :: x(:), y(:)
    integer,           intent(inout) :: len
    
    call destructProf(x, y, len)
    
    len = max(len,1) ! at least length 1
    allocate( x(len), y(len) )
    
  end subroutine constructProf
  
  
  
  ! Change the size of a profile.
  recursive subroutine resizeProf(x, y, len, newLen, keep)

    real, allocatable, intent(inout) :: x(:), y(:)
    integer,           intent(inout) :: len
    integer,           intent(in)    :: newLen
    logical,           intent(in)    :: keep ! Keep the data in the x and y arrays?
    
    real, ALLOCATABLE            :: px(:) ! copy of x (in case keep=.TRUE.)
    real, ALLOCATABLE            :: py(:) ! copy of y (in case keep=.TRUE.)
    
    ! --- Recursive call with newLen=1 if newLen < 1.
    if ( newLen < 1 ) then
      call resizeProf(x, y, len, 1, keep)
      return
    end if
    
    ! --- Backup data from profile if keep=.true.
    if ( keep ) then
      allocate( px(len), py(len) )
      if ( allocated(x) ) then
        px(1:len) = x(1:len)
      else
        px = 0.
      end if
      if ( allocated(y) ) then
        py(1:len) = y(1:len)
      else
        py = 0.
      end if
    end if
    
    ! --- Resize x and y.
    if ( allocated(x) ) deallocate(x)
    if ( allocated(y) ) deallocate(y)
    allocate( x(newLen), y(newLen) )
    
    ! --- Restore data to profile if keep=.true.
    if ( keep ) then
      x(1:min(len,newLen)) = px(1:min(len,newLen))
      y(1:min(len,newLen)) = py(1:min(len,newLen))
      deallocate( px, py )
    end if
    len = newLen
    
  end subroutine resizeProf
  
  
  
  ! Destroy a profile.
  recursive subroutine destructProf(x, y, len)

    real, allocatable, intent(inout) :: x(:), y(:)
    integer,           intent(inout) :: len
    
    if ( allocated(x) ) deallocate(x)
    if ( allocated(y) ) deallocate(y)
    len = 0
    
  end subroutine destructprof
  
  
  
  ! Read a profile from a file.
  recursive subroutine readProf(x, y, len, file)
    
    real, allocatable, intent(inout) :: x(:), y(:)
    integer,           intent(inout) :: len
    CHARACTER(LEN=*),  intent(in)    :: file    ! Filename.
    
    integer :: err
    integer :: usedLen
    real    :: xx, yy
    
    call destructProf(x, y, len)
    
    ! --- Open the file.
    OPEN(UNIT=42, FILE=file, FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
    if ( err /= 0 ) then
      write(*,*) 'ERROR in readProf: Cannot open file '//TRIM(file)//'.'
      return
    end if
    
    ! --- Construct prof with an initial length of 10.
    len = 10
    call constructProf(x, y, len)
    usedLen = 0
    
    ! --- Read profile.
    do
      READ(42, *, IOSTAT=err) xx, yy
      
      if ( err /= 0 ) exit ! end of profile reached
      
      usedLen = usedLen + 1
      
      ! Double the profile length if it becomes to small.
      if ( usedLen > len ) call resizeProf(x, y, len, 2*len, .TRUE.)
      
      x(usedLen) = xx
      y(usedLen) = yy
      
    end do
    
    ! --- Crop the profile to the length that is really used.
    call resizeProf(x, y, len, usedLen, .TRUE.)
    
    ! --- Close the file.
    CLOSE(UNIT=42)

  end subroutine readProf
  
  
  
  ! Write a profile to a file.
  recursive subroutine writProf(x, y, len, file)
    
    real, allocatable, intent(in)    :: x(:), y(:)
    integer,           intent(in)    :: len
    CHARACTER(LEN=*),  intent(in)    :: file    ! Filename.
    
    integer :: err
    integer :: i
    
    if ( (.not. allocated(x)) .or. (.not. allocated(y)) ) return
    
    ! --- Open the output file.
    OPEN(UNIT=42, FILE=file, FORM='FORMATTED', STATUS='REPLACE', ACTION='write', IOSTAT=err)
    if ( err /= 0 ) return
    
    ! --- Write data.
    do i = 1, len
      write(42, '(ES22.15,1X,ES22.15)') x(i), y(i)
    end do
    
    ! --- Close the output file.
    CLOSE(UNIT=42)

  end subroutine writProf
  
  
  
  
  ! Determine the derivative of a profile.
  recursive subroutine derivProf(x, y, len, yd)
  
    real, allocatable, intent(in)    :: x(:), y(:)
    real, allocatable, intent(inout) :: yd(:)
    integer,           intent(inout) :: len
    
    integer :: i       ! Point index of profile.
    real    :: d(-2:4) ! Distances between points.
    real    :: f(-2:4) ! Function values.
    real    :: c(-2:4) ! Coefficients for function values.
    
    ! The derivatives will be determined at the same x-positions as the profile.
    if ( allocated(yd) ) deallocate(yd)
    allocate(yd(len)) 
    
    do i = 1, len
      
      c = 0.
      f = 0.
      f(0) = y(i)
      
      if ( (i==1) .or. (i==len) ) then
        
        if ( i==1 ) then
          f(0) = y(i)
          f(1) = y(2)
          f(2) = y(3)
          
          d(1) = x(2) - x(i)
          d(2) = x(3) - x(i)
        else
          f(0) = y(i)
          f(1) = y(len-1)
          f(2) = y(len-2)
          
          d(1) = x(len-1) - x(i)
          d(2) = x(len-2) - x(i)
        end if
        c(0) = -(1/d(1)) - 1/d(2)
        c(1) = -(d(2)/(d(1)**2 - d(1)*d(2)))
        c(2) = d(1)/((d(1) - d(2))*d(2))
        
      else
        
        f(-1) = y(i-1)
        f(0)  = y(i)
        f(1)  = y(i+1)
        
        d(-1) = x(i-1) - x(i)
        d(1)  = x(i+1) - x(i)
        
        c(-1) = -(d(1)/(d(-1)**2 - d(-1)*d(1)))
        c(0)  = -(1/d(-1)) - 1/d(1)
        c(1)  = d(-1)/((d(-1) - d(1))*d(1))

      end if
      
      yd(i) = sum( c * f )
        
    end do
    
  end subroutine derivProf


  
end module profiles
