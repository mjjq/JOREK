!> Module for writing diagnostic data to the screen or to files (ascii, vtk, ...).
module mod_diag_output
  
  
  
  
  
  use parameters
  use mod_expression
  use mod_position
  use equil_info
  
  
  
  
  
  implicit none
  
  
  
  
  
  public
  private write_ascii_header, open_ascii_file, close_file, open_vtk_file
  
  
  
  
  
  ! --- Constants
  character(len=15), parameter, private :: THIS_MOD_NAME = 'mod_diag_output'
  integer,           parameter, private :: STDOUT = 6
  integer,           parameter          :: FORM_TABLE = 0
  integer,           parameter          :: FORM_LIST  = 1
  
  
  
  
  
  contains
  
  
  
  
  
  !> Reduce the dimensionality of the result array.
  subroutine reduce_result_to_0d(ierr, result, res0d, i1, i2, i3)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME)//':reduce_result_to_0d'
    
    ! --- Routine parameters.
    integer,              intent(inout) :: ierr
    real*8, allocatable,  intent(in)    :: result(:,:,:,:)
    real*8, allocatable,  intent(inout) :: res0d(:)
    integer,              intent(in)    :: i1, i2, i3
    
    ierr = 0
    
    if ( .not. allocated(result) ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': result array not allocated.'
      ierr = 101
      return
    end if
    
    if ( allocated(res0d) ) deallocate(res0d)
    
    allocate( res0d(size(result,4)) )
    res0d = result(i1,i2,i3,:)
    
  end subroutine reduce_result_to_0d
  
  
  
  
  
  !> Reduce the dimensionality of the result array.
  subroutine reduce_result_to_1d(ierr, result, res1d, i1, i2, i3)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME)//':reduce_result_to_1d'
    
    ! --- Routine parameters.
    integer,              intent(inout) :: ierr
    real*8, allocatable,  intent(in)    :: result(:,:,:,:)
    real*8, allocatable,  intent(inout) :: res1d(:,:)
    integer, optional,    intent(in)    :: i1, i2, i3
    
    ierr = 0
    
    if ( .not. allocated(result) ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': result array not allocated.'
      ierr = 101
      return
    end if
    
    if ( allocated(res1d) ) deallocate(res1d)
    
    if ( (present(i1)) .and. (present(i2)) ) then
      allocate( res1d(size(result,3),size(result,4)) )
      res1d = result(i1,i2,:,:)
    else if ( (present(i1)) .and. (present(i3)) ) then
      allocate( res1d(size(result,2),size(result,4)) )
      res1d = result(i1,:,i3,:)
    else if ( (present(i2)) .and. (present(i3)) ) then
      allocate( res1d(size(result,1),size(result,4)) )
      res1d = result(:,i2,i3,:)
    else
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': Two of the paramters i1, i2, i3 required.'
      ierr = 102
      return
    end if
    
  end subroutine reduce_result_to_1d
  
  
  
  
  
  !> Reduce the dimensionality of the result array.
  subroutine reduce_result_to_2d(ierr, result, res2d, i1, i2, i3)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME)//':reduce_result_to_2d'
    
    ! --- Routine parameters.
    integer,              intent(inout) :: ierr
    real*8, allocatable,  intent(in)    :: result(:,:,:,:)
    real*8, allocatable,  intent(inout) :: res2d(:,:,:)
    integer, optional,    intent(in)    :: i1, i2, i3
    
    ierr = 0
    
    if ( .not. allocated(result) ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': result array not allocated.'
      ierr = 101
      return
    end if
    
    if ( allocated(res2d) ) deallocate(res2d)
    
    if ( present(i1) ) then
      allocate( res2d(size(result,2),size(result,3),size(result,4)) )
      res2d = result(i1,:,:,:)
    else if ( present(i2) ) then
      allocate( res2d(size(result,1),size(result,3),size(result,4)) )
      res2d = result(:,i2,:,:)
    else if ( present(i3) ) then
      allocate( res2d(size(result,1),size(result,2),size(result,4)) )
      res2d = result(:,:,i3,:)
    else
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': parameter i1 or i2 or i3 required.'
      ierr = 102
      return
    end if
    
  end subroutine reduce_result_to_2d
  
  
  
  
  
  !> Write diagnostic output to an ascii file.
  subroutine write_ascii_0d(ierr, eq, expr_list, res0d, format, header, filename, append, comment)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME)//':write_ascii_0d'
    
    ! --- Routine parameters.
    integer,                    intent(inout) :: ierr
    type(t_equil_state),        intent(in)    :: eq
    type(t_expr_list),          intent(in)    :: expr_list
    real*8, allocatable,        intent(in)    :: res0d(:)
    integer,                    intent(in)    :: format
    logical,          optional, intent(in)    :: header
    character(len=*), optional, intent(in)    :: filename
    logical,          optional, intent(in)    :: append
    character(len=*), optional, intent(in)    :: comment
    
    ! --- Local variables.
    integer :: i_file, i
    
    ierr = 0
    
    if ( .not. allocated(res0d) ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': array res0d is not allocated.'
      ierr = 101
      return
    end if
    
    call open_ascii_file(ierr, i_file, filename, append, comment)
    
    if ( format == FORM_TABLE ) then
      call write_ascii_header(i_file, expr_list, header)
      write(i_file,'(9999es23.15)') res0d(:)
    else if ( format == FORM_LIST ) then
      do i = 1, size(res0d)
        write(i_file,'(1x,a,"=",es23.15)') expr_list%expr(i)%name, res0d(i)
      end do
    else
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': parameter format has illegal value.'
      ierr = 103
      return
    end if
    
    call close_file(ierr, i_file)
    
  end subroutine write_ascii_0d
  
  
  
  
  
  !> Write diagnostic output to an ascii file.
  subroutine write_ascii_1d(ierr, eq, expr_list, res1d, format, header, filename, append, comment)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME)//':write_ascii_1d'
    
    ! --- Routine parameters.
    integer,                    intent(inout) :: ierr
    type(t_equil_state),        intent(in)    :: eq
    type(t_expr_list),          intent(in)    :: expr_list
    real*8, allocatable,        intent(in)    :: res1d(:,:)
    integer,                    intent(in)    :: format
    logical,          optional, intent(in)    :: header
    character(len=*), optional, intent(in)    :: filename
    logical,          optional, intent(in)    :: append
    character(len=*), optional, intent(in)    :: comment
    
    ! --- Local variables.
    integer :: i_file, i
    
    ierr = 0
    
    if ( .not. allocated(res1d) ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': array res1d is not allocated.'
      ierr = 101
      return
    end if
    
    call open_ascii_file(ierr, i_file, filename, append, comment)
    
    if ( format == FORM_TABLE ) then
      call write_ascii_header(i_file, expr_list, header)
      do i = 1, size(res1d,1)
        write(i_file,'(9999es23.15)') res1d(i,:)
      end do
    else
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': parameter format has illegal value.'
      ierr = 103
      return
    end if
    
    call close_file(ierr, i_file)
    
  end subroutine write_ascii_1d
  
  
  
  
  
  !> [Private] Auxilliary routine for write_ascii routines.
  subroutine open_ascii_file(ierr, i_file, filename, append, comment)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME)//':open_ascii_file'
    
    ! --- Routine parameters
    integer,                    intent(inout) :: ierr
    integer,                    intent(out)   :: i_file
    character(len=*), optional, intent(in)    :: filename
    logical,          optional, intent(in)    :: append
    character(len=*), optional, intent(in)    :: comment
    
    ! --- Local variables
    character(len=32) :: status, access
    
    ierr = 0
    
    if ( present(filename) ) then
      i_file = 133 !###
      status = 'replace'
      access = 'sequential'
      if ( present(append) .and. (append) ) then
        status = 'old'
        access = 'append'
      end if
      open(i_file, file=trim(filename), form='formatted', status=trim(status), access=trim(access),&
        iostat=ierr)
      if ( ierr /= 0 ) then
        write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//' opening file "'//trim(filename)//'".'
        write(*,*) ierr
        return
      end if
      if ( present(append) .and. (append) ) then
        write(i_file,*)
        write(i_file,*)
      end if
      if ( present(comment) ) write(i_file,'(2a)') '# ', trim(comment)
    else
      i_file = STDOUT
      write(i_file,*)
    end if
    
  end subroutine open_ascii_file
  
  
  
  
  
  !> [Private] Auxilliary routine for write_ascii routines.
  subroutine close_file(ierr, i_file)
    
    ! --- Routine parameters
    integer,                    intent(inout) :: ierr
    integer,                    intent(in)    :: i_file
    
    ierr = 0
    
    if ( i_file /= STDOUT ) then
      close(i_file, iostat=ierr)
    else
      write(i_file,*)
    end if
    
  end subroutine close_file
  
  
  
  
  
  
  !> [Private] Auxilliary routine for write_ascii routines.
  subroutine write_ascii_header(i_file, expr_list, header)
    
    ! --- Routine parameters
    integer,           intent(in) :: i_file
    type(t_expr_list), intent(in) :: expr_list
    logical, optional, intent(in) :: header
    
    ! --- Local variables
    character(len=23) :: s
    integer :: i
    
    if ( present(header) .and. (header) ) then
      write(i_file,'(a)',advance='no') '# '
      do i = 1, expr_list%n_expr
        s = trim(expr_list%expr(i)%name)
        write(i_file,'(a)',advance='no') s
      end do
      write(i_file,'(a)')
    end if
    
  end subroutine write_ascii_header
  
  
  
  
  
  !> Write diagnostic output to a vtk file.
  subroutine write_vtk_2d(ierr, expr_list, res2d, filename, i_coord, close1, close2)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME)//':write_vtk_2d'
    
    ! --- Routine parameters.
    integer,                    intent(inout) :: ierr
    type(t_expr_list),          intent(in)    :: expr_list
    real*8, allocatable,        intent(in)    :: res2d(:,:,:)
    character(len=*),           intent(in)    :: filename
    integer,                    intent(in)    :: i_coord(2) !< Expression numbers containing coords
    logical, optional,          intent(in)    :: close1
    logical, optional,          intent(in)    :: close2
    
    ! --- Local variables.
    integer :: i_file, n(3), i, j, i_var, n_pts, n_cells, n1, n2
    character(len=80) :: buffer
    character(len=12) :: str1, str2
    character(len=1), parameter :: lf = char(10)
    
    ierr = 0
    
    if ( .not. allocated(res2d) ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': array res2d not allocated.'
      ierr = 103
      return
    end if
    n(:) = (/ size(res2d,1), size(res2d,2), size(res2d,3) /)
    
    call open_vtk_file(ierr, i_file, filename)
    
    buffer = '# vtk DataFile Version 3.0'//lf    ; write(i_file) trim(buffer)
    buffer = 'vtk output'//lf                    ; write(i_file) trim(buffer)
    buffer = 'BINARY'//lf                        ; write(i_file) trim(buffer)
    buffer = 'DATASET UNSTRUCTURED_GRID'//lf     ; write(i_file) trim(buffer)
    
    ! POINTS SECTION
    n_pts = n(1)*n(2)
    write(str1,'(i12)') n_pts
    buffer = 'POINTS '//str1//'  float'//lf      ; write(i_file) trim(buffer)
    do j = 1, n(2)
      do i = 1, n(1)
        write(i_file) real(res2d(i,j,i_coord(1)),4), real(res2d(i,j,i_coord(2)),4), real(0.d0,4)
      end do
    end do
    
    ! CELLS SECTION
    n1 = n(1)-1
    n2 = n(2)-1
    if ( present(close1) .and. close1 ) n1 = n(1)
    if ( present(close2) .and. close2 ) n2 = n(2)
    n_cells = n1*n2
    write(str1(1:12),'(i12)') n_cells
    write(str2(1:12),'(i12)') 5*(n_cells)
    buffer = lf//'CELLS '//str1//' '//str2//lf  ; write(i_file) trim(buffer)
    do j = 1, n2
      do i = 1, n1
        write(i_file) int(4,4), int(i-1+(j-1)*n(1),4), int(i-1+(mod(j-0,n(2)-1))*n(1),4),          &
          int(mod(i-0,n(1)-1)+(mod(j-0,n(2)-1))*n(1),4), int(mod(i-0,n(1)-1)+(j-1)*n(1),4)
      end do
    end do
    
    ! CELL_TYPES SECTION
    write(str1(1:12),'(i12)') n_cells
    buffer = lf//'CELL_TYPES'//str1//lf         ; write(i_file) trim(buffer)
    write(i_file) (int(9,4),i=1,n_cells)
    
    ! POINT_DATA SECTION
    write(str1(1:12),'(i12)') n_pts
    buffer = lf//'POINT_DATA '//str1            ; write(i_file) trim(buffer)
    
    do i_var = 1, n(3)
      buffer = lf//'SCALARS '//expr_list%expr(i_var)%name//' float'//lf ; write(i_file) trim(buffer)
      buffer = 'LOOKUP_TABLE default'//lf;                                write(i_file) trim(buffer)
      do j = 1, n(2)
        do i = 1, n(1)
          write(i_file) real(res2d(i,j,i_var),4)
        end do
      end do      
    enddo
    
    call close_file(ierr, i_file)
    
  end subroutine write_vtk_2d
  
  
  
  
  
  !> [Private] Auxilliary routine for write_vtk_2d routine.
  subroutine open_vtk_file(ierr, i_file, filename)
    
    ! --- Routine parameters
    integer,                    intent(inout) :: ierr
    integer,                    intent(out)   :: i_file
    character(len=*),           intent(in)    :: filename
    
    ierr = 0
    i_file = 133 !###
    
#ifdef IBM_MACHINE
    open(i_file, file=trim(filename), form='unformatted', access='stream', status='replace',       &
      iostat=ierr)
#else
    open(i_file, file=trim(filename), form='binary', convert='BIG_ENDIAN', status='replace',       &
      iostat=ierr)
#endif
    
  end subroutine open_vtk_file
  
  
  
  
  
end module mod_diag_output
