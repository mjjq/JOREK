! module dedicated to memory tracing and
! debug tracing
module tr_module

  implicit none

  interface tr_register_mem
     module procedure tr_register_mem_int4, tr_register_mem_int8
  end interface

  interface tr_unregister_mem
     module procedure tr_unregister_mem_int4, tr_unregister_mem_int8
  end interface

  !*** surdefinition for allocation ***
  interface tr_allocate
     module procedure tr_allocate1d_i, tr_allocate1d_d, &
          tr_allocate1d_c, tr_allocate2d_i, &
          tr_allocate2d_d, tr_allocate2d_c, &
          tr_allocate3d_d, tr_allocate3d_c, &
          tr_allocate3d_i
  end interface

  !*** surdefinition for deallocation ***
  interface tr_deallocate
     module procedure tr_deallocate1d_i, &
          tr_deallocate1d_d, tr_deallocate1d_c, &
          tr_deallocate2d_i, tr_deallocate2d_d, &
          tr_deallocate2d_c, tr_deallocate3d_d, &
          tr_deallocate3d_i, tr_deallocate3d_c
  end interface

  !*** surdefinition for allocation ***
  interface tr_allocatep
     module procedure tr_allocatep1d_i, tr_allocatep1d_d, & 
          tr_allocatep2d_d, tr_allocatep2d_i, &
          tr_allocatep3d_d, tr_allocatep3d_i, tr_allocatep4d_d
  end interface

  !*** surdefinition for deallocation ***
  interface tr_deallocatep
     module procedure tr_deallocatep1d_i, &
          tr_deallocatep1d_d, tr_deallocatep2d_i, &
          tr_deallocatep3d_d, tr_deallocatep2d_d, &
          tr_deallocatep3d_i, tr_deallocatep4d_d
  end interface

  !*** surdefinition for deallocation ***
  interface tr_debug_write
     module procedure tr_debug_writes, tr_debug_writei, tr_debug_writef
  end interface

  real*8 :: precond_mem = 0
  ! size of real numbers used in jorek
  integer, private, parameter :: RKIND = 8
  ! processor identity
  integer, public :: gmy_id
  ! nb of processors
  integer, private :: nbprocs
  ! used for memory size calculation
  integer*8, private :: max_allocate, nb_allocate

  character(LEN=13), private :: &
       trace_file = "trace    .out"
  integer, private :: uout_mem = 30

  real(RKIND) :: myreal
  integer     :: myint
  complex     :: mycomp
  !******************************
contains
  !******************************

  !---------------------------------------- 
  ! Init target file in each processor
  !----------------------------------------
  subroutine tr_meminit(pmy_id, pnbprocs)
    integer, intent(in) :: pmy_id, pnbprocs
    logical fexist
    gmy_id = pmy_id
    nbprocs = pnbprocs
    write(trace_file(6:9),'(I4.4)') gmy_id
    inquire(file=trace_file,exist=fexist)
    if (fexist) then
       open(uout_mem, file = trace_file, status = 'OLD', position = 'APPEND', &
            form = 'FORMATTED')
    else
       open(uout_mem, file = trace_file, status = 'NEW', &
            form = 'FORMATTED')
    end if
    write(uout_mem,*) '### meminit ### '
    call flush(uout_mem)
    close(uout_mem)
  end subroutine tr_meminit

  !---------------------------------------- 
  ! Reset file
  !----------------------------------------
  subroutine tr_resetfile()
    open(uout_mem, file = trace_file, status = 'REPLACE', &
         form = 'FORMATTED')
    write(uout_mem,*) '### meminit ### '
    call flush(uout_mem)
    close(uout_mem)
  end subroutine tr_resetfile


  !---------------------------------------- 
  ! Write special string in file trace_file
  !----------------------------------------
  subroutine tr_write(string)
    character*(*) string
    open(uout_mem, file = trace_file, status = 'OLD', &
         position = 'APPEND', form = 'FORMATTED')
    write(uout_mem,'(A)') string
    call flush(uout_mem)
    close(uout_mem)
  end subroutine tr_write

  !---------------------------------------- 
  ! Write debug remark in file trace_file
  !----------------------------------------
  subroutine tr_debug_writes(string)
    character*(*)           :: string
    character(len=1024)     :: bufstring
    write(bufstring,'(A)')string
    call tr_write("### "//trim(adjustl(bufstring))//" ###")
  end subroutine tr_debug_writes


  !---------------------------------------- 
  ! Write debug remark in file trace_file
  !----------------------------------------
  subroutine tr_debug_writei(string, int_var)
    character*(*)           :: string
    integer                 :: int_var
    character(len=1024)     :: bufstring
    write(bufstring,'(A,I20)')string,int_var
    call tr_write("### "//trim(adjustl(bufstring))//" ###")
  end subroutine tr_debug_writei


  !---------------------------------------- 
  ! Write debug remark in file trace_file
  !----------------------------------------
  subroutine tr_debug_writef(string, float_var)
    character*(*)           :: string
    real*8                  :: float_var
    character(len=1024)     :: bufstring
    write(bufstring,'(A,E20.11)')string,float_var
    call tr_write("### "//trim(adjustl(bufstring))//" ###")
  end subroutine tr_debug_writef


  !-------------------------------------------
  ! Write memory in the file trace_file (allocate)
  !-------------------------------------------
  subroutine tr_memwriteadd(size_array,type_name,var_name)
    integer*8    , intent(in)  :: size_array
    character*(*), intent(in)  :: type_name
    character*(*), intent(in)  :: var_name

#ifdef MEMTRACE
    open(uout_mem, file = trace_file, status = 'OLD', &
         position = 'APPEND', form = 'FORMATTED')
    write(uout_mem,'(A10,I15,A10,A15,A20,5X,I20)') &
         'Add', &
         size_array,' Bytes ',type_name,var_name, nb_allocate
    close(uout_mem)    
#endif
  end subroutine tr_memwriteadd

  !-------------------------------------------
  ! Write memory in the file trace_file (deallocate)
  !-------------------------------------------
  subroutine tr_memwritedel(var_name)
    character*(*), intent(in) :: var_name

#ifdef MEMTRACE
    open(uout_mem, file = trace_file, status = 'OLD', &
         position = 'APPEND', form = 'FORMATTED')
    write(uout_mem,'(50X,A20,A5,I20)') var_name, ' Supp', nb_allocate
    close(uout_mem)    
#endif
  end subroutine tr_memwritedel


  !-------------------------------------------
  ! Get memory usage in looking at /proc/self/status
  ! Return an integer representing the consumption in KB
  !-------------------------------------------
  integer*8 function get_memory_inkb(string_pattern)
    character(len=*), intent(in) :: string_pattern
    character(len=1024) :: buffer, string, adj
    integer :: pos, ierr
    integer, parameter :: fh = 15
    integer :: ios 
    integer :: line 
    ios = 0
    line = 0
    get_memory_inkb = 0
    open(fh, file='/proc/self/status', iostat=ierr)
    if (ierr .ne. 0) get_memory_inkb = -1
    ! ios is negative if an end of record condition is encountered or if                                                                                                      
    ! an endfile condition was detected.  It is positive if an error was                                                                                                      
    ! detected.  ios is zero otherwise.                                                                                                                                       
    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1
          ! Find the first instance of whitespace.                                                                                                                            
          ! Split first field and end of line                                                                                                                                 
          pos = scan(buffer,':')
          string = buffer(1:pos-1)
          if (trim(adjustl(string)) .eq. trim(adjustl(string_pattern))) then
             buffer = buffer(pos+1:)
             pos = scan(buffer,'0123456789')
             adj = buffer(pos:)
             pos = scan(adj,'k')
             string = adj(1:pos-2)
             read(string,'(I20)')get_memory_inkb
          end if
       end if
    end do
    close(fh)
    return 
  end function get_memory_inkb

  subroutine tr_set_precondmem(memused)
    real*8 memused
    precond_mem = memused
  end subroutine tr_set_precondmem


  subroutine tr_register_mem_int4(mem_in_bytes,var_name)
    integer*4          , intent(in) :: mem_in_bytes
    character*(*)      , intent(in) :: var_name

    call tr_memwriteadd(int8(mem_in_bytes),'unknown type',var_name)
    nb_allocate = nb_allocate + mem_in_bytes
  end subroutine tr_register_mem_int4

  subroutine tr_register_mem_int8(mem_in_bytes,var_name)
    integer*8          , intent(in) :: mem_in_bytes
    character*(*)      , intent(in) :: var_name

    call tr_memwriteadd(mem_in_bytes,'unknown type',var_name)
    nb_allocate = nb_allocate + mem_in_bytes
  end subroutine tr_register_mem_int8

  subroutine tr_unregister_mem_int4(mem_in_bytes,var_name)
    integer*4          , intent(in) :: mem_in_bytes
    character*(*)      , intent(in) :: var_name

    nb_allocate = nb_allocate - mem_in_bytes
    call tr_memwritedel(var_name)
  end subroutine tr_unregister_mem_int4

  subroutine tr_unregister_mem_int8(mem_in_bytes,var_name)
    integer*8          , intent(in) :: mem_in_bytes
    character*(*)      , intent(in) :: var_name

    nb_allocate = nb_allocate - mem_in_bytes
    call tr_memwritedel(var_name)
  end subroutine tr_unregister_mem_int8

  !---------------------------------------- 
  ! memory allocation for a 1D array
  !----------------------------------------
  subroutine tr_allocatep1d_i(array1d,begin_dim1,end_dim1,var_name)
    integer, dimension(:)    , pointer    :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name

    integer   :: err
    integer   :: i1
    integer*8 :: size_array 

    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array1D',var_name)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = 0
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocatep1d_i

  subroutine tr_allocatep1d_d(array1d,begin_dim1,end_dim1,var_name)
    real(RKIND), dimension(:), pointer    :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name

    integer   :: err
    integer   :: i1
    integer*8 :: size_array 

    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myreal)
    call tr_memwriteadd(size_array,'double array1D',var_name)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = 0._RKIND
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if

    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocatep1d_d

  !---------------------------------------- 
  ! memory allocation for a 2D array
  !----------------------------------------
  subroutine tr_allocatep2d_i(array2d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,var_name)
    integer    , dimension(:,:), pointer    :: array2d
    integer                    , intent(in) :: begin_dim1
    integer                    , intent(in) :: end_dim1
    integer                    , intent(in) :: begin_dim2
    integer                    , intent(in) :: end_dim2
    character*(*)  , intent(in)             :: var_name

    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 

    allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
         stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array2D',var_name)
    if (err.eq.0) then
       do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
             array2d(i1,i2) = 0
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocatep2d_i

  subroutine tr_allocatep2d_d(array2d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,var_name)
    real(RKIND)  , dimension(:,:), pointer    :: array2d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1  
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    character*(*)    , intent(in)             :: var_name

    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 

    allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
         stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1) * sizeof(myreal)
    call tr_memwriteadd(size_array,'double array2D',var_name)
    if (err.eq.0) then
       do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
             array2d(i1,i2) = 0._RKIND
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocatep2d_d

  !---------------------------------------- 
  ! memory allocation for a 3D array
  !----------------------------------------
  subroutine tr_allocatep3d_d(array3d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,var_name)
    real(RKIND), dimension(:,:,:), pointer    :: array3d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3
    character*(*)    , intent(in)             :: var_name

    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 

    allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)* sizeof(myreal)
    call tr_memwriteadd(size_array,'double array3D',var_name)
    if (err.eq.0) then
       do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
             do i1 = begin_dim1,end_dim1
                array3d(i1,i2,i3) = 0._RKIND
             end do
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocatep3d_d

  !---------------------------------------- 
  ! memory allocation for a 4D array
  !----------------------------------------
  subroutine tr_allocatep4d_d(array4d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,begin_dim4,end_dim4,var_name)
    real(RKIND), dimension(:,:,:,:), pointer    :: array4d
    integer                        , intent(in) :: begin_dim1
    integer                        , intent(in) :: end_dim1
    integer                        , intent(in) :: begin_dim2
    integer                        , intent(in) :: end_dim2
    integer                        , intent(in) :: begin_dim3
    integer                        , intent(in) :: end_dim3
    integer                        , intent(in) :: begin_dim4
    integer                        , intent(in) :: end_dim4
    character*(*)                  , intent(in) :: var_name

    integer   :: err
    integer   :: i1, i2, i3, i4
    integer*8 :: size_array 

    allocate(array4d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3,begin_dim4:end_dim4),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)*(end_dim4-begin_dim4+1)* sizeof(myreal)
    call tr_memwriteadd(size_array,'double array4D',var_name)
    if (err.eq.0) then
       do i4 = begin_dim4,end_dim4
          do i3 = begin_dim3,end_dim3
             do i2 = begin_dim2,end_dim2
                do i1 = begin_dim1,end_dim1
                   array4d(i1,i2,i3,i4) = 0._RKIND
                end do
             end do
          end do
       end do
    else
       print *,'problem in allocating ',var_name
       print *,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocatep4d_d

  subroutine tr_allocatep3d_i(array3d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,var_name)
    integer, dimension(:,:,:)    , pointer    :: array3d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3
    character*(*)    , intent(in)             :: var_name

    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 

    allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)* sizeof(myint)
    call tr_memwriteadd(size_array,'double array3D',var_name)
    if (err.eq.0) then
       do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
             do i1 = begin_dim1,end_dim1
                array3d(i1,i2,i3) = 0._RKIND
             end do
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocatep3d_i

  !---------------------------------------- 
  ! memory allocation for a 1D array
  !----------------------------------------
  subroutine tr_allocate1d_i(array1d,begin_dim1,end_dim1,var_name)
    integer, dimension(:)    , allocatable :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name

    integer   :: err
    integer   :: i1
    integer*8 :: size_array 
    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array1D',var_name)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = 0
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocate1d_i

  subroutine tr_allocate1d_d(array1d,begin_dim1,end_dim1,var_name)
    real(RKIND), dimension(:), allocatable :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name

    integer   :: err
    integer   :: i1
    integer*8 :: size_array 

    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myreal)
    call tr_memwriteadd(size_array,'double array1D',var_name)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = 0._RKIND
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if

    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocate1d_d

  subroutine tr_allocate1d_c(array1d,begin_dim1,end_dim1,var_name)
    complex      , dimension(:), allocatable :: array1d
    integer                     , intent(in) :: begin_dim1
    integer                     , intent(in) :: end_dim1
    character*(*)   , intent(in)             :: var_name

    integer   :: err
    integer   :: i1
    integer*8 :: size_array 

    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(mycomp)
    call tr_memwriteadd(size_array,'complex array1D',var_name)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = cmplx(0,0)
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocate1d_c


  !---------------------------------------- 
  ! memory allocation for a 2D array
  !----------------------------------------
  subroutine tr_allocate2d_i(array2d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,var_name)
    integer    , dimension(:,:), allocatable :: array2d
    integer                    , intent(in) :: begin_dim1
    integer                    , intent(in) :: end_dim1
    integer                    , intent(in) :: begin_dim2
    integer                    , intent(in) :: end_dim2
    character*(*)  , intent(in)             :: var_name

    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 

    allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
         stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array2D',var_name)
    if (err.eq.0) then
       do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
             array2d(i1,i2) = 0
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocate2d_i

  subroutine tr_allocate2d_d(array2d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,var_name)
    real(RKIND)  , dimension(:,:), allocatable :: array2d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1  
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    character*(*)    , intent(in)             :: var_name

    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 

    allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
         stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1) * sizeof(myreal)
    call tr_memwriteadd(size_array,'double array2D',var_name)
    if (err.eq.0) then
       do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
             array2d(i1,i2) = 0._RKIND
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocate2d_d

  subroutine tr_allocate2d_c(array2d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,var_name)
    complex      , dimension(:,:), allocatable :: array2d
    integer                       , intent(in) :: begin_dim1
    integer                       , intent(in) :: end_dim1
    integer                       , intent(in) :: begin_dim2
    integer                       , intent(in) :: end_dim2
    character*(*)     , intent(in)             :: var_name

    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 

    allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
         stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1) * sizeof(mycomp)
    call tr_memwriteadd(size_array,'complex array2D',var_name)
    if (err.eq.0) then
       do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
             array2d(i1,i2) = cmplx(0,0)
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)    
  end subroutine tr_allocate2d_c


  !---------------------------------------- 
  ! memory allocation for a 3D array
  !----------------------------------------
  subroutine tr_allocate3d_i(array3d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,var_name)
    integer    , dimension(:,:,:), allocatable :: array3d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3
    character*(*)    , intent(in)             :: var_name

    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 

    allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)* SIZEOF(myint)
    call tr_memwriteadd(size_array,'integer array3D',var_name)
    if (err.eq.0) then
       do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
             do i1 = begin_dim1,end_dim1
                array3d(i1,i2,i3) = 0
             end do
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocate3d_i

  !---------------------------------------- 
  ! memory allocation for a 3D array
  !----------------------------------------
  subroutine tr_allocate3d_c(array3d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,var_name)
    complex   , dimension(:,:,:), allocatable :: array3d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3
    character*(*)    , intent(in)             :: var_name

    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 

    allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)* SIZEOF(mycomp)
    call tr_memwriteadd(size_array,'integer array3D',var_name)
    if (err.eq.0) then
       do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
             do i1 = begin_dim1,end_dim1
                array3d(i1,i2,i3) = cmplx(0,0)
             end do
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocate3d_c

  !---------------------------------------- 
  ! memory allocation for a 3D array
  !----------------------------------------
  subroutine tr_allocate3d_d(array3d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,var_name)
    real(RKIND), dimension(:,:,:), allocatable :: array3d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3
    character*(*)    , intent(in)             :: var_name

    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 

    allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)* sizeof(myreal)
    call tr_memwriteadd(size_array,'double array3D',var_name)
    if (err.eq.0) then
       do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
             do i1 = begin_dim1,end_dim1
                array3d(i1,i2,i3) = 0._RKIND
             end do
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate  = nb_allocate + size_array
    max_allocate = max(max_allocate,nb_allocate)
  end subroutine tr_allocate3d_d

  !---------------------------------------- 
  ! memory deallocation of array 1D
  !----------------------------------------
  subroutine tr_deallocatep1d_i(array1d,var_name)
    integer, dimension(:), pointer     :: array1d
    character*(*)        , intent(in)  :: var_name

    if (associated(array1d)) then
       nb_allocate = nb_allocate - sizeof(array1d)
       call tr_memwritedel(var_name)
       deallocate(array1d)
       array1d => NULL()
    end if
  end subroutine tr_deallocatep1d_i

  subroutine tr_deallocatep1d_d(array1d,var_name)
    real(RKIND), dimension(:), pointer :: array1d
    character*(*)        , intent(in)  :: var_name

    if (associated(array1d)) then
       nb_allocate = nb_allocate - sizeof(array1d)
       call tr_memwritedel(var_name)
       deallocate(array1d)
       array1d => NULL()
    end if
  end subroutine tr_deallocatep1d_d

  subroutine tr_deallocatep2d_i(array2d,var_name)
    integer, dimension(:,:) , pointer  :: array2d
    character*(*)           , intent(in)  :: var_name

    if (associated(array2d)) then
       nb_allocate = nb_allocate - sizeof(array2d)
       call tr_memwritedel(var_name)
       deallocate(array2d)
       array2d => null()
    end if
  end subroutine tr_deallocatep2d_i

  subroutine tr_deallocatep2d_d(array2d,var_name)
    real(RKIND), dimension(:,:), pointer :: array2d
    character*(*)              , intent(in) :: var_name

    if (associated(array2d)) then
       nb_allocate = nb_allocate - sizeof(array2d)
       call tr_memwritedel(var_name)
       deallocate(array2d)
       array2d => null()
    end if
  end subroutine tr_deallocatep2d_d

  subroutine tr_deallocatep3d_d(array3d,var_name)
    real(RKIND), dimension(:,:,:) , pointer :: array3d
    character*(*)                 , intent(in) :: var_name

    if (associated(array3d)) then
       nb_allocate = nb_allocate - sizeof(array3d)
       call tr_memwritedel(var_name)
       deallocate(array3d)
       array3d => null()
    end if
  end subroutine tr_deallocatep3d_d

  subroutine tr_deallocatep4d_d(array4d,var_name)
    real(RKIND), dimension(:,:,:,:) , pointer :: array4d
    character*(*)                   , intent(in) :: var_name

    if (associated(array4d)) then
       nb_allocate = nb_allocate - sizeof(array4d)
       call tr_memwritedel(var_name)
       deallocate(array4d)
       array4d => null()
    end if
  end subroutine tr_deallocatep4d_d


  subroutine tr_deallocatep3d_i(array3d,var_name)
    integer, dimension(:,:,:) , pointer :: array3d
    character*(*)                 , intent(in) :: var_name

    if (associated(array3d)) then
       nb_allocate = nb_allocate - sizeof(array3d)
       call tr_memwritedel(var_name)
       deallocate(array3d)
       array3d => null()
    end if
  end subroutine tr_deallocatep3d_i

  !---------------------------------------- 
  ! memory deallocation of array 1D
  !----------------------------------------
  subroutine tr_deallocate1d_i(array1d,var_name)
    integer, dimension(:), allocatable  :: array1d
    character*(*)        , intent(in)  :: var_name

    if (allocated(array1d)) then
       nb_allocate = nb_allocate - sizeof(array1d)
       call tr_memwritedel(var_name)
       deallocate(array1d)
    end if
  end subroutine tr_deallocate1d_i

  subroutine tr_deallocate1d_d(array1d,var_name)
    real(RKIND), dimension(:), allocatable :: array1d
    character*(*)        , intent(in)  :: var_name

    if (allocated(array1d)) then
       nb_allocate = nb_allocate - sizeof(array1d)
       call tr_memwritedel(var_name)
       deallocate(array1d)
    end if
  end subroutine tr_deallocate1d_d

  subroutine tr_deallocate1d_c(array1d,var_name)
    complex      , dimension(:),  allocatable :: array1d
    character*(*)               , intent(in)  :: var_name

    if (allocated(array1d)) then
       nb_allocate = nb_allocate - sizeof(array1d)
       call tr_memwritedel(var_name)
       deallocate(array1d)
    end if
  end subroutine tr_deallocate1d_c

  !---------------------------------------- 
  ! memory deallocation of array 2D
  !----------------------------------------
  subroutine tr_deallocate2d_i(array2d,var_name)
    integer, dimension(:,:) , allocatable  :: array2d
    character*(*)           , intent(in)  :: var_name

    if (allocated(array2d)) then
       nb_allocate = nb_allocate - sizeof(array2d)
       call tr_memwritedel(var_name)
       deallocate(array2d)
    end if
  end subroutine tr_deallocate2d_i

  subroutine tr_deallocate2d_d(array2d,var_name)
    real(RKIND), dimension(:,:), allocatable :: array2d
    character*(*)              , intent(in) :: var_name

    if (allocated(array2d)) then
       nb_allocate = nb_allocate - sizeof(array2d)
       call tr_memwritedel(var_name)
       deallocate(array2d)
    end if
  end subroutine tr_deallocate2d_d

  subroutine tr_deallocate2d_c(array2d,var_name)
    complex      , dimension(:,:), allocatable :: array2d
    character*(*)                 , intent(in) :: var_name

    if (allocated(array2d)) then
       nb_allocate = nb_allocate - sizeof(array2d)
       call tr_memwritedel(var_name)
       deallocate(array2d)
    end if
  end subroutine tr_deallocate2d_c

  !---------------------------------------- 
  ! memory deallocation of array 3D
  !----------------------------------------
  subroutine tr_deallocate3d_d(array3d,var_name)
    real(RKIND), dimension(:,:,:) , allocatable :: array3d
    character*(*)                 , intent(in) :: var_name

    if (allocated(array3d)) then
       nb_allocate = nb_allocate - sizeof(array3d)
       call tr_memwritedel(var_name)
       deallocate(array3d)
    end if
  end subroutine tr_deallocate3d_d

  subroutine tr_deallocate3d_c(array3d,var_name)
    complex      , dimension(:,:,:), allocatable :: array3d
    character*(*)                   , intent(in) :: var_name

    if (allocated(array3d)) then
       nb_allocate = nb_allocate - sizeof(array3d)
       call tr_memwritedel(var_name)
       deallocate(array3d)
    end if
  end subroutine tr_deallocate3d_c

  subroutine tr_deallocate3d_i(array3d,var_name)
    integer      , dimension(:,:,:), allocatable :: array3d
    character*(*)                   , intent(in) :: var_name

    if (allocated(array3d)) then
       nb_allocate = nb_allocate - sizeof(array3d)
       call tr_memwritedel(var_name)
       deallocate(array3d)
    end if
  end subroutine tr_deallocate3d_i


  !***********************************************
  !  function for program analysis
  !***********************************************
  subroutine tr_print_memsize(label)
    implicit none
    character(*), intent(in) :: label
    integer*8, parameter :: GBconst = 1024_8*1024_8*1024_8
    integer*8, parameter :: MBconst = 1024_8*1024_8
    integer*8, parameter :: KBconst = 1024_8
    integer :: uout
    integer*8 :: scount, dcount, rcount, lcount
    
    rcount = KBconst * get_memory_inkb("VmRSS")
    open(uout_mem, file = trace_file, status = 'OLD', &
         position = 'APPEND', form = 'FORMATTED')
    if (nb_allocate.gt.GBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)'), label,&
            'memsize allocated within Jorek (tr_module) = ', &
            nb_allocate/dfloat(GBconst), ' GBytes'
    else 
       write(uout_mem,'(A20,A50,1f10.3,A)'), label, &
            'memsize allocated within Jorek (tr_module) = ', &
            nb_allocate/dfloat(MBconst), ' MBytes'
    end if
    lcount = rcount - nb_allocate
    if (lcount.gt.GBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)'), label, &
            'memsize occupied by libraries/others = ', &
            lcount/dfloat(GBconst), ' GBytes'
    else if (lcount.gt.MBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)'), label, &
            'memsize occupied by libraries/others = ', &
            lcount/dfloat(MBconst), ' MBytes'
    else
       write(uout_mem,'(A20,A50,1f10.3,A)'), label, &
            'memsize occupied by libraries/others = ', &
            lcount/dfloat(KBconst), ' KBytes'
    end if
    if (rcount.gt.GBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)'), label, &
            'memsize total   (RSS) = ', &
            rcount/dfloat(GBconst), ' GBytes'
    else if (rcount.gt.MBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)'), label, &
            'memsize total   (RSS) = ', &
            rcount/dfloat(MBconst), ' MBytes'
    else
       write(uout_mem,'(A20,A50,1f10.3,A)'), label, &
            'memsize total   (RSS) = ', &
            rcount/dfloat(KBconst), ' KBytes'
    end if
    close(uout_mem)
  end subroutine tr_print_memsize
end module tr_module
