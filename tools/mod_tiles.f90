! The module mod_tiles implement the tile types required by the
! mod_tile_arrays data structure. A tiled array is basically 
! an array of arrays, called tiles, such that each tile can have 
! an independet size and indexing. The idea is to split an input array 
! in multiple tiles and resize only a small amount of tiles when needed. 
! In addition, it is expected that tiled_array would allow some degree 
! of manual memory management if needed.
! WARNING: the mod_tiled_array module tries to be as less polymorphic 
!          as possible for avoiding incompatibilities with OpenMP and
!          offloading
module mod_tiles
implicit none
private 
public :: tile_int_1d, tile_int_2d
public :: tile_real8_1d, tile_real8_2d

! Types definitions -----------------------------------------------

! definition of tile types integer_1D, integer_2D, 
! double_1D and double_2D tile types
! WARNING: data in the tile must be stored with increasing index
! Attributes are:
!   data_array: (integer/real8)(N_rows/N_rows*N_columns) tile array
!               in which data are stored 
!> define 1D intger tile array
type :: tile_int_1d
  integer,dimension(:),allocatable :: data_array
contains
  !procedure(alloc_int_tile_1d),deferred :: allocate_tile => allocate_int_tile_1d
  procedure,pass(tile) :: allocate_tile   => allocate_int_tile_1d
  procedure,pass(tile) :: deallocate_tile => deallocate_int_tile_1d
  procedure,pass(tile) :: resize_tile     => resize_int_tile_1d 
end type tile_int_1d

!> define 2D integer tile array
type :: tile_int_2d
  integer,dimension(:,:),allocatable :: data_array
contains
  procedure,pass(tile) :: allocate_tile   => allocate_int_tile_2d
  procedure,pass(tile) :: deallocate_tile => deallocate_int_tile_2d
  procedure,pass(tile) :: resize_tile     => resize_int_tile_2d
end type tile_int_2d

!> define 1D double array
type :: tile_real8_1d 
  real*8,dimension(:),allocatable :: data_array
contains
  procedure,pass(tile) :: allocate_tile   => allocate_real8_tile_1d
  procedure,pass(tile) :: deallocate_tile => deallocate_real8_tile_1d
  procedure,pass(tile) :: resize_tile     => resize_real8_tile_1d
end type tile_real8_1d

!> define 2D double array
type :: tile_real8_2d
  real*8,dimension(:,:),allocatable :: data_array
contains
  procedure,pass(tile) :: allocate_tile   => allocate_real8_tile_2d
  procedure,pass(tile) :: deallocate_tile => deallocate_real8_tile_2d
  procedure,pass(tile) :: resize_tile     => resize_real8_tile_2d
end type tile_real8_2d

contains

! Allocate tile ---------------------------------------------------

! The allocate_int_tile_1d procedure allocate an integer 1d tile
! If present, tile array is initilized to dat, to 0 otherwise
! inputs:
!   tile:   (tile_int_1d) the tile to be allocated
!   N_rows: (integer) number of tile rows
!   dat:    (double)(N_rows,N_cols)(optional) data array for initialisation
! outputs:
!   tile:   (tile_int_1d) the allocated tile
subroutine allocate_int_tile_1d(tile,N_rows,dat)
  implicit none
  ! inputs-outputs
  class(tile_int_1d),intent(inout) :: tile
  ! inputs
  integer,intent(in) :: N_rows
  integer,dimension(N_rows),intent(in),optional :: dat

  if(.not.allocated(tile%data_array)) allocate(tile%data_array(N_rows))
  if(present(dat)) then
    tile%data_array = dat
  else
    tile%data_array = 0
  endif
end subroutine allocate_int_tile_1d

! The allocate_int_tile_2d procedure allocate an integer 1d tile
! If present, tile array is initilized to dat, to 0 otherwise
! inputs:
!   tile:   (tile_int_2d) the tile to be allocated
!   N_rows: (integer) number of tile rows
!   N_cols: (integer) number of tile columns
!   dat:    (double)(N_rows,N_cols)(optional) data array for initialisation
! outputs:
!   tile:   (tile_int_2d) the allocated tile
subroutine allocate_int_tile_2d(tile,N_rows,N_cols,dat)
!  implicit none
  ! inputs-outputs
  class(tile_int_2d),intent(inout) :: tile
  ! inputs
  integer,intent(in) :: N_rows,N_cols
  integer,dimension(N_rows,N_cols),intent(in),optional :: dat

  if(.not.allocated(tile%data_array)) allocate(tile%data_array(N_rows,N_cols))
  if(present(dat)) then
    tile%data_array = dat
  else
    tile%data_array = 0
  endif

end subroutine allocate_int_tile_2d

! The allocate_real8_tile_1d procedure allocate an integer 1d tile
! If present, tile array is initilized to dat, to 0 otherwise
! inputs:
!   tile:   (tile_real8_1d) the tile to be allocated
!   N_rows: (integer) number of tile rows
!   dat:    (double)(N_rows,N_cols)(optional) data array for initialisation
! outputs:
!   tile:   (tile_real8_1d) the allocated tile
subroutine allocate_real8_tile_1d(tile,N_rows,dat)
  implicit none
  ! inputs-outputs
  class(tile_real8_1d),intent(inout) :: tile
  ! inputs
  integer,intent(in) :: N_rows
  real*8,dimension(N_rows),intent(in),optional :: dat

  if(.not.allocated(tile%data_array)) allocate(tile%data_array(N_rows))
  if(present(dat)) then
    tile%data_array = dat
  else
    tile%data_array = 0.d0
  endif
end subroutine allocate_real8_tile_1d

! The allocate_real8_tile_2d procedure allocate an integer 1d tile
! If present, tile array is initilized to dat, to 0 otherwise
! inputs:
!   tile:   (tile_real8_2d) the tile to be allocated
!   N_rows: (integer) number of tile rows
!   N_cols: (integer) number of tile columns
!   dat:    (double)(N_rows,N_cols)(optional) data array for initialisation
! outputs:
!   tile:   (tile_real8_2d) the allocated tile
subroutine allocate_real8_tile_2d(tile,N_rows,N_cols,dat)
  implicit none
  ! inputs-outputs
  class(tile_real8_2d),intent(inout) :: tile
  ! inputs
  integer,intent(in) :: N_rows,N_cols
  real*8,dimension(N_rows,N_cols),intent(in),optional :: dat

  if(.not.allocated(tile%data_array)) allocate(tile%data_array(N_rows,N_cols))
  if(present(dat)) then
     tile%data_array = dat
  else
     tile%data_array = 0.d0
  endif

end subroutine allocate_real8_tile_2d

! Deallocate tile -------------------------------------------------

! Deallocate tile_int_1d data array
! inputs:
!   tile:   (tile_int_1d) the allocated tile
! outputs:
!   tile:   (tile_int_1d) the deallocated tile
subroutine deallocate_int_tile_1d(tile)
  implicit none
  ! inputs-outputs
  class(tile_int_1d),intent(inout) :: tile
  
  if(allocated(tile%data_array)) deallocate(tile%data_array)

end subroutine deallocate_int_tile_1d

! Deallocate tile_int_2d data array
! inputs:
!   tile:   (tile_int_2d) the allocated tile
! outputs:
!   tile:   (tile_int_2d) the deallocated tile
subroutine deallocate_int_tile_2d(tile)
  implicit none
  ! inputs-outputs
  class(tile_int_2d),intent(inout) :: tile
  
  if(allocated(tile%data_array)) deallocate(tile%data_array)

end subroutine deallocate_int_tile_2d

! Deallocate tile_int_1d data array
! inputs:
!   tile:   (tile_real8_1d) the allocated tile
! outputs:
!   tile:   (tile_real8_1d) the deallocated tile
subroutine deallocate_real8_tile_1d(tile)
  implicit none
  ! inputs-outputs
  class(tile_real8_1d),intent(inout) :: tile
  
  if(allocated(tile%data_array)) deallocate(tile%data_array)

end subroutine deallocate_real8_tile_1d

! Deallocate tile_int_1d data array
! inputs:
!   tile:   (tile_real8_2d) the allocated tile
! outputs:
!   tile:   (tile_real8_2d) the deallocated tile
subroutine deallocate_real8_tile_2d(tile)
  implicit none
  ! inputs-outputs
  class(tile_real8_2d),intent(inout) :: tile
  
  if(allocated(tile%data_array)) deallocate(tile%data_array)

end subroutine deallocate_real8_tile_2d

! Resize tile -----------------------------------------------------

! resize_int_tile_1d resizes a integer 1d data_array of a tile
! preserving the data
! inputs:
!   tile:       (tile_int_1d) the tile to resize
!   N_rows_new: (integer) desired size tile
!   N_data:     (integer) number of stored data
! outputs:
!  tile: (tile_int_1d) the resized tile
subroutine resize_int_tile_1d(tile,N_rows_new,N_data)
  implicit none
  ! inputs-outputs:
  class(tile_int_1d),intent(inout) :: tile
  ! inputs:
  integer,intent(in) :: N_rows_new,N_data
  ! variables
  integer,dimension(N_data) :: tmp_data !< temporary data array

  ! initialization
  if(N_data.gt.N_rows_new) then
    write(*,*) "Warning tile: number of data larger than number of rows: skip resizing"
    return
  endif

  ! resizing
  tmp_data = tile%data_array(1:N_data)
  call tile%deallocate_tile()
  call tile%allocate_tile(N_rows_new)
  tile%data_array(1:N_data) = tmp_data
  
end subroutine resize_int_tile_1d

! resize_int_tile_2d resizes a integer 2d data_array of a tile
! preserving the data
! inputs:
!   tile:        (tile_int_2d) the tile to resize
!   N_rows_new:  (integer) desired number of rows
!   N_cols_new:  (integer) desired number of columns
!   N_data_rows: (integer) number of stored data rows
!   N_data_cols: (integer) number of stored data columns
! outputs:
!  tile: (tile_int_2d) the resized tile
subroutine resize_int_tile_2d(tile,N_rows_new,N_cols_new,N_data_rows,N_data_cols)
  implicit none
  ! inputs-outputs:
  class(tile_int_2d),intent(inout) :: tile
  ! inputs:
  integer,intent(in) :: N_rows_new,N_cols_new
  integer,intent(in) :: N_data_rows,N_data_cols
  ! variables
  integer,dimension(N_data_rows,N_data_cols) :: tmp_data !< temporary data array

  ! initialization
  if(N_data_rows.gt.N_rows_new) then
    write(*,*) "Warning tile: number of data rows larger than number of rows: skip resizing"
    return
  endif
  if(N_data_cols.gt.N_cols_new) then
    write(*,*) "Warning tile: number of data columns larger than number of rows: skip resizing"
    return
  endif

  ! resizing
  tmp_data = tile%data_array(1:N_data_rows,1:N_data_cols)
  call tile%deallocate_tile()
  call tile%allocate_tile(N_rows_new,N_cols_new)
  tile%data_array(1:N_data_rows,1:N_data_cols) = tmp_data
  
end subroutine resize_int_tile_2d

! resize_real8_1d resizes a double 1d data_array of a tile
! preserving the data
! inputs:
!   tile:       (tile_real8_1d) the tile to resize
!   N_rows_new: (integer) desired size tile
!   N_data:     (integer) number of stored data
! outputs:
!  tile: (tile_real8_1d) the resized tile
subroutine resize_real8_tile_1d(tile,N_rows_new,N_data)
  implicit none
  ! inputs-outputs:
  class(tile_real8_1d),intent(inout) :: tile
  ! inputs:
  integer,intent(in) :: N_rows_new,N_data
  ! variables
  real*8,dimension(N_data) :: tmp_data !< temporary data array

  ! initialization
  if(N_data.gt.N_rows_new) then
    write(*,*) "Warning tile: number of data larger than number of rows: skip resizing"
    return
  endif

  ! resizing
  tmp_data = tile%data_array(1:N_data)
  call tile%deallocate_tile()
  call tile%allocate_tile(N_rows_new)
  tile%data_array(1:N_data) = tmp_data
  
end subroutine resize_real8_tile_1d

! resize_real8_2d resizes a double 2d data_array of a tile
! preserving the data
! inputs:
!   tile:        (tile_real8_2d) the tile to resize
!   N_rows_new:  (integer) desired number of rows
!   N_cols_new:  (integer) desired number of columns
!   N_data_rows: (integer) number of stored data rows
!   N_data_cols: (integer) number of stored data columns
! outputs:
!  tile: (tile_real8_2d) the resized tile
subroutine resize_real8_tile_2d(tile,N_rows_new,N_cols_new,N_data_rows,N_data_cols)
  implicit none
  ! inputs-outputs:
  class(tile_real8_2d),intent(inout) :: tile
  ! inputs:
  integer,intent(in) :: N_rows_new,N_cols_new
  integer,intent(in) :: N_data_rows,N_data_cols
  ! variables
  real*8,dimension(N_data_rows,N_data_cols) :: tmp_data !< temporary data array

  ! initialization
  if(N_data_rows.gt.N_rows_new) then
    write(*,*) "Warning tile: number of data rows larger than number of rows: skip resizing"
    return
  endif
  if(N_data_cols.gt.N_cols_new) then
    write(*,*) "Warning tile: number of data columns larger than number of rows: skip resizing"
    return
  endif

  ! resizing
  tmp_data = tile%data_array(1:N_data_rows,1:N_data_cols)
  call tile%deallocate_tile()
  call tile%allocate_tile(N_rows_new,N_cols_new)
  tile%data_array(1:N_data_rows,1:N_data_cols) = tmp_data
  
end subroutine resize_real8_tile_2d

!------------------------------------------------------------------

end module mod_tiles
