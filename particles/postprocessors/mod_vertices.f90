!> the module mod_vertices contains all datatypes
!> and procedures common to light and gather points
module mod_vertices
use mod_tiles, only: tile_real8_2d
implicit none

private
public :: vertices

!> Variable and type definitions ----------------------
!> vertices: abstract class containing the basic types
!> and procedures defining light or gather points
type,abstract :: vertices
  integer :: n_vertices    !< total number of vertices
  integer :: n_data_vertex !< total number of data per vertex
  integer :: n_tiles       !< total number of tiles
  !> position of each vertex in cartesian coordinates (xyz)
  type(tile_real8_2d) :: x
  !> the number of tiles .eq. number of spectra
  !> the number of rows of one tile is n_(spectral)_lines*n_vertices
  !> the number of columns are 1) values and 2) the pdf
  !> if the vertex is not samplable the pdf is se to to 1
  type(tile_real8_2d),dimension(:),allocatable :: intensity_pdf
  type(tile_real8_2d),dimension(:),allocatable :: strength_pdf

  contains

end type vertices

contains

!>-----------------------------------------------------
!>-----------------------------------------------------
!>-----------------------------------------------------
end module mod_vertices
