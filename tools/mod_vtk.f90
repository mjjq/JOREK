module mod_vtk
contains

!> Write a vtk file containing points, cells and point data (scalars and vectors)
subroutine write_vtk(filename,nnos,xyz,nel,nnoel,ien,cell_type,n_scalars,scalar_names,scalars,n_vectors,vector_names,vectors)
!> Input arguments
character*(*), intent(in) :: filename !< Output file name
integer,       intent(in) :: nnos !< Number of points
real*4,        intent(in) :: xyz(:,:) !< Point positions
integer,       intent(in) :: nel, nnoel !< Number of elements and size of element list
integer,       intent(in) :: ien(:,:) !< Element list
integer,       intent(in) :: cell_type !< Type of interpolation (vtk param)
integer,       intent(in) :: n_scalars, n_vectors !< Number of scalars/vectors to write
character*12,  intent(in) :: scalar_names(:), vector_names(:)
real*4,        intent(in) :: scalars(:,:), vectors(:,:,:)

!> Parameters
integer, parameter    :: ivtk = 22 ! an arbitrary unit number for the VTK output file

!> Internal variables
character             :: buffer*80, lf*1, str1*12, str2*12
integer :: i, j, i_var

lf = char(10) ! line feed character

#ifdef IBM_MACHINE
open(unit=ivtk,file=filename,form='unformatted',access='stream')
#else
open(unit=ivtk,file=filename,form='unformatted',access='stream',convert='BIG_ENDIAN')
#endif

buffer = '# vtk DataFile Version 3.0'//lf    ; write(ivtk) trim(buffer)
buffer = 'vtk output'//lf                    ; write(ivtk) trim(buffer)
buffer = 'BINARY'//lf                        ; write(ivtk) trim(buffer)
buffer = 'DATASET UNSTRUCTURED_GRID'//lf     ; write(ivtk) trim(buffer)

! POINTS SECTION
write(str1(1:12),'(i12)') nnos
buffer = 'POINTS '//str1//'  float'//lf      ; write(ivtk) trim(buffer)
write(ivtk) ((real(xyz(i,j),4),i=1,3),j=1,nnos)

! CELLS SECTION
write(str1(1:12),'(i12)') nel            ! number of elements (cells)
write(str2(1:12),'(i12)') nel*(1+nnoel)  ! size of the following element list (nel*(nnoel+1))
buffer = lf//'CELLS '//str1//' '//str2//lf  ; write(ivtk) trim(buffer)
write(ivtk) (int(nnoel,4),(int(ien(i,j),4),i=1,nnoel),j=1,nel)

! CELL_TYPES SECTION
write(str1(1:12),'(i12)') nel   ! number of elements (cells)
buffer = lf//'CELL_TYPES'//str1//lf         ; write(ivtk) trim(buffer)
write(ivtk) (int(cell_type,4),i=1,nel)

! POINT_DATA SECTION
write(str1(1:12),'(i12)') nnos
buffer = lf//'POINT_DATA '//str1            ; write(ivtk) trim(buffer)

do i_var =1, n_scalars
  buffer = lf//'SCALARS '//scalar_names(i_var)//' float'//lf ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf
  write(ivtk) trim(buffer)
  write(ivtk) (real(scalars(i,i_var),4),i=1,nnos)
enddo

do i_var =1, n_vectors
  buffer = lf//lf//'VECTORS '//vector_names(i_var)//' float'//lf ; write(ivtk) trim(buffer)
  write(ivtk) ((real(vectors(j,i,i_var),4),i=1,3),j=1,nnos)
enddo

close(ivtk)
end subroutine write_vtk
end module mod_vtk
