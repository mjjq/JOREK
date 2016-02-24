module mod_vtk
contains

!> Write a vtk file containing points, cells and point data (scalars and vectors)
subroutine write_vtk(filename,xyz,ien,cell_type,scalar_names,scalars,vector_names,vectors)
!> Input arguments
character*(*), intent(in) :: filename !< Output file name
real*4,        intent(in) :: xyz(:,:) !< Point positions
integer,       intent(in), optional :: ien(:,:) !< Element list ien(number of basis functions, element index)
integer,       intent(in), optional :: cell_type !< Type of interpolation (vtk param)
character*12,  intent(in), optional :: scalar_names(:), vector_names(:)
real*4,        intent(in), optional :: scalars(:,:), vectors(:,:,:) !< scalars(nnos, num_scalars)

!> Parameters
integer, parameter    :: ivtk = 22 ! an arbitrary unit number for the VTK output file TODO get automatically

!> Internal variables
character             :: buffer*80, lf*1, str1*12, str2*12
integer :: i, j, i_var
integer :: nnos, nel, nnoel

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
nnos = size(xyz,2)
write(str1(1:12),'(i12)') nnos
buffer = 'POINTS '//str1//'  float'//lf      ; write(ivtk) trim(buffer)
write(ivtk) ((real(xyz(i,j),4),i=1,3),j=1,nnos)

! CELLS SECTION
if (present(ien)) then
  nel   = size(ien,2)
  nnoel = size(ien,1)
  write(str1(1:12),'(i12)') nel            ! number of elements (cells)
  write(str2(1:12),'(i12)') nel*(1+nnoel)  ! size of the following element list (nel*(nnoel+1))
  buffer = lf//'CELLS '//str1//' '//str2//lf  ; write(ivtk) trim(buffer)
  write(ivtk) (int(nnoel,4),(int(ien(i,j),4),i=1,nnoel),j=1,nel)

  ! CELL_TYPES SECTION
  if (present(cell_type)) then
    write(str1(1:12),'(i12)') nel   ! number of elements (cells)
    buffer = lf//'CELL_TYPES'//str1//lf         ; write(ivtk) trim(buffer)
    write(ivtk) (int(cell_type,4),i=1,nel)
  endif
endif

! POINT_DATA SECTION
write(str1(1:12),'(i12)') nnos
buffer = lf//'POINT_DATA '//str1            ; write(ivtk) trim(buffer)

if (present(scalars)) then
  do i_var = 1, size(scalars,2)
    if (present(scalar_names)) then
      buffer = lf//'SCALARS '//scalar_names(i_var)//' float'//lf ; write(ivtk) trim(buffer)
    else
      write(str1(1:12),'(i12)') i_var
      buffer = lf//'SCALARS '//str1//' float'//lf ; write(ivtk) trim(buffer)
    endif
    buffer = 'LOOKUP_TABLE default'//lf
    write(ivtk) trim(buffer)
    write(ivtk) (real(scalars(i,i_var),4),i=1,size(scalars,1))
  enddo
endif

if (present(vectors)) then
  do i_var = 1, size(vectors,3)
    if (present(vector_names)) then
      buffer = lf//lf//'VECTORS '//vector_names(i_var)//' float'//lf ; write(ivtk) trim(buffer)
    else
      write(str1(1:12),'(i12)') i_var
      buffer = lf//lf//'VECTORS '//str1//' float'//lf ; write(ivtk) trim(buffer)
    endif
    write(ivtk) ((real(vectors(j,i,i_var),4),i=1,3),j=1,size(scalars,1))
  enddo
endif

close(ivtk)
end subroutine write_vtk
end module mod_vtk
