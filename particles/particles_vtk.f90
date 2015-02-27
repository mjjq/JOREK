subroutine particles_vtk(particle_list,particle_file)
use parameters
use data_structure
use mod_particles

implicit none

type (type_particle_list) :: particle_list
type (type_particle)      :: particle
character*(*),intent(in)  :: particle_file

integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors, i, j, i_var
real*4,allocatable    :: xyz (:,:), scalars(:,:), vectors(:,:,:)
character*12, allocatable :: scalar_names(:), vector_names(:)
integer, parameter    :: ivtk = 22 ! an arbitrary unit number for the VTK output file
character             :: buffer*80, lf*1, str1*12, str2*12

write(*,*) '***********************************'
write(*,*) '*    export particles to VTK      *'
write(*,*) '***********************************'
write(*,*) ' number of particles : ',particle_list%n_particles

open(unit=ivtk,file=particle_file,form='binary',convert='BIG_ENDIAN')

lf = char(10)

buffer = '# vtk DataFile Version 3.0'//lf    ; write(ivtk) trim(buffer)
buffer = 'vtk output'//lf                    ; write(ivtk) trim(buffer)
buffer = 'BINARY'//lf                        ; write(ivtk) trim(buffer)
buffer = 'DATASET UNSTRUCTURED_GRID'//lf     ; write(ivtk) trim(buffer)


nnos = particle_list%n_particles
n_vectors = 1
n_scalars = 1

allocate(scalar_names(n_scalars), vector_names(n_vectors))
allocate(xyz(3,nnos),scalars(nnos,1:n_scalars),vectors(nnos,3,1:n_vectors))

scalar_names(1)='charge      '
vector_names(1)='velocity    '

xyz = 0.

do j=1, nnos

  particle = particle_list%particle(j)

  xyz(1,j) = particle%x(1) * cos(particle%x(3))
  xyz(3,j) = particle%x(1) * sin(particle%x(3))
  xyz(2,j) = particle%x(2)

  scalars(j,1) = particle%q

  vectors(j,1,1) = particle%v(1) * cos(particle%x(3)) -  particle%v(3) * sin(particle%x(3))
  vectors(j,3,1) = particle%v(1) * sin(particle%x(3)) +  particle%v(3) * cos(particle%x(3))
  vectors(j,2,1) = particle%v(2)

!  write(*,'(A,3e14.6)') 'check velocity : ', vectors(j,:,1)

enddo

! POINTS SECTION
write(str1(1:12),'(i12)') nnos
buffer = 'POINTS '//str1//'  float'//lf      ; write(ivtk) trim(buffer)
write(ivtk) ((real(xyz(i,j),4),i=1,3),j=1,nnos)

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

write(*,*) '* done export particles to VTK    *'
write(*,*) '***********************************'

return
end