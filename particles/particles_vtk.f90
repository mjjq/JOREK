!> This subroutine gathers and exports all particles to a VTK file
!! The VTK file contains
!! - XYZ coordinates
!!   - charge (scalar)
!!   - velocity (vector at t-dt/2 NB!)
subroutine particles_vtk(particle_list,particle_file)
use parameters
use data_structure
use mod_particles
use mpi_mod

implicit none

type (type_particle_list), intent(in) :: particle_list
character*(*),             intent(in) :: particle_file

type (type_particle)      :: particle

integer               :: nnoel, nnos, nnos_local, nnos_max, nnos_zero, nel
integer               :: nsub, inode, ielm, n_scalars, n_vectors, i, j, k, i_var, ierr
integer               :: nsend, nrecv
integer*4             :: my_id, n_cpu
real,allocatable      :: xyz (:,:), scalars(:,:), vectors(:,:,:)
character*12, allocatable :: scalar_names(:), vector_names(:)
integer, parameter    :: ivtk = 22 ! an arbitrary unit number for the VTK output file
character             :: buffer*80, lf*1, str1*12, str2*12
integer               :: status(MPI_STATUS_SIZE)


call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs

if (my_id .eq. 0) then
  write(*,*) '***********************************'
  write(*,*) '*    export particles to VTK      *'
  write(*,*) '***********************************'

  open(unit=ivtk,file=particle_file,access='stream',form='unformatted',convert='BIG_ENDIAN')

  lf = char(10)

  buffer = '# vtk DataFile Version 3.0'//lf    ; write(ivtk) trim(buffer)
  buffer = 'vtk output'//lf                    ; write(ivtk) trim(buffer)
  buffer = 'BINARY'//lf                        ; write(ivtk) trim(buffer)
  buffer = 'DATASET UNSTRUCTURED_GRID'//lf     ; write(ivtk) trim(buffer)
endif


nnos_local = particle_list%n_particles
call MPI_AllReduce(nnos_local,nnos,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(nnos_local,nnos_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
if (my_id .eq.0) write(*,*) ' total number of particles : ',nnos,nnos_max

n_vectors = 1
n_scalars = 1

allocate(scalar_names(n_scalars), vector_names(n_vectors))
allocate(xyz(3,nnos_max),scalars(nnos_max,1:n_scalars),vectors(nnos_max,3,1:n_vectors))

scalar_names(1)='charge      '
vector_names(1)='velocity    '

xyz = 0.

! Process local particles
nnos_zero = nnos_local
do j=1, nnos_local
  particle = particle_list%particle(j)

  xyz(1,j) = particle%x(1) * cos(particle%x(3))
  xyz(2,j) = particle%x(1) * sin(particle%x(3))
  xyz(3,j) = particle%x(2)

  scalars(j,1) = particle%q

  vectors(j,1,1) = particle%v(1) * cos(particle%x(3)) -  particle%v(3) * sin(particle%x(3))
  vectors(j,2,1) = particle%v(1) * sin(particle%x(3)) +  particle%v(3) * cos(particle%x(3))
  vectors(j,3,1) = particle%v(2)
enddo


! Write XYZ coordinates
if (my_id .eq.0) then
  ! POINTS SECTION
  write(str1(1:12),'(i12)') nnos
  buffer = 'POINTS '//str1//'  float'//lf      ; write(ivtk) trim(buffer)
  write(ivtk) ((real(xyz(i,k),4),i=1,3),k=1,nnos_zero)

  ! Gather from other processors
  do j=1,n_cpu-1
    call mpi_recv(nnos_local,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)

    if (nnos_local .gt. 0) then
      nrecv = 3*nnos_local ! 3 doubles (xyz array)
      call mpi_recv(xyz,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
      write(ivtk) ((real(xyz(i,k),4),i=1,3),k=1,nnos_local)
    endif
  enddo
else
  call mpi_send(nnos_local, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
  if (nnos_local .gt. 0) then
    nsend = 3*nnos_local
    call mpi_send(xyz, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
  endif
endif


! POINT_DATA SECTION
if (my_id .eq.0) then
  write(str1(1:12),'(i12)') nnos
  buffer = lf//'POINT_DATA '//str1            ; write(ivtk) trim(buffer)
endif

do i_var =1, n_scalars
  if (my_id .eq. 0) then
    buffer = 'SCALARS '//scalar_names(i_var)//' float'//lf                              ; write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf                                                 ; write(ivtk) trim(buffer)
    write(ivtk) (real(scalars(i,i_var),4),i=1,nnos_zero)
  endif

  if (my_id .eq. 0) then
    do j=1,n_cpu-1
      call mpi_recv(nnos_local,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)

      if (nnos_local .gt. 0) then
        nrecv = nnos_local
        call mpi_recv(scalars(1:nnos_local,i_var),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        write(ivtk) (real(scalars(i,i_var),4),i=1,nnos_local)
      endif
    enddo
  else
    call mpi_send(nnos_local, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    if (nnos_local .gt. 0) then
      nsend = nnos_local
      call mpi_send(scalars(1:nnos_local,i_var), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
    endif
  endif
enddo


! VECTORS
do i_var =1, n_vectors
  if (my_id .eq. 0) then
    buffer = lf//lf//'VECTORS '//vector_names(i_var)//' float'//lf ; write(ivtk) trim(buffer)
    write(ivtk) ((real(vectors(k,i,i_var),4),i=1,3),k=1,nnos_zero)
  endif

  if (my_id .eq. 0) then
    do j=1,n_cpu-1
      call mpi_recv(nnos_local,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      if (nnos_local .gt. 0) then
        nrecv = 3*nnos_local
        call mpi_recv(vectors(1:nnos_local,:,i_var),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        write(ivtk) ((real(vectors(k,i,i_var),4),i=1,3),k=1,nnos_local)
      endif
    enddo
  else
    call mpi_send(nnos_local, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    if (nnos_local .gt. 0) then
      nsend = 3*nnos_local
      call mpi_send(vectors(1:nnos_local,:,i_var), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
    endif
  endif
enddo


if (my_id .eq.0) then
  close(ivtk)

  write(*,*) '* done export particles to VTK    *'
  write(*,*) '***********************************'
endif

end subroutine particles_vtk
