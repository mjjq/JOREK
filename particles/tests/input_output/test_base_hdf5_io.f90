program test_boris_hdf5_io
use mod_particle_io_boris
use mpi
implicit none

class(particle_base), allocatable, dimension(:) :: particles
class(particle_base), allocatable, dimension(:) :: particles_read
class(particle_hdf5_io), allocatable :: io
integer :: i, ierr, provided

call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)

! Prepare HDF5 IO
allocate(particle_hdf5_io_boris::io)
allocate(particle_boris::io%particle_type(2))
call io%set_data_type

! Prepare some particles
allocate(particle_boris::particles(10))
select type (p => particles)
type is (particle_boris)
  do i=1,size(p,1)
    p(i)%x = real((/i,i+1,i+2/),8)
    p(i)%mass = real(i,4)/3
    p(i)%weight = real(i+1,4)/4
    p(i)%st = real((/i,i+1/),8)/(2.d0*size(p,1))
    p(i)%i_elm = 100*i
    p(i)%q = 3*i
    p(i)%label = 4*i
    p(i)%lost = .true.
    p(i)%v = real((/2*i,2*i+1,2*i+2/),8)
  end do
end select

! Test write
write(*,*) "Start writing"
call io%write('test_base_hdf5_io.h5', particles)
write(*,*) "End writing"


! Test read
write(*,*) "Start reading"
call io%read('test_base_hdf5_io.h5', particles_read)
write(*,*) "End reading"

select type (p1 => particles)
type is (particle_boris)
  select type (p2 => particles_read)
  type is (particle_boris)
    if (size(p1,1) .ne. size(p2,1)) write(*,*) "lengths wrong"

    do i=1,size(p1,1)
      if (norm2(p1(i)%x-p2(i)%x) .gt. 1d-30) write(*,*) "error in x"
      if (abs(p1(i)%mass-p2(i)%mass) .gt. 1d-30) write(*,*) "error in mass"
      if (abs(p1(i)%weight-p2(i)%weight) .gt. 1d-30) write(*,*) "error in weight"
      if (norm2(p1(i)%st-p2(i)%st) .gt. 1d-30) write(*,*) "error in st"
      if (p1(i)%i_elm .ne. p2(i)%i_elm) write(*,*) "error in i_elm"
      if (p1(i)%q .ne. p2(i)%q) write(*,*) "error in q"
      if (p1(i)%label .ne. p2(i)%label) write(*,*) "error in label"
      if (p1(i)%lost .ne. p2(i)%lost) write(*,*) "error in lost"
      if (norm2(p1(i)%v-p2(i)%v) .gt. 1d-30) write(*,*) "error in v"
    end do
  end select
end select

call MPI_Finalize(ierr)
end program test_boris_hdf5_io
