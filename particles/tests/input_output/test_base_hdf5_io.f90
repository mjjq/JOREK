program test_boris_hdf5_io
use mod_particle_io_boris
implicit none

class(particle_boris), allocatable, dimension(10) :: particles
class(particle_hdf5_io), allocatable :: io
integer :: i

! Prepare HDF5 IO
allocate(particle_hdf5_io_boris::io); allocate(particle_boris::io%particle_type)
call io%set_data_type

! Prepare some particles
allocate(particle_boris::particles)
do i=1,size(particles,1)
  particles(i)%x = real((/i,i+1,i+2/),8)
  particles(i)%v = real((/2*i,2*i+1,2*i+2/),8)
end do

! Test write
call io%write('test_base_hdf5_io.h5', particles)

deallocate(particles)
! Test read
call io%read('test_base_hdf5_io.h5', particles)

do i=1,size(particles,1)
  if (particles(i)%x .ne. real((/i,i+1,i+2/),8)) then
    write(*,*) "error in read x"
  end if
  if (particles(i)%v .ne. real((/2*i,2*i+1,2*i+2/),8)) then
    write(*,*) "error in read y"
  end if
end do

end program test_boris_hdf5_io
