module base_hdf5_io_test
use mod_particle_io_boris
use hdf5
implicit none

contains

!> Helper function for testing
subroutine write_read_test(particles, filename)
  class(particle_base), allocatable, dimension(:), intent(in) :: particles
  character(len=*) :: filename

  class(particle_hdf5_io), allocatable :: io
  class(particle_base), allocatable, dimension(:) :: particles_read
  integer :: i

  allocate(particle_hdf5_io_boris::io)
  ! Prepare HDF5 IO
  allocate(particle_boris::io%particle_type(2))
  call io%set_data_type

  ! Test write
  call io%read(filename, particles_read)

  ! Test read
  call io%read(filename, particles_read)

  select type (p1 => particles)
  type is (particle_boris)
    select type (p2 => particles_read)
    type is (particle_boris)
      if (size(p1,1) .ne. size(p2,1)) write(*,*) "lengths wrong", size(p2,1)

      do i=1,size(p2,1)
        if (norm2(p1(i)%x-p2(i)%x) .gt. 1d-30) write(*,*) "error in x"
        if (abs(p1(i)%mass-p2(i)%mass) .gt. 1d-30) write(*,*) "error in mass"
        if (abs(p1(i)%weight-p2(i)%weight) .gt. 1d-30) write(*,*) "error in weight"
        if (norm2(p1(i)%st-p2(i)%st) .gt. 1d-30) write(*,*) "error in st"
        if (p1(i)%i_elm .ne. p2(i)%i_elm) write(*,*) "error in i_elm"
        if (p1(i)%q .ne. p2(i)%q) write(*,*) "error in q"
        if (p1(i)%label .ne. p2(i)%label) write(*,*) "error in label"
        if (p1(i)%lost .neqv. p2(i)%lost) write(*,*) "error in lost"
        if (norm2(p1(i)%v-p2(i)%v) .gt. 1d-30) write(*,*) "error in v"
      end do
    end select
  end select
end subroutine

!> Helper function for allocating particles
subroutine allocate_particles(particles, n)
  class(particle_base), allocatable, dimension(:), intent(out) :: particles
  integer, intent(in) :: n
  integer :: i
  ! Prepare some particles
  allocate(particle_boris::particles(n))
  select type (p => particles)
  type is (particle_boris)
    do i=1,size(p,1)
      p(i)%x = real((/i,i+1,i+2/),8)
      p(i)%mass = real(i,4)/3
      p(i)%weight = real(i+1,4)/4
      p(i)%st = real((/i,i+1/),8)/(2.d0*size(p,1))
      p(i)%i_elm = 100*i
      p(i)%q = int(3*i,1)
      p(i)%label = int(4*i,1)
      p(i)%lost = .true.
      p(i)%v = real((/2*i,2*i+1,2*i+2/),8)
    end do
  end select
end subroutine

!> Setup routine for testing
subroutine base_hdf5_io_test_setup
  integer :: ierr
  call h5open_f(ierr)
end subroutine base_hdf5_io_test_setup





!> Test writing a single particle
subroutine test_write_single
  class(particle_base), allocatable, dimension(:) :: particles
  call allocate_particles(particles, 1)
  call write_read_test(particles, 'test_base_hdf5_io.h5')
end subroutine test_write_single


!> Test writing many particles
subroutine test_write_many
  class(particle_base), allocatable, dimension(:) :: particles
  call allocate_particles(particles, 100000)
  call write_read_test(particles, 'test_base_hdf5_io.h5')
end subroutine test_write_many


end module base_hdf5_io_test
