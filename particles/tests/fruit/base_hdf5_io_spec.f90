module base_hdf5_io_spec
!< This module contains some testcases for hdf5 io
!< Cases that should be added:
!< - Read with different type as written
!< - Write different type from io%particle_type
!< - How to read from unknown type?
use mod_particle_io_boris
use fruit
implicit none

contains



subroutine test_write_single_particle_new_file
  class(particle_base), allocatable, dimension(:) :: particles
  character(len=*), parameter :: filename = '/tmp/base_hdf5_io_test.h5'
  integer :: u, stat

  call allocate_particles(particles, 1)
  ! Delete file if it exists
  open(newunit=u, iostat=stat, file=filename, status='old')
  if (stat == 0) close(u, status='delete')

  call write_read_particles_test(particles, filename)

  ! Delete the file again
  open(newunit=u, iostat=stat, file=filename, status='old')
  if (stat == 0) close(u, status='delete')
end subroutine test_write_single_particle_new_file


subroutine test_write_single_particle_existing_file
  class(particle_base), allocatable, dimension(:) :: particles
  character(len=*), parameter :: filename = '/tmp/base_hdf5_io_test.h5'
  integer :: u, stat

  call allocate_particles(particles, 1)
  ! Create file
  open(newunit=u, iostat=stat, file=filename, status='replace')
  write(u, "(A)") "garbage"
  close(u)

  call write_read_particles_test(particles, filename)
  
  ! Delete the file again
  open(newunit=u, iostat=stat, file=filename, status='old')
  if (stat == 0) close(u, status='delete')
end subroutine test_write_single_particle_existing_file


subroutine test_write_many_particles_new_file
  class(particle_base), allocatable, dimension(:) :: particles
  character(len=*), parameter :: filename = '/tmp/base_hdf5_io_test.h5'
  integer :: u, stat

  call allocate_particles(particles, 100000)
  call write_read_particles_test(particles, filename)

  ! Delete the file again
  open(newunit=u, iostat=stat, file=filename, status='old')
  if (stat == 0) close(u, status='delete')
end subroutine test_write_many_particles_new_file



!> Helper function for testing
subroutine write_read_particles_test(particles, filename)
  class(particle_base), allocatable, dimension(:), intent(in) :: particles
  character(len=*) :: filename

  class(particle_hdf5_io), allocatable :: io
  class(particle_base), allocatable, dimension(:) :: particles_read

  allocate(particle_hdf5_io_boris::io)
  ! Prepare HDF5 IO
  allocate(particle_boris::io%particle_type(2))
  call io%set_data_type

  call io%write(filename, particles)
  call io%read(filename, particles_read)

  select type (p1 => particles)
  type is (particle_boris)
    select type (p2 => particles_read)
    type is (particle_boris)
      call assert_equals(size(p1,1), size(p2,1), "read particle file length")

      ! Test the first and last particle
      call assert_true(particles_same(p1(lbound(p1,1)), p2(lbound(p2,1))), "First particle comparison")
      call assert_true(particles_same(p1(ubound(p1,1)), p2(ubound(p2,1))), "Last particle comparison")
    end select
  end select
end subroutine write_read_particles_test

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
      call random_number(p(i)%st)
      p(i)%i_elm = 100*i
      p(i)%q = int(3*i,1)
      p(i)%label = int(4*i,1)
      p(i)%lost = .true.
      p(i)%v = real((/2*i,2*i+1,2*i+2/),8)
    end do
  end select
end subroutine allocate_particles

!> Helper function for comparing particles
function particles_same(p1, p2) result(same)
  class(particle_base), intent(in) :: p1, p2
  real*8, parameter :: tolerance = 1d-30
  logical :: same
  same = .true.

  ! Base attributes testing
  if (norm2(p1%x-p2%x)         .gt. tolerance) same = .false.
  if (abs(p1%mass-p2%mass)     .gt. tolerance) same = .false.
  if (abs(p1%weight-p2%weight) .gt. tolerance) same = .false.
  if (norm2(p1%st-p2%st)       .gt. tolerance) same = .false.
  if (p1%i_elm .ne. p2%i_elm)  same = .false.
  if (p1%q     .ne. p2%q)      same = .false.
  if (p1%label .ne. p2%label)  same = .false.
  if (p1%lost  .neqv. p2%lost) same = .false.

  select type(p1 => p1)
    type is (particle_boris)
      select type (p2 => p2)
        type is (particle_boris)
          if (norm2(p1%v-p2%v) .gt. tolerance) same = .false.
        class default
          same = .false.
      end select
  end select
end function particles_same
end module base_hdf5_io_spec
