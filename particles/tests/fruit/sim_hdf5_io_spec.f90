!> This module contains some testcases for hdf5 io
!> Cases that should be added:
!> - Read with different type as written
!> - Write different type from io%particle_type
!> - How to read from unknown type?
module sim_hdf5_io_spec
use mod_io_actions
use mod_particle_sim
use mod_boris
use fruit
implicit none

contains

!> Test the filename generating routine for times > 0
subroutine test_get_filename
  class(write_action), allocatable :: writer
  allocate(writer)

  call assert_equals('part000.00000000.h5', trim(writer%get_filename(0.d0)), 'test default settings')
  call assert_equals('part001.10000000.h5', trim(writer%get_filename(1.1d0)), 'test default settings 2')
  writer%decimal_digits = 2; writer%fractional_digits = 0
  call assert_equals('part21.h5', trim(writer%get_filename(21d0)), 'decimal point should be removed if fractional_digits = 0')
  writer%decimal_digits = 0; writer%fractional_digits = 0
  call assert_equals('part.h5', trim(writer%get_filename(12.d0)), 'test without numbers')
  writer%decimal_digits = 1; writer%fractional_digits = 3
  call assert_equals('part1.123.h5', trim(writer%get_filename(1.123d0)), 'test manual format')
  call assert_equals('part*.123.h5', trim(writer%get_filename(13.123d0)), 'test manual format overflow')
  writer%basename = 'testing'; writer%extension = '.rst'
  call assert_equals('testing2.123.rst', trim(writer%get_filename(2.123d0)), 'test manual format with extension and basename')
end subroutine test_get_filename


subroutine test_write_read_sim_time
  type(particle_sim) :: sim_to_write, sim_to_read
  class(write_action), allocatable :: writer
  class(read_action), allocatable  :: reader
  logical :: file_exists
  integer :: i, u, stat
  allocate(writer, reader)

  sim_to_write%time = 21.19d0
  writer%decimal_digits = 2; writer%fractional_digits = 0

  call writer%run(sim_to_write)
  ! test if a file with the right name was created
  inquire(file='part21.h5', exist=file_exists)
  call assert_true(file_exists, 'file with the right name should be created')

  reader%filename = 'part21.h5'
  call reader%run(sim_to_read)

  ! Test that the right time was read
  call assert_equals(sim_to_write%time, sim_to_read%time, "time should be read from the file")
  ! Delete the file
  open(newunit=u, iostat=stat, file='part21.h5', status='old')
  if (stat .eq. 0) close(u, status='delete')
end subroutine test_write_read_sim_time

subroutine test_write_sim_one_particle_kinetic_leapfrog
  type(particle_sim) :: sim_to_write, sim_to_read
  class(write_action), allocatable :: writer
  class(read_action), allocatable  :: reader
  logical :: file_exists
  integer :: i, u, stat, n_groups, n_particles
  allocate(writer, reader)

  allocate(sim_to_write%groups(1))
  call allocate_particles(sim_to_write%groups(1)%particles, 1)
  sim_to_write%time = 21.d0

  call writer%run(sim_to_write)
  ! test if a file with the right name was created
  inquire(file='part021.00000000.h5', exist=file_exists)
  call assert_true(file_exists, 'file with the right name should be created')

  reader%time = sim_to_write%time
  call reader%run(sim_to_read)
  ! Test that we have the right stuff in sim_to_read now
  n_groups = size(sim_to_read%groups, 1)
  call assert_equals(1, n_groups, 'should have one group exactly')
  if (n_groups .eq. 1 .and. lbound(sim_to_read%groups,1) .eq. 1) then
    n_particles = size(sim_to_read%groups(1)%particles, 1)
    call assert_equals(1, n_particles, 'should have one particle exactly')
    if (n_particles .eq. 1 .and. lbound(sim_to_read%groups(1)%particles,1) .eq. 1) then
      call assert_true(particles_same(sim_to_write%groups(1)%particles(1), sim_to_read%groups(1)%particles(1)), &
          'particle i must be as written')
    end if
  end if

  ! Delete the file
  open(newunit=u, iostat=stat, file='part021.00000000.h5', status='old')
  if (stat .eq. 0) close(u, status='delete')
end subroutine test_write_sim_one_particle_kinetic_leapfrog

subroutine test_write_sim_one_group_boris
  type(particle_sim) :: sim_to_write, sim_to_read
  class(write_action), allocatable :: writer
  class(read_action), allocatable  :: reader
  logical :: file_exists
  integer :: i, u, stat, n_groups, n_particles
  allocate(writer, reader)

  allocate(sim_to_write%groups(1))
  call allocate_particles(sim_to_write%groups(1)%particles, 2)
  sim_to_write%time = 21.d0

  call writer%run(sim_to_write)
  ! test if a file with the right name was created
  inquire(file='part021.00000000.h5', exist=file_exists)
  call assert_true(file_exists, 'file with the right name should be created')

  reader%time = sim_to_write%time
  call reader%run(sim_to_read)
  ! Test that we have the right stuff in sim_to_read now
  n_groups = size(sim_to_read%groups, 1)
  call assert_equals(1, n_groups, 'should have one group exactly')
  if (n_groups .eq. 1 .and. lbound(sim_to_read%groups,1) .eq. 1) then
    n_particles = size(sim_to_read%groups(1)%particles, 1)
    call assert_equals(2, n_particles, 'should have two particles exactly')
    if (n_particles .eq. 2 .and. lbound(sim_to_read%groups(1)%particles,1) .eq. 1) then
      do i=1,2
        call assert_true(particles_same(sim_to_write%groups(1)%particles(i), sim_to_read%groups(1)%particles(i)), &
            'particle i must be as written')
      end do
    end if
  end if

  ! Delete the file
  open(newunit=u, iostat=stat, file='part021.00000000.h5', status='old')
  if (stat .eq. 0) close(u, status='delete')
end subroutine test_write_sim_one_group_boris

subroutine test_write_sim_two_groups_boris
  type(particle_sim) :: sim_to_write, sim_to_read
  class(write_action), allocatable :: writer
  class(read_action), allocatable  :: reader
  logical :: file_exists
  integer :: i, j, u, stat, n_groups, n_particles
  allocate(writer, reader)

  allocate(sim_to_write%groups(2))
  call allocate_particles(sim_to_write%groups(1)%particles, 2)
  call allocate_particles(sim_to_write%groups(2)%particles, 2)
  sim_to_write%time = 21.d0

  call writer%run(sim_to_write)
  ! test if a file with the right name was created
  inquire(file='part021.00000000.h5', exist=file_exists)
  call assert_true(file_exists, 'file with the right name should be created')

  reader%time = sim_to_write%time
  call reader%run(sim_to_read)
  ! Test that we have the right stuff in sim_to_read now
  n_groups = size(sim_to_read%groups, 1)
  call assert_equals(2, n_groups, 'should have one group exactly')
  if (n_groups .eq. 2 .and. lbound(sim_to_read%groups,1) .eq. 1) then
    do i=1,2
      n_particles = size(sim_to_read%groups(i)%particles, 1)
      call assert_equals(2, n_particles, 'should have two particles exactly')
      if (n_particles .eq. 2 .and. lbound(sim_to_read%groups(i)%particles,1) .eq. 1) then
        do j=1,2
          call assert_true(particles_same(sim_to_write%groups(i)%particles(j), sim_to_read%groups(i)%particles(j)), &
              'particle j must be as written')
        end do
      end if
    end do
  end if

  ! Delete the file
  open(newunit=u, iostat=stat, file='part021.00000000.h5', status='old')
  if (stat .eq. 0) close(u, status='delete')
end subroutine test_write_sim_two_groups_boris


!> Helper function for allocating particles
subroutine allocate_particles(particles, n)
  class(particle_base), allocatable, dimension(:), intent(out) :: particles
  integer, intent(in) :: n
  integer :: i
  ! Prepare some particles
  allocate(particle_kinetic_leapfrog::particles(n))
  select type (p => particles)
  type is (particle_kinetic_leapfrog)
    do i=1,size(p,1)
      p(i)%x = real((/i,i+1,i+2/),8)
      p(i)%m = real(i,4)/3
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
  if (abs(p1%m-p2%m)     .gt. tolerance) same = .false.
  if (abs(p1%weight-p2%weight) .gt. tolerance) same = .false.
  if (norm2(p1%st-p2%st)       .gt. tolerance) same = .false.
  if (p1%i_elm .ne. p2%i_elm)  same = .false.
  if (p1%q     .ne. p2%q)      same = .false.
  if (p1%label .ne. p2%label)  same = .false.
  if (p1%lost  .neqv. p2%lost) same = .false.

  select type(p1 => p1)
    type is (particle_kinetic_leapfrog)
      select type (p2 => p2)
        type is (particle_kinetic_leapfrog)
          if (norm2(p1%v-p2%v) .gt. tolerance) same = .false.
        class default
          same = .false.
      end select
  end select
end function particles_same
end module sim_hdf5_io_spec
