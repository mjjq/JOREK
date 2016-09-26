!> Redistribute particles over CPUs to do load-balancing
module mod_redistribute_particles
  use mod_particle_types
contains

!> This function creates a derived MPI type for the particle and returns it
!> If it already exists the old handle is returned
!> TODO: alter for polymorphic particles
function get_particle_derived_type() result(dtype_out)
  use mpi
  use parameters

  implicit none

  integer               :: ierr, dtype_out
  integer, save         :: dtype
  logical, save         :: dtype_set = .false.

  integer :: len(9) = (/n_dim,3,3,1,1,1,1,1,1/), t(9) = (/ &
    MPI_REAL8,MPI_REAL8,MPI_REAL8,MPI_REAL4,MPI_REAL4, &
    MPI_INTEGER,MPI_INTEGER1,MPI_INTEGER1,MPI_INTEGER1/) ! MPI_INTEGER1 == MPI_LOGICAL1

  integer(kind=MPI_ADDRESS_KIND) :: base, disp(9)
  type(type_particle) :: particle

  dtype_out = dtype
  if (dtype_set) return

  ! Get memory addresses in the type
  call MPI_Get_address(particle,        base,    ierr)
  call MPI_Get_address(particle%st,     disp(1), ierr)
  call MPI_Get_address(particle%x,      disp(2), ierr)
  call MPI_Get_address(particle%v,      disp(3), ierr)
  call MPI_Get_address(particle%mass,   disp(4), ierr)
  call MPI_Get_address(particle%weight, disp(5), ierr)
  call MPI_Get_address(particle%i_elm,  disp(6), ierr)
  call MPI_Get_address(particle%q,      disp(7), ierr)
  call MPI_Get_address(particle%label,  disp(8), ierr)
  call MPI_Get_address(particle%lost,   disp(9), ierr)

  ! Rebase to particle memory beginning
  disp = disp - base

  ! Commit the structured type
  call MPI_Type_create_struct(9, len, disp, t, dtype, ierr)
  if (ierr .ne. 0) write(*,*) "Error creating particle datatype: ", ierr
  call MPI_Type_commit(dtype, ierr)
  if (ierr .ne. 0) write(*,*) "Error committing particle datatype: ", ierr

  ! Set the save bit
  dtype_set = .true.
  dtype_out = dtype
  return
end function get_particle_derived_type

!> Append a single particle to the list and grow it if needed
pure subroutine append_particle_to_list(particle_list, particle)
  implicit none

  type(type_particle_list), intent(inout) :: particle_list
  type(type_particle), intent(in) :: particle

  if (size(particle_list%particle,1) .lt. particle_list%n_particles+1) then
    call grow_particle_list(particle_list)
  endif
  particle_list%n_particles = particle_list%n_particles + 1
  particle_list%particle(particle_list%n_particles) = particle
end subroutine append_particle_to_list


!> Append a list of particles to the list and grow it if necessary
pure subroutine append_particles_to_list(particle_list, particles)
  implicit none

  type(type_particle_list), intent(inout) :: particle_list
  type(type_particle), dimension(:), intent(in) :: particles

  ! Grow it until it fits
  do while (size(particle_list%particle,1) .lt. particle_list%n_particles+size(particles,1))
    call grow_particle_list(particle_list)
  enddo
  particle_list%particle(particle_list%n_particles+1:particle_list%n_particles+size(particles,1)) = particles
  particle_list%n_particles = particle_list%n_particles + size(particles,1)
end subroutine append_particles_to_list


!> Grow the particle list by a specific factor
pure subroutine grow_particle_list(particle_list)
  implicit none
  type(type_particle_list), intent(inout) :: particle_list
  type(type_particle), dimension(:), allocatable :: temp

  real*8, parameter :: growth_factor = 1.2d0

  allocate(temp(lbound(particle_list%particle,1):ubound(particle_list%particle,1)+ &
      int((growth_factor - 1.d0) * real(size(particle_list%particle,1), 8))))
  temp(lbound(particle_list%particle,1):ubound(particle_list%particle,1)) = particle_list%particle
  call move_alloc(from=temp,to=particle_list%particle) ! deallocates temp as well
end subroutine grow_particle_list

!> Redistribute particles over all participating processes
!> Calculates the speed of calculation on this processor
!> And sends particles to try to get the calculations times equal
!> NB. This sends lost particles as well, even though these do not contribute
!> to the computation time. If there are many lost particles this could be
!> inefficient. Does not take into account different integrators for different species.
subroutine redistribute_particles(particles, wtime)
  use mpi
  implicit none

  type(type_particle_list), intent(inout) :: particles
  real*8, intent(in) :: wtime !< Particle pusher calculation time on this CPU

  type(type_particle), allocatable, dimension(:) :: particles_recv
  integer :: n_particles_recv
  integer :: my_id, n_cpu, ierr

  integer, allocatable, dimension(:) :: wants, nums
  real*8, allocatable, dimension(:) :: cis
  real*8 :: speed
  integer :: num_particles
  integer :: i, send_begin, j

  integer :: dtype, status(MPI_STATUS_SIZE)
  integer :: lost
  integer :: from, to, n, cond, nonlost

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs

  allocate(wants(0:n_cpu-1),cis(0:n_cpu-1),nums(0:n_cpu-1))
  dtype = get_particle_derived_type()

  ! Calculate the particle processing time of this proc
  ! First get the number of lost particles (this could have also been saved,
  ! stupid step!)
  lost=0
  !$omp parallel do default(none) shared(particle_list) reduction(+:lost) private(i)
  do i=1,particle_list%n_particles
    if (particle_list%particle(i)%lost) lost = lost + 1
  enddo
  !$omp end parallel do

  ! Get number of non-lost particles and particle processing time
  num_particles = particle_list%n_particles - lost
  speed = num_particles / wtime

  ! Communicate this to every node
  call MPI_AllGather(num_particles, 1, MPI_INTEGER, nums, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
  call MPI_AllGather(speed, 1, MPI_REAL8, cis, 1, MPI_REAL8, MPI_COMM_WORLD, ierr)

  ! Only run if the imbalance is greater than 1 percent
  if (maxval(nums/cis,1)/minval(nums/cis,1) .le. 1.01) return

  ! Calculate how many non-lost particles eachs processor wants to have or lose
  wants = sum(nums) * (cis/sum(cis)) - nums

  ! Debug output
  if (my_id .eq. 0) write(*,"(2000i9)") wants

  ! Wants now contains the requested change per cpu, which should sum to
  ! something close to zero

  ! Use a heuristic algorithm to set up these moves
  ! Algorithm: take largest absolute value, fill it with largest value of
  ! opposite sign. Repeat until a specific threshold.
  i = 0
  do
    i = i + 2 ! 2 messages per iteration
    from = minloc(wants,1)-1 ! 1-based result
    to = maxloc(wants,1)-1 ! 1-based result
    n = min(wants(to),-wants(from)) ! number of nonlost particles to send
    cond = wants(to)-wants(from)

    ! Stop if the difference between wants and has is < this
    if (cond .lt. 10 .or. n .le. 0) exit ! do not set this <= 1

    ! Update list of wants
    wants(from) = wants(from) + n
    wants(to)   = wants(to)   - n

    ! Set up mpi data transfer on the procs involved
    if (my_id .eq. from) then
      ! find out how many particles to send to include n nonlost
      nonlost = 0
      do j=particle_list%n_particles,1,-1
        if (.not. particle_list%particle(j)%lost) nonlost = nonlost + 1
        if (nonlost .eq. n) then
          send_begin = j
          exit
        endif
      enddo
      write(*,"(A,i7,A,i4,A,i4,A,i7,A)") "Send ", particle_list%n_particles-send_begin+1, " particles from ", from, "->", to, &
        " including ", particle_list%n_particles-send_begin+1-n, " lost"
      ! Send inclusive, from send_begin to n_particles
      call MPI_Send(particle_list%n_particles-send_begin+1, 1, MPI_INTEGER, to, i, MPI_COMM_WORLD, ierr)
      call MPI_Send(particle_list%particle(send_begin:particle_list%n_particles), &
        particle_list%n_particles-send_begin+1, dtype, to, i+1, MPI_COMM_WORLD, ierr)
      ! Reset the end of the particle list
      particle_list%n_particles = send_begin-1

    else if (my_id .eq. to) then
      call MPI_Recv(n_particles_recv, 1, MPI_INTEGER, from, i, MPI_COMM_WORLD, status, ierr)
      allocate(particles_recv(n_particles_recv))
      call MPI_Recv(particles_recv, n_particles_recv, dtype, from, i+1, MPI_COMM_WORLD, status, ierr)

      call append_particles_to_list(particle_list, particles_recv)
      deallocate(particles_recv)
      !write(*,"(A,i7,A,i4,A,i4)") "Recv ", n_particles_recv, " particles from ", from, "->", to
    endif
  enddo

  ! Debug output
  if (my_id .eq. 0) write(*,"(2000i9)") wants

  deallocate(wants,nums,cis)

end subroutine redistribute_particles
end module mod_redistribute_particles
