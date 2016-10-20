!> Module for the initialization of particles in configuration space (6D)
!> by rejection sampling
module mod_initialise_particles
use mod_rng
use data_structure
use mod_particle_types
use constants
implicit none
private
public seed_positions, no_transform, adjust_particle_weights
public set_velocity_from_T
contains
!> Set positions for particles by rejection sampling from geometric and mhd
!> variables after collecting with transform, within Rbound, Zbound and Phibound
!> if present. See [[test_rejection_sampling]] for examples.
subroutine seed_positions(particles, node_list, element_list, &
        rng, variables, transform, f, Rbound, Zbound, Phibound)
  use mpi
  use mod_sampling
  use mod_random_seed
  use mod_interp4
  !$ use omp_lib
  implicit none

  class(particle_base), dimension(:), intent(inout) :: particles
  type(type_node_list), intent(in)                  :: node_list
  type(type_element_list), intent(in)               :: element_list
  class(type_rng), intent(in)                       :: rng
  integer, dimension(:), intent(in), optional       :: variables !< Which variables from JOREK to use. If absent, sample uniformly.
  real*8, external, optional                        :: transform !< Merge variables into a single criterium between 0 and 1 for rej.  sampling
  !< Special values: 0 = 1, -1 = R, -2 = Z, -3 = Phi. Must be in ascending order!
  real*8, intent(in), optional                      :: f !< Weighting factor: f=0 indicates uniform weights, f=1 indicates uniform distribution
  !< (particle weight proportional to transform(P) at that point.) If omitted take f=0.
  real*8, dimension(2), intent(in), optional        :: Rbound, Zbound, Phibound !< Between which coordinates to sample (RZPhi).
  !< if omitted, determine automatically from node_list

  ! Internal variables
  real*8  :: R, Z, phi, s, t, DUMMY_REAL
  real*8  :: Rbox(2), Zbox(2), Phibox(2)
  integer :: i, j, k, ifail
  real*8  :: ran(4)
  integer :: i_elm
  real*8  :: t0, t1, ostart, oend
  integer :: seq, n_streams, n_threads, i_thread
  integer :: n_geom, n_mhd
  integer :: my_id, n_cpu
  integer :: seed
  real*8, dimension(:), allocatable :: P
  class(type_rng), allocatable, dimension(:) :: rngs ! The RNGs for all the threads
  integer, dimension(:), allocatable :: i_to_find
  logical, dimension(:), allocatable :: not_found

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ifail)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ifail)
  if (present(variables)) then
    if (.not. present(transform)) then
      write(*,*) "ERROR: if variables are present in set_particle_position_rejection_sampling transform must also be present"
      call MPI_ABORT(MPI_COMM_WORLD, 10, ifail)
    end if
    ! Get the number of mhd variables to use
    allocate(P(size(variables,1)))
    n_mhd = count(variables .gt. 0)
    n_geom = size(variables, 1) - n_mhd
  else
    n_mhd = 0
    n_geom = 0
  end if

  ! Setup bounding boxes
  call domain_bounding_box(node_list, element_list, Rbox(1), Rbox(2), Zbox(1), Zbox(2))
  Phibox = [0.d0, TWOPI]
  ! Check if the requested bounding boxes have any overlap with the domain
  ! Does not check for combinations of R and Z
  if (present(Rbound)) then
    if (maxval(Rbound) .gt. minval(Rbox) .and. maxval(Rbox) .gt. minval(Rbound)) then
      Rbox = Rbound
    else
      write(*,*) "ERROR: no overlap between domain and requested bounding box in R, domain=", Rbox, ", box=", Rbound
      write(*,*) "Sampling from whole domain in R"
    end if
  end if
  if (present(Zbound)) then
    if (maxval(Zbound) .gt. minval(Zbox) .and. maxval(Zbox) .gt. minval(Zbound)) then
      Zbox = Zbound
    else
      write(*,*) "ERROR: no overlap between domain and requested bounding box in Z, domain=", Zbox, ", box=", Zbound
      write(*,*) "Sampling from whole domain in Z"
    end if
  end if
  if (present(Phibound)) PhiBox = Phibound

  ! Calculate a single random seed and communicate it over MPI
  if (my_id .eq. 0) seed = random_seed()
  call MPI_Bcast(seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ifail)

  ! Prepare list of particles to seed
  allocate(i_to_find(size(particles,1)),not_found(size(particles,1)))
  i_to_find = [(i, i=1,size(particles,1))] ! which particles still to do
  not_found = .true. ! whether this one has been sampled succesfully

  ! Setup (Q)RNGs, one per thread
  n_threads = 1
  !$ n_threads = omp_get_max_threads()
  allocate(rngs(0:n_threads-1), source=rng)
  n_streams = n_cpu*n_threads ! Works only for homogeneous environments!

  do i_thread=0,n_threads-1
    seq = my_id*n_threads + i_thread + 1
    call rngs(i_thread)%initialize(4, seed, n_streams, seq, ifail)
    if (ifail .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, -1, ifail)
  end do

  call cpu_time(t0)
  !$ ostart = omp_get_wtime()

  ! Filter over all particles to sample them, repeat for rejected positions until
  ! empty. This is required if the distribution has some correlation with the
  ! samples mod something (as is the case for the Sobol sequence).
  ! Even then, an inbalance in openmp scheduling or number of particles per node
  ! could cause a slight correlation
  ! in the output. This could be worse if the distribution is very narrow.
  ! In that case the Sobol series should also be implemented in 64-bits as the total
  ! number of values is 2^31 now.
  do while (any(not_found))
  !$omp parallel default(none) &
  !$omp   shared(particles, node_list, element_list, Rbox, Zbox, PhiBox, variables, &
  !$omp          rngs, n_threads, n_streams, seed, my_id, n_mhd, n_geom, i_to_find, not_found) &
  !$omp   private(j, i, R, Z, phi, i_elm, s, t, ifail, seq, ran, i_thread, P, DUMMY_REAL)
  i_thread = 0
  !$ i_thread=omp_get_thread_num()
  !$omp do
  do i=1,size(i_to_find,1)
    j = i_to_find(i)
    ! Generate a random position to put this particle
    call rngs(i_thread)%next(ran)
    call transform_uniform_cylindrical(ran(1:3), Rbox, Zbox, PhiBox, R, Z, phi)

    call find_RZ(node_list,element_list,R,Z,DUMMY_REAL,DUMMY_REAL,i_elm,s,t,ifail)
    if (ifail .eq. 0) then
      if (present(variables)) then
        ! Select the mhd variables requested
        if (n_mhd .ge. 1) then
          call interp4(node_list,element_list,i_elm,variables(n_geom:n_geom+n_mhd),n_mhd,s,t,phi,P(n_geom:n_geom+n_mhd))
        end if
        do k=1,n_geom
          select case (variables(k))
          case (0);  P(k) = 1.d0
          case (-1); P(k) = R
          case (-2); P(k) = Z
          case (-3); P(k) = phi
          end select
        end do

        if (present(transform)) then
          if (ran(4) .lt. transform(P)) then
            particles(j)%x = [R, Z, phi]
            particles(j)%i_elm = i_elm
            particles(j)%st = [s, t]
            not_found(i) = .false.
          end if
        end if
      else
        particles(j)%x = [R, Z, phi]
        particles(j)%i_elm = i_elm
        particles(j)%st = [s, t]
        not_found(i) = .false.
      end if
    end if
  enddo
  !$omp end do
  !$omp end parallel
  ! Now pack only the indices of particles we still need to do
  i_to_find = pack(i_to_find, not_found) ! implicitly allocates
  deallocate(not_found); allocate(not_found(size(i_to_find,1)))
  not_found = .true.
  end do

  call cpu_time(t1)
  !$ oend = omp_get_wtime()
  write(*,'(i5,A,2f12.4)') my_id, ' Time particle initialize cpu/wall :',t1-t0, oend-ostart
  if (my_id .eq. 0) then
    write(*,*) '* done initialising particles    *'
    write(*,*) '**********************************'
  endif
end subroutine seed_positions


!> Calculate the size of a box around the domain (in RZ)
subroutine domain_bounding_box(node_list, element_list, Rmin, Rmax, Zmin, Zmax)
  type(type_node_list), intent(in)                  :: node_list
  type(type_element_list), intent(in)               :: element_list
  real*8, intent(out)                               :: Rmin, Rmax, Zmin, Zmax
  real*8 :: el_Rmin, el_Rmax, el_Zmin, el_Zmax
  integer :: i
  ! Initial setting
  call RZ_minmax(node_list, element_list, 1, Rmin, Rmax, Zmin, Zmax)
  do i=2,element_list%n_elements
    call RZ_minmax(node_list, element_list, i, el_Rmin, el_Rmax, el_Zmin, el_Zmax)
    Rmin = min(Rmin, el_Rmin)
    Rmax = max(Rmax, el_Rmax)
    Zmin = min(Zmin, el_Zmin)
    Zmax = max(Zmax, el_Zmax)
  enddo
end subroutine domain_bounding_box


!> Dummy function to use when no transform is desired. Copies the first parameter
!> into the output or sets out to 1.
pure function no_transform(in) result(out)
  real*8, dimension(:), intent(in) :: in
  real*8 :: out
  if (size(in,1) .gt. 0) then
    out = in(1)
  else
    out = 1.d0
  end if
end function no_transform


!> Adjust weights on all particles to have the correct number of atoms in total
subroutine adjust_particle_weights(particles, num_atoms_total)
use mpi
class(particle_base), intent(inout), dimension(:) :: particles
real*8, intent(in)                                :: num_atoms_total !< What the sum of the weights should be
real*8 :: local_weights, sum_weights
integer :: ifail
local_weights = sum(particles(:)%weight)
call MPI_AllReduce(local_weights,sum_weights,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ifail)
! Divide all weights by the sum of weights and multiply by the requested number of atoms
particles(:)%weight = particles(:)%weight / sum_weights * num_atoms_total
end subroutine adjust_particle_weights

!> Set v of a particle for use with kinetic codes
subroutine set_velocity_from_T(particles, mass, node_list, element_list, rng, cor, v_par, grad_T)
use constants
use data_structure
use phys_module, only: central_density, central_mass, F0
use mod_sampling
use mpi
use mod_random_seed
use mod_coordinate_transforms
use mod_coronal
implicit none

class(particle_base), intent(inout), dimension(:) :: particles !< Particle to initialize
real*8, intent(in)                                :: mass
type(type_node_list), intent(in)                  :: node_list
type(type_element_list), intent(in)               :: element_list
class(type_rng), intent(in)                       :: rng
type(coronal), intent(in), optional               :: cor !< Coronal equilibrium datatype for this particle. If unset, do not alter q
logical, intent(in), optional                     :: v_par !< Include the parallel velocity if present and true
logical, intent(in), optional                     :: grad_T !< Include the temperature gradient if present and true

class(type_rng), allocatable :: my_rng
integer :: i, ifail, seed, my_id, n_cpu
real*8, dimension(4) :: P, P_s, P_t, P_phi
real*8 :: v_out(4), ran4(4), Psi_R, Psi_Z, B(3)
real*8 :: R, R_s, R_t, Z, Z_s, Z_t
real*8 :: background_kbT, background_Kelvin, background_density, V_thermal
real*8 :: DUMMY_REAL, Z_coronal

allocate(my_rng, source=rng)
! Calculate a single random seed and communicate it over MPI
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ifail)
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ifail)
if (my_id .eq. 0) seed = random_seed()
call MPI_Bcast(seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ifail)
call my_rng%initialize(4, seed, n_cpu, my_id, ifail)

#if (JOREK_MODEL < 300)
if (present(v_par) .and. v_par) then
  write(*,*) "ERROR: initialization with v// not possible with this model"
  call MPI_ABORT(-1, MPI_COMM_WORLD, ifail)
end if
#endif
if (present(grad_T) .and. grad_T) then
  write(*,*) "ERROR: initialization with gradT not implemented"
  call MPI_ABORT(-1, MPI_COMM_WORLD, ifail)
end if

do i=1,size(particles)
#if (JOREK_MODEL == 400)
  call interp_PRZ(node_list,element_list,particles(i)%i_elm,[1,5,8,7],4,particles(i)%st(1),particles(i)%st(2),particles(i)%x(3),&
      P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
#else
  call interp_PRZ(node_list,element_list,particles(i)%i_elm,[1,5,6,7],4,particles(i)%st(1),particles(i)%st(2),particles(i)%x(3),&
      P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
#endif

  background_density = P(2) * 1d20                           ! plasma density [1/m^3]
  ! Assume that the particles have the same temperature as the electrons
#if (JOREK_MODEL == 400)
  background_kbT = P(3)/(MU_ZERO*central_density*1.d20)      ! P(1) contains the electron temperature
#else
  background_kbT = P(3)/(2.d0*MU_ZERO*central_density*1.d20) ! P(1) contains the total plasma temperature in J/kB = T = Te + Ti
#endif
  V_thermal = sqrt(background_kbT / (2.d0*mass*ATOMIC_MASS_UNIT))      ! [m/s]

  call my_rng%next(ran4)
  v_out = boxmueller_transform(ran4) * V_thermal ! [m/s], vx, vy, vz, dummy
  if (present(v_par) .and. v_par) then
    psi_R = (  P_s(4) * Z_t - P_t(4) * Z_s )/(R_s * Z_t - R_t * Z_s)
    psi_Z = (- P_s(4) * R_t + P_t(4) * R_s )/(R_s * Z_t - R_t * Z_s)
    B     = [psi_Z, -psi_R, F0]/R
    v_out(1:3) = v_out(1:3) + P(4)*B/norm2(B)
  end if

  background_kelvin  = background_kbT / K_BOLTZ              ! electron temperature [K]
  if (background_density .le. 0.d0 .or. background_kelvin .le. 0.d0) then
    Z_coronal = 0.d0
  else
    call cor%interp(log10(background_density),log10(background_kelvin),Z_coronal,DUMMY_REAL)
  endif

  select type (p => particles(i))
  type is (particle_kinetic_leapfrog)
    p%v = vector_rotation(v_out(1:3), p%x(3)) ! rotate the velocity vector to the correct angle
    p%q = nint(Z_coronal,1)                   !< charge
  end select
end do
end subroutine set_velocity_from_T
end module mod_initialise_particles
