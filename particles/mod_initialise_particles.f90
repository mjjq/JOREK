!> Module for the initialization of particles in configuration space (6D)
!> by rejection sampling
module mod_initialise_particles
use mod_rng
use data_structure
use mod_particle_types
use constants
use mod_interp_PRZ
implicit none
private
public initialise_particles, no_transform, adjust_particle_weights
public set_velocity_from_T, domain_bounding_box, initialise_particles_H_mu_psi
public set_particle_weights_canonical_maxwellian, normalize_with_projection
contains
!> Set positions for particles by rejection sampling from geometric and mhd
!> variables after collecting with transform, within Rbound, Zbound and Phibound
!> if present. See [[test_rejection_sampling]] for examples.
subroutine initialise_particles(particles, node_list, element_list, &
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
  class(type_rng), intent(in)                       :: rng !< What type of random number generator to use. Is re-seeded in the subroutine.
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
  real*8  :: ran(7)
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
    call rngs(i_thread)%initialize(7, seed, n_streams, seq, ifail)
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
  ! TODO fix also for MPI or broadcast to nodes
  do while (any(not_found))
  !$omp parallel default(none) &
  !$omp   shared(particles, node_list, element_list, Rbox, Zbox, PhiBox, variables, &
  !$omp          rngs, n_threads, n_streams, seed, my_id, n_mhd, n_geom, i_to_find, not_found) &
  !$omp   private(j, i, R, Z, phi, i_elm, s, t, ifail, seq, ran, i_thread, P, DUMMY_REAL)
  i_thread = 0
  !$ i_thread=omp_get_thread_num()
  !$omp do schedule(static)
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
          if (ran(4) .lt. transform(p)) then
            particles(j)%x = [r, z, phi]
            particles(j)%i_elm = i_elm
            particles(j)%st = [s, t]
            select type (pa => particles(j))
            type is (particle_kinetic_leapfrog)
              pa%v = ran(5:7) ! save other components of this point for velocity init in a later routine
            end select
            not_found(i) = .false.
          end if
        end if
      else
        particles(j)%x = [r, z, phi]
        particles(j)%i_elm = i_elm
        particles(j)%st = [s, t]
        select type (pa => particles(j))
        type is (particle_kinetic_leapfrog)
          pa%v = ran(5:7) ! save other components of this point for velocity init in a later routine
        end select
        not_found(i) = .false.
      end if
    end if
  enddo
  !$omp end do
  !$omp end parallel
  ! now pack only the indices of particles we still need to do
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
end subroutine initialise_particles

!> Initialise particle positions in H, mu, psi, theta, phi, gamma (gyrophase) space.
!> Set H_transform to transform from [0,1] to your desired range, Psi_transform to do the same (optional)
!>
!> Does not do weighting of the particles.
!>
!> **This subroutine does not support MPI or openMP yet!**
subroutine initialise_particles_H_mu_psi(particles, node_list, element_list, rng_base, mass, &
        H_transform, Theta_transform, Psi_transform, cor)
  use data_structure
  use mod_rng
  use mod_random_seed
  use constants
  use phys_module, only: F0
  use mod_coronal
  use mod_boris, only: left_handed_cross_product
  implicit none
  class(particle_base), dimension(:), intent(inout) :: particles
  type(type_node_list), intent(in)                  :: node_list
  type(type_element_list), intent(in)               :: element_list
  class(type_rng), intent(in)                       :: rng_base !< What type of random number generator to use (will be reseeded here)
  real*8, intent(in)                                :: mass
  real*8, external                                  :: H_transform !< Function to transform 0-1 to the H-domain (eV)
  real*8, external, optional                        :: Theta_transform !< Function to transform 0-1 to the theta-domain
  real*8, external, optional                        :: Psi_transform !< Function to transform 0-1 to the Psi-domain
  !< if omitted, determine automatically from node_list
  type(coronal), intent(in), optional               :: cor !< Coronal equilibrium datatype for this particle. If unset, do not alter q

  ! Internal variables
  class(type_rng), allocatable :: rng
  real*8  :: ran(6), Rbox(2), Zbox(2)
  real*8  :: H, muB, v_perp, v_par, gamma
  real*8  :: psi, psimin, psimax, theta, phi
  real*8  :: R, Z, inv_st_jac, psi_r, psi_z, B(3), B_hat(3)
  real*8, dimension(1) :: P, P_s, P_t, P_phi
  real*8  :: R_s, R_t, Z_s, Z_t
  real*8  :: s, t, u_init_max
  real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis
  integer :: i_elm, i, ifail, n_try
  real*8, dimension(element_list%n_elements,2) :: psi_minmax_list

  psimin= 1d10
  psimax=-1d10
  ! Preparatory work: determine psi_min,max
  !$omp parallel do default(shared) &
  !$omp private(i_elm) reduction(min:psimin) &
  !$omp reduction(max:psimax)
  do i_elm=1,element_list%n_elements
    call psi_minmax(node_list,element_list,i_elm,psi_minmax_list(i_elm,1),psi_minmax_list(i_elm,2))
    psimin = min(psi_minmax_list(i_elm,1),psimin)
    psimax = max(psi_minmax_list(i_elm,2),psimax)
  end do
  !$omp end parallel do
  ! Preparatory work: setup RNG
  allocate(rng,source=rng_base)
  call rng%initialize(6, random_seed(), 1, 1, ifail)
  ! Preparatory work: get R_axis, Z_axis
  call find_axis(0, node_list, element_list, psi_axis, R_axis, Z_axis, s_axis, t_axis, i_elm, ifail)


  ! Loop over all points in the series until we have enough particles
  n_try = 0
  do i=1,size(particles)
    i_elm = 0
    do while (i_elm .eq. 0)
      call rng%next(ran)

      if (present(Psi_transform)) then
        psi = Psi_transform(ran(1))
      else
        psi = (psimax-psimin)*ran(1)+psimin
      end if
      ! Try to find this position
      if (present(Theta_transform)) then
        theta = Theta_transform(ran(5))
      else
	theta = TWOPI*ran(5)
      end if
      phi = TWOPI*ran(4)
      ! 1. Find R, Z corresponding to psi, theta
      call find_theta_psi(node_list,element_list,psi_minmax_list,theta,psi,phi,R_axis,Z_axis,i_elm,s,t,R,Z)
      n_try = n_try + 1
    end do
    ! If we are here a suitable position has been found
    ! For now set the guiding center position of the particle
    particles(i)%i_elm = i_elm
    particles(i)%st = [s,t]
    particles(i)%x = [R,Z,phi]

    ! 1. Get B at this position
    call       interp_PRZ(node_list,element_list,i_elm,[1],1,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
    inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
    psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
    psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
    ! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
    B        = [+psi_Z, -psi_R, F0] / R
    B_hat = B/norm2(B)

    ! 2. Calculate v_perp and v_par
    H  = H_transform(ran(2)) ! [eV]
    muB = H*(ran(3)-0.5d0) ! uniformly distributed between -H and H [eV]
    v_perp = sqrt(2*abs(muB*EL_CHG)/(mass*ATOMIC_MASS_UNIT)) ! [m/s]
    v_par  = sign(sqrt(2*(H-muB)*EL_CHG/(mass*ATOMIC_MASS_UNIT)),muB)

    ! 3. Calculate charge (if cor is present)
    if (present(cor)) then
      select type (p => particles(i))
      type is (particle_kinetic_leapfrog)
        p%q = q_coronal(node_list, element_list, s, t, phi, i_elm, cor)
      end select
    end if

    ! 4. Output to particles (dependent on type of particle)
    gamma = TWOPI*ran(6)
    select type(p => particles(i))
    type is (particle_kinetic_leapfrog)
      p%v = v_par * B_hat + v_perp * &
      ((cos(gamma) * [0.d0, B_hat(3), -B_hat(2)]) + &
        sin(gamma) * (B_hat(1) * B_hat - [1.d0, 0.d0, 0.d0]))
      if (p%q .gt. 0) p%x = p%x + (mass*ATOMIC_MASS_UNIT*left_handed_cross_product(p%v,B_hat))/(real(p%q,8)*EL_CHG*norm2(B))
    end select
  end do
end subroutine initialise_particles_H_mu_psi

!> Calculate the particle weights according to the canonical maxwellian distribution function
!> (no electric fields):
!> \[ F_{MC}(P_\phi,H,\mu) = \frac{n(\bar\psi)}{\left[2\pi \bar T(\bar\psi)/m\right]^3/2}
!>                          exp\left{-\frac{H}{\bar T(\bar\psi)}\right} \]
!> where \(\bar\psi = P_\phi/q\), \(P_\phi = q\psi + m R v_\phi\),
!> \( H = m/2 v_\parallel^2 + \mu B \) and \(\mu = \frac{m v_\perp^2}{2B}\)
!>
!> The particle weight is set to the value of this distribution function at the specific point.
!> \(\bar T(\bar\psi)\) is approximated by \(T(\psi)\) if it is missing.
!> \(\bar n(\bar\psi)\) is approximated by 1 if it is missing.
subroutine set_particle_weights_canonical_maxwellian(particles, node_list, element_list, mass, n_psibar, T_psibar)
  use data_structure
  use constants
  use phys_module, only: central_density, central_mass
  implicit none
  class(particle_base), dimension(:), intent(inout) :: particles
  type(type_node_list), intent(in)                  :: node_list
  type(type_element_list), intent(in)               :: element_list
  real*8, intent(in)                                :: mass
  real*8, external, optional                        :: n_psibar
  real*8, external, optional                        :: T_psibar

  interface
    function n_psibar(psibar)
      real*8, intent(in) :: psibar
    end function n_psibar
    function T_psibar(psibar)
      real*8, intent(in) :: psibar
    end function T_psibar
  end interface

  integer :: i
  real*8  :: t_norm, psibar, H, n, T
  real*8, dimension(1) :: P, P_s, P_t, P_phi
  real*8  :: R, R_s, R_t, Z, Z_s, Z_t

  t_norm = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)

  !$omp parallel do default(none) private(i, psibar, H, n, T, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t) &
  !$omp shared(particles, node_list, element_list, mass, central_density)
  do i=1,size(particles,1)
    call       interp_PRZ(node_list,element_list,particles(i)%i_elm,[1],1, &
        particles(i)%st(1),particles(i)%st(2),particles(i)%x(3),P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
    select type (pa => particles(i))
    type is (particle_kinetic_leapfrog)
      psibar = real(pa%q,8) * P(1) * EL_CHG + mass * ATOMIC_MASS_UNIT * R * pa%v(3)
      H = mass*ATOMIC_MASS_UNIT*0.5d0 * norm2(pa%v)
    class default
      write(*,*) "ERROR: add code for your type here"
    end select

    if (present(n_psibar)) then
      n = n_psibar(psibar) ! [m^-3]
    else
      n = 1 ! units irrelevant if normalized later to a total number of particles
    end if
    if (present(T_psibar)) then
      T = T_psibar(psibar) ! [K]
    else
      ! Calculate the local temperature and use this instead
      call       interp_PRZ(node_list,element_list,particles(i)%i_elm,[6],1, &
        particles(i)%st(1),particles(i)%st(2),particles(i)%x(3),P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
      T = P(1)/(2.d0*MU_ZERO*central_density*1.d20*K_BOLTZ) ! [K]
#if (JOREK_MODEL == 400)
      T = T*2d0 ! P(1) contains the ion temperature in this model, reverse previous correction
#endif
      ! Workaround for low-temperature regions
      ! Because we give weights based on the energy and temperature areas with lower temperature are getting
      ! too high weights. Work around this by defining a minimum temperature to stop the outer regions from 
      ! dominating the projection.
      T = max(1d7,T)
    end if
    particles(i)%weight = n/(TWOPI*T/(mass*ATOMIC_MASS_UNIT)) * exp(-H/T)
  end do
  !$omp end parallel do
end subroutine set_particle_weights_canonical_maxwellian

!> Normalize particles with the result of the projection of the first group
subroutine normalize_with_projection(proj, particles, i_group)
  use mod_project_particles
  type(project_to_vtk), intent(in) :: proj
  class(particle_base), dimension(:), intent(inout) :: particles
  integer, intent(in), optional :: i_group

  integer :: group = 1
  integer :: i
  real*8, dimension(1) :: P, P_s, P_t, P_phi
  real*8 :: R, R_s, R_t, Z, Z_s, Z_t
  if (present(i_group)) group = i_group

  do i=1,size(particles,1)
    if (particles(i)%i_elm .ne. 0) then
      call interp_PRZ(proj%node_list,proj%element_list,particles(i)%i_elm,[group],1, &
        particles(i)%st(1),particles(i)%st(2),particles(i)%x(3),P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
      particles(i)%weight = particles(i)%weight/P(1)
    end if
  end do
end subroutine normalize_with_projection


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

function q_coronal(node_list, element_list, s, t, phi, i_elm, cor)
use data_structure
use phys_module, only: central_density
use mod_coronal
type(type_node_list), intent(in)                  :: node_list
type(type_element_list), intent(in)               :: element_list
real*8, intent(in)                                :: s, t, phi
integer, intent(in)                               :: i_elm
type(coronal), intent(in)                         :: cor
integer                                           :: q_coronal

real*8, dimension(2) :: P, P_s, P_t, P_phi
real*8               :: R, R_s, R_t, Z, Z_s, Z_t, q
real*8 :: local_Te, local_Ne, DUMMY_REAL
call interp_PRZ(node_list,element_list,i_elm,&
#if (JOREK_MODEL == 400)
      [5,8],& ! electron temperature
#else
      [5,6],& ! electron temperature + ion temperature (assumed equal)
#endif
          2,s,t,phi,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)

local_Ne = P(1) * 1d20                           ! plasma density [1/m^3]
local_Te = P(2)/(2.d0*MU_ZERO*central_density*1.d20)/K_BOLTZ
#if (JOREK_MODEL == 400)
local_Te = local_T_e*2d0 ! P(1) contains the electron temperature, reverse previous correction
#endif

if (local_Ne .le. 0.d0 .or. local_Te .le. 0.d0) then
  q_coronal = 0
else
  call cor%interp(log10(local_Ne),log10(local_Te),q,DUMMY_REAL)
  q_coronal = nint(q,1)
endif
end function

!> Set v of a particle for use with kinetic codes
subroutine set_velocity_from_T(particles, mass, node_list, element_list, cor, v_par)
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
type(coronal), intent(in), optional               :: cor !< Coronal equilibrium datatype for this particle. If unset, do not alter q
logical, intent(in), optional                     :: v_par !< Include the parallel velocity if present and true

class(type_rng), allocatable :: my_rng
integer :: i, ifail, seed, my_id, n_cpu
real*8, dimension(4) :: P, P_s, P_t, P_phi
real*8 :: v_out(3)
real*8 :: R, R_s, R_t, Z, Z_s, Z_t, Psi, Psi_R, Psi_Z, B_hat(3)
real*8, parameter :: r_hat(3) = [1.d0, 0.d0, 0.d0]
real*8 :: background_kbT, background_Kelvin, background_density, V_thermal
real*8 :: DUMMY_REAL, Z_coronal, t_norm

t_norm = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)

! Calculate a single random seed and communicate it over MPI
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ifail)
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ifail)
if (my_id .eq. 0) seed = random_seed()
call MPI_Bcast(seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ifail)

#if (JOREK_MODEL < 300)
if (present(v_par) .and. v_par) then
  write(*,*) "ERROR: initialization with v// not possible with this model"
  call MPI_ABORT(-1, MPI_COMM_WORLD, ifail)
end if
#endif

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

  ! Only an implementation for particle_kinetic_leapfrog now
  select type (pa => particles(i))
  type is (particle_kinetic_leapfrog)
    ! v_out now contains parallel and perpendicular velocities and the gyrophase
    v_out(1:2) = boxmueller_transform(pa%v(1:2))*V_thermal ! 2 gaussian distributed random numbers
    v_out(2)   = v_out(2)*sqrt(2.d0) ! equipartition of energy in parallel and 2 perpendicular directions
    v_out(3)   = pa%v(3)*TWOPI ! gyrophase, between 0 and 2PI
    if (present(v_par) .and. v_par) then
      v_out(1) = v_out(1) + P(4)/t_norm
    end if

    background_kelvin  = background_kbT / K_BOLTZ              ! electron temperature [K]
    if (background_density .le. 0.d0 .or. background_kelvin .le. 0.d0) then
      Z_coronal = 0.d0
    else
      call cor%interp(log10(background_density),log10(background_kelvin),Z_coronal,DUMMY_REAL)
    endif

    ! Calculate b^ (unit vector in direction of B)
    psi_R = (  P_s(4) * Z_t - P_t(4) * Z_s )/(R_s * Z_t - R_t * Z_s)
    psi_Z = (- P_s(4) * R_t + P_t(4) * R_s )/(R_s * Z_t - R_t * Z_s)
    B_hat = [psi_Z, -psi_R, F0]/(R)
    B_hat = B_hat/norm2(B_hat)

    ! Transform parallel and perpendicular velocities to R, Z, Phi
    ! To get the perpendicular vector, get a single vector perpendicular to b (b x r)
    ! and rotate it by another vector perpendicular to b.
    ! use the vector triple product to simplify.
    ! I'm not sure if this formula is the same in a right-handed coordinate system...
    ! this might change the direction of the rotation, but that is not important.
    pa%v = v_out(1) * B_hat + v_out(2) * &
      ((cos(v_out(3)) * [0.d0, B_hat(3), -B_hat(2)]) + &
        sin(v_out(3)) * (B_hat(1) * B_hat - r_hat))
    pa%q = nint(Z_coronal,1)                   !< charge
  class default
    write(*,*) "set_velocity_from_T not implemented for this particle type"
  end select
end do
end subroutine set_velocity_from_T
end module mod_initialise_particles
