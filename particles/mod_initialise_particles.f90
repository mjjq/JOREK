!> Module for the initialization of particles in configuration space (6D)
!> by rejection sampling
module mod_initialise_particles
use mod_rng
use data_structure
use mod_particle_types
use constants
implicit none
private
public particle_init_params
public set_particle_position_rejection_sampling, no_transform

type particle_init_params
  real*8 :: f !< Weighting factor: f=0 indicates uniform weights, f=1 indicates uniform distribution
end type particle_init_params

contains

!> Set positions for particles by rejection sampling from geometric and mhd
!> variables after collecting with transform, within Rbound, Zbound and Phibound
!> if present.
subroutine set_particle_position_rejection_sampling(particles, node_list, element_list, &
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
  interface
    pure function transform(in)
      real*8, intent(in), dimension(:) :: in !< Values of variables(:) at the particle position
    end function transform
  end interface

  ! Internal variables
  real*8  :: R, Z, phi, s, t, DUMMY_REAL
  real*8  :: Rbox(2), Zbox(2), Phibox(2)
  integer :: i, j, ifail
  real*8  :: ran(4)
  integer :: i_elm
  real*8  :: t0, t1, ostart, oend
  integer :: seq, n_streams, n_threads, i_thread
  integer :: n_geom, n_mhd
  integer :: my_id, n_cpu
  integer :: seed
  real*8, dimension(:), allocatable :: P
  class(type_rng), allocatable, dimension(:) :: rngs ! The RNGs for all the threads

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
  if (present(Rbound)) Rbox = Rbound
  if (present(Zbound)) Zbox = Zbound
  if (present(Phibound)) PhiBox = Phibound

  ! Setup (Q)RNGs, one per thread
  n_threads = 1
  !$ n_threads = omp_get_max_threads()
  allocate(rngs(0:n_threads-1), source=rng)
  n_streams = n_cpu*n_threads ! Works only for homogeneous environments!

  ! Calculate a single random seed and communicate it over MPI
  if (my_id .eq. 0) seed = random_seed()
  call MPI_Bcast(seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ifail)

  call cpu_time(t0)
  !$ ostart = omp_get_wtime()

  !$omp parallel default(none) &
  !$omp   shared(node_list, element_list, Rbox, Zbox, PhiBox, variables, &
  !$omp          n_cpu, rngs, n_threads, n_streams, seed, my_id, n_mhd, n_geom) &
  !$omp   private(j, i, R, Z, phi, i_elm, s, t, ifail, seq, ran, i_thread, P, DUMMY_REAL)

  i_thread = 0
  !$ i_thread = omp_get_thread_num()
  seq = my_id*n_threads + i_thread + 1
  ! generate random numbers in 4 dimensions
  call rngs(i_thread)%initialize(4, seed, n_streams, seq, ifail)
  if (ifail .ne. 0) then
    write(*,*) "Error seeding rng: ", ifail
    call MPI_ABORT(MPI_COMM_WORLD, -1, ifail)
  endif

  !$omp do
  do j=1,size(particles,1)
    ifail = 1
    do while (ifail .ne. 0)
      ! Generate a random position to put this particle
      call rngs(i_thread)%next(ran)
      call transform_uniform_cylindrical(ran(1:3), Rbox, Zbox, PhiBox, R, Z, phi)

      call find_RZ(node_list,element_list,R,Z,DUMMY_REAL,DUMMY_REAL,i_elm,s,t,ifail)
      if (ifail .ne. 0) cycle ! out of domain

      ! Select the mhd variables requested
      if (n_mhd .ge. 1) then
        call interp4(node_list,element_list,i_elm,variables(n_geom:n_geom+n_mhd),n_mhd,s,t,phi,P(n_geom:n_geom+n_mhd))
      end if
      do i=1,n_geom
        select case (variables(i))
        case (0);  P(i) = 1.d0
        case (-1); P(i) = R
        case (-2); P(i) = Z
        case (-3); P(i) = phi
        end select
      end do

      if (present(variables) .and. present(transform)) then
        if (ran(4) .lt. transform(P)) then
          particles(j)%x = [R, Z, phi]
          particles(j)%i_elm = i_elm
          particles(j)%st = [s, t]
          ifail = 0
        else
          ifail = 1
        endif
      else
        particles(j)%x = [R, Z, phi]
        particles(j)%i_elm = i_elm
        particles(j)%st = [s, t]
        ifail = 0
      end if
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call cpu_time(t1)
  !$ oend = omp_get_wtime()
  write(*,'(i5,A,2f12.4)') my_id, ' Time particle initialize cpu/wall :',t1-t0, oend-ostart
  if (my_id .eq. 0) then
    write(*,*) '* done initialising particles    *'
    write(*,*) '**********************************'
  endif
end subroutine set_particle_position_rejection_sampling


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


!> Set v and q of a particle for use with kinetic codes
subroutine set_charge_from_coronal_eq(particles, cor, ran4, ifail)
use constants
use data_structure
use nodes_elements
use phys_module, only : central_density, central_mass
use mod_openadas
use mod_sampling
use mod_coronal
implicit none

class(particle_base), intent(inout), dimension(:) :: particles !< Particle to initialize
type(type_coronal), intent(in)     :: cor !< Coronal equilibrium datatype for this particle
real*8, dimension(4), intent(in)   :: ran4 !< Four uniform random numbers (or sobol subset)
integer, intent(out) :: ifail

integer :: i_var(4), i
real*8, dimension(4) :: P, P_s, P_t, P_phi
real*8 :: v_out(4)
real*8 :: R, R_s, R_t, Z, Z_s, Z_t
real*8 :: background_density, background_kbT, background_kelvin, V_thermal
real*8 :: v_norm
real*8 :: Z_coronal, radiation_coronal
real*8 :: mass_main_ion

select type (particles)
type is (particle_kinetic_leapfrog)
  do i=1,size(particles)
    i_var = (/ 1, 5, 6, 7 /)
    call interp_PRZ(node_list,element_list,particles(i)%i_elm,i_var,size(i_var),particles(i)%st(1),particles(i)%st(2),particles(i)%x(3),&
        P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)

    background_density = P(2) * 1d20                                    ! plasma density [1/m^3]
    background_kbT     = P(3) /(MU_ZERO*central_density*1.d20)          ! Total plasma temperature in J/kB = T = Te + Ti
    background_kelvin  = background_kbT / K_BOLTZ / 2.d0                ! electron temperature [K]
    ! This is not valid for model400 (for that, remove the factor 2 above and below)

    V_thermal = sqrt(background_kbT / (2.d0*particles(i)%m*ATOMIC_MASS_UNIT))      ! [m/s]

    v_out = boxmueller_transform(ran4) * V_thermal ! [m/s], vx, vy, vz, dummy

    particles(i)%v(1) =   v_out(1) * cos(particles(i)%x(3)) + v_out(2) * sin(particles(i)%x(3))   ! V_R
    particles(i)%v(3) = - v_out(1) * sin(particles(i)%x(3)) + v_out(2) * cos(particles(i)%x(3))   ! V_phi [physical component]
    particles(i)%v(2) =   v_out(3)                                    ! V_Z

    mass_main_ion        = mass_proton * central_mass

    v_norm = sqrt(mu_zero * mass_main_ion * central_density * 1.d20)      ! JOREK normalisation for velocity
    particles(i)%v = particles(i)%v * v_norm

    if (background_density .le. 0.d0 .or. background_kelvin .le. 0.d0) then
      ifail = 1
    else
      ifail = 0
      call interpolate_coronal(cor, log10(background_density),log10(background_kelvin),Z_coronal,radiation_coronal)
      particles(i)%q       = nint(Z_coronal,1)                     !< charge (initialised with the coronal equilibrium value)
    endif
  end do
end select

end subroutine set_charge_from_coronal_eq
end module mod_initialise_particles
