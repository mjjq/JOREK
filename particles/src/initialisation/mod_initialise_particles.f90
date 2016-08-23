!> Module for the initialization of particles in configuration space (6D)
!> by rejection sampling
module mod_initialise_particles
use mod_coronal
type particle_init_params
  integer*1 :: species = 0     !< Atomic number Z of the particles (-1) for electrons
  real*4    :: atomic_mass     !< Atomic mass in a.m.u.

  !> TODO make the below object-oriented
  character(len=80) :: location_accept_function = 'location_accept_any' !< Which function to use for particle position rejection sampling
  real*4  :: location_accept_parameters(1:9) = 0 !< Extra arguments for this function
  integer :: particle_seed = 0 !< Seed for PCG random sequence used for particle init
  character(len=6) :: particle_initializer = 'pcg32' !< Method to use for seeding particles (options: pcg32, sobol)
end type particle_init_params
contains

!> Initialize_particles creates n_particles*(species>0), divided over all processors
subroutine initialise_particles(my_id,n_cpu,coronal,particle_list,particle_list_GC)

!$ use omp_lib
use constants
use data_structure
use mod_particles
use nodes_elements
use mod_sampling
use mod_rng
use mod_pcg32_rng
use mod_sobseq_rng
use_mpi

implicit none

! Routine parameters
type (type_coronal), intent(in) :: coronal(1:n_species) !< Coronal equilibrium data for all species
type (type_particle_list), intent(out) :: particle_list !< Output particle list
type (type_particle_list), intent(out) :: particle_list_GC !< Output guiding centre particle list
integer, intent(in) :: my_id, n_cpu !< MPI ID of this proc (0..n_cpi-1)

! Internal variables
integer, parameter :: N_d = 8 ! number of dimensions of the RNG
integer :: n_p(N_species) ! Number of particles on this cpu
type (type_particle)      :: particle
real*8  :: R, Z, phi
integer :: i, j, ifail
real*8 :: ran(N_d), t0, t1, ostart, oend
real*8 :: Rbox(2), Zbox(2), Rmin, Rmax, Zmin, Zmax
integer :: seq, n_streams, n_threads, i_thread
class(type_rng), allocatable, dimension(:) :: rngs ! The RNGs for all the threads

if (my_id .eq. 0) then
  write(*,*) '**********************************'
  write(*,*) '*      initialising particles    *'
  write(*,*) '**********************************'
endif

! Divide the particles over all processors, give the first processor the remainder
n_p = 0
do i=1,N_species
  if (species(i) .le. 0) cycle
  n_p(i) = n_particles(i)/n_cpu
  if (my_id .eq. 0) then ! Add the remainder to the first cpu
    n_p(i) = n_p(i) + modulo(n_particles(i), n_cpu)
  endif
enddo

! Allocate space for normal and GC particles
particle_list%n_particles = sum(n_p, .not. particle_GC)
allocate(particle_list%particle(particle_list%n_particles), stat=ifail)
if (ifail .gt. 0) write(*,"(i3,a,i8,a)") my_id, "unable to allocate particle_list%particle(", particle_list%n_particles, ")"
particle_list_GC%n_particles = sum(n_p,particle_GC)
allocate(particle_list_GC%particle(particle_list_GC%n_particles), stat=ifail)
if (ifail .gt. 0) write(*,"(i3,a,i8,a)") my_id, "unable to allocate particle_list_GC%particle(", particle_list_GC%n_particles, ")"


call cpu_time(t0)
!$ ostart = omp_get_wtime()

Rbox = (/1d9,-1d9/) ! Large and inversed bounds so that they are always updated
Zbox = (/1d9,-1d9/)
! Get domain size: min and max R and Z (phi is always between 0 and TWOPI)
do i=1,element_list%n_elements
  call RZ_minmax(node_list, element_list, i, Rmin, Rmax, Zmin, Zmax)
  Rbox(1) = min(Rbox(1), Rmin)
  Rbox(2) = max(Rbox(2), Rmax)
  Zbox(1) = min(Zbox(1), Zmin)
  Zbox(2) = max(Zbox(2), Zmax)
enddo

do i=1,N_species
  if (n_p(i) .le. 0) cycle
  if (atomic_mass(i) .eq. 0) then
    write(*,*) "Error: atomic mass of species ", i, " is zero! exiting."
    call MPI_ABORT(MPI_COMM_WORLD, -1, ifail)
    call exit(1)
  endif

  ! Setup (Q)RNGs, one per thread
  n_threads = 1
  !$ n_threads = omp_get_max_threads()
  select case (particle_initializer(i))
    case ("sobol")
      allocate(sobseq_rng::rngs(0:n_threads-1)) ! C indexing
    case ("pcg32")
      allocate(pcg32_rng::rngs(0:n_threads-1))
    case default
      write(*,*) "ERROR: Unrecognized value for particle_initializer(", i, "), exiting"
      call MPI_ABORT(MPI_COMM_WORLD, -1, ifail)
  end select

  !$omp parallel default(none) &
  !$omp   shared(particle_list, particle_list_GC, node_list, element_list, &
  !$omp          species, atomic_mass, Rbox, Zbox, particle_GC, i, n_p, coronal, &
  !$omp          my_id, particle_seed, particle_initializer, n_cpu, rngs, n_threads) &
  !$omp   private(j, R, Z, phi, ifail, particle, seq, n_streams, ran, i_thread)

  ! Seed RNGs
  n_streams = n_cpu*n_threads
  i_thread = 0
  !$ i_thread = omp_get_thread_num()
  seq = my_id*n_threads + i_thread + 1
  call rngs(i_thread)%initialize(n_dims=N_d, seed=particle_seed(i), n_streams=n_streams, i_stream=seq, ifail=ifail)
  if (ifail .ne. 0) then
    write(*,*) "Error seeding rng: ", ifail
    call MPI_ABORT(MPI_COMM_WORLD, -1, ifail)
  endif

  !$omp do
  do j=1,n_p(i)
    ifail = 1
    do while (ifail .ne. 0)
      ! Generate a random position to put this particle
      call rngs(i_thread)%next(ran)
      call transform_uniform_cylindrical(ran(1:3), Rbox, Zbox, (/0.d0, TWOPI/), R, Z, Phi)

      ! Deny or accept this location
      if (ran(4) .lt. accept_location(i, R, Z, phi)) then
        call particle_init_default(i, R, Z, phi, particle, ifail)
        if (ifail .eq. 0) then
          call particle_init(particle, coronal(i), ran(5:8), ifail)
        endif
      else
        ifail = 1
      endif
    enddo
    ! Save the result
    if (particle_GC(i)) then
      particle_list_GC%particle(j) = particle
    else
      particle_list%particle(j) = particle
    endif
  enddo
  !$omp end do
  !$omp end parallel
enddo

call cpu_time(t1)
!$ oend = omp_get_wtime()
write(*,'(i5,A,2f12.4)') my_id, ' Time particle initialize cpu/wall :',t1-t0, oend-ostart
if (my_id .eq. 0) then
  write(*,*) '* done initialising particles    *'
  write(*,*) '**********************************'
endif

end subroutine initialise_particles


!> Returns the probability of accepting a location with location_accept_function(i)
function accept_location(i,R,Z,phi)
use mod_particles, only: location_accept_function, location_accept_parameters
implicit none

integer, intent(in) :: i !< Species number
real*8,  intent(in) :: R,Z,phi !< Location coordinates
real*8 :: accept_location !< Acceptance probability

select case (location_accept_function(i))
  case ('joined_gaussian')
    accept_location = joined_gaussian(R,Z,phi,location_accept_parameters(:,i))
  case ('location_accept_any')
    accept_location = 1.d0
  case DEFAULT
    write(*,*) "WARNING: invalid value for location_accept_function: ", location_accept_function(i)
    accept_location = -1.d0
    call exit(1)
end select
end function accept_location


!> A joined gaussian function for accepting or rejecting locations
!! See http://jorek.eu/wiki/doku.php?id=particle_init#joined_gaussian
function joined_gaussian(R,Z,phi,p)
use nodes_elements

real*8, intent(in) :: R,Z,phi !< Coordinates
real*4, intent(in) :: p(1:9) !< Location_accept_parameters
real*8 :: joined_gaussian !< Acceptance probability

real*8 :: R_out, Z_out, s, t
integer :: i_elm, ifail, variable(1)
real*8 :: R_s, R_t, Z_s, Z_t ! Dummy variables
real*8, dimension(1) :: P_s, P_t, P_phi, Px
real*8 :: x

! Find our variable
select case(int(p(1)))
  case (-1)
    x = R
  case (-2)
    x = Z
  case (-3)
    x = phi
  case DEFAULT
    call find_RZ(node_list,element_list,R,Z,R_out,Z_out,i_elm,s,t,ifail)
    ! Select the mhd variable
    variable = (/int(p(1))/) ! Because the intel compiler is a nazi
    if (variable(1) .le. 0 .or. variable(1) .gt. n_var) then
      write(*,*) "ERROR: invalid MHD variable selected: ", variable
      call exit(1)
    endif
    if (ifail .eq. 0) then
      call interp_PRZ(node_list,element_list,i_elm,variable,1,s,t,phi,Px,P_s,P_t,P_phi,R_out,R_s,R_t,Z_out,Z_s,Z_t)
      x = Px(1)
    else
      ! Outside of domain
      joined_gaussian = 0.d0
      return
    endif
end select

! Calculate which section of the domain we are on
if (x .gt. p(2) .and. x .lt. p(4)) then ! In the joined part
  joined_gaussian = 1.d0
else ! In one of the two gaussians or between them if p(2) .gt. p(4)
  if (x .gt. (p(3)*p(4)+p(2)*p(5))/(p(3)+p(5))) then ! if x larger than the midpoint
    ! Use the one which is larger
    if (p(2) .ge. p(4)) then
      joined_gaussian = exp(-(x-p(2))**2/p(3)**2)
    else
      joined_gaussian = exp(-(x-p(4))**2/p(5)**2)
    endif
  else
    ! else the smaller one
    if (p(2) .lt. p(4)) then
      joined_gaussian = exp(-(x-p(2))**2/p(3)**2)
    else
      joined_gaussian = exp(-(x-p(4))**2/p(5)**2)
    endif
  endif
endif
end function joined_gaussian

!> Initialize a single particle at R_in, Z_in, phi_in
!> Does not initialize v and q, but only x, st, i_elm, mass, label
subroutine particle_init_default(i, R_in, Z_in, phi_in, particle, ifail)
use constants
use data_structure
use mod_particles
use nodes_elements
implicit none

real*8, intent(in)   :: R_in, Z_in, phi_in !< Particle coordinates
integer, intent(in)  :: i !< Particle species
type(type_particle), intent(out) :: particle !< Output particle
integer, intent(out) :: ifail

real*8 :: R_out, Z_out, s, t
integer :: i_elm

! Find s and t at this position and ielm
call find_RZ(node_list,element_list,R_in,Z_in,R_out,Z_out,i_elm,s,t,ifail)

! NB: this routine does not initialise v and q!
if (ifail .eq. 0) then
  particle%x       = (/ R_out, Z_out, phi_in /)    !< particle position in real space
  particle%st      = (/ s, t /)                    !< particle position in the finite element (i_elm)
  particle%i_elm   = i_elm                         !< element containing this particle
  particle%mass    = atomic_mass(i)                !< mass in AMU
  particle%label   = int(i,1)                      !< Type of particle
  particle%weight  = 1       !< weight (i.e. number of particles) TODO calculate and set number density here
  particle%lost    = .false.                       !< active particle
endif
end subroutine particle_init_default


!> Set v and q of a particle for use with kinetic codes
pure subroutine particle_init(particle, cor, ran4, ifail)
use constants
use data_structure
use mod_particles
use nodes_elements
use phys_module, only : central_density, central_mass
use mod_openadas
use mod_sampling
implicit none

interface
  pure subroutine interp_PRZ(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
    import :: type_node_list, type_element_list
    type (type_node_list),    intent(in)  :: node_list
    type (type_element_list), intent(in)  :: element_list
    integer,                  intent(in)  :: i_elm
    integer,                  intent(in)  :: n_v, i_v(n_v)
    real*8,                   intent(in)  :: s, t, phi
    real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v)
    real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
    real*8,                   intent(out) :: P_phi(n_v)
  end subroutine interp_PRZ
end interface

type(type_particle), intent(inout) :: particle !< Particle to initialize
type(type_coronal), intent(in)     :: cor !< Coronal equilibrium datatype for this particle
real*8, dimension(4), intent(in)   :: ran4 !< Four uniform random numbers (or sobol subset)
integer, intent(out) :: ifail

integer :: i_var(4)
real*8, dimension(4) :: P, P_s, P_t, P_phi
real*8 :: v_out(4)
real*8 :: R, R_s, R_t, Z, Z_s, Z_t
real*8 :: background_density, background_kbT, background_kelvin, V_thermal
real*8 :: v_norm
real*8 :: Z_coronal, radiation_coronal
real*8 :: mass_main_ion

i_var = (/ 1, 5, 6, 7 /)
call interp_PRZ(node_list,element_list,particle%i_elm,i_var,size(i_var),particle%st(1),particle%st(2),particle%x(3),&
    P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)

background_density = P(2) * 1d20                                    ! plasma density [1/m^3]
background_kbT     = P(3) /(MU_ZERO*central_density*1.d20)          ! Total plasma temperature in J/kB = T = Te + Ti
background_kelvin  = background_kbT / K_BOLTZ / 2.d0                ! electron temperature [K]
! This is not valid for model400 (for that, remove the factor 2 above and below)

V_thermal = sqrt(background_kbT / (2.d0*particle%mass*ATOMIC_MASS_UNIT))      ! [m/s]

v_out = boxmueller_transform(ran4) * V_thermal ! [m/s], vx, vy, vz, dummy

particle%v(1) =   v_out(1) * cos(particle%x(3)) + v_out(2) * sin(particle%x(3))   ! V_R
particle%v(3) = - v_out(1) * sin(particle%x(3)) + v_out(2) * cos(particle%x(3))   ! V_phi [physical component]
particle%v(2) =   v_out(3)                                    ! V_Z

mass_main_ion        = mass_proton * central_mass

v_norm = sqrt(mu_zero * mass_main_ion * central_density * 1.d20)      ! JOREK normalisation for velocity
particle%v = particle%v * v_norm

if (background_density .le. 0.d0 .or. background_kelvin .le. 0.d0) then
  ifail = 1
else
  ifail = 0
  call interpolate_coronal(cor, log10(background_density),log10(background_kelvin),Z_coronal,radiation_coronal)
  particle%q       = nint(Z_coronal,1)                     !< charge (initialised with the coronal equilibrium value)
endif

end subroutine particle_init
end module mod_initialise_particles
