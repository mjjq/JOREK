module mod_initialise_particles
contains

!> Initialize_particles creates n_particles*(species>0), divided over all processors
subroutine initialise_particles(my_id,n_cpu,particle_list,particle_list_GC)

!$ use omp_lib
use constants
use data_structure
use mod_particles
use nodes_elements

implicit none

! Routine parameters
type (type_particle_list), intent(out) :: particle_list
type (type_particle_list), intent(out) :: particle_list_GC
integer, intent(in) :: my_id, n_cpu

! Internal variables
integer :: n_p(N_species) ! Number of particles on this cpu
type (type_particle)      :: particle
real*8  :: R, Z, phi
integer :: i, j, ifail
real*8 :: ran3(3), t0, t1, ostart, oend
real*8 :: Rbox(2), Zbox(2)

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

! Call find_RZ once to initialise elements_minmax before going OMP (parameters are random-ish and irrelevant)
call find_RZ(node_list,element_list,0.d0,0.d0,R,Z,i,ran3(1),ran3(2),ifail)

call cpu_time(t0)
!$ ostart = omp_get_wtime()
call random_seed()

! Get domain size: min and max R and Z (phi is always between 0 and TWOPI)
Rbox = (/0.77d0, 4.89d0/)
Zbox = (/-3.8d0, 3.8d0/)
! TODO get these from the elements instead of hardcoding

! Create particles of all kinds given in the input file
do i=1,N_species
  if (n_p(i) .le. 0) cycle
  if (atomic_mass(i) .eq. 0) then
    write(*,*) "Error: atomic mass of species ", i, " is zero! exiting."
    call exit(1)
  endif

  !$omp parallel do default(none) &
  !$omp   shared(particle_list, particle_list_GC, node_list, element_list, &
  !$omp          species, atomic_mass, Rbox, Zbox, particle_GC, i) &
  !$omp   private(j, R, Z, phi, ifail, particle, ran3)
  do j=1,n_p(i)
    ifail = 1
    do while (ifail .ne. 0)
      ! Generate a random position to put this particle
      call random_number(ran3)
      ! Use inversion sampling to correct for cylindrical coordinates
      ! r = sqrt(rand() (B^2 - A^2) + A^2) for min and max radius A and B
      R   = sqrt(ran3(1) * (Rbox(2)**2-Rbox(1)**2) + Rbox(1)**2)
      Z   = (Zbox(2)-Zbox(1))*(ran3(2)-0.5d0)
      phi = TWOPI*(ran3(3)-0.5d0)

      ! Deny or accept this location
      if (accept_location(i, R,Z,phi)) then
        call particle_init_default(i, R, Z, phi, particle, ifail)
        if (ifail .eq. 0) then
          call particle_init(particle, ifail)
        endif
      else
        ifail = 1
      endif
    enddo
    ! Save the result
    if (particle_GC(i)) then
      particle_list_GC%particle(i) = particle
    else
      particle_list%particle(i) = particle
    endif
  enddo
  !$omp end parallel do
enddo

call cpu_time(t1)
!$ oend = omp_get_wtime()
write(*,'(i5,A,2f12.4)') my_id, ' Time particle initialize cpu/wall :',t1-t0, oend-ostart
if (my_id .eq. 0) then
  write(*,*) '* done initialising particles    *'
  write(*,*) '**********************************'
endif

end subroutine initialise_particles

!> Returns whether or not a location is accepted according to location_accept_function(i)
function accept_location(i,R,Z,phi)
use mod_particles, only: location_accept_function, location_accept_parameters
implicit none

integer, intent(in) :: i !< Species number
real*8,  intent(in) :: R,Z,phi
logical :: accept_location

real*8 :: ran

call random_number(ran)
select case (location_accept_function(i))
  case ('joined_gaussian')
    accept_location = (ran .le. joined_gaussian(R,Z,phi,location_accept_parameters(:,i)))
  case ('location_accept_any')
    accept_location = .true.
  case DEFAULT
    write(*,*) "WARNING: invalid value for location_accept_function: ", location_accept_function(i)
    accept_location = .false.
    call exit(1)
end select
end function accept_location

!> A joined gaussian function for accepting or rejecting locations
!! See http://jorek.eu/wiki/doku.php?id=particle_init#joined_gaussian
function joined_gaussian(R,Z,phi,p)
use nodes_elements

real*8, intent(in) :: R,Z,phi
real*4, intent(in) :: p(1:9)
real*8 :: joined_gaussian

real*8 :: R_out, Z_out, s, t
integer :: i_elm, ifail
real*8 :: P_s, P_t, P_phi, R_s, R_t, Z_s, Z_t ! Useless variables
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
    if (ifail .eq. 0) then
      call interp_PRZ(node_list,element_list,i_elm,p(1),1,s,t,phi,x,P_s,P_t,P_phi,R_out,R_s,R_t,Z_out,Z_s,Z_t)
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
subroutine particle_init_default(i, R_in, Z_in, phi_in, particle, ifail)
use constants
use data_structure
use mod_particles
use nodes_elements
implicit none

real*8, intent(in)   :: R_in, Z_in, phi_in
integer, intent(in)  :: i 
type(type_particle), intent(out) :: particle
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
  particle%species = int(species(i),1)             !< Type of particle
  particle%weight  = 1       !< weight (i.e. number of particles) TODO calculate and set number density here
  particle%lost    = .false.                       !< active particle
endif
end subroutine particle_init_default

!> Set v and q of a particle for use with kinetic codes
subroutine particle_init(particle, ifail)
use constants
use data_structure
use mod_particles
use nodes_elements
use phys_module, only : central_density, central_mass
use openadas
implicit none

type(type_particle), intent(inout) :: particle
integer, intent(out) :: ifail

integer :: i_var(4), i_elm
real*8, dimension(4) :: P, P_s, P_t, P_phi, ran4
real*8 :: R, R_s, R_t, Z, Z_s, Z_t
real*8 :: background_density, background_kbT, background_kelvin, V_thermal
real*8 :: vx, vy, vz, v_norm
real*8 :: Z_coronal, radiation_coronal
real*8 :: mass_main_ion

i_var = (/ 1, 5, 6, 7 /)
call interp_PRZ(node_list,element_list,i_elm,i_var,size(i_var),particle%st(1),particle%st(2),particle%x(3),&
    P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)

background_density = P(2) * 1d20                             ! plasma density [1/m^3]
background_kbT     = P(3) /(MU_ZERO*central_density*1.d20)   ! plasma temperature [J]
background_kelvin  = background_kbT / EL_CHG                 ! plasma temperature [K]

V_thermal = sqrt(background_kbT / (atomic_mass(particle%species)*ATOMIC_MASS_UNIT))      ! [m/s]

call random_number(ran4)

vx = V_thermal * sqrt(-2*Log(ran4(1))) * cos(TWOPI*ran4(2))  ! [m/s]  uses box-mueller transform
vy = V_thermal * sqrt(-2*Log(ran4(1))) * sin(TWOPI*ran4(2))
vz = V_thermal * sqrt(-2*Log(ran4(3))) * sin(TWOPI*ran4(4))

! TODO Why 1d-6?
call interpolate_coronal(dlog10(background_density*1d-6),dlog10(background_kelvin),Z_coronal,radiation_coronal)

particle%v(1) =   vx * cos(particle%x(3)) + vy * sin(particle%x(3))   ! V_R
particle%v(3) = - vx * sin(particle%x(3)) + vy * cos(particle%x(3))   ! V_phi [physical component]
particle%v(2) =   vz                                    ! V_Z

mass_main_ion        = mass_proton * central_mass
v_norm = sqrt(mu_zero * mass_main_ion * central_density * 1.d20)      ! JOREK normalisation for velocity
particle%v = particle%v * v_norm

particle%q       = int(Z_coronal,1)                     !< charge (initialised with the coronal equilibrium value)

ifail = 0

end subroutine particle_init
end module mod_initialise_particles
