module mod_update_particles
contains
!> Update_particles performs n_step steps of the Boris method with size t_step
!!
!! The Boris method is adjusted for cylindrical coordinates.
!! Some care must be taken when element boundaries are crossed.
!! It is parallelized with OMP.
!! See G.L. Delzanno, E. Camporeale / JCP 253 (2013) 259-277 for details
subroutine update_particles(my_id, i_species, particle_list, t_step, n_step, adf11, &
        energy_list, momentum_list, toroidal_field_factor, field_interp_time)

!$ use omp_lib
use parameters
use data_structure
use nodes_elements
use constants
use phys_module, only : F0, central_mass, central_density
use mod_particles
use openadas
use mod_ionisation_recombination
use mod_particle_diagnostics

implicit none

! -- Routine parameters
integer, intent(in)       :: my_id              !< Id of the current process
integer, intent(in)       :: i_species          !< The number of the current particle species
type (type_particle_list), intent(inout) :: particle_list      !< The particles we will march forward in time
real*8,  intent(in)       :: t_step             !< The size of each timestep
integer, intent(in)       :: n_step             !< The number of timesteps we will perform
type (type_adf11_all), intent(in) :: adf11
real*8,  intent(out), dimension(:), optional :: energy_list !< Energy of the particles at the next-to(!) final timestep
real*8,  intent(out), dimension(:), optional :: momentum_list !< Generalized toroidal momentum of the particles at the next-to(!) final timestep
real*8,  intent(in), optional :: toroidal_field_factor !< Multiply B_phi with this WARNING: use only for testing!
logical, intent(in),  optional :: field_interp_time !< Interpolate the fields linearly in time as if the first step was in the previous fields (almost) and the last in the current

! -- Local variables
type (type_particle)      :: particle
real*8                    :: B0(3), E0(3) ! Local B and E field at particle position
real*8                    :: x(3), st(2), v(3), x_prev(3), v_prev(3) ! (Previous) values of position and velocity
real*8                    :: v_tmp(3), R_update, RPhi_update ! Temporary values for the coordinate system transformation
real*8                    :: qom, B02, B_phi_factor, q, m, eom
real*8                    :: R, Z, psi, psid, U, Ud, energy_lost_particles, R_inv
real*8                    :: fE, fB, t_norm
real*8                    :: R_out, Z_out, s_out, t_out
integer                   :: i, j, i_elm, ifail, ielm_out
logical                   :: do_substep

! -- Output variables
integer                   :: n_lost, pir(DOMAIN_PLASMA:DOMAIN_LOWER_PRIVATE)
real*8                    :: t0, t1, ostart, oend, delta_fraction
integer                   :: find_RZ_count
real*8                    :: stats(5), total_lost
real*8, save              :: total_energy_lost_particles = 0.d0

interface
  subroutine calc_EB(i_elm,st,phi,E,B,psi,U,delta_fraction)
    integer, intent(in) :: i_elm
    real*8, intent(in)  :: st(2), phi
    real*8, intent(in), optional :: delta_fraction

    real*8, intent(out) :: E(3), B(3), psi, U
  end subroutine calc_EB
end interface

if (present(toroidal_field_factor)) then
  B_phi_factor = toroidal_field_factor
  write(*,*) 'INFO: setting toroidal magnetic field factor for particle propagator testcase to', B_phi_factor
else
  B_phi_factor = 1.d0
endif

do_substep = .false. ! Might need better name
if (present(field_interp_time)) then
  if (field_interp_time) do_substep = .true.
endif

if (present(energy_list))   energy_list   = 0.d0
if (present(momentum_list)) momentum_list = 0.d0
energy_lost_particles = 0.d0
n_lost = 0

find_RZ_count = 0
call cpu_time(t0)
!$ t0 = omp_get_wtime()
if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '* JOREK2 : update particles           *'
  write(*,*) '***************************************'
  write(*,*) 'integrating trajectories for ', n_step, ' steps'
endif
t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20)    ! jorek time normalisation

!$omp parallel default(none) &
!$omp   shared(particle_list, node_list, element_list, t_step,n_step, F0, t_norm, B_phi_factor, central_mass, energy_list, &
!$omp   momentum_list, do_substep, adf11, particle_ion_rec, i_species) &
!$omp   reduction(+:find_RZ_count,energy_lost_particles,n_lost) &
!$omp   private(i, j, particle, x, v, st, i_elm, R, Z, &
!$omp           qom, B0, B02, E0, v_tmp, fE, fB,             &
!$omp           x_prev, v_prev, R_update, RPhi_update, R_out ,Z_out, ielm_out, s_out, &
!$omp           t_out, ifail, psi, U, psid, Ud, q, m, eom, R_inv, delta_fraction)

!$omp do
do i = 1, particle_list%n_particles

  particle = particle_list%particle(i)

  ! Skip this particle if it left the domain
  if (particle%lost) cycle

  ! Get s,t,phi,velocity and element index of the particle
  x     = particle%x
  st    = particle%st
  v     = particle%v
  i_elm = particle%i_elm

  ! TODO better initialization here
  psi = 0
  U = 0
  eom = EL_CHG / (particle%mass * ATOMIC_MASS_UNIT) * t_norm

  do j = 1, n_step
    R = x(1)
    Z = x(2)

    ! q/m in JOREK units, including correction for E and B having semi-SI units
    qom = particle%q * eom

    if (do_substep) then
      delta_fraction = 1.d0 - real(j)/real(n_step)
      call calc_EB(i_elm,st,x(3),E0,B0,psi,U,delta_fraction)
    else
      call calc_EB(i_elm,st,x(3),E0,B0,psi,U)
    endif
    B0(3) = B0(3)*B_phi_factor
    B02   = dot_product(B0,B0)

    ! Calculate the geometric factor f = tan(q/m delta_t/2 |B|)/|B|
    fB = tan(qom * t_step * 0.5d0 * sqrt(B02)) / sqrt(B02)
    fE = qom * t_step * 0.5d0

    ! Calculate the electric field update (v^n-1/2 -> v-) with the Boris method
    v = v + fE * E0
    ! Calculate the rotation
    v = (v + 2.d0*fB/(1.d0+fB*fB*B02)*( cross_product(v,B0) &
      - fB * v * B02 &
      + fB * B0 * dot_product(v,B0)))
    ! Calculate the next electric field update (v+ -> v^n+1/2)
    v = v + fE * E0

    R_update = R + v(1) * t_step
    RPhi_update = v(3) * t_step

    ! Calculate the new R and Z
    x(1) = sqrt(R_update**2 + RPhi_update**2)
    x(2) = Z + t_step * v(2)
    ! Calculate the new phi
    x(3) = x(3) + asin(RPhi_update / x(1))
    R_inv = 1.d0/x(1)

    ! Adjust R and Phi velocities to the new reference frame (z stays the same)
    v_tmp = v
    v(1)  = (R_update     * v_tmp(1) + RPhi_update * v_tmp(3))*R_inv
    v(3)  = (-RPhi_update * v_tmp(1) + R_update    * v_tmp(3))*R_inv

    ! Find new st coordinates and new element
    call find_RZ_nearby(node_list,element_list,x(1:2),(/R,Z/),st,st,i_elm,i_elm,ifail)
    ! If ifail in 2..5 we used find_RZ
    if (ifail .ge. 2 .and. ifail .le. 5) then
      find_RZ_count = find_RZ_count + 1
    endif

    if (ifail .eq. -1 .or. i_elm .eq. 0) then
      particle%lost = .true.
      n_lost = n_lost + 1
      ! Save lost energy (using inaccurate velocity at t-dt/2, use variables
      ! because omp atomic does not work with particle%q)
      q = real(particle%q) * EL_CHG
      m = particle%mass * ATOMIC_MASS_UNIT
      energy_lost_particles = energy_lost_particles + 0.5 * m * dot_product(v,v) &
                            + q * F0 * U * t_norm
      exit
    endif

    ! Calculate ionisation and recombination odds to find new charge
    if (particle_ion_rec(i_species)) then
      call update_particle_charge(node_list, element_list, particle, adf11, t_step*t_norm)
    endif
  enddo

  ! Save the new values for this particle
  particle_list%particle(i)%st    = st
  particle_list%particle(i)%v     = v
  particle_list%particle(i)%x     = x
  particle_list%particle(i)%i_elm = i_elm
  particle_list%particle(i)%lost  = particle%lost

  if (.not. particle%lost .and. (present(momentum_list) .or. present(energy_list))) then
    ! Calculate the fields at the new position
    call calc_EB(i_elm,st,x(3),E0,B0,psi,U)

    ! Perform a final update (half-step) to get the velocity at this time
    ! (TODO split out code below into boris method function)
    B0(3) = B0(3)*B_phi_factor
    B02    = dot_product(B0,B0)

    ! Calculate the geometric factor f = tan(q/m delta_t/2 |B|)/|B| !half-step!
    fB = 0.5d0 * tan(qom * t_step * 0.5d0 * sqrt(B02)) / sqrt(B02)
    fE = 0.5d0 * qom * t_step * 0.5d0

    ! Calculate the electric field update (v^n-1/2 -> v-) with the Boris method
    v = v + fE * E0
    ! Calculate the rotation
    v = (v + 2.d0*fB/(1.d0+fB*fB*B02)*( cross_product(v,B0) &
      - fB * v * B02 &
      + fB * B0 * dot_product(v,B0)))
    ! Calculate the next electric field update (v+ -> v^n+1/2)
    v = v + fE * E0

    R_update = R + v(1) * t_step
    RPhi_update = v(3) * t_step

    ! Calculate the new R and Z
    x(1) = sqrt(R_update**2 + RPhi_update**2)
    x(2) = Z + t_step * v(2)
    ! Calculate the new phi
    x(3) = x(3) + asin(RPhi_update / x(1))

    ! Adjust R and Phi velocities to the new reference frame (z stays the same)
    v_tmp = v
    v(1) =  R_update/x(1)    * v_tmp(1) + RPhi_update/x(1) * v_tmp(3)
    v(3) = -RPhi_update/x(1) * v_tmp(1) + R_update/x(1)    * v_tmp(3)

    if (present(energy_list)) then
      energy_list(i) = 0.5 * particle_list%particle(i)%mass * ATOMIC_MASS_UNIT * dot_product(v,v) &
                     + EL_CHG * particle%q * F0 * U * t_norm
    endif
    if (present(momentum_list)) then
      momentum_list(i) = particle_list%particle(i)%x(1) * v(3) + qom * psi
    endif
  endif

enddo
!$omp end do
!$omp end parallel

call cpu_time(t1)
!$ t1 = omp_get_wtime()

!! Output timing and diagnostics values
!! min, mean, max, stddev
stats = mpi_stats(t1-t0)
if (my_id .eq. 0) write(*,'(A,4f9.5)') 'Time particle update (min/mean/max/stddev/[total]): ', stats(1:4)
stats = mpi_stats(real(find_RZ_count,8)*100.d0/real(particle_list%n_particles*n_step,8))
if (my_id .eq. 0) write(*,'(A,4f9.5,A)') '   Find_RZ used in ', stats(1:4), ' % of the runs'
stats = mpi_stats(real(n_lost,8))
if (my_id .eq. 0) write(*,'(A,5f9.5)')   '   number of lost particles: ', stats
stats = mpi_stats(energy_lost_particles)
if (my_id .eq. 0) write(*,'(A,5g12.4)')  '   particle energy left domain: ', stats
pir = particles_in_regions(node_list, element_list, particle_list)
if (my_id .eq. 0) write(*,'(A,5i8)') '   particle locations (plasma, sol, out, up, low): ', pir


total_energy_lost_particles = total_energy_lost_particles + energy_lost_particles
stats = mpi_stats(total_energy_lost_particles)
if (.not. present(energy_list)) then
  if (my_id .eq. 0) write(*,'(A,5g18.10)') '   lost particle energy: ', stats
else
  total_lost = stats(5)
  call statistics_no_zero_MPI(energy_list, stats, n_lost)
  if (my_id .eq. 0) write(*,'(A,6g12.4)') '   energy min/mean/max/stddev/total/total_lost :',&
    stats, total_lost 
endif
if (present(momentum_list)) then
  call statistics_no_zero_MPI(momentum_list, stats, n_lost)
  if (my_id .eq. 0) write(*,'(A,6g12.4)') '   momentum min/mean/max/stddev/total/n_lost :',&
    stats,n_lost
endif
end subroutine update_particles


!> Calculate statistics over MPI
!! Return stats(1) = min, stats(2) = mean, stats(3) = max, stats(4) = stddev
function mpi_stats(var)
  use mpi_mod
  implicit none
  real*8, intent(in) :: var
  real*8, dimension(5) :: mpi_stats
  integer :: n_cpu, ierr
  real*8, allocatable, dimension(:) :: vars

  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)
  allocate(vars(n_cpu))
  vars = 0.d0

  call MPI_Gather(var, 1, MPI_REAL8, vars, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  mpi_stats(1) = minval(vars)
  mpi_stats(2) = sum(vars)/size(vars, 1)
  mpi_stats(3) = maxval(vars)
  mpi_stats(4) = sqrt(sum((vars-mpi_stats(2))**2)/size(vars,1))
  mpi_stats(5) = sum(vars)
end function mpi_stats


!> Calculate statistics on sets of numbers ignoring zeros
!! Performs MPI comunication for the mean
!! Only returns usable values on the root node
subroutine statistics_no_zero_MPI(list,stats,num_zeros)
  use mpi_mod
  implicit none
  real*8, intent(in), dimension(:) :: list
  real*8, intent(out), dimension(5) :: stats !=(/minv, mean, maxv, sd, total/)
  integer, intent(out) :: num_zeros
  integer :: num_values, ierr
  logical, dimension(:), allocatable :: mask

  allocate(mask(size(list)))
  ! Everything > 0 is true and will be used
  mask = abs(list) > 0.d0
  num_values = count(mask)
  num_zeros = size(list) - num_values

  stats = 0.d0
  call MPI_Reduce(minval(list,mask), stats(1), 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(maxval(list,mask), stats(3), 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sum(list,mask),    stats(5), 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(count(mask), num_values, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  stats(2) = stats(5)/real(num_values,8)
  call MPI_Bcast(stats(2), 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sum((list-stats(2))**2,mask), stats(4), 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  stats(4) = sqrt(stats(4)/real(num_values,8))
end subroutine statistics_no_zero_MPI
end module mod_update_particles
