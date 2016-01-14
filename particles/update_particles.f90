!> Update_particles performs n_step steps of the Boris method with size t_step
!!
!! The Boris method is adjusted for cylindrical coordinates.
!! Some care must be taken when element boundaries are crossed.
!! It is parallelized with OMP.
!! See G.L. Delzanno, E. Camporeale / JCP 253 (2013) 259-277 for details
subroutine update_particles(my_id, particle_list, t_step, n_step, toroidal_field_factor)

use parameters
use data_structure
use nodes_elements
use constants
use phys_module, only : F0, central_mass, central_density
use mod_particles

implicit none

! -- Routine parameters
type (type_particle_list) :: particle_list      !< The particles we will march forward in time
real*8,  intent(in)       :: t_step             !< The size of each timestep
integer, intent(in)       :: n_step             !< The number of timesteps we will perform
integer, intent(in)       :: my_id              !< Id of the current process
real*8,  intent(in), optional :: toroidal_field_factor !< Multiply B_phi with this WARNING: use only for testing!

! -- Local variables
type (type_particle)      :: particle
real*8                    :: B0(3), E0(3) ! Local B and E field at particle position
real*8                    :: x(3), st(2), v(3), x_prev(3), v_prev(3) ! (Previous) values of position and velocity
real*8                    :: v_tmp(3), R_update, RPhi_update ! Temporary values for the coordinate system transformation
real*8                    :: qom, B02, B_phi_factor
real*8                    :: R, Z, psi, psi_prev, U, U_prev, particle_energy
real*8                    :: fE, fB, t_norm
real*8                    :: R_out, Z_out, s_out, t_out
integer                   :: i, j, k, m, i_elm, n_done, ifail, ielm_out, n_lost
logical                   :: changed, lost, search
real*8                    :: t0, t1
integer                   :: find_RZ_count

real*8, allocatable       :: rp(:,:), zp(:,:), tp(:), wp(:,:), mp(:,:), pp(:,:), E_particles(:), E_particles_lost(:)


if (present(toroidal_field_factor)) then
  B_phi_factor = toroidal_field_factor
  write(*,*) 'INFO: setting toroidal magnetic field factor for particle propagator testcase to', B_phi_factor
else
  B_phi_factor = 1.d0
endif

!allocate(rp(n_step,particle_list%n_particles),zp(n_step,particle_list%n_particles),pp(n_step,particle_list%n_particles))
!allocate(Wp(n_step,particle_list%n_particles),tp(n_step),mp(n_step,particle_list%n_particles))
allocate(E_particles(particle_list%n_particles),E_particles_lost(particle_list%n_particles))
E_particles = 0.d0
E_particles_lost = 0.d0 ! XXX too much storage
find_RZ_count = 0

call cpu_time(t0)

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '* JOREK2 : update particles           *'
  write(*,*) '***************************************'
  write(*,*) 'integrating trajectories for ', n_step, ' steps'
endif

t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20)    ! jorek time normalisation

!$omp parallel default(none) &
!$omp   shared(particle_list, node_list, element_list, t_step,n_step, F0, t_norm, B_phi_factor, central_mass, mp, wp, find_RZ_count, E_particles, E_particles_lost) &
!$omp   private(i, j, k, particle, x, v, st, i_elm, R, Z, &
!$omp           qom, B0, B02, E0, v_tmp, fE, fB, changed, lost, search,            &
!$omp           x_prev, v_prev, R_update, RPhi_update, R_out ,Z_out, ielm_out, s_out, &
!$omp           t_out, ifail, psi, U, psi_prev, U_prev, particle_energy)

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

  do j = 1, n_step
    R = x(1)
    Z = x(2)

    ! q/m in JOREK units, including correction for E and B having semi-SI units
    qom = real(particle%q) * EL_CHG / (particle%mass * ATOMIC_MASS_UNIT) * t_norm

    psi_prev = psi
    U_prev = U
    call calc_EB(i_elm,st,x(3),E0,B0,psi,U)
    B0(3) = B0(3)*B_phi_factor
    B02    = dot_product(B0,B0)

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

    ! Adjust R and Phi velocities to the new reference frame (z stays the same)
    v_tmp = v
    v(1) =  R_update/x(1)    * v_tmp(1) + RPhi_update/x(1) * v_tmp(3)
    v(3) = -RPhi_update/x(1) * v_tmp(1) + R_update/x(1)    * v_tmp(3)

    ! Find new st coordinates and new element
    call find_RZ_nearby(node_list,element_list,x(1:2),(/R,Z/),st,st,i_elm,i_elm,ifail)
    ! If ifail in 2..5 we used find_RZ
    if (ifail .ge. 2 .and. ifail .le. 5) find_RZ_count = find_RZ_count + 1 ! XXX not threadsafe! works if find_RZ is not often used, which we want

    if (ifail .eq. -1 .or. i_elm .eq. 0) then
      particle%lost = .true.
      exit
    endif

    ! Save energy statistics
    !mp(j,i) = x(1) * v(3) + 0.5 * qom * (psi + psi_prev)
    !Wp(j,i) = 0.5 * particle%mass * (v(1)*v(1) + v(2)*v(2) + v(3)*v(3)) &
    !    + 0.5 * EL_CHG / MASS_PROTON * t_norm * particle%q * F0 * (U + U_prev)

  enddo

  ! Save the new values for this particle
  particle_list%particle(i)%st    = st
  particle_list%particle(i)%v     = v
  particle_list%particle(i)%x     = x
  particle_list%particle(i)%i_elm = i_elm
  particle_list%particle(i)%lost  = particle%lost

  ! Calculate the fields (or at least U) at the new position
  if (.not. particle%lost) then
    U = U_prev
    call calc_EB(i_elm,st,x(3),E0,B0,psi,U)
  else
    ! Do not calculate the fields at the new position as this is outside of the
    ! domain
  endif

  particle_energy = 0.5 * particle%mass * ATOMIC_MASS_UNIT * dot_product(v,v) &
                  + 0.5 * EL_CHG / MASS_PROTON * t_norm * particle%q * F0 * (U + U_prev)
  if (.not. particle%lost) then
    E_particles(i) = particle_energy
  else ! Only happens if the particle is lost in this loop
    E_particles_lost(i) = particle_energy
  endif

enddo
!$omp end do
!$omp end parallel

call cpu_time(t1)

write(*,'(i5,A,f12.4)') my_id, ' Elapsed time particle update :',t1-t0
write(*,'(i5,A,f7.3,A)') my_id, '   Find_RZ used in ', &
  real(find_RZ_count)*100.d0/real(particle_list%n_particles*n_step), ' % of the runs'
write(*,'(i5,A,g18.10)') my_id, ' Total active particle energy:',sum(E_particles)
write(*,'(i5,A,g18.10)') my_id, '  particle energy left domain:',sum(E_particles_lost)

deallocate(E_particles)


n_done = n_step
!write(*,'(A,4e16.8)') 'mean: ',sum(abs(Wp(1:n_done,1)-Wp(1,1)))/wp(1,1)/n_done,sum(abs(mp(1:n_done,1)-mp(1,1)))/mp(1,1)/n_done
!write(*,'(A,4e16.8)') 'max : ',maxval(abs(Wp(1:n_done,1)-Wp(1,1)))/wp(1,1),maxval(abs(mp(1:n_done,1)-mp(1,1)))/mp(1,1)

!call lplot6(21,1,rp(1:n_done,1),zp(1:n_done,1),n_done,' ')
!call lincol(0)
!call lplot6(1,2,tp(1:n_done),Wp(1:n_done,1)-Wp(1,1),n_done,'Energy')
!call lplot6(1,3,tp(1:n_done),mp(1:n_done,1)-mp(1,1),n_done,'moment')
!call lplot6(1,1,tp(1:n_done),pp(1:n_done,1),n_done,'toroidal angle')

n_lost = 0
do j=1, particle_list%n_particles
  if (particle_list%particle(j)%lost) n_lost = n_lost + 1
enddo
write(*,*) my_id,'lost particles : ',n_lost

!open(21,file='traject.txt')
!write(21,*) '          X               Y               Z                R              Wp             Mp         time'
!do i=1,n_done
!  write(21,'(8e16.8)') rp(i,1)*cos(pp(i,1)),zp(i,1),rp(i,1)*sin(pp(i,1)),rp(i,1),wp(i,1),mp(i,1),tp(i)
!enddo
!close(21)

return
end
