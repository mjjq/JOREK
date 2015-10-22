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
real*8                    :: x(3), v(3), x_prev(3), v_prev(3) ! (Previous) values of position and velocity
real*8                    :: P(2), P_s(2), P_t(2) ! Placeholder for evaluating variables and derivatives locally
real*8                    :: v_tmp(3), R_update, RPhi_update ! Temporary values for the coordinate system transformation
real*8                    :: qom, B02, psi_R, psi_Z, U_R, U_Z, U_phi, U, B_phi_factor
real*8                    :: psi, psi_s, psi_t, u_s, u_t, psi_prev, U_prev
real*8                    :: R, R_s, R_t, Z, Z_s, Z_t, st_jac
real*8                    :: R_step, Z_step, fE, fB, t_norm
real*8                    :: R_out, Z_out, s_out, t_out
integer                   :: i, j, k, m, i_elm, j_elm, i_var(2), n_done, ifail, ielm_out, n_lost
logical                   :: changed, lost, search
real*8                    :: t0, t1

real*8, allocatable       :: rp(:,:), zp(:,:), tp(:), wp(:,:), mp(:,:), pp(:,:)

if (present(toroidal_field_factor)) then
  B_phi_factor = toroidal_field_factor
  write(*,*) 'WARNING: disabling toroidal magnetic field for particle propagator'
else
  B_phi_factor = 1.d0
endif

allocate(rp(n_step,particle_list%n_particles),zp(n_step,particle_list%n_particles),pp(n_step,particle_list%n_particles))
allocate(Wp(n_step,particle_list%n_particles),tp(n_step),mp(n_step,particle_list%n_particles))

call cpu_time(t0)

! Select the first and second physics variables
i_var = (/1, 2/)

t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20)    ! jorek time normalisation

!$omp parallel default(none) &
!$omp   shared(particle_list, node_list, element_list, t_step,n_step, F0, t_norm, i_var, B_phi_factor, central_mass, mp, wp) &
!$omp   private(i, j, k, particle, x, v, i_elm, j_elm, psi, psi_s, psi_t, psi_R, psi_Z, psi_prev, U_prev, R, R_s, R_t, Z, Z_s, Z_t,st_jac,              &
!$omp           qom, B0, B02, E0, v_tmp, fE, fB, R_step, Z_step, changed, lost, search,            &
!$omp           U, U_s, U_t, U_R, U_Z, U_phi, x_prev, v_prev, R_update, RPhi_update, P, P_s, P_t, R_out ,Z_out, ielm_out, s_out, t_out, ifail)

!$omp do
do i = 1, particle_list%n_particles

  particle = particle_list%particle(i)

  ! Skip this particle if it left the domain
  if (particle%lost) cycle

  ! Get s,t,phi,velocity and element index of the particle
  x(1:2) = particle%st
  x(3)   = particle%x(3)
  v      = particle%v
  i_elm  = particle%i_elm

  ! Interpolate the fields to get psi and U at the current position
  call interp_PRZ(node_list,element_list,i_elm,i_var,2,x(1),x(2),x(3),P,P_s,P_t,R,R_s,R_t,Z,Z_s,Z_t)

  !write(*,*) 'initial velocity : ',v
  !write(*,*) 'initial position : ',R,Z,x(3)
  psi = P(1)
  U   = P(2)

  ! Perform n_step Boris method steps here
  do j = 1, n_step

    ! write(*,'(A,2i5,3f12.8,A,3f12.8)') 'STEP : ',j,i_elm,x,' RZPhi: ',R,Z,x(3)

    ! This is the jacobian of the transformation from RZ to st
    st_jac = R_s * Z_t - R_t * Z_s

    ! And these are the local derivatives to s and t
    psi_s = P_s(1); psi_t = P_t(1)
    u_s   = P_s(2); u_t   = P_t(2)

    psi_prev   = psi
    U_prev     = U

    ! Calculate the derivatives to R and Z
    psi_R    = (  psi_s * Z_t - psi_t * Z_s ) / st_jac
    psi_Z    = (- psi_s * R_t + psi_t * R_s ) / st_jac
    U_R      = (  u_s   * Z_t - u_t   * Z_s ) / st_jac
    U_Z      = (- u_s   * R_t + u_t   * R_s ) / st_jac
    ! And assume for now no electric field in the phi-direction
    U_phi    = 0.d0

    ! Calculate the normalized fraction q/m in JOREK units, including correction for E and B having semi-SI units
    qom = real(particle%q) * EL_CHG / (particle%mass * ATOMIC_MASS_UNIT) * t_norm

    ! Calculate the magnetic field multiplied by q/m Delta_t / 2 and its magnitude
    B0     = (/ + psi_Z, - psi_R, F0*B_phi_factor /) / R
    B02    = dot_product(B0,B0)

    ! Calculate the local electric field, obtained from E=-Grad (u F0), multiplied by q/m Delta_t/2
    E0     = (/ - F0 * U_R, - F0 * U_Z, - F0 * U_phi / R /)

    ! Calculate the geometric factor f = tan(q/m delta_t/2 |B|)/|B|
    ! The t_norm is to normalize for B being in SI units instead of time-nondimensionalized
    fB = tan(qom * t_step * 0.5d0 * sqrt(B02)) / sqrt(B02)
    !fB = qom * t_step * 0.5d0
    fE = qom * t_step * 0.5d0

    ! Calculate the electric field update (v^n-1/2 -> v-)
    v = v + fE * E0
    ! Calculate the rotation
    v = (v + 2.d0*fB/(1.d0+fB*fB*B02)*( cross_product(v,B0) &
      - fB * v * B02 &
      + fB * B0 * dot_product(v,B0)))
    ! Calculate the next electric field update
    v = v + fE * E0

    ! Calculate the correction step (x = change in R, y is change in phi)
    R_update = R + v(1) * t_step
    RPhi_update = v(3) * t_step

    ! Calculate the new R and Z
    R_step = sqrt(R_update**2 + RPhi_update**2)
    Z_step = Z + t_step * v(2)
    ! Calculate the new phi
    x(3) = x(3) + asin(RPhi_update / R_step)

    ! Adjust R and Phi velocities to the new reference frame (z stays the same)
    v_tmp = v
    v(1) =  R_update/R_step    * v_tmp(1) + RPhi_update/R_step * v_tmp(3)
    v(3) = -RPhi_update/R_step * v_tmp(1) + R_update/R_step    * v_tmp(3)

    !write(*,*) R_step, Z_step, x(3)


    ! Perform at most 3 newton iteration steps to find the new values of s and t in this element
    do k = 1, 3
      x_prev = x
      x(1) = x(1) + ( Z_t * (R_step-R) - R_t * (Z_step-Z)) / st_jac
      x(2) = x(2) + (-Z_s * (R_step-R) + R_s * (Z_step-Z)) / st_jac

      ! Get the local derivatives at x
      call interp3_RZ(node_list,element_list,i_elm,x(1),x(2),R,R_s,R_t,Z,Z_s,Z_t)

      ! Update local jacobian
      st_jac = R_s * Z_t - R_t * Z_s

      !write(*,'(A,2i3,4f12.8,2e12.4)') 'first : ',j,k,x(1:2),R_step,Z_step,sqrt((R_step - R)**2 + (Z_step - Z)**2)

      if ((R_step - R)**2 + (Z_step - Z)**2 < 1.d-16) then
        exit ! Converged enough, we're done
      else if (k == 3) then
        write (*,*) "WARNING: Newton iteration failed, change = ", sqrt((R_step-R)**2 + (Z_step-Z)**2)
      endif
    enddo

    ! Check if the particle left this element in any of the coordinates. Keep going until we find it
    do m = 1, 4

      ! See if the particle is lost, in another element nearby or we need to search for it
      call check_element_boundary(element_list,i_elm,x(1:2),j_elm,x(1:2),changed,lost,search)

      !write(*,*) 'changed,lost,search : ',m,i_elm,changed,lost,search

      ! If we cannot find the particle
      if (search) then
        ! Use the nuclear option
        call find_RZ(node_list,element_list,R_step,Z_step,R_out,Z_out,ielm_out,s_out,t_out,ifail)
        x(1)  = s_out
        x(2)  = t_out
        i_elm = ielm_out
        exit ! Stop this loop, we're done
      endif

      if (lost) then
        particle%lost = .true.
        exit
      endif

      if (changed) then
        ! Perform newton iteration to find the correct position in the new element i_elm
        do k = 1, 3
          x_prev = x
          call interp3_RZ(node_list,element_list,i_elm,x(1),x(2),R,R_s,R_t,Z,Z_s,Z_t)
          st_jac = R_s * Z_t - R_t * Z_s
          x(1) = x(1) + ( Z_t * (R_step-R) - R_t * (Z_step-Z)) / st_jac
          x(2) = x(2) + (-Z_s * (R_step-R) + R_s * (Z_step-Z)) / st_jac

          !write(*,'(A,2i3,24f12.8,2e12.4)') 'second : ',j,k,x(1:2),R_step,Z_step,sqrt((R_step - R)**2 + (Z_step - Z)**2)

          if ((R_step - R)**2 + (Z_step - Z)**2 < 1.d-16) then
            exit ! Converged enough, we're done
          else if (k == 3) then
            write (*,*) "Newton iteration after element change failed, change = ", sqrt((R_step-R)**2 + (Z_step-Z)**2)
          endif
        enddo
      else
        exit
      endif

    enddo

    if (maxval(abs(x(1:2)-0.5)) .gt. 0.5) then
      write(*,*) "Error finding the right element"
      write(*,*) "Search for particle again, first routine did not work"
      call find_RZ(node_list,element_list,R_step,Z_step,R_out,Z_out,ielm_out,s_out,t_out,ifail)
      x(1)  = s_out
      x(2)  = t_out
      i_elm = ielm_out
    endif

    ! Update variables for the next iteration
    call interp_PRZ(node_list,element_list,i_elm,i_var,2,x(1),x(2),x(3),P,P_s,P_t,R,R_s,R_t,Z,Z_s,Z_t)

    ! Debug output
    !write(*,*) real(j)*t_step, R_step, Z_step, x(3), i_elm

    R = R_step
    Z = Z_step

    psi = P(1)
    U   = P(2)

    mp(j,i) = R_step * v(3) + 0.5 * qom * (psi + psi_prev)
    Wp(j,i) = 0.5 * particle%mass * (v(1)*v(1) + v(2)*v(2) + v(3)*v(3)) + 0.5 * EL_CHG / MASS_PROTON * t_norm * particle%q * F0 * (U + U_prev)

  enddo

  ! Save the new values for this particle
  particle_list%particle(i)%st    = x(1:2)
  particle_list%particle(i)%v     = v
  particle_list%particle(i)%x     = (/ R, Z, x(3) /)
  particle_list%particle(i)%i_elm = i_elm
  particle_list%particle(i)%lost  = particle%lost

enddo
!$omp end do
!$omp end parallel

call cpu_time(t1)

write(*,'(i5,A,f12.4)') my_id, ' Elapsed time particle update :',t1-t0

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
