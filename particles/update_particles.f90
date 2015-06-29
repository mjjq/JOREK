!> Update_particles performs n_step steps of the Boris method with size t_step
!!
!! The Boris method is adjusted for cylindrical coordinates.
!! Some care must be taken when element boundaries are crossed.
!! It is parallelized with OMP.
subroutine update_particles(my_id,particle_list, t_step, n_step)

use parameters
use data_structure
use nodes_elements
use constants
use clock_module
use phys_module, only : F0, central_mass
use mod_particles

implicit none

! -- Routine parameters
type (type_particle_list) :: particle_list      !< The particles we will march forward in time
real*8,  intent(in)       :: t_step             !< The size of each timestep
integer, intent(in)       :: n_step             !< The number of timesteps we will perform
integer, intent(in)       :: my_id              !< Id of the current process

! -- Local variables
type (type_particle)      :: particle
real*8                    :: B(3), B0(3), E0(3), x_prev(3), v_prev(3), P(2), P_s(2), P_t(2)
real*8                    :: x(3), v(3), v0(3), v_up(3), delta_v(3), delta_phi, f
real*8                    :: qom, B02, psi_R, psi_Z, psi_prev, U_R, U_Z, U_phi, U, U_prev
real*8                    :: psi, psi_s, psi_t, u_s, u_t, tsecond, omega_norm
real*8                    :: R, R_s, R_t, Z, Z_s, Z_t, st_jac, delta_x(2)
real*8                    :: R_prev, Z_prev, R_step, Z_step, error_RZ, R_axis, Z_axis, psi_axis, CR, CZ, r_crit2
real*8                    :: R_out, Z_out, s_out, t_out, value_out
integer                   :: i, j, k, m, i_elm, j_elm, i_var(2), n_done, ifail, ielm_out, n_lost
logical                   :: changed, lost, search
real*8                    :: t0, t1

real*8, allocatable       :: rp(:,:), zp(:,:), tp(:), Wp(:,:), mp(:,:), pp(:,:)

!allocate(rp(n_step,particle_list%n_particles),zp(n_step,particle_list%n_particles),pp(n_step,particle_list%n_particles))
!allocate(Wp(n_step,particle_list%n_particles),tp(n_step),mp(n_step,particle_list%n_particles))

write(*,'(A)') '*********************************************'
write(*,'(A,i12,f12.8)') 'updating particles : ',n_step,t_step
write(*,'(A)') '*********************************************'
write(*,'(i5,A,i12)') my_id,'  number of particles : ',particle_list%n_particles

call cpu_time(t0)

! Select the first and second physics variables
i_var = (/1, 2/)

! Omega_norm is the normalization to JOREK units of an angular frequency, per B
! It is Omega = qB/m sqrt(mu_0 rho) with q = e, m = m_p*central_mass and n = 10^20, divided by B
omega_norm = EL_CHG / (MASS_PROTON * central_mass) * SQRT(MU_ZERO * MASS_PROTON * central_mass * 1.D20)

! This is a measure for the size of the central element
r_crit2 = 0.01 * (node_list%node(element_list%element(1)%vertex(1))%x(1,1) - node_list%node(element_list%element(1)%vertex(2))%x(1,1))**2

!$omp parallel default(none) &
!$omp   shared(particle_list, node_list, element_list, t_step,n_step, F0, omega_norm, i_var, tp, rp, zp, pp, wp, mp, n_done, r_crit2) &
!$omp   private(i, j, k, particle, x, v, i_elm, j_elm, psi, psi_s, psi_t, psi_R, psi_Z, R, R_s, R_t, Z, Z_s, Z_t,st_jac,              &
!$omp           psi_prev, qom, B0, B02, E0, v0, delta_v, R_prev, Z_prev, delta_phi, R_step, Z_step, changed, lost, search,            &
!$omp           R_axis, Z_axis, psi_axis, CR, CZ, U_prev,                                                                             &
!$omp           U, U_s, U_t, U_R, U_Z, U_phi, delta_x, x_prev, v_prev, P, P_s, P_t, error_RZ, R_out ,Z_out, ielm_out, s_out, t_out, ifail)

!$omp do
do i = 1, particle_list%n_particles

  particle = particle_list%particle(i)

  if (particle%lost) cycle

  ! Get s,t,phi
  x(1:2) = particle%st
  x(3)   = particle%x(3)
  v      = particle%v
  i_elm  = particle%i_elm

  ! Interpolate the fields to get psi and U at the current position
  call interp_PRZ(node_list,element_list,i_elm,i_var,2,x(1),x(2),x(3),P,P_s,P_t,R,R_s,R_t,Z,Z_s,Z_t)

  ! Save the curent values. This implies a zero previous step n-1/2 in the Boris method.
  psi = P(1)
  U   = P(2)

  ! Perform n_step Boris method steps here
  do j = 1, n_step

    ! This is the jacobian of the transformation from RZ to st
    st_jac = R_s * Z_t - R_t * Z_s

    ! And these are the local derivatives
    psi_s = P_s(1); psi_t = P_t(1)
    u_s   = P_s(2); u_t   = P_t(2)

    ! Save the previous values
    psi_prev   = psi
    U_prev     = U

    ! Calculate the derivatives to R and Z
    psi_R    = (  psi_s * Z_t - psi_t * Z_s ) / st_jac
    psi_Z    = (- psi_s * R_t + psi_t * R_s ) / st_jac
    U_R      = (  u_s   * Z_t - u_t   * Z_s ) / st_jac
    U_Z      = (- u_s   * R_t + u_t   * R_s ) / st_jac
    ! And assume for now no change in the phi-direction
    U_phi    = 0.d0

    ! Calculate the normalized fraction q/m
    qom = particle%q / particle%mass * omega_norm

    ! Calculate the magnetic field multiplied by q/m Delta_t / 2 and its magnitude
    B0     = (/ + psi_Z, - psi_R, F0 /) / R * qom * t_step / 2.d0
    B02    = dot_product(B0,B0)

    ! Calculate the local electric field, obtained from E=-Grad (u F0), multiplied by q/m Delta_t/2
    E0     = (/ - F0 * U_R, - F0 * U_Z, - F0 * U_phi / R /) * qom * t_step / 2.

    ! Calculate the geometric factor f = tan(q/m delta_t/2 |B|)/|B|
    f = tan(qom * t_step / 2.d0 * sqrt(B02)) / sqrt(B02)

    ! Perform the first half step update (v0 = v^-)
    v0       = v + E0

    ! Perform the full step rotation
    v0       = (v0 + 2.d0 f / (1.d0 + f**2 * B02) * (+ cross_product(v0,B0) &
                                                     + f * B0 * dot_product(v0,B0) &
                                                     - f * v0 * B02)
    
    ! Perform the second half step update (calculate v^n+1/2 from v+)
    v        = v0 + E0

    ! Update the position
    R_prev = R
    Z_prev = Z

    R_step = R_prev + t_step * v(1)
    Z_step = Z_prev + t_step * v(2)

    delta_phi = t_step * v(3) / R_step

    x(3) = x(3) + delta_phi

    do k = 1, 3

      x(1) = x(1) + ( Z_t * (R_step-R) - R_t * (Z_step-Z)) / st_jac
      x(2) = x(2) + (-Z_s * (R_step-R) + R_s * (Z_step-Z)) / st_jac

      call interp3_RZ(node_list,element_list,i_elm,x(1),x(2),R,R_s,R_t,Z,Z_s,Z_t)

      st_jac = R_s * Z_t - R_t * Z_s

    enddo

    v_prev = v

    v(1) =   cos(delta_phi) * v_prev(1) + sin(delta_phi) * v_prev(3)
    v(3) = - sin(delta_phi) * v_prev(1) + cos(delta_phi) * v_prev(3)

    R_step = 0.5d0 * ( R_step + R_prev + t_step * v(1) )
    Z_step = 0.5d0 * ( Z_step + Z_prev + t_step * v(2) )

    x(1) = x(1) + ( Z_t * (R_step-R) - R_t * (Z_step-Z)) / st_jac
    x(2) = x(2) + (-Z_s * (R_step-R) + R_s * (Z_step-Z)) / st_jac

    do m = 1, 2

      call check_element_boundary(element_list,i_elm,x(1:2),x_prev(1:2),j_elm,x(1:2),delta_x(1:2),changed, lost, search)

      if (search) then
        call find_RZ(node_list,element_list,R_step,Z_step,R_out,Z_out,ielm_out,s_out,t_out,ifail)
        x(1)  = s_out
        x(2)  = t_out
        i_elm = ielm_out
        exit
      endif

      if (lost) exit

      if (changed) then

        do k = 1, 3

          call interp3_RZ(node_list,element_list,i_elm,x(1),x(2),R,R_s,R_t,Z,Z_s,Z_t)

          st_jac = R_s * Z_t - R_t * Z_s

          x(1) = x(1) + ( Z_t * (R_step-R) - R_t * (Z_step-Z)) / st_jac
          x(2) = x(2) + (-Z_s * (R_step-R) + R_s * (Z_step-Z)) / st_jac

!          error_RZ = sqrt((R-R_step)**2 + (Z - Z_step)**2)

        enddo

      else
        exit
      endif

    enddo

    if (maxval(abs(x(1:2)-0.5)) .gt. 0.5) then
      call find_RZ(node_list,element_list,R_step,Z_step,R_out,Z_out,ielm_out,s_out,t_out,ifail)
      x(1)  = s_out
      x(2)  = t_out
      i_elm = ielm_out
    endif

    if (lost) then
      particle%lost = lost
      exit
    endif

    call interp_PRZ(node_list,element_list,i_elm,i_var,2,x(1),x(2),x(3),P,P_s,P_t,R,R_s,R_t,Z,Z_s,Z_t)

    psi = P(1)
    U   = P(2)

    R = R_step
    Z = Z_step

 !   mp(j,i) = R_step * v(3) + 0.5 * qom * (psi + psi_prev)
 !   Wp(j,i) = 0.5 * particle%mass * (v(1)*v(1) + v(2)*v(2) + v(3)*v(3)) + 0.5 * omega_norm * particle%q * F0 * (U + U_prev)

 !   rp(j,i) = R
 !   zp(j,i) = Z
 !   pp(j,i) = x(3)
 !   tp(j)   = j * t_step

    n_done = j

  enddo

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

write(*,*) 'Done particles update'

!open(21,file='traject.txt')
!write(21,*) '          X               Y               Z                R              Wp             Mp         time'
!do i=1,n_done
!  write(21,'(8e16.8)') rp(i,1)*cos(pp(i,1)),zp(i,1),rp(i,1)*sin(pp(i,1)),rp(i,1),wp(i,1),mp(i,1),tp(i)
!enddo
!close(21)

return
end
