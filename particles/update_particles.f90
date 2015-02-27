subroutine boris_update(particle_list, t_step, n_step)

use parameters
use data_structure
use mod_particles
use nodes_elements
implicit none

type (type_particle_list) :: particle_list
real*8,  intent(in)       :: t_step
integer, intent(in)       :: n_step

type (type_particle)      :: particle
real*8                    :: B(3), E(3), B0(3), E0(3), x_prev(3), v_prev(3)
real*8                    :: v0(3), v_up(3), delta_v(3), delta_phi
real*8                    :: qom, t_step, twoB02, B02, psi_R, psi_Z, psi_prev
real*8                    :: psi, psi_s, psi_t, t1, t2
real*8                    :: R, R_s, R_t, Z, Z_s, Z_t, st_jac, delta_x(2), R1, Z1, R_prev, Z_prev, R_step, Z_step
integer                   :: i, j, n_step, i_var, i_elm, j_elm, format_rst, error
logical                   :: changed

real*8, allocatable :: rp(:,:), zp(:,:), tp(:), Wp(:,:), mp(:,:)
allocate(rp(n_step,particle_list%n_particles),zp(n_step,particle_list%n_particles),pp(n_step,particle_list%n_particles))
allocate(Wp(n_step,,particle_list%n_particles),tp(n_step),mp(n_step,,particle_list%n_particles))


do  i = 1,  particle_list%n_particles

  particle = particle_list%particle(j)

  call interp_PRZ(node_list,element_list,particle%i_elm,1,1,particle%x(1),particle%x(2),psi,psi_s,psi_t,R,R_s,R_t,Z,Z_s,Z_t)

  do j = 1 ,n_step

    psi_prev = psi
    psi_R    = (  psi_s * Z_t - psi_t * Z_s ) / st_jac
    psi_Z    = (- psi_s * R_t + psi_t * R_s ) / st_jac

    st_jac = R_s * Z_t - R_t * Z_s

    q0m = particle%q / particle%mass

    B0     = (/ - psi_Z, + psi_R, F0 /) / R * qom * t_step / 2.  ! physical components
    B02    = inner_product(B0,B0)
    twoB02 = 2.d0 / (1.d0 + B02)

    E      = (/ 0.0, 0.0, 0.0 /)
    E0     = E * qom * t_step / 2.

    v0       = v + E0                                 ! physical components

    delta_v  = cross(v0,B0) + dot_product(v0,B0) * B0 - B02 * v0

    v        = v0 + twoB02 * delta_v + E0              ! physical components

    R_prev = R
    Z_prev = Z

    R_step = R_prev + t_step * v(1)
    Z_step = Z_prev + t_step * v(2)

    delta_phi = t_step * v(3) / R_step

    x(3) = x(3) + delta_phi

    do j=1,3

      st_jac = R_s * Z_t - R_t * Z_s

      x(1) = x(1) + ( Z_t * (R_step-R) - R_t * (Z_step-Z)) / st_jac
      x(2) = x(2) + (-Z_s * (R_step-R) + R_s * (Z_step-Z)) / st_jac

      call interp_RZ(node_list,element_list,i_elm,x(1),x(2),R,R_s,R_t,Z,Z_s,Z_t)

    enddo

    v_prev = v

    v(1) =   cos(delta_phi) * v_prev(1) + sin(delta_phi) * v_prev(3)
    v(3) = - sin(delta_phi) * v_prev(1) + cos(delta_phi) * v_prev(3)

    R_step = 0.5 * ( R + R_prev + t_step * v(1) )
    Z_step = 0.5 * ( Z + Z_prev + t_step * v(2) )

    x(1) = x(1) + ( Z_t * (R_step-R) - R_t * (Z_step-Z)) / st_jac
    x(2) = x(2) + (-Z_s * (R_step-R) + R_s * (Z_step-Z)) / st_jac

    changed = .false.

    call check_element_boundary(element_list,i_elm,x,x_prev,j_elm,x,delta_x,changed)

    if (changed) then

      do j=1,3

        call interp_RZ(node_list,element_list,i_elm,x(1),x(2),R,R_s,R_t,Z,Z_s,Z_t)

        st_jac = R_s * Z_t - R_t * Z_s

        x(1) = x(1) + ( Z_t * (R_step-R) - R_t * (Z_step-Z)) / st_jac
        x(2) = x(2) + (-Z_s * (R_step-R) + R_s * (Z_step-Z)) / st_jac

      enddo

    endif

    call interp_PRZ(node_list,element_list,i_elm,i_var,1,x(1),x(2),psi,psi_s,psi_t,R,R_s,R_t,Z,Z_s,Z_t)

  enddo

  particle_list%particle(i)%st    = x
  particle_list%particle(i)%v     = v
  particle_list%particle(i)%x     = (/R, Z, phi /)
  particle_list%particle(i)%i_elm = i_elm

  mp(j,i) = R * v(3) - 0.5 * qom * (psi + psi_prev)
  Wp(j,i) = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)

  rp(j,i) = R
  zp(j,i) = Z
  pp(j,i) = x(3)
  tp(j)   = j * t_step

enddo

return
end