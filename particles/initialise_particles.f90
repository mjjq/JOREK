subroutine initialise_particles(my_id,n_cpu,node_list,element_list,particle_list,boxsize)

use constants
use data_structure
use phys_module, only : F0, central_density, central_mass
use mod_particles
use openadas

implicit none

type (type_node_list)     :: node_list
type (type_element_list)  :: element_list
type (type_particle_list) :: particle_list
type (type_particle)      :: particle

integer :: my_id
integer :: n_cpu
real*8  :: boxsize(3)

real*8  :: R_in, Z_in, R_out, Z_out, phi_in, vx, vy, vz
real*8  :: R, R_s, R_t, Z, Z_s, Z_t, B_field(3), B_0, st_jac
real*8  :: psi_s, psi_t, psi_R, psi_Z, s_elm, t_elm, mass_ion, ran3(3), ran4(4), P(4), P_s(4), P_t(4), P_phi(4)
real*8  :: V_thermal, V_norm, mass_main_ion, XYZ_in(3), background_density, background_kbT, background_kelvin
real*8  :: Z_coronal, radiation_coronal, mass_impurity, atomic_mass_impurity
integer :: i_var(4), i_elm, ifail, i, i_done, n_done, i_max, n_max, n_threads, omp_tid, n_particles_thread

integer, external :: omp_get_num_threads, omp_get_thread_num

particle_list%n_particles = 100000

if (particle_list%n_particles .gt. n_particles_max) then
  STOP 'particle_list%n_particles .gt. n_particles_max'
endif

allocate(particle_list%particle(particle_list%n_particles))

write(*,*) '**********************************'
write(*,*) '*      initialising particles    *'
write(*,*) '**********************************'
write(*,*) ' number of particles : ',particle_list%n_particles

!$omp parallel default(none) &
!$omp   shared(particle_list, node_list, element_list, n_max, boxsize, i_var, mass_main_ion, atomic_mass_impurity, &
!$omp          mass_impurity, central_density, central_mass, n_threads)                                            &
!$omp   private(i_max, ran3, ran4, XYZ_in, R_in, Z_in, phi_in, R_out, Z_out, i_elm, s_elm, t_elm, ifail,           &
!$omp           P, P_s, P_t, R, R_s, R_t, Z, Z_s, Z_t, background_density, background_kbT, background_kelvin,      &
!$omp           V_thermal, v_norm, vx, vy, vz, Z_coronal, radiation_coronal, particle, n_done,                     &
!$omp           omp_tid, n_particles_thread, P_phi)

n_threads = omp_get_num_threads()
omp_tid   = omp_get_thread_num()

n_done = 0
i_max  = 0
n_particles_thread = particle_list%n_particles / n_threads

n_max  = 10 * particle_list%n_particles

call random_seed()

mass_main_ion        = mass_proton * central_mass
atomic_mass_impurity = 183.84
mass_impurity        = mass_proton * atomic_mass_impurity ! Tungsten

!$omp do
do i_max = 1, n_max

  if (n_done .ge. n_particles_thread) cycle

  if (mod(i_max,10000) .eq. 0) write(*,*) ' progress : ',i_max,n_done

  call random_number(ran3)

  XYZ_in = - boxsize + 2.* ran3 * boxsize
  R_in   = sqrt(XYZ_in(1)**2 + XYZ_in(2)**2)
  Z_in   = XYZ_in(3)
  phi_in = atan2(XYZ_in(2),XYZ_in(1))

  call find_RZ(node_list,element_list,R_in,Z_in,R_out,Z_out,i_elm,s_elm,t_elm,ifail)

!  if ( sqrt((R_in-3.3)**2 + Z_in**2) .gt. 0.4 ) cycle

  if (ifail .eq. 0) then

    i_var = (/ 1, 5, 6, 7 /)
    call interp_PRZ(node_list,element_list,i_elm,i_var,size(i_var),s_elm,t_elm,phi_in,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)

    background_density = P(2) * 1d20                             ! plasma density [1/m^3]
    background_kbT     = P(3) /(MU_ZERO*central_density*1.d20)   ! plasma temperature [J]
    background_kelvin  = background_kbT / EL_CHG                 ! plasma temperature [K]

    V_thermal = sqrt(background_kbT / mass_impurity)      ! [m/s]

    v_norm = sqrt(mu_zero * mass_main_ion * central_density * 1.d20)      ! JOREK normalisation for velocity

    call random_number(ran4)

    vx = V_thermal * sqrt(-2*Log(ran4(1))) * cos(TWOPI*ran4(2))  ! [m/s]  uses box-mueller transform
    vy = V_thermal * sqrt(-2*Log(ran4(1))) * sin(TWOPI*ran4(2))
    vz = V_thermal * sqrt(-2*Log(ran4(3))) * sin(TWOPI*ran4(4))

    call interpolate_coronal(alog10(background_density*1d-6),alog10(background_kelvin),Z_coronal,radiation_coronal)

    particle%v(1) =   vx * cos(phi_in) + vy * sin(phi_in)   ! V_R
    particle%v(3) = - vx * sin(phi_in) + vy * cos(phi_in)   ! V_phi [physical component]
    particle%v(2) =   vz                                    ! V_Z

    particle%v = particle%v * v_norm

    particle%x       = (/ R_out, Z_out, phi_in /)    !< particle position in real space
    particle%st      = (/ s_elm, t_elm /)            !< particle position in the finite element (i_elm)
    particle%i_elm   = i_elm                         !< element containing this particle
    particle%q       = Z_coronal                     !< charge (initialised with the coronal equilibrium value)
    particle%mass    = atomic_mass_impurity          !< mass (Tungsten)
    particle%weight  = 1.                            !< weight (i.e. number of particles)
    particle%lost    = .false.                       !< active particle

    particle_list%particle(omp_tid*n_particles_thread+n_done+1) = particle

!    write(*,'(4i5,3f14.6)') omp_tid,i_max,n_done,omp_tid*n_particles_thread+n_done+1,ran3

    n_done = n_done + 1

    if (i_max .eq. n_max) write(*,'(A,3i6)') 'WARNING : ',omp_tid,i_max,n_max

  endif

enddo
!$omp end do

particle_list%n_particles = n_particles_thread * n_threads

!$omp end parallel

!do i=1, particle_list%n_particles
!  write(*,'(i5,3f14.6)') particle_list%particle(i)%i_elm,particle_list%particle(i)%x
!enddo

write(*,*) ' number of particles : ',particle_list%n_particles
write(*,*) '* done initialising particles    *'
write(*,*) '**********************************'

return
end

subroutine initialise_particles_simon(my_id,n_cpu,node_list,element_list,particle_list)

use mod_particles
use constants
use data_structure
use phys_module, only : F0, central_density, central_mass

implicit none

type (type_node_list)     :: node_list
type (type_element_list)  :: element_list
type (type_particle_list) :: particle_list
type (type_particle)      :: particle

integer :: my_id
integer :: n_cpu

real*8  :: particle_energy(3), particle_energy_perp(3), R_in, Z_in, phi_in, R_out, Z_out, v_R, v_phi
real*8  :: P(1), P_s(1), P_t(1), P_phi(1), R, R_s, R_t, Z, Z_s, Z_t, B_field(3), B_0, st_jac, v_norm, v_perp2, v_par
real*8  :: psi_s, psi_t, psi_R, psi_Z, s_elm, t_elm, mass_ion
integer :: i_var(1), i_elm, ifail, i_part

particle_list%n_particles = 1


!------------ simons test (see thesis http://www2.ipp.mpg.de/~Simon.Pinches/thesis/node57.html)

R_in   = 3.025
Z_in   = 0.0
phi_in = 0.0

call find_RZ(node_list,element_list,R_in,Z_in,R_out,Z_out,i_elm,s_elm,t_elm,ifail)

particle_energy      = (/ 170., 50.,    164.    /)      ! [eV]
particle_energy_perp = (/  40., 49.585, 161.832 /)      ! [eV]

mass_ion = mass_proton * central_mass

v_norm = sqrt(mu_zero * mass_ion * central_density * 1.d20)

i_var = (/ 1 /)
call interp_PRZ(node_list,element_list,i_elm,i_var,1,s_elm,t_elm,phi_in,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)

st_jac = R_s * Z_t - R_t * Z_s
psi_s  = P_s(1); psi_t = P_t(1);
psi_R  = (  psi_s * Z_t - psi_t * Z_s ) / st_jac
psi_Z  = (- psi_s * R_t + psi_t * R_s ) / st_jac

B_field     = (/ + psi_Z, - psi_R, F0 /) / R
B_0         = sqrt(dot_product(B_field,B_field))

i_part = 3

v_par = sqrt((particle_energy(i_part) - particle_energy_perp(i_part))* 2. / mass_ion * el_chg)

v_phi = v_par * B_0 / B_field(3)

v_perp2 = particle_energy_perp(i_part)* 2. / mass_ion * el_chg

write(*,'(A,3e16.8)') ' perpendicular velocity : ',sqrt(v_perp2),v_norm,sqrt(v_perp2) * v_norm
write(*,'(A,3e16.8)') ' gyro radius            : ',mass_ion * sqrt(v_perp2) / (el_chg * B_0)
write(*,'(A,3e16.8)') ' gyro frequency         : ',el_chg * B_0 / mass_ion, el_chg * B_0 / mass_ion * v_norm

v_R = sqrt(v_perp2 - v_phi**2*(1. - (B_field(3)/B_0)**2)**2 - v_phi * B_field(3)**2 * B_field(2)**2 / B_0**4)

write(*,'(A,3e16.8)') 'CHECK energy : ', v_R,v_phi,(v_R**2+v_phi**2) * 0.5 * mass_ion / el_chg

particle%v(1) = - v_R * v_norm
particle%v(2) = 0.
particle%v(3) = - v_phi * v_norm

particle%x       = (/ R_out, Z_out, 0. /)        !< particle position in real space
particle%st      = (/ s_elm, t_elm /)            !< particle position in the finite element (i_elm)
particle%i_elm   = i_elm
particle%q       = 1.                            !< charge
particle%mass    = central_mass                  !< mass
particle%weight  = 1.                            !< weight (i.e. number of particles)

particle_list%particle(1) = particle

write(*,'(A,3e16.8)') 'Velocity : ',particle%v

return
end
