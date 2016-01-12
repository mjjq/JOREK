!> Initialize_particles creates n_particles, divided over all processors
!!
!! n_particles/n_cpu particles are created per processor.
!! The total number may therefore be slightly lower than expected
!! It is parallelized with OMP.
subroutine initialise_particles(my_id,n_cpu,node_list,element_list,particle_list,boxcenter,boxwidth,n_particles)

use constants
use data_structure
use phys_module, only : F0, central_density, central_mass
use mod_particles
use openadas

implicit none

! Routine parameters
type (type_node_list)    , intent(in)  :: node_list
type (type_element_list) , intent(in)  :: element_list
type (type_particle_list), intent(out) :: particle_list
integer, intent(in) :: my_id, n_cpu, n_particles
real*8 , intent(in) :: boxcenter(3), boxwidth(3)

! Internal variables
type (type_particle)      :: particle
real*8  :: R_in, Z_in, R_out, Z_out, phi_in, vx, vy, vz
real*8  :: R, R_s, R_t, Z, Z_s, Z_t, B_field(3), B_0, st_jac
real*8  :: psi_s, psi_t, psi_R, psi_Z, s_elm, t_elm, mass_ion, ran3(3), ran4(4), P(4), P_s(4), P_t(4), P_phi(4)
real*8  :: V_thermal, V_norm, mass_main_ion, XYZ_in(3), background_density, background_kbT, background_kelvin
real*8  :: Z_coronal, radiation_coronal, mass_impurity, atomic_mass_impurity
integer :: i_var(4), i_elm, ifail, i, i_done, n_done, i_max, n_max, n_threads, omp_tid, n_particles_thread

integer, external :: omp_get_num_threads, omp_get_thread_num

if (n_particles .gt. n_particles_max) then
  STOP 'n_particles .gt. n_particles_max'
endif

! Divide the particles over all processors
particle_list%n_particles = n_particles/n_cpu
allocate(particle_list%particle(particle_list%n_particles))


if (my_id .eq. 0) then
  write(*,*) '**********************************'
  write(*,*) '*      initialising particles    *'
  write(*,*) '**********************************'
  write(*,*) ' number of particles : ',n_particles
endif

!$omp parallel default(none) &
!$omp   shared(particle_list, node_list, element_list, n_max, boxcenter, boxwidth, i_var, mass_main_ion, atomic_mass_impurity, &
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
mass_impurity        = mass_proton * atomic_mass_impurity ! Tungsten

v_norm = sqrt(mu_zero * mass_main_ion * central_density * 1.d20)      ! JOREK normalisation for velocity
write(*,*) "sqrt(mu_0 rho_0) = ", v_norm


!$omp do
do i_max = 1, n_max

  if (n_done .ge. n_particles_thread) cycle

  if (mod(i_max,10000) .eq. 0) write(*,*) ' progress : ',i_max,n_done

  call random_number(ran3)

  R_in   = boxcenter(1) + boxwidth(1)*(ran3(1)-0.5d0)
  Z_in   = boxcenter(2) + boxwidth(2)*(ran3(2)-0.5d0)
  phi_in = boxcenter(3) + boxwidth(3)*(ran3(3)-0.5d0)

  call find_RZ(node_list,element_list,R_in,Z_in,R_out,Z_out,i_elm,s_elm,t_elm,ifail)

!  if ( sqrt((R_in-3.3)**2 + Z_in**2) .gt. 0.4 ) cycle

  if (ifail .eq. 0) then

    i_var = (/ 1, 5, 6, 7 /)
    call interp_PRZ(node_list,element_list,i_elm,i_var,size(i_var),s_elm,t_elm,phi_in,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)

    background_density = P(2) * 1d20                             ! plasma density [1/m^3]
    background_kbT     = P(3) /(MU_ZERO*central_density*1.d20)   ! plasma temperature [J]
    background_kelvin  = background_kbT / EL_CHG                 ! plasma temperature [K]

    V_thermal = sqrt(background_kbT / mass_impurity)      ! [m/s]

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

write(*,*) ' number of particles : ',particle_list%n_particles
if (my_id .eq. 0) then
  write(*,*) '* done initialising particles    *'
  write(*,*) '**********************************'
endif

end
