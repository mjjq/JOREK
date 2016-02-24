!> Initialize_particles creates n_particles, divided over all processors
!!
!! n_particles/n_cpu particles are created per processor (first processor has
!! some more).
!! It is parallelized with OMP.
subroutine initialise_particles(my_id,n_cpu,node_list,element_list,particle_list,boxcenter,boxwidth,n_particles)

!$ use omp_lib
use constants
use data_structure
use phys_module, only : F0, central_density, central_mass, atomic_mass_impurity
use mod_particles

implicit none

! Routine parameters
type (type_node_list)    , intent(in)  :: node_list
type (type_element_list) , intent(in)  :: element_list
type (type_particle_list), intent(out) :: particle_list
integer, intent(in) :: my_id, n_cpu, n_particles
real*8 , intent(in) :: boxcenter(3), boxwidth(3)

! Internal variables
type (type_particle)      :: particle
real*8  :: R_in, Z_in, phi_in
integer :: i, ifail, n_threads
real*8 :: mass_impurity, ran3(3), t0, t1, ostart, oend

if (my_id .eq. 0) then
  write(*,*) '**********************************'
  write(*,*) '*      initialising particles    *'
  write(*,*) '**********************************'
  write(*,*) ' number of particles : ',n_particles
endif

! Divide the particles over all processors, give the first processor the remainder
if (my_id .eq. 0) then
  particle_list%n_particles = n_particles/n_cpu + modulo(n_particles, n_cpu)
else
  particle_list%n_particles = n_particles/n_cpu
endif
write(*,*) my_id, particle_list%n_particles, n_particles, n_cpu
allocate(particle_list%particle(particle_list%n_particles), stat=ifail)
if (ifail .gt. 0) write(*,"(i3,a,i8,a)") my_id, "unable to allocate particle_list%particle(", particle_list%n_particles, ")"


call cpu_time(t0)
!$ ostart = omp_get_wtime()

mass_impurity = mass_proton * atomic_mass_impurity ! Tungsten

! Call find_RZ once to initialise elements_minmax before going OMP (parameters are random-ish and irrelevant)
call find_RZ(node_list,element_list,boxcenter(1),boxcenter(2),R_in,Z_in,i,ran3(1),ran3(2),ifail)

call random_seed()
!$omp parallel do default(none) &
!$omp   shared(particle_list, node_list, element_list, boxcenter, boxwidth, &
!$omp          mass_impurity, central_density, central_mass, n_threads) &
!$omp   private(i, R_in, Z_in, phi_in, ifail, particle, ran3)
do i=1,particle_list%n_particles
  ifail = 1
  do while (ifail .ne. 0)
    ! Generate a random position to put this particle
    call random_number(ran3)
    ! Use inversion sampling to correct for cylindrical coordinates
    R_in   = boxcenter(1) + boxwidth(1)*(sqrt(ran3(1))-0.5d0)
    Z_in   = boxcenter(2) + boxwidth(2)*(ran3(2)-0.5d0)
    phi_in = boxcenter(3) + boxwidth(3)*(ran3(3)-0.5d0)

    call initialize_particle_coronal(node_list, element_list, R_in, Z_in, phi_in, mass_impurity, particle, ifail)
  enddo
  ! Save the result
  particle_list%particle(i) = particle
enddo
!$omp end parallel do

call cpu_time(t1)
!$ oend = omp_get_wtime()
write(*,'(i5,A,2f12.4)') my_id, ' Time particle initialize cpu/wall :',t1-t0, oend-ostart
if (my_id .eq. 0) then
  write(*,*) '* done initialising particles    *'
  write(*,*) '**********************************'
endif

end

!> Function to return a single particle at R_in, Z_in, phi_in at the coronal
!! equilibrium ionization and at the local thermal velocity
subroutine initialize_particle_coronal(node_list, element_list, R_in, Z_in, phi_in, particle_mass, particle, ifail)
use constants
use data_structure
use mod_particles
use phys_module, only : F0, central_density, central_mass, atomic_mass_impurity
use openadas
implicit none

real*8, intent(in) :: R_in, Z_in, phi_in, particle_mass
type (type_node_list)    , intent(in)  :: node_list
type (type_element_list) , intent(in)  :: element_list
type(type_particle), intent(out)       :: particle
integer, intent(out) :: ifail

real*8 :: R_out, Z_out, s_elm, t_elm
integer :: i_var(4), i_elm
real*8, dimension(4) :: P, P_s, P_t, P_phi, ran4
real*8 :: R, R_s, R_t, Z, Z_s, Z_t
real*8 :: background_density, background_kbT, background_kelvin, V_thermal
real*8 :: vx, vy, vz, v_norm
real*8 :: Z_coronal, radiation_coronal
real*8 :: mass_main_ion

! Find s and t at this position and ielm
call find_RZ(node_list,element_list,R_in,Z_in,R_out,Z_out,i_elm,s_elm,t_elm,ifail)

if (ifail .eq. 0) then

  i_var = (/ 1, 5, 6, 7 /)
  call interp_PRZ(node_list,element_list,i_elm,i_var,size(i_var),s_elm,t_elm,phi_in,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)

  background_density = P(2) * 1d20                             ! plasma density [1/m^3]
  background_kbT     = P(3) /(MU_ZERO*central_density*1.d20)   ! plasma temperature [J]
  background_kelvin  = background_kbT / EL_CHG                 ! plasma temperature [K]

  V_thermal = sqrt(background_kbT / particle_mass)      ! [m/s]

  call random_number(ran4)

  vx = V_thermal * sqrt(-2*Log(ran4(1))) * cos(TWOPI*ran4(2))  ! [m/s]  uses box-mueller transform
  vy = V_thermal * sqrt(-2*Log(ran4(1))) * sin(TWOPI*ran4(2))
  vz = V_thermal * sqrt(-2*Log(ran4(3))) * sin(TWOPI*ran4(4))

  call interpolate_coronal(alog10(background_density*1d-6),alog10(background_kelvin),Z_coronal,radiation_coronal)

  particle%v(1) =   vx * cos(phi_in) + vy * sin(phi_in)   ! V_R
  particle%v(3) = - vx * sin(phi_in) + vy * cos(phi_in)   ! V_phi [physical component]
  particle%v(2) =   vz                                    ! V_Z

  mass_main_ion        = mass_proton * central_mass
  v_norm = sqrt(mu_zero * mass_main_ion * central_density * 1.d20)      ! JOREK normalisation for velocity
  particle%v = particle%v * v_norm

  particle%x       = (/ R_out, Z_out, phi_in /)    !< particle position in real space
  particle%st      = (/ s_elm, t_elm /)            !< particle position in the finite element (i_elm)
  particle%i_elm   = i_elm                         !< element containing this particle
  particle%q       = Z_coronal                     !< charge (initialised with the coronal equilibrium value)
  particle%mass    = atomic_mass_impurity          !< mass (Tungsten)
  particle%weight  = 1.                            !< weight (i.e. number of particles)
  particle%lost    = .false.                       !< active particle

  !particle_list%particle(omp_tid*n_particles_thread+n_done+1) = particle

  !n_done = n_done + 1

  !if (i_max .eq. n_max) write(*,'(A,3i6)') 'WARNING : ',omp_tid,i_max,n_max
endif
! if ifail != 0 particle is empty

end subroutine initialize_particle_coronal
