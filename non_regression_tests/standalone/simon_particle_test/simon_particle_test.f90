!> This program tests the reproduction of particle trajectories in a
!! large aspect ratio tokamak with circular flux surfaces.
!!
!! See http://www2.ipp.mpg.de/~Simon.Pinches/thesis/node57.html for details.
!! This program reads a namelist input file containing the JOREK settings, and uses
!! tstep and nstep to determine the integration parameters.
program simon_particle_test

use data_structure
use phys_module
use basis_at_gaussian
use nodes_elements
use mod_particles
use clock_module
use parameters
use constants
use tr_module
use mumps_module

implicit none

interface
  subroutine equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint2,xcase2, nice_q)
    use data_structure
    integer(kind=4),             intent(in)    :: my_id
    integer(kind=4),             intent(in)    :: xcase2
    type (type_node_list),       intent(inout) :: node_list
    type (type_element_list),    intent(inout) :: element_list
    type (type_bnd_node_list)   ,intent(inout) :: bnd_node_list
    type (type_bnd_element_list),intent(inout) :: bnd_elm_list
    logical(kind=4),             intent(in)    :: xpoint2
    logical(kind=4),             intent(in)    :: nice_q
  end subroutine equilibrium

  function guiding_center_position(particle, dt) result(x_gc)
    use mod_particles
    type (type_particle), intent(in) :: particle
    real*8, intent(in) :: dt
    real*8 :: x_gc(3)
  end function guiding_center_position
end interface

! MPI parameters
integer :: required, provided, StatInfo, my_id, n_cpu
integer :: ierr
integer*4  :: rank, comm_size

! Private variables
logical, parameter :: save_vtk = .false.
integer :: i,j
integer :: MPI_GROUP_WORLD
character(len=18) :: filename, fileout_bin
real*8 :: particle_energies(3)

type (type_particle_list) :: particle_list


!! Initialize MPI
required = MPI_THREAD_MULTIPLE

call MPI_Init_thread(required, provided, StatInfo)

call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
my_id = rank
call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
n_cpu = comm_size

call tr_meminit(my_id, n_cpu)


if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '* JOREK2 : Simon Particle test        *'
  write(*,*) '***************************************'
endif


!! Read input parameters
if (my_id .eq. 0) then
  call initialise_parameters(my_id, "__NO_FILENAME__")
endif


!! Initialize grid and solution
call initialise_basis()
element_list%n_elements      = 0
bnd_elm_list%n_bnd_elements  = 0
node_list%n_nodes            = 0

call define_boundary()
if (n_radial <= 0 .or. n_pol <= 0) then
  write(*,*) "Not enough grid cells specified, exiting"
  call exit(1)
endif
call grid_polar_bezier(R_geo, Z_geo, amin, 0.d0, 0.d0, fbnd, fpsi, mf, n_radial, n_pol,    &
  node_list, element_list)

! --- Determine boundary information from the grid
!call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

call MPI_Barrier(MPI_COMM_WORLD,ierr)

call update_neighbours(element_list,node_list)
call broadcast_elements(my_id, element_list)       ! elements
call broadcast_nodes(my_id, node_list)             ! nodes

! Initialize mumps solver
call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
call MPI_GROUP_INCL(MPI_GROUP_WORLD,1,0,MPI_GROUP_MUMPS_EQUIL,ierr)
call MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_GROUP_MUMPS_EQUIL,MPI_COMM_MUMPS_EQUIL,ierr)
if (my_id == 0) call initialise_mumps(MPI_COMM_MUMPS_EQUIL)


!! Calculate equilibrium field
call equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint,xcase, .true.)

! Clean up mumps
mumps_par%JOB = -2
if (my_id == 0) call DMUMPS(mumps_par)


!! Write equilibrium to file
if (my_id .eq. 0) then
  fileout_bin = "jorek_restart.rst"
  write (*,*) " =============>, jorek2, filename = ", fileout_bin
  call export_binary_restart(node_list, element_list, fileout_bin, 1)
endif


!! Initialize particles
call initialise_particles_simon(node_list,element_list,particle_list, t_step_particles)


if (save_vtk) then
  write(filename,"(A,I0.5,A)") "particles", 0, ".vtk"
  call particles_vtk(particle_list,filename)
endif

!! Open output file
open(file="positions_RZ.dat",status="replace",unit=21)
write(21,"(A3,100A16)") "t  ", "R1", "Z1", "R2", "Z2", "R3", "Z3"
open(file="gc_RZ.dat",status="replace",unit=22)
write(22,"(A3,100A16)") "t  ", "R1", "Z1", "R2", "Z2", "R3", "Z3"

! Calculate energies in advance
do j=1,particle_list%n_particles
  particle_energies(j) = sum(particle_list%particle(j)%v**2)
enddo

!! Perform time-stepping
do i=1,n_step_particles/nout_particles
  call update_particles(my_id, particle_list, t_step_particles, nout_particles, 1.d0) ! TODO fix bug with present of toroidal field factor
  write(21,"(100g16.8)") i*nout_particles*t_step_particles, (particle_list%particle(j)%x(1:2), j=1,particle_list%n_particles)
  write(22,"(100g16.8)") i*nout_particles*t_step_particles, (guiding_center_position(particle_list%particle(j), t_step_particles), j=1,particle_list%n_particles)
  ! Save particle position to file (or stdout in this case)
  if (save_vtk) then
    write(filename,"(A,I0.5,A)") "particles", i*nout_particles, ".vtk"
    call particles_vtk(particle_list,filename)
  endif
  write(*,*) "current iteration:", i*nout_particles
enddo
write(*,*) "Energy errors:"
write(*,"(100g16.8)") t_step_particles, (sum(particle_list%particle(j)%v**2)-particle_energies(j), j=1,particle_list%n_particles)

call MPI_FINALIZE(ierr)

contains





!> Setup the 3 particles we will use and test some invariants
subroutine initialise_particles_simon(node_list,element_list,particle_list,dt)

use mod_particles
use constants
use data_structure
use phys_module, only: F0, central_density, central_mass

implicit none

type (type_node_list), intent(in)        :: node_list
type (type_element_list), intent(in)     :: element_list
type (type_particle_list), intent(inout) :: particle_list
real*8, intent(in) :: dt

type (type_particle)      :: particle

real*8  :: particle_energy(3), particle_energy_perp(3), R_in, Z_in, phi_in, R_out, Z_out
real*8  :: E(3), B(3), B_norm, psi, U, t_norm, v_perp, v_par
real*8  :: s_elm, t_elm, mass_ion, v(3), f, qom
integer :: i_var(1), i_elm, ifail, i_part


particle_list%n_particles = 3
allocate(particle_list%particle(particle_list%n_particles))

! Start position of GC of all 3 particles
! TODO calculate guiding center position of particle with negative q to calculate particle position
R_in   = 3.025
Z_in   = 0.0
phi_in = 0.0

call find_RZ(node_list,element_list,R_in,Z_in,R_out,Z_out,i_elm,s_elm,t_elm,ifail)

! These are probably fast alpha particles
particle_energy      = (/ 170., 50.,    164.    /)      ! [eV]
particle_energy_perp = (/  40., 49.585, 161.832 /)      ! [eV]

mass_ion = 4.d0 * mass_proton ! XXX impurity_atomic_mass not used
t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds


call calc_EB(i_elm,(/s_elm,t_elm/),phi_in,E,B,psi,U)
B_norm         = norm2(B)

do i_part=1,3
  ! Assumes B_field only in toroidal direction!
  v_par  = sqrt(2.d0*(particle_energy(i_part) - particle_energy_perp(i_part)) * el_chg / mass_ion)
  v_perp = sqrt(particle_energy_perp(i_part) * 2.d0 * el_chg / mass_ion)

  write(*,'(A,3e16.8)') ' perpendicular velocity : ',v_perp,t_norm,v_perp * t_norm
  write(*,'(A,3e16.8)') ' gyro radius            : ',mass_ion * v_perp / (el_chg * B_norm)
  write(*,'(A,3e16.8)') ' gyro frequency         : ',el_chg * B_norm / mass_ion, el_chg * B_norm / mass_ion * t_norm

  write(*,'(A,4e16.8)') 'CHECK energy : ', v_par,v_perp,(v_par**2+v_perp**2) * 0.5 * mass_ion / el_chg, particle_energy(i_part)

  ! Velocity at t=0
  particle%v = -t_norm/B_norm * (v_par * B + v_perp * cross_product((/1.d0,0.d0,0.d0/),B))

  particle%x       = (/ R_out, Z_out, 0.d0 /)        !< particle position in real space
  particle%st      = (/ s_elm, t_elm /)            !< particle position in the finite element (i_elm)
  particle%i_elm   = i_elm
  particle%q       = 2                             !< charge
  particle%mass    = 4.d0                          !< mass in amu
  particle%weight  = 1.                            !< weight (i.e. number of particles)
  particle%lost    = .false.
  
  qom     = real(particle%q) * el_chg / (particle%mass * atomic_mass_unit)
  f = -qom*dt*0.25d0*t_norm ! Perform the initial half-step without tan correction

  v = particle%v + f*E
  v = (v + 2.d0*f/(1.d0+f**2*dot_product(B,B)) * (cross_product(v,B) - f*v*dot_product(B,B) + f*B*dot_product(v,B)))
  v = v + f*E

  ! particle velocity at t=-1/2 dt
  particle%v = v

  particle_list%particle(i_part) = particle

  write(*,'(A,3e16.8)') 'Velocity : ',particle%v
enddo
end subroutine initialise_particles_simon
end program simon_particle_test
