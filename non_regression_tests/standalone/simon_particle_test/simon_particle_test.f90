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
call initialise_particles_simon(node_list,element_list,particle_list)


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
do i=1,nstep/nout
  call update_particles(my_id, particle_list, tstep, nout, 1.d0) ! TODO fix bug with present of toroidal field factor
  write(21,"(100g16.8)") i*nout*tstep, (particle_list%particle(j)%x(1:2), j=1,particle_list%n_particles)
  write(22,"(100g16.8)") i*nout*tstep, (guiding_center_position(particle_list%particle(j), tstep), j=1,particle_list%n_particles)
  ! Save particle position to file (or stdout in this case)
  if (save_vtk) then
    write(filename,"(A,I0.5,A)") "particles", i*nout, ".vtk"
    call particles_vtk(particle_list,filename)
  endif
enddo
write(*,*) "Energy errors:"
write(*,"(100g16.8)") tstep, (sum(particle_list%particle(j)%v**2)-particle_energies(j), j=1,particle_list%n_particles)

call MPI_FINALIZE(ierr)

contains

!> Setup the 3 particles we will use and test some invariants
subroutine initialise_particles_simon(node_list,element_list,particle_list)

use mod_particles
use constants
use data_structure
use phys_module, only: F0, central_density, central_mass

implicit none

type (type_node_list), intent(in)        :: node_list
type (type_element_list), intent(in)     :: element_list
type (type_particle_list), intent(inout) :: particle_list

type (type_particle)      :: particle

real*8  :: particle_energy(3), particle_energy_perp(3), R_in, Z_in, phi_in, R_out, Z_out
real*8  :: P(1), P_s(1), P_t(1), R, R_s, R_t, Z, Z_s, Z_t, B_field(3), B_0, st_jac, v_norm, v_perp, v_par
real*8  :: psi_s, psi_t, psi_R, psi_Z, s_elm, t_elm, mass_ion
integer :: i_var(1), i_elm, ifail, i_part


particle_list%n_particles = 3
allocate(particle_list%particle(particle_list%n_particles))

! Start position of all 3 particles
R_in   = 3.025
Z_in   = 0.0
phi_in = 0.0

call find_RZ(node_list,element_list,R_in,Z_in,R_out,Z_out,i_elm,s_elm,t_elm,ifail)

particle_energy      = (/ 170., 50.,    164.    /)      ! [eV]
particle_energy_perp = (/  40., 49.585, 161.832 /)      ! [eV]

mass_ion = 4.d0 * mass_proton
v_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

i_var = (/ 1 /)
call interp_PRZ(node_list,element_list,i_elm,i_var,1,s_elm,t_elm,phi_in,P,P_s,P_t,R,R_s,R_t,Z,Z_s,Z_t)

st_jac = R_s * Z_t - R_t * Z_s
psi_s  = P_s(1); psi_t = P_t(1);
psi_R  = (  psi_s * Z_t - psi_t * Z_s ) / st_jac
psi_Z  = (- psi_s * R_t + psi_t * R_s ) / st_jac

B_field     = (/ + psi_Z, - psi_R, F0 /) / R
B_0         = sqrt(dot_product(B_field,B_field))

do i_part=1,3
  ! Assumes B_field only in toroidal direction!
  v_par  = sqrt(2.d0*(particle_energy(i_part) - particle_energy_perp(i_part)) * el_chg / mass_ion)
  v_perp = sqrt(particle_energy_perp(i_part) * 2.d0 * el_chg / mass_ion)

  write(*,'(A,3e16.8)') ' perpendicular velocity : ',v_perp,v_norm,v_perp * v_norm
  write(*,'(A,3e16.8)') ' gyro radius            : ',mass_ion * v_perp / (el_chg * B_0)
  write(*,'(A,3e16.8)') ' gyro frequency         : ',el_chg * B_0 / mass_ion, el_chg * B_0 / mass_ion * v_norm

  write(*,'(A,4e16.8)') 'CHECK energy : ', v_par,v_perp,(v_par**2+v_perp**2) * 0.5 * mass_ion / el_chg, particle_energy(i_part)

  particle%v = -v_norm/B_0 * (v_par * B_field + v_perp * cross_product((/1.d0,0.d0,0.d0/),B_field))

  particle%x       = (/ R_out, Z_out, 0.d0 /)        !< particle position in real space
  particle%st      = (/ s_elm, t_elm /)            !< particle position in the finite element (i_elm)
  particle%i_elm   = i_elm
  particle%q       = 2                             !< charge
  particle%mass    = 4.d0                          !< mass in amu
  particle%weight  = 1.                            !< weight (i.e. number of particles)
  particle%lost    = .false.

  particle_list%particle(i_part) = particle

  write(*,'(A,3e16.8)') 'Velocity : ',particle%v
enddo
end subroutine initialise_particles_simon
end program simon_particle_test
