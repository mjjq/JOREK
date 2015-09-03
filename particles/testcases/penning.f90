!> This program tests the reproduction of the trajectory of a single particle in a penning trap
!!
!! The electric and magnetic field are expressed in terms of the local coordinates s and t on the JOREK mesh
!! A single particle is released at position x0 with velocity v0 and tracked until t_end
!! The output is a hdf5 file containing the error in particle position for different timestep sizes dt
program penning

use data_structure
use phys_module
use basis_at_gaussian
use nodes_elements
use mod_particles
use clock_module
use coordinate_transforms ! For solution of penning trap trajectory

implicit none

! Define our particle list
type (type_particle_list) :: particle_list !< This contains all particles used in the simulations

! Penning trap parameters
real*8, parameter :: omega_e = 4.9d0
real*8, parameter :: omega_b = 25.d0
real*8, parameter :: qom     = 1.d0
real*8, parameter :: epsilon = -1.d0
real*8, parameter :: x0(3)   = [10.d0,0.d0,0.d0] ! in RZPhi
real*8, parameter :: v0(3)   = [50.d0,20.d0,0.d0] ! in RZPhi = [100,0,20] in xyz
real*8, parameter :: B0      = omega_b*qom
real*8, parameter :: Phi0    = epsilon*omega_e**2*qom*0.5d0
! The particle remains between +- 3 in Z and 8 and 13 in R for these parameters. set this in an input file

! Local variables
integer :: i,ielm_out,ifail
real*8  :: R,R_s,R_t,R_st,Z,Z_s,Z_t,Z_st
real*8  :: s,t
real*8  :: x_a(3), x_e(3), v(3)

! Fake MPI presence
integer, parameter :: my_id = 0


write(*,*) '***************************************'
write(*,*) '* JOREK2 : Penning trap test          *'
write(*,*) '***************************************'

call initialise_parameters(my_id, "__NO_FILENAME__")

! Only 1 toroidal mode (n=0)
mode(1) = 0

! Create a simple grid (is good enough for these fields)
call grid_bezier_square(n_R, n_Z, R_begin, R_end, Z_begin, Z_end, .false., node_list, element_list)
call update_neighbours(element_list,node_list)


write(*,*) 'Initializing fields'
! Set the fields in the first two parameters
do i=1,node_list%n_nodes
  ! Get value and derivatives
  R    = node_list%node(i)%x(1,1)
  R_s  = node_list%node(i)%x(2,1)
  R_t  = node_list%node(i)%x(3,1)
  R_st = node_list%node(i)%x(4,1)
  Z    = node_list%node(i)%x(1,2)
  Z_s  = node_list%node(i)%x(2,2)
  Z_t  = node_list%node(i)%x(3,2)
  Z_st = node_list%node(i)%x(4,2)

  ! B has only a z component, so F = 0, Psi(R) = -B0/2 R^2
  node_list%node(i)%values(1,1,1) = -0.5d0*B0*R**2
  node_list%node(i)%values(1,2,1) = -B0*R*R_s
  node_list%node(i)%values(1,3,1) = -B0*R*R_t
  node_list%node(i)%values(1,4,1) = -B0*R_s*R_t - B0*R*R_st

  ! E = -Grad F0*U, U = E0(R^2 - 2 Z^2) (see model001/initial_conditions.f90 for reference)
  node_list%node(i)%values(1,1,2) = Phi0*(R**2-2.d0*Z**2)
  node_list%node(i)%values(1,2,2) = Phi0*(2.d0*R*R_s - 4.d0*Z*Z_s)
  node_list%node(i)%values(1,3,2) = Phi0*(2.d0*R*R_t - 4.d0*Z*Z_t)
  node_list%node(i)%values(1,4,2) = Phi0*(2.d0*R_s*R_t - 4.d0*Z_s*Z_t &
                                        + 2.d0*R*R_st  - 4.d0*Z*Z_st) ! U_RZ = 0

enddo
! The electric potential is F0*U, so set F0 to 1
F0 = 1.d0

write(*,*) 'Initializing particle'
! Initialize the particle
particle_list%n_particles = 1
allocate(particle_list%particle(particle_list%n_particles))

call find_RZ(node_list,element_list,x0(1),x0(2),R,Z,ielm_out,s,t,ifail)
if (ifail .ne. 0) write(*,*) "Error finding initial particle in grid"
particle_list%particle(1)%st = [s,t]
particle_list%particle(1)%i_elm = ielm_out
particle_list%particle(1)%x = [R,Z,x0(3)]
particle_list%particle(1)%v = v0
particle_list%particle(1)%q = 1.d0
particle_list%particle(1)%mass = 1.d0
particle_list%particle(1)%weight = 1.d0
particle_list%particle(1)%lost = .false.

! First test a constant magnetic field in z-direction
F0 = 0.d0 ! By disabling the magnetic field
call update_particles(1, particle_list, tstep, nstep)
! And test the conservation of energy in the xy-plane
v = particle_list%particle(1)%v
write(*,*) "V: ", v(1)**2+v(3)**2, "Err: ", v(1)**2+v(3)**2-50.d0**2

! Reset the particle position
particle_list%particle(1)%x = [R,Z,x0(3)]
particle_list%particle(1)%v = v0
particle_list%particle(1)%st = [s,t]
particle_list%particle(1)%i_elm = ielm_out
! Enable the electric field
F0 = 1.d0

! Run for n_steps with timestep dt
call update_particles(1, particle_list, tstep, nstep, 0.d0)

! Check position against analytical result
write(*,*) "======== x"
x_e = particle_list%particle(1)%x
x_a = analytical_trajectory(tstep*real(nstep))
write(*,*) x_e
write(*,*) x_a
! Difference in xy coordinate
write(*,*) sqrt(x_e(1)**2+x_a(1)**2 - 2.d0*x_e(1)*x_a(1)*(cos(x_e(3))*cos(x_a(3))+sin(x_e(3))*sin(x_a(3))))

contains
  !> This function calculates the analytical trajectory of a particle in a penning trap
  !!
  !! The parameters epsilon, omega_e, omega_b, x0, v0 are used and the solution is
  !! returned in RZPhi coordinates.
  pure function analytical_trajectory(t) result(x)
    real*8, intent(in) :: t !< The time at which to calculate the solution value
    real*8             :: x(3) !< The position of the particle at time t (in RZPhi)

    ! Internal variables
    real*8 :: omega_plus, omega_minus
    real*8 :: R_plus, R_minus
    real*8 :: T_plus, T_minus
    real*8 :: omega
    complex(kind=8) :: w
    real*8 :: x0_xyz(3), v0_xyz(3)

    ! Initialize in xyz coordinates
    x0_xyz = RZPhiToXYZ(x0)
    v0_xyz = RZPhiToXYZ(v0)

    ! Some initialization
    omega = sqrt(-2*epsilon)*omega_e
    omega_plus  = 0.5*(omega_b + sqrt(omega_b**2 + 4*epsilon*omega_e**2))
    omega_minus = 0.5*(omega_b - sqrt(omega_b**2 + 4*epsilon*omega_e**2))
    R_minus = (omega_plus * x0_xyz(1) + v0_xyz(2))/(omega_plus - omega_minus)
    R_plus  = x0_xyz(1) - R_minus
    T_minus = (omega_plus * x0_xyz(2) - v0_xyz(1))/(omega_plus - omega_minus)
    T_plus  = x0_xyz(2) - T_minus

    ! Calculate the result in the x-y plane in terms of w = x + iy
    w = cmplx(R_plus ,T_plus , 8)*exp(cmplx(0.d0,-omega_plus*t, 8)) + &
        cmplx(R_minus,T_minus, 8)*exp(cmplx(0.d0,-omega_minus*t, 8))

    ! Return result in RZPhi coordinates
    x(1) = abs(w) ! R = |w|
    x(2) = x0_xyz(3)*cos(omega*t)+v0_xyz(3)*sin(omega*t)/omega
    x(2) = atan2(dimag(w), dble(w)) ! Phi = atan2(y,x)
  end function analytical_trajectory
end
