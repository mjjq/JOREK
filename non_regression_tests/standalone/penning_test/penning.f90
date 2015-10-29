!> This program tests the reproduction of the trajectory of in a penning trap
!!
!! The electric and magnetic field are expressed in terms of
!! the local coordinates s and t on the JOREK mesh.
!! A single particle is released at position x0 with velocity v0
!! and tracked until t_end. The output code signifies whether the
!! tests passed the preset accuracy.
!!
!! The code reads a namelist input file to determine the grid used,
!! the number of grid cells to use and the time step sizes to use.
program penning

use data_structure
use phys_module
use basis_at_gaussian
use nodes_elements
use mod_particles
use clock_module
use mod_coordinate_transforms ! For solution of penning trap trajectory
use parameters
use constants

implicit none

! Define our particle list
type (type_particle_list) :: particle_list !< This contains all particles used in the simulations

! Penning trap parameters (in SI units)
real*8, parameter :: omega_e = 4.9d0
real*8, parameter :: omega_b = 25.d0
real*8, parameter :: epsilon = -1.d0
real*8, parameter :: x0(3)   = [10.d0,0.d0,0.d0] ! in RZPhi
real*8, parameter :: v0_SI(3)   = [50.d0,20.d0,0.d0] ! in RZPhi, = [50,0,20] in xyz
! Corresponding JOREK units
real*8 :: v0(3)
! Penning trap parameters (for fields)
integer, parameter :: charge = 1 ! in electron charges
real*8, parameter :: mass    = 1.d0 ! in atomic mass units
real*8 :: t_norm
real*8 :: qom ! In SI units
real*8 :: B0
real*8 :: Phi0
! The particle remains between +- 3 in Z and 8 and 13 in R for these parameters. set this in an input file

! Local variables
integer :: i,ielm_out,ifail
real*8  :: R,R_s,R_t,R_st,Z,Z_s,Z_t,Z_st
real*8  :: s,t
real*8  :: x_a(3), x_e(3), v(3), err_norm, err_ref
! For the initial half-step
real*8  :: E(3), B(3), B2, f

! Fake MPI presence
integer, parameter :: my_id = 0


write(*,*) '***************************************'
write(*,*) '* JOREK2 : Penning trap test          *'
write(*,*) '***************************************'

call initialise_parameters(my_id, "__NO_FILENAME__")

! Initialize fields in JOREK units
t_norm  = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds
qom     = real(charge) * el_chg / (mass * atomic_mass_unit)
B0      = omega_b/qom ! In T
Phi0    = epsilon*omega_e**2/qom/2.d0*t_norm ! In JOREK units: E_SI*t_norm
v0      = v0_SI * t_norm

! Only 1 toroidal mode (n=0)
mode(1) = 0

! Create a grid according to the variables present in the input file
if ((n_R > 0) .and. (n_Z > 0) .and. (n_radial > 0)) then
  call grid_bezier_square_polar(n_R, n_Z, n_radial, R_begin, R_end, Z_begin, Z_end, R_geo,   &
    Z_geo, amin, fbnd, fpsi, mf, .true., node_list, element_list)
else if ((n_R > 0) .and. (n_Z > 0) ) then
  call grid_bezier_square(n_R, n_Z, R_begin, R_end, Z_begin, Z_end, .true., node_list,       &
    element_list)
else if ((n_radial > 0) .and. (n_pol > 0) ) then
  call grid_polar_bezier(R_geo, Z_geo, amin, 0.d0, 0.d0, fbnd, fpsi, mf, n_radial, n_pol,    &
    node_list, element_list)
else
  write(*,*) 'FATAL : no valid combination of grid-sizes specified'
  stop
end if
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

  ! E = -Grad F0*U, U = Phi0(R^2 - 2 Z^2) (see model001/initial_conditions.f90 for reference)
  node_list%node(i)%values(1,1,2) = Phi0*(R**2-2.d0*Z**2)
  node_list%node(i)%values(1,2,2) = Phi0*(2.d0*R*R_s - 4.d0*Z*Z_s)
  node_list%node(i)%values(1,3,2) = Phi0*(2.d0*R*R_t - 4.d0*Z*Z_t)
  node_list%node(i)%values(1,4,2) = Phi0*(2.d0*R_s*R_t - 4.d0*Z_s*Z_t &
                                        + 2.d0*R*R_st  - 4.d0*Z*Z_st) ! U_RZ = 0

enddo
! The electric potential is F0*U, so set F0 to 1
F0 = 1.d0

! Write a restart file containing the grid
write(*,*) "INFO: Exporting grid to jorek_restart.rst"
call export_binary_restart(node_list,element_list,'jorek_restart.rst',0)

! Initialize the particle list
particle_list%n_particles = 1
allocate(particle_list%particle(particle_list%n_particles))

do i=1,size(tstep_n)
if (tstep_n(i) .le. 1.1) cycle ! Skip anything below 1 (the default value if not entered)
  write(*,*) "WARNING: overriding nstep_n specified"
  nstep_n(i) = floor(1.6e8 / tstep_n(i)) ! Override nstep_n
  call reset_particle
  ! Calculate the correct velocity for a half-step backwards to obtain second-order convergence
  E = [-2.d0*Phi0*x0(1),0.d0,0.d0]
  B = [0.d0, B0, 0.d0]
  B2 = B0**2
  f = -qom*tstep_n(i)*0.25d0*t_norm ! Perform the initial half-step without tan correction

  v = v0 + f*E
  v = (v + 2.d0*f/(1.d0+f**2*B2) * (cross_product(v,B) - f*v*B2 + f*B*dot_product(v,B)))
  v = v + f*E

  particle_list%particle(1)%v = v
  call update_particles(my_id, particle_list, tstep_n(i), nstep_n(i), 0.d0)

  ! Check position against analytical result
  x_e = particle_list%particle(1)%x
  x_a = analytical_trajectory(tstep_n(i)*real(nstep_n(i)))
  err_norm = sqrt(x_e(1)**2+x_a(1)**2 - 2.d0*x_e(1)*x_a(1)*(cos(x_e(3))*cos(x_a(3))+sin(x_e(3))*sin(x_a(3))))

  ! Exit the test if the error is too large
  ! Norm error scales as dt^2
  err_ref = 3.14084d-8*tstep_n(i)**2
  write(*,*) "RESULT: n_pol= ", n_radial, " dt= ", tstep_n(i), " error= ", err_norm, " reference= ", err_ref
  if (isnan(err_norm) .or. err_norm .gt. err_ref*1.2d0) then
    write(*,*) "Penning test failed"
    stop 1
  endif
enddo
write(*,*) "Tests successfull at all timestep sizes"

contains
  !> This function calculates the analytical trajectory of a particle in a penning trap
  !!
  !! The parameters epsilon, omega_e, omega_b, x0, v0 are used and the solution is
  !! returned in RZPhi coordinates.
  pure function analytical_trajectory(t_jorek) result(x)
    real*8, intent(in) :: t_jorek !< The time (in JOREK units) at which to calculate the solution value
    real*8             :: x(3) !< The position of the particle at time t (in RZPhi)

    ! Internal variables
    real*8 :: t !< The time (in SI units) at which to calculate the solution value
    real*8 :: omega_plus, omega_minus
    real*8 :: R_plus, R_minus
    real*8 :: T_plus, T_minus
    real*8 :: omega
    complex(kind=8) :: w
    real*8 :: x0_xyz(3), v0_xyz(3)

    t = t_jorek * t_norm ! t_norm is defined in the program above
    ! Initialize in xyz coordinates
    x0_xyz = RZPhiToXYZ(x0)
    v0_xyz = RZPhiToXYZ(v0_SI)

    ! Some initialization
    omega = sqrt(-2.d0*epsilon)*omega_e
    omega_plus  = 0.5d0*(omega_b + sqrt(omega_b**2 + 4.d0*epsilon*omega_e**2))
    omega_minus = 0.5d0*(omega_b - sqrt(omega_b**2 + 4.d0*epsilon*omega_e**2))
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
    x(3) = -atan2(dimag(w), dble(w)) ! Phi = -atan2(y,x)
  end function analytical_trajectory

  !> Reset a particle to the initial conditions
  subroutine reset_particle() ! Uses the global variables
    call find_RZ(node_list,element_list,x0(1),x0(2),R,Z,ielm_out,s,t,ifail)
    if (ifail .ne. 0) then
      write(*,*) "CRITICAL: could not find initial particle in grid", &
          "Particle location: ", x0(1), x0(2), " R,z_out= ", R, Z
      stop 1
    endif
    particle_list%particle(1)%st = [s,t]
    particle_list%particle(1)%i_elm = ielm_out
    particle_list%particle(1)%x = [R,Z,x0(3)]
    particle_list%particle(1)%v = v0
    particle_list%particle(1)%q = charge
    particle_list%particle(1)%mass = mass
    particle_list%particle(1)%weight = 1.d0
    particle_list%particle(1)%lost = .false.
  end subroutine reset_particle
end program penning
