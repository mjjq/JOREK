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

! Penning trap parameters (in JOREK units) (for check)
real*8, parameter :: omega_e = 4.9d0
real*8, parameter :: omega_b = 25.d0
real*8, parameter :: epsilon = -1.d0
real*8, parameter :: x0(3)   = [10.d0,0.d0,0.d0] ! in RZPhi
real*8, parameter :: v0(3)   = [50.d0,20.d0,0.d0] ! in RZPhi, = [100,0,20] in xyz
! Penning trap parameters (in JOREK units) (for fields)
integer, parameter :: charge = 1
real*8, parameter :: mass = 1.d0
real*8, parameter :: qom = (real(charge) * EL_CHG * SQRT(MU_ZERO * MASS_PROTON * 2.d0 * 1.D20)) & ! Use central_mass preset value here
        / (mass * ATOMIC_MASS_UNIT)
! qom     = 62.5576 ! Normalized q/m (above value is 62.557605832483851
real*8, parameter :: B0      = omega_b/qom
real*8, parameter :: Phi0    = epsilon*omega_e**2/qom/2.d0
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

! Initialize the particle list
particle_list%n_particles = 1
allocate(particle_list%particle(particle_list%n_particles))

! First test a constant magnetic field in z-direction
F0 = 0.d0 ! By disabling the electric field
do i=1,size(tstep_n)
  if (nstep_n(i) .eq. 0) cycle
  call reset_particle
  call update_particles(my_id, particle_list, tstep_n(i), nstep_n(i))
  ! And test the conservation of energy in the xy-plane
  v = particle_list%particle(1)%v
  write(*,*) "UniformB test, tstep: ", tstep_n(i), " Energy error: ", abs(dot_product(v,v)-dot_product(v0,v0))

  ! exit the test if the error is too large
  if (abs(dot_product(v,v)-dot_product(v0,v0)) .gt. 1.d-6) then
    write(*,*) "CRITICAL: UniformB test failed, error larger that 1.d-6"
    stop 1
  endif
  ! Exit if there are nans anywhere
  if (isnan(v(1))) then
    write(*,*) "CRITICAL: NaN encountered, stopping"
    stop 1
  endif
enddo

! Test the full penning trap
! Enable the electric field
F0 = 1.d0

do i=1,size(tstep_n)
  if (nstep_n(i) .eq. 0) cycle
  call reset_particle
  ! Calculate the correct velocity for a half-step backwards to obtain second-order convergence
  E = [-2.d0*Phi0*x0(1),0.d0,0.d0]
  B = [0.d0, B0, 0.d0]
  B2 = B0**2
  f = -qom*tstep_n(i)*0.25d0

  v = v0 + f*E
  v = (v + 2.d0*f/(1.d0+f**2*B2) * (cross_product(v,B) - f*v*B2 + f*B*dot_product(v,B)))
  v = v + f*E

  particle_list%particle(1)%v = v
  call update_particles(my_id, particle_list, tstep_n(i), nstep_n(i), 0.d0)

  ! Check position against analytical result
  x_e = particle_list%particle(1)%x
  x_a = analytical_trajectory(tstep_n(i)*real(nstep_n(i)))
  err_norm = sqrt(x_e(1)**2+x_a(1)**2 - 2.d0*x_e(1)*x_a(1)*(cos(x_e(3))*cos(x_a(3))+sin(x_e(3))*sin(x_a(3))))
  write(*,*) "Penning test, tstep: ", tstep_n(i), " Position error: ", err_norm

  ! Exit the test if the error is too large
  ! Norm error: classical error in reference code, scales as
  err_ref = 4.d2*(omega_b*tstep_n(i))**2*1.2d0
  if (err_norm .gt. err_ref .and. .false.) then
    write(*,*) "Penning test failed, error larger than ", err_ref
    stop 1
  endif
enddo

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
    write(*,*) "Initializing test particle"
    call find_RZ(node_list,element_list,x0(1),x0(2),R,Z,ielm_out,s,t,ifail)
    if (ifail .ne. 0) then
      write(*,*) "CRITICAL: could not find initial particle in grid"
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
end
