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

! Testing parameters
integer, parameter :: n_dt = 16
integer            :: n_steps

! Penning trap parameters
real*8, parameter :: omega_e = 4.9
real*8, parameter :: omega_b = 25.0
real*8, parameter :: qom     = 1.0
real*8, parameter :: epsilon = -1.0
real*8, parameter :: x0(3)   = [10,0,0] ! in RZPhi
real*8, parameter :: v0(3)   = [100,20,0] ! in RZPhi = [100,0,100] in xyz
real*8, parameter :: b_z     = omega_b*qom
! The particle remains between +- 3 in Z and 8 and 13 in R for these parameters

! Local variables
integer :: i
real*8  :: stats(4,n_dt)
real*8  :: x_cart(3)
real*8  :: t_begin, t_end
real*8  :: dt

! Fake MPI presence
integer, parameter :: my_id = 0


write(*,*) '***************************************'
write(*,*) '* JOREK2 : Penning trap test          *'
write(*,*) '***************************************'

call initialise_parameters(my_id, "__NO_FILENAME__")

! Only 1 toroidal mode (n=0)
mode(1) = 0

! TODO replace this with a better grid
call tr_resetfile()
element_list%n_elements      = 0
bnd_elm_list%n_bnd_elements  = 0
node_list%n_nodes            = 0

call define_boundary
call grid_polar_bezier(R_geo, Z_geo, amin, 0.d0, 0.d0, fbnd, fpsi, mf, n_radial, n_pol, node_list, element_list)

!call import_binary_restart(node_list,element_list, 'jorek_restart.rst', rst_format, ierr)

!call initialise_basis                              ! define the basis functions at the Gaussian points

!call update_neighbours(element_list,node_list)     ! update neighbour information in the element_list

!$omp parallel do default(private) &
!$omp   shared(stats) &
!$omp   schedule(dynamic)
! Run in parallel for all different timesteps sizes dt
dt_loop: do i=n_dt,1,-1 ! ordering is reverse for better thread scheduling
  dt = 10.d0**(-i*0.25)
  n_steps = int(t_end/dt) ! Floors, so we are less than dt away from t_end
  call cpu_time(t_begin)
  
  ! Initialize the particle

  ! Run for n_steps with timestep dt
  ! Run the extra (on average) half-step needed to arrive at t=t_end

  call cpu_time(t_end)
  

  stats(1,i) = dt
  stats(2,i) = t_end-t_begin
  ! Convert to cartesian coordinates for analysis
  !call RZPhiToXYZ(x, x_cart)
  ! Compare the result to the analytical value
  !x_cart = x_cart - analytical_trajectory(t_end)
  stats(3,i) = norm2(x_cart(1:2))
  stats(4,i) = abs(x_cart(3))
end do dt_loop
!$omp end parallel do

contains
  !> This function calculates the analytical trajectory of a particle in a penning trap
  !!
  !! The parameters epsilon, omega_e, omega_b, x0, v0 are used and the solution is
  !! returned in cartesian coordinates.
  pure function analytical_trajectory(t) result(x)
    real*8, intent(in) :: t !< The time at which to calculate the solution value
    real*8             :: x(3) !< The (cartesian) position of the particle at time t

    ! Internal variables
    real*8    :: omega_plus, omega_minus
    real*8    :: R_plus, R_minus
    real*8    :: T_plus, T_minus
    real*8    :: omega
    complex(kind=8) :: w

    ! Some initialization
    omega = sqrt(-2*epsilon)*omega_e
    omega_plus  = 0.5*(omega_b + sqrt(omega_b**2 + 4*epsilon*omega_e**2))
    omega_minus = 0.5*(omega_b - sqrt(omega_b**2 + 4*epsilon*omega_e**2))
    R_minus = (omega_plus * x0(1) + v0(2))/(omega_plus - omega_minus)
    R_plus  = x0(1) - R_minus
    T_minus = (omega_plus * x0(2) - v0(1))/(omega_plus - omega_minus)
    T_plus  = x0(2) - T_minus

    ! Calculate the result in the x-y plane in terms of w = x + iy
    w = cmplx(R_plus ,T_plus , 8)*exp(cmplx(0.d0,-omega_plus*t, 8)) + &
        cmplx(R_minus,T_minus, 8)*exp(cmplx(0.d0,-omega_minus*t, 8))

    x(1) = real(w)
    x(2) = aimag(w)
    x(3) = x0(3)*cos(omega*t)+v0(3)*sin(omega*t)/omega
  end function analytical_trajectory

  !> This function returns the quadrupole field in cylindrical coordinates
  pure function E_quadrupole(RZPhi) result(E)
    real*8, intent(in) :: RZPhi(3)
    real*8             :: E(3)
    E(1) = RZPhi(1)
    E(2) = -2.d0*RZPhi(2)
    E(3) = 0
    E = -epsilon*omega_e**2*qom*E ! epsilon and omega_e are obtained from the program
  end function E_quadrupole
end
