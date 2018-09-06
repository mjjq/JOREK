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
use particle_tracer
use mod_coordinate_transforms ! For solution of penning trap trajectory
use mod_parameters
use constants
use mod_fields_linear
use mod_export_restart
use mod_neighbours
use mod_find_rz_nearby
use mod_penning_case
use mod_penning_case_jorek

implicit none

type(jorek_fields_interp_linear) :: fields

! Define our particle list
type(particle_kinetic_leapfrog) :: particle

! Penning trap parameters (in SI units)
! The particle remains between +- 3 in Z and 8 and 13 in R for these parameters. set this in an input file

! Local variables
integer :: i,ielm_out,ifail,j,i_elm_old
real*8  :: R,R_s,R_t,R_st,Z,Z_s,Z_t,Z_st
real*8  :: s,t
real*8  :: x_a(3), x_e(3), err_norm, err_ref
real*8  :: rz_old(2), st_old(2)
! For the initial half-step
real*8  :: E(3), B(3), psi, U
real*8  :: t_norm, qom, B0, Phi0

! Fake MPI presence
integer, parameter :: my_id = 0


write(*,*) '***************************************'
write(*,*) '* JOREK2 : Penning trap test          *'
write(*,*) '***************************************'

call initialise_parameters(my_id, "__NO_FILENAME__")

t_norm  = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds
qom     = real(charge) * el_chg / (mass * atomic_mass_unit)
B0      = omega_b/qom ! In T
Phi0    = epsilon*omega_e**2/qom/2.d0*t_norm ! In JOREK units: E_SI*t_norm

fields%static = .true.
rst_hdf5 = 1

! Only 1 toroidal mode (n=0)
mode(1) = 0

! Create a grid according to the variables present in the input file
if ((n_R > 0) .and. (n_Z > 0) .and. (n_radial > 0)) then
  call grid_bezier_square_polar(n_R, n_Z, n_radial, R_begin, R_end, Z_begin, Z_end, R_geo,   &
    Z_geo, amin, fbnd, fpsi, mf, .true., fields%node_list, fields%element_list)
else if ((n_R > 0) .and. (n_Z > 0) ) then
  call grid_bezier_square(n_R, n_Z, R_begin, R_end, Z_begin, Z_end, .true., fields%node_list,       &
    fields%element_list)
else if ((n_radial > 0) .and. (n_pol > 0) ) then
  call grid_polar_bezier(R_geo, Z_geo, amin, 0.d0, 0.d0, fbnd, fpsi, mf, n_radial, n_pol,    &
    fields%node_list, fields%element_list)
else
  write(*,*) 'FATAL : no valid combination of grid-sizes specified'
  stop
end if
call update_neighbours(fields%node_list, fields%element_list)

call jorek_penning_fields(fields%node_list, fields%element_list)


! Write a restart file containing the grid
write(*,*) "INFO: Exporting grid to jorek_restart.h5"
call export_restart(fields%node_list,fields%element_list,'jorek_restart')

! Initialize the particle list

do i=1,size(tstep_n)
if (tstep_n(i) .le. 1.1) cycle ! Skip anything below 1 (the default value if not entered)
  write(*,*) "WARNING: overriding nstep_n specified"
  nstep_n(i) = floor(1.6e8 / tstep_n(i)) ! Override nstep_n
  call reset_particle
  ! Calculate the correct velocity for a half-step backwards to obtain second-order convergence
  E = [-2.d0*Phi0*x0(1),0.d0,0.d0]/t_norm
  B = [0.d0, B0, 0.d0]
  call boris_initial_half_step_backwards_RZPhi(particle, real(mass,8), E, B, tstep_n(i)*t_norm)

  do j=1,nstep_n(i)
    call fields%calc_EBpsiU(0.d0, particle%i_elm, particle%st, particle%x(3), E, B, psi, U)
    rz_old    = particle%x(1:2)
    st_old    = particle%st
    i_elm_old = particle%i_elm
    call boris_push_cylindrical(particle, real(mass,8), E, B, tstep_n(i)*t_norm)
    call find_RZ_nearby(fields%node_list, fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
        particle%x(1), particle%x(2), particle%st(1), particle%st(2), particle%i_elm, ifail)
  end do

  ! Check position against analytical result
  x_e = particle%x
  x_a = penning_trajectory(tstep_n(i)*real(nstep_n(i)))
  err_norm = sqrt(x_e(1)**2+x_a(1)**2 - 2.d0*x_e(1)*x_a(1)*(cos(x_e(3))*cos(x_a(3))+sin(x_e(3))*sin(x_a(3))))

  ! Exit the test if the error is too large
  ! Norm error scales as dt^2
  err_ref = 3.14084d-8*tstep_n(i)**2
  if (n_pol .gt. 0) err_ref = max(err_ref, 2.d0*3.675d3*real(n_pol)**(-4)) ! only use this condition for polar grids, with huge margin (2x)
  write(*,"(A,i3,A,g12.4,A,g16.8,a,g16.8)") "RESULT: n_radial= ", n_radial, " dt= ", tstep_n(i), " error= ", err_norm, " reference= ", err_ref
  if (isnan(err_norm) .or. err_norm .gt. err_ref*1.2d0) then
    write(*,*) "Penning test failed"
    stop 1
  endif
enddo
write(*,*) "Tests successfull at all timestep sizes"

contains
  !> Reset a particle to the initial conditions
  subroutine reset_particle() ! Uses the global variables
    implicit none
    call find_RZ(fields%node_list,fields%element_list,x0(1),x0(2),R,Z,ielm_out,s,t,ifail)
    if (ifail .ne. 0) then
      write(*,*) "CRITICAL: could not find initial particle in grid", &
          "Particle location: ", x0(1), x0(2), " R,z_out= ", R, Z
      stop 1
    endif
    particle%st = [s,t]
    particle%i_elm = ielm_out
    particle%x = [R,Z,x0(3)]
    particle%v = v0
    particle%q = charge
    particle%weight = 1.d0
  end subroutine reset_particle
end program penning
