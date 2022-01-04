!> mod_particle_types_test contains variables and 
!> procedure for testing the methods coded in
!> mod_particle_types
module mod_particle_types_test
use fruit
use mod_particle_types
implicit none

private
public :: run_fruit_particle_types

!> Variables ------------------------------------
integer,parameter :: n_particles=1000
integer,dimension(2),parameter :: q_interval=(/1,100/)
integer,dimension(2),parameter :: i_elm_interval=(/1,10000000/)
integer,dimension(2),parameter :: i_life_interval=(/1,10000000/)
real*8,dimension(2),parameter  :: t_birth_interval=(/0.d0,3.45d4/)
real*8,dimension(2),parameter  :: st_interval=(/0.d0,1.d0/)
real*8,dimension(2),parameter  :: mass_interval=(/5.485d-4,124.d0/)
real*8,dimension(2),parameter  :: v_interval=(/-6.75d3,8.45d3/)
real*8,dimension(2),parameter  :: Ekin_interval=(/0.d0,1.d7/)
real*8,dimension(2),parameter  :: mu_interval=(/0.d0,1.d-5/)
real*8,dimension(2),parameter  :: Bnorm_interval=(/0.d0,1.4d1/)
real*8,dimension(2),parameter  :: weight_interval=(/0.d0,1.d3/)
real*8,dimension(3),parameter  :: x_lowbnd=(/-5.d2,-1.d2,1.d0/)
real*8,dimension(3),parameter  :: x_uppbnd=(/7.d2,2.d2,4.d2/)
real*8,dimension(3),parameter  :: vp3d_lowbnd=(/-1.25d3,-7.5d2,-8.d1/)
real*8,dimension(3),parameter  :: vp3d_uppbnd=(/7.5d1,2.35d2,4.85d3/)
real*8,dimension(3),parameter  :: ABE_lowbnd=(/-2.67d0,-9.85d0,3.5d-1/)
real*8,dimension(3),parameter  :: ABE_uppbnd=(/7.8d-1,2.35d0,5.67d0/)
class(particle_base),dimension(:),allocatable :: particle_fieldline_list
class(particle_base),dimension(:),allocatable :: particle_gc_list
class(particle_base),dimension(:),allocatable :: particle_gc_vpar_list
class(particle_base),dimension(:),allocatable :: particle_gc_Qin_list
class(particle_base),dimension(:),allocatable :: particle_kinetic_list
class(particle_base),dimension(:),allocatable :: particle_kinetic_leapfrog_list
class(particle_base),dimension(:),allocatable :: particle_kinetic_relativistic_list
class(particle_base),dimension(:),allocatable :: particle_gc_relativistic_list
!> Interfaces -----------------------------------

contains

!> Fruit basket ---------------------------------
!> fruit basket containing all set-up, test and
!> tear-down procedures
subroutine run_fruit_particle_types()
  implicit none
  
  write(*,'(/A)') "  ... setting-up: particle types tests"
  call setup
  write(*,'(/A)') "  ... running: particle types tests"
  write(*,'(/A)') "  ... tearing-up: particle types tests"
  call teardown
end subroutine run_fruit_particle_types

!> Set-up and tear-down -------------------------
!> set-up particle types test features
subroutine setup()
  implicit none

  !> allocate particle arrays
  allocate(particle_fieldline::particle_fieldline_list(n_particles))
  allocate(particle_gc::particle_gc_list(n_particles))
  allocate(particle_gc_vpar::particle_gc_vpar_list(n_particles))
  allocate(particle_gc_Qin::particle_gc_Qin_list(n_particles))
  allocate(particle_kinetic::particle_kinetic_list(n_particles))
  allocate(particle_kinetic_leapfrog::particle_kinetic_leapfrog_list(n_particles))
  allocate(particle_kinetic_relativistic::particle_kinetic_relativistic_list(n_particles))
  allocate(particle_gc_relativistic::particle_gc_relativistic_list(n_particles))

end subroutine setup

!> tear-down particle types test features
subroutine teardown()
  implicit none
  !> deallocate particle arrays
  deallocate(particle_fieldline_list)
  deallocate(particle_gc_list)
  deallocate(particle_gc_vpar_list)
  deallocate(particle_gc_Qin_list)
  deallocate(particle_kinetic_list)
  deallocate(particle_kinetic_leapfrog_list)
  deallocate(particle_kinetic_relativistic_list)
  deallocate(particle_gc_relativistic_list)
end subroutine teardown

!> Tests ----------------------------------------
!> Tools ----------------------------------------
!>-----------------------------------------------
end module mod_particle_types_test
