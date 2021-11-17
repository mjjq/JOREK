!> the mod_particle_type_openmp module tests the
!> the compatibility of coupling strategies
!> between the template-like particle type and
!> the openmp directives
module mod_particle_type_openmp
use fruit
use mod_particle_types, only: particle_base
use mod_particle_types, only: particle_kinetic_relativistic
use mod_particle_types, only: particle_gc
use mod_particle_sim, only: particle_sim
implicit none

private
public :: run_fruit_particle_type_openmp

! Parameters ----------------------------------------
!> variables *_rel_fo: relativistic full orbits
!> variables *_gc: guiding center
integer :: n_fields_rel_fo=8
integer :: n_fields_gc=7
integer*1,parameter :: q_rel_fo=1
integer*1,parameter :: q_gc=3
integer,parameter :: i_elm_rel_fo=512
integer,parameter :: i_elm_gc=1024
real*8,parameter :: E_gc=1.d2
real*8,parameter :: mu_gc=1.d-10
real*8,dimension(3),parameter :: x_rel_fo=(/3.2d0,1.d0,7.d-1/)
real*8,dimension(3),parameter :: x_gc=(/1.3d0,-5.d-2,3.d0/)
real*8,dimension(3),parameter :: p_rel_fo=(/2.d0,1.d-5,4.d1/)
real*8,dimension(2),parameter :: st_rel_fo=(/1.d-1,9.d-1/)
real*8,dimension(2),parameter :: st_gc=(/5.d-2,5.5d-1/)

! Variables -----------------------------------------
type(particle_sim) :: sim !< constant p-type per group
integer :: n_groups,n_particles
!> solutions for particle lists
integer*1,dimension(:),allocatable  :: int1_rel_fo_sol,int1_gc_sol
integer,dimension(:),allocatable  :: int_rel_fo_sol,int_gc_sol
real*8,dimension(:,:),allocatable :: real8_rel_fo_sol,real8_gc_sol

contains
! Test baskets --------------------------------------
subroutine run_fruit_particle_type_openmp()
  implicit none

  write(*,"(/A)") "  ... set-up particle type openmp tests"
  call setup
  write(*,"(/A)") "  ... run particle type openmp tests"
  write(*,"(/A)") "  ... tear-down particle type openmp tests"
  call teardown
  
end subroutine run_fruit_particle_type_openmp

! Set-up tear-down ----------------------------------
!> initialise tests
subroutine setup()
  implicit none

  !> variables
  integer :: ii

  !> allocate groups
  allocate(sim%groups(n_groups))
  !> groups are initialised one by one because the particle
  !> type must be specified
  allocate(particle_kinetic_relativistic::sim%groups(1)%particles(n_particles))
  allocate(particle_gc::sim%groups(2)%particles(n_particles))
  allocate(int1_rel_fo_sol(n_particles)); 
  allocate(int1_gc_sol(n_particles));
  allocate(int_rel_fo_sol(n_particles)); 
  allocate(int_gc_sol(n_particles));
  allocate(real8_rel_fo_sol(n_fields_rel_fo,n_particles)); 
  allocate(real8_gc_sol(n_fields_gc,n_particles));

  !> initialisation
  !> fill the constant type particle arrays
  int1_rel_fo_sol = q_rel_fo; int1_gc_sol = q_gc
  int_rel_fo_sol = i_elm_rel_fo; int_gc_sol = i_elm_gc;
  do ii=1,n_particles
    !> copy relativistic full orbit double variables
    real8_rel_fo_sol(1:3,ii) = x_rel_fo
    real8_rel_fo_sol(4:5,ii) = st_rel_fo
    real8_rel_fo_sol(6:8,ii) = p_rel_fo
    select type (p=>sim%groups(1)%particles(ii))
    type is (particle_kinetic_relativistic)
      p%x = x_rel_fo; p%st = st_rel_fo; p%i_elm = i_elm_rel_fo;
      p%p = p_rel_fo; p%q = q_rel_fo;
    end select
  enddo
  do ii=1,n_particles
    !> copy gc double varibales
    real8_gc_sol(1:3,ii) = x_gc
    real8_gc_sol(4:5,ii) = st_gc
    real8_gc_sol(6:7,ii) = (/E_gc,mu_gc/)   
    select type (p=>sim%groups(2)%particles(ii))
    type is (particle_gc)
      p%x = x_gc; p%st = st_gc; p%i_elm = i_elm_gc;
      p%E = E_gc; p%mu = mu_gc; p%q = q_gc;
    end select
  enddo

end subroutine setup

!> clean up simulation varibales
subroutine teardown()
  implicit none
  !> deallocate simulations and all their allocatables
  deallocate(sim%groups)
  deallocate(int1_rel_fo_sol); deallocate(int1_gc_sol);
  deallocate(int_rel_fo_sol);  deallocate(int_gc_sol);
  deallocate(real8_rel_fo_sol); 
  deallocate(real8_gc_sol);
end subroutine teardown

! Tests ---------------------------------------------

end module mod_particle_type_openmp
