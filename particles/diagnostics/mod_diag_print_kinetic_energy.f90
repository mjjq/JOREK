!> Module containing action to print diagnostics on the kinetic energy
!> of all particles in the simulation.
module mod_diag_print_kinetic_energy
use mod_particle_sim
use mod_particle_types
use mod_event
implicit none

type, extends(action) :: diag_print_kinetic_energy
contains
  procedure :: do => do_print_kinetic_energy
end type diag_print_kinetic_energy

contains

!> Print some statistics on the kinetic energy of all particles in the simulation (with MPI and openmp support).
!> Writing analysis scripts in this way aids reusability, as they can also be called inline in a simulation then.
subroutine do_print_kinetic_energy(this, sim, ev)
  use mod_mpi_stats
  class(diag_print_kinetic_energy), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), optional :: ev
  integer :: i, j
  real*8, dimension(:), allocatable :: tmp

  ! For each of the groups, calculate the kinetic energy with the
  ! appropriate method
  do i=1,size(sim%groups)
    select type (p => sim%groups(i)%particles)
    type is (particle_kinetic_leapfrog)
      tmp = boris_kinetic_energy(p)
    type is (particle_kinetic_relativistic)
      tmp = particle_relativistic_kinetic_energy(p,sim%groups(i)%mass)
    class default
      write(*,*) "do_print_kinetic_energy not implemented for this particle type"
      return
    end select

    write(*,"(f14.7,A,i3,A,5g14.7)") sim%time, "Group", i, " kinetic energy min/mean/max/stddev/sum ", mpi_stats_list(tmp)
    deallocate(tmp)
  end do
end subroutine do_print_kinetic_energy


!> The impure keyword is used for openmp support in fortran 2008.
!> It requires ifort 16.0 or up. Remove the simd declaration if needed
!> The impure keyword should be [unnecessary](https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/591902)
!> starting with implementation of the openmp 4.1 spec
!>
!> if there are any problems, it can be removed along with the simd instruction
!> The speed improvements of this still have to be tested
!>
!> Probably this is a bug because if a function has to be used
!> within a SIMD region it must be declared to be SIMD compatible 
!> with the pragma: !$omp declare simd (options)
impure elemental function boris_kinetic_energy(particle) result(energy)
  class(particle_kinetic_leapfrog), intent(in) :: particle
  real*8 :: energy
  energy = dot_product(particle%v, particle%v)
end function boris_kinetic_energy

!> This function computes the kinetic energy of a realtivistic
!> particle in eV which is equal to:
!> E_{kin} = c\[\sqrt{\(mc\)^2 + p^2}-mc\]/q_proton
!>
!> The impure elemental declaration is preserved for consistency with
!> boris_kinetic_energy() function despite I hardly belive this will make any
!> difference in terms of performances due to the lack of data continuity.
!>
!> In addition, if it has to be used within a SIMD region it must be declare as
!> SIMD compatible function using the pragama: !$omp declare simd (options)
!> inputs:
!>   particle: (particle_kinetic_relativistic) a relativistic kinetic particle type
!>   mass:     (real8) the particle mass
!> outputs:
!>   energy: (real8) the particle kinetic energy in eV
impure elemental function particle_relativistic_kinetic_energy(particle,mass) result(energy)
  use constants, only: SPEED_OF_LIGTH,ATOMIC_MASS_UNIT,EL_CHG
  ! declare input variables
  class(particle_kinetic_relativistic),intent(in) :: particle !< relativistic particle
  real(kind=8),intent(in) :: mass !< particle mass
  ! declare output variables
  real(kind=8) :: energy !< particle kinetic energy

  ! compute the relativistic particle kinetic energy
  energy = ATOMIC_MASS_UNIT*SPEED_OF_LIGTH*(sqrt((mass*SPEED_OF_LIGTH)&
  *(mass*SPEED_OF_LIGTH)+dot_product(particle%p,particle%p))-&
  mass*SPEED_OF_LIGTH)/EL_CHG  

end function

end module mod_diag_print_kinetic_energy

