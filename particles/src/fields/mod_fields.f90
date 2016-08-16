module mod_fields
  type, abstract :: fields
    contains
    procedure(at_particle), deferred :: at_particle
  end type fields

  interface
    pure subroutine at_particle(this, particle, t, E, B, psi, U)
      use mod_particle_base
      import :: fields
      implicit none
      class(fields), intent(in)         :: this
      class(particle_base), intent(in)  :: particle
      real*8, intent(in)                :: t !< The current time
      real*8, dimension(3), intent(out) :: E, B
      real*8, intent(out)               :: psi, U
    end subroutine at_particle
  end interface
end module mod_fields
