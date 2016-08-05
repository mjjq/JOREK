module mod_fields
  type, abstract :: type_fields
    contains
    procedure(at_particle), deferred :: at_particle
  end type type_fields

  interface
    pure subroutine at_particle(this, particle, E, B, psi, U)
      use mod_particle_type
      import :: type_fields
      implicit none
      class(type_fields), intent(in)        :: this
      class(type_particle), intent(in)      :: particle
      real*8, dimension(3), intent(out)     :: E, B
      real*8, intent(out)                   :: psi, U
    end subroutine at_particle
  end interface
end module mod_fields
