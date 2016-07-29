!> Mod_particle_type contains the default particle type, type_particle
module mod_particle_base
  type, abstract :: particle_base
    real*8    :: x(3)             !< particle position in real space
    real*4    :: mass             !< mass [atomic mass units]
    real*4    :: weight           !< weight (i.e. number of particles)
    real*8    :: st(2)            !< JOREK integration: particle position in the element
    integer*4 :: i_elm            !< JOREK integration: index in element_list
    integer*1 :: q                !< charge [e]
    integer*1 :: label            !< Particle type number (i in species(i))
    logical*1 :: lost             !< particle is active or lost
  end type particle_base
end module mod_particle_base
