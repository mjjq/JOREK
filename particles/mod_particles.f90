module mod_particles

  use parameters

  type type_particle

    real*8  :: x(3)             !< particle position in real space (R,Z,phi)
    real*8  :: st(n_dim)        !< particle position in the finite element (i_elm)
    real*8  :: v(3)             !< particle velocity in (R,Z,phi)
    real*8  :: q                !< charge
    real*8  :: mass             !< mass
    real*8  :: weight           !< weight (i.e. number of particles)
    integer :: i_elm            !< the index of the element containing the particle in the element_list
    logical :: lost             !< particle is active or lost

  end type type_particle

  type type_particle_list

    integer                                         :: n_particles                 !< the number of particles in the list
    type (type_particle), dimension(:), allocatable :: particle   !< an allocatable list of particles

  end type type_particle_list

contains

function cross_product(a,b)
!----------------------------------------
! input  :  (a_R,a_Z,A_phi), (b_R,b_Z,b_phi)
! output :  (a x b)_(R,Z,phi)
!----------------------------------------
implicit none
real*8 :: cross_product(3)
real*8 :: a(3), b(3)

cross_product(1) = a(2)*b(3) - a(3)*b(2)
cross_product(2) = a(3)*b(1) - a(1)*b(3)
cross_product(3) = a(1)*b(2) - a(2)*b(1)

return
end function cross_product

subroutine export_particles(particle_list,particle_file)

implicit none

type (type_particle_list) :: particle_list
type (type_particle)      :: particle
character*(*),intent(in)  :: particle_file
integer                   :: i

write(*,*) '***********************************'
write(*,*) '*       export particles          *'

open(22,file=particle_file,form='unformatted', status='replace', action='write')

write(22) particle_list%n_particles

do i=1, particle_list%n_particles
  write(22) particle_list%particle(i)
enddo

close(22)

write(*,*) '*       particles exported        *'
write(*,*) '***********************************'

return
endsubroutine export_particles

subroutine import_particles(particle_list)

implicit none

type (type_particle_list) :: particle_list
type (type_particle)      :: particle
integer :: i

write(*,*) '***********************************'
write(*,*) '*       import particles          *'

open(22,file='particles.rst')

read(22) particle_list%n_particles

do i=1, particle_list%n_particles
  read(22) particle_list%particle(i)
enddo

close(22)

write(*,*) '*       particles imported        *'
write(*,*) '***********************************'

return
endsubroutine import_particles

end module mod_particles
