!> Testing conversion of particles from one type into another
!>
!> Compile with `make test_convert_particles`
!> Run with `./test_convert_particles < JOREK_namelist`

program test_convert_particles

use particle_tracer
use mod_particle_io
use mod_particle_diagnostics
use mod_fields_linear
use mod_gc_relativistic                    
implicit none

! Set up the simulation variables
integer(kind=4)                     :: n_part, ifail
type(diag_print_kinetic_energy)     :: print_kinetic_energy
type(write_particle_diagnostics)    :: diag
real*8                              :: gyro_angle
type(particle_gc_relativistic)      :: particle_in
!type(particle_kinetic_relativistic) :: particle_out
type(particle_gc)                   :: particle_out
real(kind=8)                        :: psi, U
real(kind=8),dimension(3)           :: E, B

call sim%initialize(num_groups=1)

! Set up the diagnostics output
diag = write_particle_diagnostics(filename='diag.h5',only=[1,2,6,12,13,14]) ! store total and kinetic energies, p_phi, phi, R, Z

sim%time = 1.d-7  ! start time 

! Allocate a group and particle(s) of type particle_gc_relativistic
n_part = 1
allocate(particle_gc_relativistic::sim%groups(1)%particles(n_part))
sim%groups(1)%mass = 5.4857990907016d-4 !< particle mass in AMU

! Set events to write output data and stop the simulation
events = [event(read_jorek_fields_interp_linear(i=-1)), & 
          event(diag,start=sim%time,step=1d-8),         &
          event(stop_action(),start=sim%time+2.d-8)]

! Run first event to read the JOREK fields
call with(sim, events, at=0.d0)

select type (p=>sim%groups(1)%particles(1))
type is (particle_gc_relativistic)
  p%x = [3.d0,0.d0,0.d0]
  call find_RZ(sim%fields%node_list, sim%fields%element_list, &
               p%x(1), p%x(2), & ! inputs
               p%x(1), p%x(2), p%i_elm, p%st(1), p%st(2), ifail) ! outputs
  p%p = [1.6d6,1.d11]
  p%q = -1
end select

particle_in = sim%groups(1)%particles(1)

! Print particle information before conversion
write(*,*) 'BEFORE CONVERSION'
write(*,*) 'x, y, z: ', particle_in%x(1), particle_in%x(2), particle_in%x(3)
write(*,*) 'p(1), p(2): ', particle_in%p(1), particle_in%p(2)
call print_kinetic_energy%do(sim,events(1))  !< print particle kinetic energy

! Convert particle
gyro_angle = 0.
call sim%fields%calc_EBpsiU(sim%time, sim%groups(1)%particles(1)%i_elm, sim%groups(1)%particles(1)%st, &              
                            sim%groups(1)%particles(1)%x(3), E, B, psi, U)
call relativistic_gc_to_particle(sim%fields%node_list, sim%fields%element_list, &
                                 particle_in, particle_out, sim%groups(1)%mass, &
				 B, gyro_angle)

! Print particle information after conversion
write(*,*) 'AFTER CONVERSION'
write(*,*) 'x, y, z: ', particle_out%x(1), particle_out%x(2), particle_out%x(3)
!write(*,*) 'p(1), p(2), p(3): ', particle_out%p(1), particle_out%p(2), particle_out%p(3)
write(*,*) 'E, mu: ', particle_out%E, particle_out%mu

! Convert back
!call relativistic_kinetic_to_particle(sim%fields%node_list, sim%fields%element_list,    &
!                                      particle_out, particle_in, sim%groups(1)%mass, B)
particle_in = gc_to_relativistic_gc(particle_out,sim%groups(1)%mass,B)

! Print particle information after converting back
write(*,*) 'AFTER CONVERTING BACK'
write(*,*) 'x, y, z: ', particle_in%x(1), particle_in%x(2), particle_in%x(3)
write(*,*) 'p(1), p(2): ', particle_in%p(1), particle_in%p(2)

! Finalize the simulation
call sim%finalize

end program test_convert_particles

