program project_particle_restart
use particle_tracer
use mod_particle_io
implicit none

type(projection)                 :: proj
type(event)                      :: fieldreader

! Start up MPI, jorek
call sim%initialize(num_groups=1)

! Read a sim
call read_simulation_hdf5(sim, 'part_restart.h5')

! Set up the field reader
fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
call with(sim, fieldreader)

proj = new_projection(sim%fields%node_list, sim%fields%element_list, smoothing=0d-3, smoothing2=0d0, &
                      to_h5=.true.,basename='proj')

call with(sim, proj)
end program project_particle_restart
