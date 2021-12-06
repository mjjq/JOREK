!> The mod_light_vertices module contains variables
!> and procedures used for defining and defining actions
!> of the light points
module mod_light_vertices
use mod_vertices, only: vertices
implcit none

private
public :: light_vertices

!> Variables --------------------------------------------
type,abstract,extends(vertices) :: light_vertices
  contains
  procedure,nopass :: initialise_jorek_sim
  procedure(init_lights_parts),deferred,pass(light_vertices)  :: init_lights_from_particles
  procedure(direct_funct),deferred,pass(light_vertices)       :: directionality_funct
end type light_vertices

!> Interfaces -------------------------------------------
interface
  !> computes and store the coordinates and properties of 
  !> lights from particle simulations
  !> inputs:
  !>   light_vert:   (light_vertices) empty light vertices
  !>   sim_particle: (particle_sim) initialised particle simulation
  !> outputs:
  !>   light_vert: (light_vertices) filled light vertices
  subroutine init_lights_parts(light_vert,sim_particle)
    use mod_particle_sim, only: particle_sim
    implicit none
    !> inputs-outputs
    class(light_vertices),intent(inout) :: light_vert
    !> inputs
    class(particle_sim),intent(in) :: sim_particle
  end subroutine init_lights_parts

  !> computes the directionality function for a given point
  !> in space (cartesian coordinate) and a given light for
  !> all spectra wavelengths
  !> inputs:
  !>   light_vert: (light_vertices) initialised light vertices
  !>   spectra:    (spectrum_base) initilises spectra
  !>   light_id:   (integer) id of the light to use
  !>   x_shaded:   (real8)(n_x) point illuminated by the light
  !> outputs:
  !>   light_vert: (light_vertices) initialised light vertices
  !>   light_dstb: (real8)(n_points,n_spectra) light intensity distribution
  !>               from the light_id light to the x_shaded point for all
  !>               spectra points and all spectra
  subroutine direct_funct(light_vert,spectra,light_id,x_shaded,light_dstb)
    use mod_spectra,      only: spectrum_base
    use mod_particle_sim, only: particle_sim
    implicit none
    !> inputs-outpus
    class(light_vertices),intent(inout) :: light_vert
    !> inputs
    class(spectrym_base),intent(in) :: spectra
    integer,intent(in) :: light_id
    real*8,dimension(light_vert%n_x),intent(in) :: x_shaded
    !> outputs
    real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_dstb
  end subroutine direct_funct
end interface

contains

!> Procedures -------------------------------------------
!> initialise and load a particle simulations.:
!> inputs:
!>   sim_particle:      (particle_sim) particle simulation to be loaded
!>   n_groups:          (integer) number of particle sim groups
!>   rank:              (integer) MPI task rank
!>   n_tasks:           (integer) number of MPI tasks
!>   particle_filename: (character) name of the particle file to be loaded
!>   use_fields_hb_in:  (logical)(optional) if present and true the 
!>                      the hermite birkhoff interpolant is used.
!>                      Default: .false.
!>   jorek_id_in:       (integer) id of the jorek restart file to read
!>                      set to -1 to not include. Default: 0
!>   rst_format_in:     (integer) format of the restart file is .rst or .hdf5
!>                      Default: 0
!>   basename_in:       (charechter)(80) basename of the jorek restart file
!>                      Default: jorek
!> outputs:
!>   sim_particle: (particle_sim) loaded particle simulation
subroutine initialise_jorek_sim(sim_particle,n_groups,&
rank,n_tasks,particle_filename,use_fields_hb_in,jorek_id_in,&
rst_format_in,basename_in)
  use mod_fields,                 only: fields_base
  use mod_fields_linear,          only: jorek_fields_interp_linear
  use mod_fields_linear,          only: read_jorek_fields_interp_linear
  use mod_fields_hermite_birkhoff,only: jorek_fields_interp_hermite_birkhoff
  use mod_fields_hermite_birkhoff,only: read_jorek_fields_interp_hermite_birkhoff
  use mod_particle_sim,           only: particle_sim
  use mod_particle_io,            only: read_simulation_hdf5
  implicit none
  !> inputs-outputs
  type(particle_sim),intent(inout) :: sim_particle
  class(fields_base),intent(inout) :: jorek_fields
  !> inputs
  integer,intent(in)                    :: n_groups,rank,n_tasks
  character(len=*),intent(in)           :: particle_filename
  logical,intent(in),optional           :: use_fields_hb_in
  integer,intent(in),optional           :: jorek_id_in,rst_format_in
  character(len=80),intent(in),optional :: basename_in
  !> variables
  type(read_jorek_fields_interp_linear)           :: read_jorek_interp_lin
  type(read_jorek_fields_interp_hermite_birkhoff) :: read_jorek_interp_hb
  logical :: use_fields_hb
  integer :: jorek_id,rst_format
  character(len=80) :: basename

  !> se optional parameters
  use_fields_hb = .false.; jorek_id = 0;
  rst_format = 0; basename = "jorek";
  if(present(use_fields_hb_in)) use_fields_hb = use_fields_hb_in
  if(present(jorek_id_in)) jorek_id = jorek_id_in
  if(present(rst_format_in)) rst_format = rst_format_in
  if(present(basename_in)) basename = trim(basename_in)

  !> initialise and load particle simulation
  sim_particle%initialize(n_groups,.false.,rank,n_tasks)
  call read_simulation_hdf5(sim_particle,trim(particle_filename))

  !> load jorek fields
  if(use_fields_hb) then
    allocate(jorek_fields_interp_hermite_birkhoff::particle_sim%fields)
    read_jorek_interp_hb = read_jorek_fields_interp_hermite_birkhoff(&
    basename,jorek_id,rst_format)
    call read_jorek_interp_hb%do(sim_particle)
  else
    allocate(jorek_fields_interp_linear::particle_sim%fields)
    read_jorek_interp_lin = read_jorek_fields_interp_linear(&
    basename,jorek_id,rst_format)
    call read_jorek_interp_lin%do(sim_particle)
  endif
  
end subroutine initialise_jorek_sim

!> Tools ------------------------------------------------
!>-------------------------------------------------------
end module mod_light_vertices
