program test_initialisation_phase_space
!> program used for testing the initialisation method: 
!>   initialise_particles_in_phase_space 
!> contained in mod_initialise_particles
!> the idea is to provide a reject_uniform function
!> which reject a uniformly distributed random sample
!> between [0 1] as a function of a user defined 
!> 6D distribution. The distribution function as to be
!> in R,Z,phi coordinates for the physical space, 
!> momentum, pitch and gyro angles coordinates for the
!> momentum space and the distribution of electric charges
!> be aware that the electric charge should be a double 
!use mod_particle_sim
!use mod_particle_types,       only: particle_base,particle_kinetic_relativistic 
!use mod_fields_linear,        only: jorek_field_interp_linear,read_jorek_fields_interp_linear
!use mod_initialise_particles, only: initialise_particles_in_phase_space
use constants, only: TWOPI,PI
use mod_random_seed
use particle_tracer
implicit none

!> Variable declarations --------------------------------------------------------------------
integer,parameter            :: n_groups=1
type(sobseq_rng)             :: sob_rng
type(event)                  :: field_reader
integer                      :: ii,n_particles
real*8                       :: start_time,mass
real*8,dimension(2)          :: Rbox,Zbox
real*8,dimension(2)          :: Rbound,Zbound,Phibound,Ekinbound,Pitchbound,Chibound,Chargebound
character(len=:),allocatable :: jorek_filename

!> MPI and groups initialisation ------------------------------------------------------------
call sim%initialize(num_groups=n_groups)

!>-------------------------------------------------------------------------------------------
!> Define inputs ----------------------------------------------------------------------------
n_particles = 1000
start_time  = 0.d0
mass        = 5.48579909065d-4 !< electron mass in AMU
Rbound      = [0.d0,9.99d2]
Zbound      = [-9.99d2,9.99d2]
Phibound    = 5.d-1*[PI,3.d0*PI]
Ekinbound   = 2d7
Pitchbound  = [PI-2.95d-1,PI]
Chibound    = [0.d0,TWOPI]
Chargebound = -1.d0
allocate(character(len=25)::jorek_filename)
jorek_filename = 'jorek_equilibrium' 

!> Initialisation ---------------------------------------------------------------------------
!> read jorek field
field_reader = event(read_jorek_fields_interp_linear(basename=trim(jorek_filename),i=-1))
call with(sim,field_reader)
!> initailise simulation parameters
sim%time = start_time
do ii=1,n_groups
  sim%groups(ii)%mass = mass
  allocate(particle_kinetic_relativistic::sim%groups(ii)%particles(n_particles))
enddo
call domain_bounding_box(sim%fields%node_list,sim%fields%element_list,Rbox(1),Rbox(2),Zbox(1),Zbox(2))

!> Test particle initialisation -------------------------------------------------------------
do ii=1,n_groups
  call initialise_particles_in_phase_space(sim%groups(ii)%particles,sim%fields,sob_rng,&
  reject_uniform,sim%groups(ii)%mass,start_time,Ekinbound,Pitchbound,Chibound,&
  Rbound,Zbound,Phibound,chargebound)
enddo

!> Produce histogram from particles ---------------------------------------------------------
if(Rbox(1).ge.Rbound(1)) Rbound(1) = Rbox(1)
if((Rbox(2).lt.Rbound(2)).and.((Rbox(2)-Rbound(1)).gt.0.d0)) Rbound(2) = Rbox(2)
if((Zbox(1).gt.0.d0).and.(Zbound(1).ge.Zbox(1))) Zbound(1) = Zbox(1)
if((Zbox(1).lt.0.d0).and.(Zbound(1).lt.Zbox(1))) Zbound(1) = Zbox(1)
if((Zbox(2).gt.0.d0).and.(Zbound(2).ge.Zbox(2)).and.((Zbox(2)-Zbound(1)).gt.0.d0)) Zbound(2) = Zbox(2)
if((Zbox(2).lt.0.d0).and.(Zbound(2).lt.Zbox(2)).and.((Zbox(2)-Zbound(1)).gt.0.d0)) Zbound(2) = Zbox(2)


!> Clean-up ---------------------------------------------------------------------------------
call sim%finalize()

contains

!> Generate an equidistant mesh
!> inputs:
!>   n_points: (integer) number of nodes
!>   bounds:   (real8)(2) upper and lower mesh extrema
!> outputs:
!>   mesh:     (real8)(n_points) equidistant mesh nodes
subroutine compute_equidistant_mesh(mesh,n_points,bounds)
  implicit none
  !> input variables
  integer,intent(in)                     :: n_points
  real*8,dimension(2),intent(in)         :: bounds
  !> output variables
  real*8,dimension(n_points),intent(out) :: mesh
  !> variables
  real*8                                 :: delta
  
  delta = (bounds(2)-bounds(1))/(n_points-1)
  do ii=1,n_points
    mesh(ii) = bounds(1)+real(delta*(ii-1),kind=8)
  enddo
end subroutine compute_equidistant_mesh

!> Prepare a histogram given equidistant meshes
!> inputs:
!>   n_particles: (integer) number of particles
!>   particles:   (particle_kinetic_relativistic)(n_particles) particle list
!>   time:        (real8) time of the required MHD fields
!>   nR:          (integer) size of the mesh along the major radius
!>   nZ:          (integer) size of the mesh along the vertical coordinate  
!>   nphi:        (integer) size of the mesh along the toroidal angle
!>   np:          (integer) size of the mesh along the total momentum
!>   npitch:      (integer) size of the mesh along the pitch angle
!>   ngyro:       (integer) size of the mesh along the gyro angle
!>   Rmesh:       (real8)(nR) major radius equidistant mesh
!>   Zmesh:       (real8)(nZ) vertical coordinate equidistant mesh
!>   phimesh:     (real8)(nphi) toroidal angle equidistant mesh
!>   pmesh:       (real8)(np) total momentum equidistant mesh
!>   pitchmesh:   (real8)(npitch) pitch angle equidistant mesh
!>   gyromesh:    (real8)(ngyro) gyro angle equidistant mesh
!> outputs:
!>   histogram:   (real8)(nR,nZ,nphi,np,npitch,ngyro) 6D histogram
subroutine equidistant_histogram(histogram,n_particles,particles,time,nR,nZ,nphi,&
  np,npitch,ngyro,Rmesh,Zmesh,phimesh,pmesh,pitchmesh,gyromesh,fields)
  use constants,                 only: TWOPI
  use mod_coordinate_transforms, only: vector_cartesian_to_cylindrical
  use mod_fields,                only: fields_base
  use mod_particle_types,        only: particle_kinetic_relativistic
  use mod_pusher_tools,          only: get_orthonormals
  implicit none
  !> inputs:
  integer,intent(in)                                       :: n_particles
  type(particle_kinetic_relativistic),dimension(n_particles),intent(in) :: particles
  class(fields_base),intent(in)                            :: fields
  real*8,intent(in)                                        :: time
  integer,intent(in)                                       :: nR,nZ,nphi,np,npitch,ngyro
  real*8,dimension(nR),intent(in)                          :: Rmesh
  real*8,dimension(nZ),intent(in)                          :: Zmesh
  real*8,dimension(nphi),intent(in)                        :: phimesh
  real*8,dimension(np),intent(in)                          :: pmesh
  real*8,dimension(npitch),intent(in)                      :: pitchmesh
  real*8,dimension(ngyro),intent(in)                       :: gyromesh
  !> outputs:
  real*8,dimension(nR,nZ,nphi,np,npitch,ngyro),intent(out) :: histogram
  !> internal variables
  integer              :: ii
  real*8               :: dR,dZ,dphi,dp,dpitch,dgyro
  real*8               :: psi,U,ppar
  integer,dimension(6) :: ids
  real*8,dimension(3)  :: B,E,e1,e2,pcyl
  real*8,dimension(6)  :: dist,coord,mesh_init

  !> Initialisations-----------------------------------------------------------------------
  histogram = 0.d0
  !> compute interval distances
  dist = [Rmesh(2)-Rmesh(1),Zmesh(2)-Zmesh(1),phimesh(2)-phimesh(1),&
  pmesh(2)-pmesh(1),pitchmesh(2)-pitchmesh(1),gyromesh(2)-gyromesh(1)]
  mesh_init = [Rmesh(1),Zmesh(1),phimesh(1),pmesh(1),pitchmesh(1),gyromesh(1)]
  !> Compute the histogram
  !$omp parallel do default(shared) firstprivate(n_particles,time) &
  !$omp private(ii,coord,E,B,psi,U,pcyl,e1,e2,ppar,ids) &
  !$omp reduction(+:histogram)
  do ii=1,n_particles
    coord(1:3) = particles(ii)%x
    call fields%calc_EBpsiU(time,particles(ii)%i_elm,particles(ii)%st,&
    coord(3),E,B,psi,U)
    pcyl = vector_cartesian_to_cylindrical(coord(3),particles(ii)%p)
    B = B/norm2(B)
    call get_orthonormals(B,e1,e2)
    coord(4)  = norm2(pcyl)
    ppar  = dot_product(pcyl,B)
    coord(5) = acos(ppar/coord(4)) 
    coord(6)  = atan2(dot_product(pcyl-ppar*B,e2),dot_product(pcyl-ppar*B,e1))
    if(coord(6).lt.0.d0) coord(6) = TWOPI+coord(6)
    coord(6) = mod(coord(6),TWOPI)
    ids = floor((coord-mesh_init)/dist)+1
    histogram(ids(1),ids(2),ids(3),ids(4),ids(5),ids(6)) = &
    histogram(ids(1),ids(2),ids(3),ids(4),ids(5),ids(6)) + particles(ii)%weight
  enddo
  !$omp end parallel do
end subroutine equidistant_histogram

!> Phase space distributions for initialisation
!> inputs:
!>   n_x:    (integer)    number of variables
!>   x:      (real8)(n_x) variables (R,Z,phi,energy,pitch,gyro,charge)
!>   st:     (real8)(2)   particle position in local mesh coord.
!>   i_elm:  (integer)    element index containing the particle
!>   rand:   (real8)      random number
!>   fields: (field_base) MHD fields class
!> outputs:
!>   reject_uniform: (logical) if true the sample is rejected
function reject_uniform(n_x,x,st,i_elm,rand,fields)
  use mod_fields, only: fields_base
  implicit none
  !> constants
  real*8,parameter                 :: threshold=5.d-1
  !> inputs:
  integer,intent(in)               :: n_x,i_elm
  real*8,intent(in)                :: rand
  real*8,dimension(n_x),intent(in) :: x
  real*8,dimension(2),intent(in)   :: st
  class(fields_base),intent(in)    :: fields
  !> outputs:
  logical                          :: reject_uniform

  !> check if the particle is valid
  reject_uniform = .true.
  if(i_elm.le.0) return
  if((st(1).lt.0.d0).or.(st(1).gt.1.d0)) return
  if((st(2).lt.0.d0).or.(st(2).gt.1.d0)) return
  !> rject of accept solution
  if(rand.le.threshold) reject_uniform = .false.

end function reject_uniform

end program test_initialisation_phase_space
