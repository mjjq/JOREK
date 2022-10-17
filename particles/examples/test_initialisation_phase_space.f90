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

!> Interfaces -------------------------------------------------------------------------------
interface
  function pdf_f(nx,x,i_elm,time,st,fields) result(pdf)
    use mod_fields, only: fields_base
    implicit none
    !> inputs:
    integer,intent(in)              :: nx
    real*8,intent(in)               :: time
    real*8,dimension(nx),intent(in) :: x
    integer,intent(in)              :: i_elm
    real*8,dimension(2),intent(in)  :: st
    class(fields_base),intent(in)   :: fields
    !> outputs:
    real*8 :: pdf
  end function pdf_f
end interface

!> Variable declarations --------------------------------------------------------------------
type(sobseq_rng)             :: sob_rng
type(event)                  :: field_reader
integer                      :: ii,n_particles,nR,nZ,nphi,np,npitch,nchi
real*8                       :: start_time,mass,charge,error
real*8,dimension(2)          :: Rbox,Zbox,Rbound,Zbound,Phibound
real*8,dimension(2)          :: Ekinbound,Pbound,Pitchbound,Chibound,Chargebound
real*8,dimension(:),allocatable           :: Rmesh,Zmesh,phimesh,pmesh,pitchmesh,chimesh
real*8,dimension(:,:,:,:,:,:),allocatable :: histogram,number_of_particles
character(len=:),allocatable              :: jorek_filename

!> MPI and groups initialisation ------------------------------------------------------------
call sim%initialize(num_groups=1)

!>-------------------------------------------------------------------------------------------
!> Define inputs ----------------------------------------------------------------------------
n_particles = 1000
nR          = 5
nZ          = 5
nphi        = 5
np          = 5 
npitch      = 5
nchi        = 5
start_time  = 0.d0
mass        = 5.48579909065d-4 !< electron mass in AMU
Rbound      = [0.d0,9.99d2]
Zbound      = [-9.99d2,9.99d2]
Phibound    = 5.d-1*[PI,3.d0*PI]
Ekinbound   = [2d7-1,2d7+1]
Pitchbound  = [PI-2.95d-1,PI]
Chibound    = [0.d0,TWOPI]
Chargebound = -1.d0
charge      = -1.d0 
allocate(character(len=25)::jorek_filename)
jorek_filename = 'jorek_equilibrium' 

!> Initialisation ---------------------------------------------------------------------------
!> allocate arrays and initialise them to 0
write(*,*) "Test: initialise_particle_in_phase_space: started."
write(*,*) "... setting-up test features"
allocate(Rmesh(nR)); allocate(Zmesh(nZ)); allocate(phimesh(nphi)); allocate(pmesh(np));
allocate(pitchmesh(npitch)); allocate(chimesh(nchi));
allocate(histogram(nR,nZ,nphi,np,npitch,nchi));
allocate(number_of_particles(nR,nZ,nphi,np,npitch,nchi));
Rmesh = 0.d0; Zmesh = 0.d0; phimesh = 0.d0; pmesh = 0.d0; pitchmesh = 0.d0; chimesh = 0.d0;
histogram = 0.d0; number_of_particles = 0.d0;
!> read jorek field
field_reader = event(read_jorek_fields_interp_linear(basename=trim(jorek_filename),i=-1))
call with(sim,field_reader)
!> initailise simulation parameters
sim%time = start_time
sim%groups(1)%mass = mass
allocate(particle_kinetic_relativistic::sim%groups(1)%particles(n_particles))
call domain_bounding_box(sim%fields%node_list,sim%fields%element_list,Rbox(1),Rbox(2),Zbox(1),Zbox(2))

!> Test particle initialisation -------------------------------------------------------------
write(*,*) "... initialising particles in phase space"
call initialise_particles_in_phase_space(sim%groups(1)%particles,sim%fields,sob_rng,&
reject_uniform,sim%groups(1)%mass,start_time,Ekinbound,Pitchbound,Chibound,&
Rbound,Zbound,Phibound,chargebound)

!> Produce histogram from particles ---------------------------------------------------------
write(*,*) "... building particle histogram"
if(Rbox(1).ge.Rbound(1)) Rbound(1) = Rbox(1)
if((Rbox(2).lt.Rbound(2)).and.((Rbox(2)-Rbound(1)).gt.0.d0)) Rbound(2) = Rbox(2)
if((Zbox(1).gt.0.d0).and.(Zbound(1).ge.Zbox(1))) Zbound(1) = Zbox(1)
if((Zbox(1).lt.0.d0).and.(Zbound(1).lt.Zbox(1))) Zbound(1) = Zbox(1)
if((Zbox(2).gt.0.d0).and.(Zbound(2).ge.Zbox(2)).and.((Zbox(2)-Zbound(1)).gt.0.d0)) Zbound(2) = Zbox(2)
if((Zbox(2).lt.0.d0).and.(Zbound(2).lt.Zbox(2)).and.((Zbox(2)-Zbound(1)).gt.0.d0)) Zbound(2) = Zbox(2)
Pbound = mass*SPEED_OF_LIGHT*sqrt(((EL_CHG*Ekinbound/(mass*SPEED_OF_LIGHT**2))+1.d0)**2-1.d0)
call compute_equidistant_mesh(Rmesh,nR,Rbound)
call compute_equidistant_mesh(Zmesh,nZ,Zbound)
call compute_equidistant_mesh(phimesh,nphi,Phibound) 
call compute_equidistant_mesh(pmesh,np,Pbound) 
call compute_equidistant_mesh(pitchmesh,npitch,Pitchbound)
call compute_equidistant_mesh(chimesh,nchi,Chibound)
select type(plist=>sim%groups(1)%particles)
  type is (particle_kinetic_relativistic)
  call equidistant_histogram(histogram,n_particles,plist,start_time,nR,nZ,nphi,np,&
  npitch,nchi,Rmesh,Zmesh,phimesh,pmesh,pitchmesh,chimesh,sim%fields)
end select

!> Compute the expected number of particles from the PDF for each kinetic mesh element and
!> the L2 error w.r.t. the particle histogram --------------------------------------------
write(*,*) "... computing expected number of particles from the pdf (rectangle integration)"
call compute_nparticles_from_f_equidistant_RZphipthetachi_mesh(&
     number_of_particles,nR,nZ,nphi,np,npitch,nchi,Rmesh,Zmesh,phimesh,&
     pmesh,pitchmesh,chimesh,charge,start_time,pdf_uniform,sim%fields)
write(*,*) "... computing L2 error"
call compute_error_norm2_ndim6(error,nR-1,nZ-1,nphi-1,np-1,npitch-1,nchi-1,&
histogram,number_of_particles) 

!> Log test results -------------------------------------------------------------------------
write(*,*) "... logging test results"
write(*,*) "L2 error between the particle histogram and the expected number of particles: ",error
write(*,*) " "

!> Clean-up ---------------------------------------------------------------------------------
deallocate(Rmesh); deallocate(Zmesh); deallocate(phimesh); deallocate(pmesh); 
deallocate(pitchmesh); deallocate(chimesh); deallocate(histogram); 
deallocate(number_of_particles);
call sim%finalize()
write(*,*) "Test: initialise_particle_in_phase_space: completed."

contains

!> Compute the L2 error of 6D-arrays
subroutine compute_error_norm2_ndim6(error,n1,n2,n3,n4,n5,n6,array1,array2)
  implicit none
  !> inputs
  integer,intent(in)                             :: n1,n2,n3,n4,n5,n6
  real*8,dimension(n1,n2,n3,n4,n5,n6),intent(in) :: array1,array2
  !> outputs
  real*8,intent(out) :: error
  !> variables
  integer :: ii,jj,kk,pp,qq

  !> initialisation
  error = 0.d0
  !> compute norm2 error, the openmp collapse clause is used at first. Manual collapse
  !> of the indices should be tried in case of reduced performance.
  !$omp parallel do default(none) firstprivate(n2,n3,n4,n5,n6) &
  !$omp shared(array1,array2) private(ii,jj,kk,pp,qq) reduction(+:error) &
  !$omp collapse(5)
  do ii=1,n6
    do jj=1,n5
      do kk=1,n4
        do pp=1,n3
          do qq=1,n2
            error = error + dot_product((array2(:,qq,pp,kk,jj,ii)-array1(:,qq,pp,kk,jj,ii)),&
            (array2(:,qq,pp,kk,jj,ii)-array1(:,qq,pp,kk,jj,ii)))
          enddo
        enddo
      enddo
    enddo
  enddo
  !$omp end parallel do
  error = sqrt(error)
end subroutine compute_error_norm2_ndim6

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
!>   histogram:   (real8)(nR-1,nZ-1,nphi-1,np-1,npitch-1,ngyro-1) 6D histogram
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
  real*8,dimension(nR-1,nZ-1,nphi-1,np-1,npitch-1,ngyro-1),intent(out) :: histogram
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
!  !$omp parallel do default(shared) firstprivate(n_particles,time) &
!  !$omp private(ii,coord,E,B,psi,U,pcyl,e1,e2,ppar,ids) &
!  !$omp reduction(+:histogram)
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
!  !$omp end parallel do
end subroutine equidistant_histogram

!> compute the number of particles in an equidistant cyclindrial in space
!> spherical in momentum space mesh given a particle distributon function
!> using the 6D rectangle integration method
!> inputs:
!>   nR:        (integer) number of the major radius mesh nodes
!>   nZ:        (integer) number of the vertical mesh nodes
!>   nphi:      (integer) number of the toroidal angle mesh nodes
!>   np:        (integer) number of the momentum mesh nodes
!>   npitch:    (integer) number of the pitch-angle mesh nodes
!>   ngyro:     (integer) number of the gyro-angle mesh nodes
!>   Rmesh:     (real8)(nR) major radius mesh
!>   Zmesh:     (real8)(nZ) vertical coordinate mesh
!>   phimesh:   (real8)(nphi) toroidal mesh
!>   pitchmesh: (real8)(npitch) pitch angle mesh
!>   gyromesh:  (real8)(ngyro) gyro angle mesh
!>   pdf:       (procedure) particle distribution function
!>   fields:    (fields_base) jorek MHD fields
!> outputs:
!>   particle_number: (real8)(nR-1,nZ-1,nphi-1,np-1,npitch-1,ngyro-1) 
!>                    number of particles at the mesh element center 
!>                    integrated via rectangle method  
subroutine compute_nparticles_from_f_equidistant_RZphipthetachi_mesh(&
particle_number,nR,nZ,nphi,np,npitch,ngyro,Rmesh,Zmesh,phimesh,pmesh,& 
pitchmesh,gyromesh,charge,time,pdf,fields)
   use mod_fields, only: fields_base
  implicit none
  !> Parameters
  integer,parameter                   :: nx=7
  !> Inputs
  class(fields_base),intent(in)       :: fields
  integer,intent(in)                  :: nR,nZ,nphi,np,npitch,ngyro
  real*8,intent(in)                   :: charge,time
  real*8,dimension(nR),intent(in)     :: Rmesh
  real*8,dimension(nZ),intent(in)     :: Zmesh
  real*8,dimension(nphi),intent(in)   :: phimesh
  real*8,dimension(np),intent(in)     :: pmesh
  real*8,dimension(npitch),intent(in) :: pitchmesh
  real*8,dimension(ngyro),intent(in)  :: gyromesh
  procedure(pdf_f)                    :: pdf
  !> Outputs:
  real*8,dimension(nR-1,nZ-1,nphi-1,np-1,npitch-1,ngyro-1),intent(out) :: particle_number
  !> Variables:
  integer                  :: ii,jj,kk,pp,qq,rr
  integer                  :: i_elm,ifail
  real*8                   :: one_over_six,dummy_double_1,dummy_double_2
  real*8,dimension(2)      :: st
  real*8,dimension(nx)     :: x_midpoints
  real*8,dimension(nR)     :: R2mesh
  real*8,dimension(np)     :: p3mesh
  real*8,dimension(npitch) :: cospitchmesh 
  !> initialisation
  one_over_six = 1.d0/6.d0; cospitchmesh = cos(pitchmesh);
  p3mesh = pmesh**3; R2mesh = Rmesh**2; particle_number = 0.d0; 
  x_midpoints = [0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,charge];
  !> loop over all midpoints, try with openmp collapse clause first. If slow,
  !> manually collapse all loops in one
  !$omp parallel do default(shared) &
  !$omp firstprivate(ngyro,npitch,np,nphi,nZ,nR,one_over_six,time,x_midpoints) &
  !$omp private(ii,jj,kk,pp,qq,rr,i_elm,st,dummy_double_1,dummy_double_2) &
  !$omp collapse(6)
  do ii=1,ngyro-1
    do jj=1,npitch-1
      do kk=1,np-1
        do pp=1,nphi-1
          do qq=1,nZ-1
            do rr=1,nR-1
              !> compute midpoint in the local mesh cell
              x_midpoints(1:6) = [Rmesh(rr),Zmesh(qq),phimesh(pp),pmesh(kk),&
              pitchmesh(jj),gyromesh(ii)] + 5.d-1*([Rmesh(rr+1),Zmesh(qq+1),&
              phimesh(pp+1),pmesh(kk+1),pitchmesh(jj+1),gyromesh(ii+1)] - &
              [Rmesh(rr),Zmesh(qq),phimesh(pp),pmesh(kk),pitchmesh(jj),gyromesh(ii)])
              !> find local JOREK mesh coordinates 
              call find_RZ(fields%node_list,fields%element_list,x_midpoints(1),&
              x_midpoints(2),dummy_double_1,dummy_double_2,i_elm,st(1),st(2),ifail)
              !> estimate the number of particles at the kinetic mesh midpoint
              particle_number(rr,qq,pp,kk,jj,ii) = pdf(nx,x_midpoints,i_elm,time,st,fields)*&
              (p3mesh(kk+1)-p3mesh(kk))*(cospitchmesh(jj+1)-cospitchmesh(jj))*&
              (gyromesh(ii+1)-gyromesh(ii))*(R2mesh(rr+1)-R2mesh(rr))*(Zmesh(qq+1)-Zmesh(qq))*&
              (phimesh(pp+1)-phimesh(pp))*one_over_six
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
 !$omp end parallel do
end subroutine compute_nparticles_from_f_equidistant_RZphipthetachi_mesh

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
function reject_uniform(n_x,x,st,time,i_elm,rand,fields)
  use mod_fields, only: fields_base
  implicit none
  !> constants
  real*8,parameter                 :: threshold=5.d-1
  !> inputs:
  integer,intent(in)               :: n_x,i_elm
  real*8,intent(in)                :: rand,time
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

!> Phase space distribution for testing
function pdf_uniform(nx,x,i_elm,time,st,fields) result(pdf)
  use mod_fields, only: fields_base
  implicit none
  !> constants
  real*8,parameter :: val=5.d-1
  !> Inputs:
  integer,intent(in)              :: nx,i_elm
  real*8,intent(in)               :: time
  real*8,dimension(nx),intent(in) :: x
  real*8,dimension(2),intent(in)  :: st
  class(fields_base),intent(in)   :: fields
  !> Outputs:
  real*8 :: pdf
  !> Evalutate pdf
  pdf = val
end function pdf_uniform

end program test_initialisation_phase_space
