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
type(sobseq_rng)             :: sob_rng
type(event)                  :: field_reader
integer                      :: ii,n_particles,nR,nZ,nphi,np,npitch,nchi
integer                      :: n_int_pdf_param,n_real_pdf_param
integer,dimension(:),allocatable :: int_pdf_param
real*8                       :: start_time,mass,charge,error,error_norm
real*8                       :: error_avg_norm,pdf_upper_bound
real*8,dimension(2)          :: Rbox,Zbox,Rbound,Zbound,Phibound
real*8,dimension(2)          :: Ekinbound,Pbound,Pitchbound,Chibound,Chargebound
real*8,dimension(6)          :: var_min,var_max   
real*8,dimension(:),allocatable           :: Rmesh,Zmesh,phimesh,pmesh,pitchmesh,chimesh
real*8,dimension(:),allocatable           :: real_pdf_param
real*8,dimension(:,:,:,:,:,:),allocatable :: expected_pdf,pdf_at_midpoints
character(len=:),allocatable              :: jorek_filename

!> MPI and groups initialisation ------------------------------------------------------------
call sim%initialize(num_groups=1)

!>-------------------------------------------------------------------------------------------
!> Define inputs ----------------------------------------------------------------------------
n_particles = 10
nR          = 10
nZ          = 10
nphi        = 10
np          = 10
npitch      = 10
nchi        = 10
start_time  = 0.d0
mass        = 5.48579909065d-4 !< electron mass in AMU
Rbound      = [0.d0,9.99d2]
Zbound      = [-9.99d2,9.99d2]
Phibound    = 5.d-1*[PI,3.d0*PI]
Ekinbound   = [2d7-1d6,2d7+1d6]
Pitchbound  = [PI-2.95d-1,PI]
Chibound    = [0.d0,TWOPI]
Chargebound = -1.d0
charge      = -1.d0
allocate(character(len=25)::jorek_filename)
jorek_filename   = 'jorek_equilibrium' 
n_int_pdf_param  = 0
n_real_pdf_param = 0

!> Initialisation ---------------------------------------------------------------------------
!> allocate arrays and initialise them to 0
write(*,*) "Test: initialise_particle_in_phase_space: started."
write(*,*) "... setting-up test features"
allocate(Rmesh(nR)); allocate(Zmesh(nZ)); allocate(phimesh(nphi)); allocate(pmesh(np));
allocate(pitchmesh(npitch)); allocate(chimesh(nchi));
allocate(expected_pdf(nR,nZ,nphi,np,npitch,nchi));
allocate(pdf_at_midpoints(nR,nZ,nphi,np,npitch,nchi));

Rmesh = 0.d0; Zmesh = 0.d0; phimesh = 0.d0; pmesh = 0.d0; pitchmesh = 0.d0; chimesh = 0.d0;
expected_pdf = 0.d0; pdf_at_midpoints = 0.d0;
!> read jorek field
field_reader = event(read_jorek_fields_interp_linear(basename=trim(jorek_filename),i=-1))
call with(sim,field_reader)
!> initailise simulation parameters
sim%time = start_time
sim%groups(1)%mass = mass
allocate(particle_kinetic_relativistic::sim%groups(1)%particles(n_particles))
call domain_bounding_box(sim%fields%node_list,sim%fields%element_list,Rbox(1),Rbox(2),Zbox(1),Zbox(2))
if(Rbox(1).ge.Rbound(1)) Rbound(1) = Rbox(1)
if((Rbox(2).lt.Rbound(2)).and.((Rbox(2)-Rbound(1)).gt.0.d0)) Rbound(2) = Rbox(2)
if((Zbox(1).gt.0.d0).and.(Zbound(1).ge.Zbox(1))) Zbound(1) = Zbox(1)
if((Zbox(1).lt.0.d0).and.(Zbound(1).lt.Zbox(1))) Zbound(1) = Zbox(1)
if((Zbox(2).gt.0.d0).and.(Zbound(2).ge.Zbox(2)).and.((Zbox(2)-Zbound(1)).gt.0.d0)) Zbound(2) = Zbox(2)
if((Zbox(2).lt.0.d0).and.(Zbound(2).lt.Zbox(2)).and.((Zbox(2)-Zbound(1)).gt.0.d0)) Zbound(2) = Zbox(2)
Pbound = mass*SPEED_OF_LIGHT*sqrt(((EL_CHG*Ekinbound/(ATOMIC_MASS_UNIT*mass*SPEED_OF_LIGHT**2))+1.d0)**2-1.d0)
var_min = [Rbound(1),Zbound(1),Phibound(1),Pbound(1),Pitchbound(1),Chibound(1)]
var_max = [Rbound(2),Zbound(2),Phibound(2),Pbound(2),Pitchbound(2),Chibound(2)]
!> compute the pdf upper bound
pdf_upper_bound = sup_pdf_uniform(6,var_min,var_max,n_real_pdf_param,&
real_pdf_param,n_int_pdf_param,int_pdf_param);

!> Test particle initialisation -------------------------------------------------------------
write(*,*) "... initialising particles in phase space"
call initialise_particles_in_phase_space(sim%groups(1)%particles,sim%fields,sob_rng,&
pdf_uniform,pdf_upper_bound,sim%groups(1)%mass,start_time,Ekinbound,&
Pitchbound,Chibound,Rbound,Zbound,Phibound,chargebound,n_real_pdf_param,&
real_pdf_param,n_int_pdf_param,int_pdf_param)

!> Produce the expected pdf fromt the particle histogram --------------------------------------
write(*,*) "... building particle histogram and computing the expected pdf"
call compute_equidistant_mesh(Rmesh,nR,Rbound)
call compute_equidistant_mesh(Zmesh,nZ,Zbound)
call compute_equidistant_mesh(phimesh,nphi,Phibound) 
call compute_equidistant_mesh(pmesh,np,Pbound) 
call compute_equidistant_mesh(pitchmesh,npitch,Pitchbound)
call compute_equidistant_mesh(chimesh,nchi,Chibound)
select type(plist=>sim%groups(1)%particles)
  type is (particle_kinetic_relativistic)
  call estimate_pdf_from_histogram(expected_pdf,n_particles,plist,start_time,nR,nZ,nphi,np,&
  npitch,nchi,Rmesh,Zmesh,phimesh,pmesh,pitchmesh,chimesh,sim%fields)
end select

!> Compute the input pdf at the mesh element midpoints and the L2 error w.r.t. the 
!> the expected pdf from the  particle histogram --------------------------------------------
write(*,*) "... computing the input pdf at the midpoints of the mesh elements"
call evaluate_pdf_at_midpoints(pdf_at_midpoints,nR,nZ,nphi,np,npitch,nchi,Rmesh,Zmesh,phimesh,&
     pmesh,pitchmesh,chimesh,charge,start_time,pdf_uniform,sim%fields,n_real_pdf_param,&
     real_pdf_param,n_int_pdf_param,int_pdf_param)
write(*,*) "... computing L2 error"
call compute_error_norm2_ndim6(error,error_norm,error_avg_norm,nR-1,nZ-1,nphi-1,np-1,&
npitch-1,nchi-1,expected_pdf,pdf_at_midpoints,pdf_upper_bound) 

!> Log test results -------------------------------------------------------------------------
write(*,*) "... logging test results"
write(*,*) "L2 error between the expected pdf from particle histogram and the input pdf at mid points: ",error
write(*,*) "L2 error normalized to the maximum of the input pdf: ",error_norm
write(*,*) "L2 error averaged and normalised to the maximum of the input pdf: ",error_avg_norm
write(*,*) " "

!> Clean-up ---------------------------------------------------------------------------------
deallocate(Rmesh); deallocate(Zmesh); deallocate(phimesh); deallocate(pmesh); 
deallocate(pitchmesh); deallocate(chimesh); deallocate(expected_pdf); 
deallocate(pdf_at_midpoints);
call sim%finalize()
write(*,*) "Test: initialise_particle_in_phase_space: completed."

contains

!> Compute the L2 error of 6D-arrays
subroutine compute_error_norm2_ndim6(error_L2,error_L2_norm,error_L2_avg_norm,&
n1,n2,n3,n4,n5,n6,array1,array2,sup_array2)
  implicit none
  !> inputs
  integer,intent(in)                             :: n1,n2,n3,n4,n5,n6
  real*8,intent(in)                              :: sup_array2
  real*8,dimension(n1,n2,n3,n4,n5,n6),intent(in) :: array1,array2
  !> outputs
  real*8,intent(out) :: error_L2,error_L2_norm,error_L2_avg_norm
  !> variables
  integer :: ii,jj,kk,pp,qq
  real*8  :: error

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
  error_L2 = sqrt(error); error_L2_norm = error_L2/abs(sup_array2);
  error_L2_norm = sqrt(error/(n1*n2*n3*n4*n5*n6))/abs(sup_array2)
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

!> Estimate the pdf from the histogram
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
!>   estimated_pdf:   (real8)(nR-1,nZ-1,nphi-1,np-1,npitch-1,ngyro-1) 6D pdf estimation
subroutine estimate_pdf_from_histogram(estimated_pdf,n_particles,particles,time,nR,&
  nZ,nphi,np,npitch,ngyro,Rmesh,Zmesh,phimesh,pmesh,pitchmesh,gyromesh,fields)
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
  real*8,dimension(nR-1,nZ-1,nphi-1,np-1,npitch-1,ngyro-1),intent(out) :: estimated_pdf
  !> internal variables
  integer              :: ii,jj,kk,pp,qq,ss
  real*8               :: dR,dZ,dphi,dp,dpitch,dgyro
  real*8               :: psi,U,ppar,one_third
  integer,dimension(6) :: ids
  real*8,dimension(3)  :: B,E,e1,e2,pcyl
  real*8,dimension(6)  :: dist,coord,mesh_init

  !> Initialisations-----------------------------------------------------------------------
  estimated_pdf = 0.d0; one_third = 1.d0/3.d0;
  !> compute interval distances
  dist = [Rmesh(2)-Rmesh(1),Zmesh(2)-Zmesh(1),phimesh(2)-phimesh(1),&
  pmesh(2)-pmesh(1),pitchmesh(2)-pitchmesh(1),gyromesh(2)-gyromesh(1)]
  mesh_init = [Rmesh(1),Zmesh(1),phimesh(1),pmesh(1),pitchmesh(1),gyromesh(1)]
  !> Compute the histogram
  !$omp parallel do default(shared) firstprivate(n_particles,time) &
  !$omp private(ii,coord,E,B,psi,U,pcyl,e1,e2,ppar,ids) &
  !$omp reduction(+:estimated_pdf)
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
    estimated_pdf(ids(1),ids(2),ids(3),ids(4),ids(5),ids(6)) = &
    estimated_pdf(ids(1),ids(2),ids(3),ids(4),ids(5),ids(6)) + particles(ii)%weight
  enddo
  !$omp end parallel do
  !> normalize for the the total weight mass
  estimated_pdf = estimated_pdf/sum(particles(:)%weight)
  !> divide for the volume of the mesh element
  !$omp parallel do default(none) &
  !$omp shared(estimated_pdf,Rmesh,Zmesh,phimesh,pmesh,pitchmesh,chimesh) &
  !$omp firstprivate(ngyro,npitch,np,nphi,nZ,nR,one_third) &
  !$omp private(ii,jj,kk,pp,qq,ss,dist) &
  !$omp collapse(6)
  do ii=1,ngyro-1
    do jj=1,npitch-1
      do kk=1,np-1
        do pp=1,nphi-1
          do qq=1,nZ-1
            do ss=1,nR-1
              dist = [5d-1*(Rmesh(ss+1)**2-Rmesh(ss)**2),Zmesh(qq+1)-Zmesh(qq),phimesh(pp+1)-phimesh(pp),&
                     one_third*(pmesh(kk+1)**3-pmesh(kk)**3),cos(pitchmesh(jj))-cos(pitchmesh(jj+1)),&
                     chimesh(ii+1)-chimesh(ii)]
              estimated_pdf(ss,qq,pp,kk,jj,ii) = estimated_pdf(ss,qq,pp,kk,jj,ii)/&
              product(dist,mask=(abs(dist).gt.0.d0))       
            enddo          
          enddo
        enddo
      enddo
    enddo
  enddo  
  !$omp end parallel do
end subroutine estimate_pdf_from_histogram

!> evaluate the pdf at the mesh element midpoints
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
!>   pdf_at_midpoints: (real8) value of the pdf at the midpoints
subroutine evaluate_pdf_at_midpoints(&
pdf_midpoints,nR,nZ,nphi,np,npitch,ngyro,Rmesh,Zmesh,phimesh,pmesh,& 
pitchmesh,gyromesh,charge,time,pdf,fields,n_real_pdf_param_in,&
real_pdf_param_in,n_int_pdf_param_in,int_pdf_param_in)
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
  integer,intent(in),optional         :: n_real_pdf_param_in,n_int_pdf_param_in
  integer,dimension(:),allocatable,intent(in),optional :: int_pdf_param_in
  real*8,dimension(:),allocatable,intent(in),optional  :: real_pdf_param_in
  procedure(pdf_f)                    :: pdf
  !> Outputs:
  real*8,dimension(nR-1,nZ-1,nphi-1,np-1,npitch-1,ngyro-1),intent(out) :: pdf_midpoints
  !> Variables:
  integer                  :: ii,jj,kk,pp,qq,rr
  integer                  :: i_elm,ifail,n_real_pdf_param,n_int_pdf_param
  integer,dimension(:),allocatable :: int_pdf_param
  real*8                   :: one_over_six,dummy_double_1,dummy_double_2
  real*8,dimension(2)      :: st
  real*8,dimension(nx)     :: x_midpoints,x_min,x_max
  real*8,dimension(nR)     :: R2mesh
  real*8,dimension(np)     :: p3mesh
  real*8,dimension(npitch) :: cospitchmesh
  real*8,dimension(:),allocatable :: real_pdf_param
  !> initialisation
  one_over_six = 1.d0/6.d0; pdf_midpoints = 0.d0; 
  x_midpoints = [0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,charge];
  x_min = [Rmesh(1),Zmesh(1),phimesh(1),pmesh(1),pitchmesh(1),gyromesh(1),charge]
  x_max = [Rmesh(nR),Zmesh(nZ),phimesh(nphi),pmesh(np),pitchmesh(npitch),gyromesh(ngyro),charge]
  n_real_pdf_param = 0; if(present(n_real_pdf_param_in)) n_real_pdf_param = n_real_pdf_param_in;
  if(present(real_pdf_param_in).and.(n_real_pdf_param_in.gt.0)) then
    allocate(real_pdf_param(n_real_pdf_param)); real_pdf_param = real_pdf_param_in;
  endif
  n_int_pdf_param = 0; if(present(n_int_pdf_param_in)) n_int_pdf_param = n_int_pdf_param_in;
  if(present(int_pdf_param_in).and.(n_int_pdf_param_in.gt.0)) then
    allocate(int_pdf_param(n_int_pdf_param)); int_pdf_param = int_pdf_param_in;
  endif
  !> loop over all midpoints, try with openmp collapse clause first. If slow,
  !> manually collapse all loops in one
  !$omp parallel do default(shared) &
  !$omp firstprivate(ngyro,npitch,np,nphi,nZ,nR,one_over_six,time,x_midpoints,x_min,x_max,&
  !$omp n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param) &
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
              pdf_midpoints(rr,qq,pp,kk,jj,ii) = pdf(nx,x_midpoints,st,time,i_elm,fields,&
              x_min,x_max,n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !$omp end parallel do
  !> Clean-up
  if(allocated(real_pdf_param)) deallocate(real_pdf_param);
  if(allocated(int_pdf_param))  deallocate(int_pdf_param);
end subroutine evaluate_pdf_at_midpoints

!> Phase space distribution for testing
!> inputs:
!>   nx:           (integer) number of variables
!>   x:            (real8)(nx) random state to accept
!>   i_elm:        (integer) jorek mesh element number
!>   st:           (real8)(2) local mesh coordinates
!>   fields:       (fields_base) jorek MHD fields
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!> outputs:
!>   pdf: (real8) value of the probability density 
function pdf_uniform(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,i_elm,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: pdf_uniform
  !> Evalutate pdf
  pdf_uniform = 6.d0/((x_max(1)**2-x_min(1)**2)*(x_max(2)-x_min(2))*&
  (x_max(3)-x_min(3))*(x_max(4)**3-x_min(4)**3)*&
  (cos(x_min(5))-cos(x_max(5)))*(x_max(6)-x_min(6)))
end function pdf_uniform

!> Upper bound of the phase space distribution for testing
!> inputs:
!>   nx:           (integer) number of variables
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!> outputs:
!>   sup_pdf:      (real8) value of the probability density upper bound
function sup_pdf_uniform(nx,x_min,x_max,n_real_param,real_param,&
n_int_param,int_param) result(sup_pdf)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,dimension(nx),intent(in)             :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: sup_pdf
  !> Evalutate pdf
  sup_pdf = 6.d0/((x_max(1)**2-x_min(1)**2)*(x_max(2)-x_min(2))*&
  (x_max(3)-x_min(3))*(x_max(4)**3-x_min(4)**3)*&
  (cos(x_min(5))-cos(x_max(5)))*(x_max(6)-x_min(6)))
end function sup_pdf_uniform

end program test_initialisation_phase_space
