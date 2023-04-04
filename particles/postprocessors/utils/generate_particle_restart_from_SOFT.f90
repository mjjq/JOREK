!> generate_particle_restart_from_SOFT generates a 
!> JOREK particle restart file from a orbits file
!> produced by the SOFT code
program generate_particle_restart_from_SOFT
use constants,       only: PI,SPEED_OF_LIGHT,ATOMIC_MASS_UNIT,EL_CHG
use mod_mpi_tools,   only: init_mpi_threads,finalize_mpi_threads
use mod_particle_io, only: write_simulation_hdf5
use particle_tracer

implicit none

!> Data-types ------------------------------------------------------------------------------
!> derived datatype describing the particle pdf used in SODT
type type_soft_pdf
  real*8,dimension(:),allocatable   :: xi,p !< cos(pitch_angle) and momentum sizes
  real*8,dimension(:,:),allocatable :: pdf  !< particle pdf
end type type_soft_pdf

!> Variables -------------------------------------------------------------------------------
type(pcg32_rng)      :: rng_type
type(event)          :: field_reader
type(type_soft_pdf),dimension(:),allocatable :: soft_pdf_list
logical              :: do_write_particles_in_hdf5
integer              :: my_id,n_cpus,ierr,n_vec,n_groups,n_phi,n_r_pdf_mesh
integer,dimension(2) :: dims
real*8               :: time,mass,charge
real*8,dimension(2)  :: phi_interval
real*8,dimension(:),allocatable   :: soft_orbit_ppar_local,soft_orbit_pperp_local
real*8,dimension(:),allocatable   :: soft_orbit_Jdtdrho_local,soft_orbit_ppar
real*8,dimension(:),allocatable   :: soft_orbit_pperp,soft_orbit_Jdtdrho,soft_pdf_r_mesh
real*8,dimension(:),allocatable   :: soft_R_mesh,soft_Z_mesh
real*8,dimension(:,:),allocatable :: soft_orbit_x,soft_orbit_x_local,soft_r_minor
character(len=17)    :: fields_filename,soft_pdf_filename
character(len=28)    :: soft_magfield_filename
character(len=125)   :: particle_filename
character(len=193)   :: soft_orbit_filename
character(len=25)    :: filename_jorek_hdf5
!> Variable definitions --------------------------------------------------------------------
do_write_particles_in_hdf5 = .false.         !< writeh particle in hdf5 if true
n_groups               = 1                   !< number of jorek particle groups
n_vec                  = 3                   !< component of a vector
n_phi                  = 6                   !< number of toroidal positions to be sampled for each particle
time                   = 0d0                 !< simulation time
mass                   = 5.48579909065d-4    !< electron mass in AMU 
charge                 = -1d0                !< electron charge
phi_interval           = [0d0,PI]            !< toroidal angle interval in which particles are sampled
fields_filename        = 'jorek_equilibrium' !< jorek restart filename
soft_magfield_filename = 'magnetic_field_jorek_to_soft' !< soft magnetic field file
particle_filename      = 'part_restart_soft_orbits.h5'  !< particle restart filename 
soft_pdf_filename      = 'pdf_jorek_to_soft'            !< soft distribution field input from jorek
soft_orbit_filename    = 'orbit_test_jorek_JET_pulse95135_t48dot54_parabolic_qprofile_q95_6dot8_press0_res1r5dot88en1m5_res2r4dot705en1m4_Ip612en1MA_Ekin20MeV_np10_theta1_2dot85_itheta2_pi_nitheta100_norbits100_a96_wall' !< soft orbit filename
filename_jorek_hdf5 = 'jorek_particles_from_soft'
!> Initialisation --------------------------------------------------------------------------
!> initialise the MPI communicator
call init_mpi_threads(my_id,n_cpus,ierr)
!> read mhd data
write(*,*) "Reading MHD data ..."
call sim%initialize(n_groups,.true.,my_id,n_cpus)
field_reader = event(read_jorek_fields_interp_linear(basename=trim(fields_filename),i=-1))
call with(sim,field_reader)
write(*,*) "Reading MHD data: completed!"
!> Read soft input and generate JOREK particle restart -------------------------------------
write(*,*) "Reading SOFT magnetic field, pdf and orbit files ..."
!> compute and broadcast the minor radii array
call read_and_broadcast_soft_minor_radii(my_id,soft_magfield_filename,soft_R_mesh,soft_Z_mesh,soft_r_minor)
if(my_id.eq.0) then
  call read_soft_orbit_file(soft_orbit_filename,n_vec,dims(1),soft_orbit_x,&
  soft_orbit_ppar,soft_orbit_pperp,soft_orbit_Jdtdrho)
  !> scatter the global arrays to each mpi process)
  dims(2) = dims(1)/n_cpus;
  if(dims(2)*n_cpus.lt.dims(1)) dims(2) = dims(2)+1
  !> read the soft_pdf_file
  call read_soft_pdf_file(soft_pdf_filename,n_r_pdf_mesh,soft_pdf_r_mesh,soft_pdf_list)
endif
call MPI_Bcast(dims,size(dims),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
allocate(soft_orbit_x_local(n_vec,dims(2)));  soft_orbit_x_local     = 0d0;
allocate(soft_orbit_ppar_local(dims(2)));     soft_orbit_ppar_local  = 0d0;
allocate(soft_orbit_pperp_local(dims(2)));    soft_orbit_pperp_local = 0d0;
allocate(soft_orbit_Jdtdrho_local(dims(2)));  soft_orbit_Jdtdrho_local = 0d0;
call scatter_2D_arrays(my_id,n_cpus,n_vec,dims(2),dims(1),soft_orbit_x,soft_orbit_x_local)
call MPI_Scatter(soft_orbit_ppar,dims(2),MPI_REAL8,soft_orbit_ppar_local,dims(2),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_Scatter(soft_orbit_pperp,dims(2),MPI_REAL8,soft_orbit_pperp_local,dims(2),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_Scatter(soft_orbit_Jdtdrho,dims(2),MPI_REAL8,soft_orbit_Jdtdrho_local,dims(2),&
MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call broadcast_soft_pdf_list(my_id,n_r_pdf_mesh,soft_pdf_r_mesh,soft_pdf_list)
if(allocated(soft_orbit_x))       deallocate(soft_orbit_x)
if(allocated(soft_orbit_ppar))    deallocate(soft_orbit_ppar)
if(allocated(soft_orbit_pperp))   deallocate(soft_orbit_pperp)
if(allocated(soft_orbit_Jdtdrho)) deallocate(soft_orbit_Jdtdrho)
write(*,*) "Reading SOFT magnetic field, pdf and orbit files: completed!"
!> Generate JOREK relativistic gc from SOFT orbits -----------------------------------------
write(*,*) 'Converting soft orbit to JOREK relativistic gc ...'
call convert_soft_orbits_in_jorek_relativistic_gcs(sim,rng_type,my_id,n_cpus,time,\
mass,charge,n_vec,dims(2),n_phi,phi_interval,soft_orbit_x_local,\
soft_orbit_ppar_local,soft_orbit_pperp_local,soft_orbit_Jdtdrho_local)
write(*,*) 'Converting soft orbit to JOREK relativistic gc: completed!'
write(*,*) 'Writing JOREK relativistic gc in ',trim(particle_filename),' ...'
call write_simulation_hdf5(sim,trim(particle_filename))
write(*,*) 'Writing JOREK relativistic gc in ',trim(particle_filename),' completed!'
!> Write data in HDF5 file -----------------------------------------------------------------
if(do_write_particles_in_hdf5) then
  write(*,*) 'Writing JOREK particles in HDF5 file ...'
  call write_particles_in_hdf5(my_id,filename_jorek_hdf5,n_cpus,n_vec,sim)
  write(*,*) 'Writing JOREK particles in HDF5 file: completed!'
endif
!> Finalisation ----------------------------------------------------------------------------
if(allocated(soft_orbit_x))             deallocate(soft_orbit_x)
if(allocated(soft_orbit_ppar))          deallocate(soft_orbit_ppar)
if(allocated(soft_orbit_pperp))         deallocate(soft_orbit_pperp)
if(allocated(soft_orbit_x_local))       deallocate(soft_orbit_x_local)
if(allocated(soft_orbit_ppar_local))    deallocate(soft_orbit_ppar_local)
if(allocated(soft_orbit_pperp_local))   deallocate(soft_orbit_pperp_local)
if(allocated(soft_orbit_Jdtdrho_local)) deallocate(soft_orbit_Jdtdrho_local)
if(allocated(soft_pdf_r_mesh))          deallocate(soft_pdf_r_mesh)
if(allocated(soft_pdf_list))            deallocate(soft_pdf_list)
if(allocated(soft_R_mesh))              deallocate(soft_R_mesh)
if(allocated(soft_Z_mesh))              deallocate(soft_Z_mesh)
if(allocated(soft_r_minor))             deallocate(soft_r_minor)
call finalize_mpi_threads(ierr)
write(*,*) 'Program terminated: good bye!'

contains

!> generate jorek particle relativistic gc from soft orbits 
!> inputs:
!>   sim:                (particle_sim) particle simulation type to be initialised
!>   rng_base:           (type_rng) type of the random number generator
!>   my_id:              (integer) MPI task rank
!>   n_cpus:             (integer) number of MPI tasks
!>   time:               (real8) simulation time
!>   mass:               (real8) particle mass in AMU
!>   charge:             (real8) particle charge (RE: -1)
!>   n_vec:              (integer) size of the position vector
!>   n_points:           (integer) number of soft valid points
!>   n_phi:              (integer) number of toroidal samples per particles
!>   phi_interval:       (real8)(2) toroidal angle interval for sampling
!>   soft_orbit_x:       (real8)(n_vec,n_points) soft positions in xyz
!>   soft_orbit_ppar:    (real8)(n_points) soft parallel momentum
!>   soft_orbit_pperp:   (real8)(n_points) soft perpendicular momentum
!>   soft_orbit_Jdtdrho: (real8)(n_points) soft jacobian*dpoloidal*dminor_radius
!> outpus:
!>   sim: (particle_sim) initialised particle simulation
subroutine convert_soft_orbits_in_jorek_relativistic_gcs(sim,rng_base,my_id,&
n_cpus,time,mass,charge,n_vec,n_points,n_phi,phi_interval,&
soft_orbit_x,soft_orbit_ppar,soft_orbit_pperp,soft_orbit_Jdtdrho)
  use mpi
  use constants,                 only: SPEED_OF_LIGHT
  use mod_coordinate_transforms, only: cartesian_to_cylindrical
  use mod_random_seed,           only: random_seed
  use mod_rng
  !$ use omp_lib
  implicit none
  !> inputs-outputs:
  type(particle_sim),intent(inout) :: sim
  !> inputs:
  class(type_rng),intent(in)       :: rng_base
  integer,intent(in)               :: my_id,n_cpus,n_vec
  integer,intent(in)               :: n_points,n_phi
  real*8,intent(in)                :: time,mass,charge
  real*8,dimension(2),intent(in)   :: phi_interval
  real*8,dimension(n_points),intent(in)       :: soft_orbit_ppar
  real*8,dimension(n_points),intent(in)       :: soft_orbit_pperp
  real*8,dimension(n_points),intent(in)       :: soft_orbit_Jdtdrho
  real*8,dimension(n_vec,n_points),intent(in) :: soft_orbit_x

  !> variables:
  class(type_rng),dimension(:),allocatable :: rngs
  integer   :: ii,jj,i_elm,ifail,n_threads,thread_id
  real*8                  :: U,psi
  real*8,dimension(2)     :: st,Rbox,Zbox
  real*8,dimension(3)     :: RZphi,B_field,E_field
  real*8,dimension(n_phi) :: phi_array
  !> initialise the random number generator
  n_threads = 1
!$ n_threads = omp_get_max_threads()
  allocate(rngs(n_threads),source=rng_base) 
  do ii=1,n_threads
    call rngs(ii)%initialize(n_phi,random_seed(),n_cpus*n_threads,my_id*n_threads+ii,ifail)
    if(ifail.ne.0) call MPI_Abort(MPI_COMM_WORLD,-1,ifail)
  enddo
  !> initialise and allocate particle simulation array
  call domain_bounding_box(sim%fields%node_list,sim%fields%element_list,&
  Rbox(1),Rbox(2),Zbox(1),Zbox(2))
  sim%time = time; sim%groups(1)%mass = mass; 
  allocate(particle_gc_relativistic::sim%groups(1)%particles(n_points*n_phi))
  !> loop on the soft orbits
  !$omp parallel default(shared) firstprivate(n_points,n_vec,time,charge,n_phi,mass,&
  !$omp Rbox,Zbox,phi_interval) private(ii,jj,RZphi,i_elm,st,B_field,E_field,U,&
  !$omp psi,ifail,phi_array,thread_id)
  thread_id = 1
  !$ thread_id = omp_get_thread_num()+1
  !$omp do
  do ii=1,n_points
    !> transform the soft orbit coordinates in jorek global/local coordinates
    RZphi = cartesian_to_cylindrical(soft_orbit_x(:,ii)); i_elm = -1;
    if(((RZphi(1).ge.Rbox(1)).and.(RZphi(1).le.Rbox(2))).and.((RZphi(2).ge.Zbox(1)).and.(RZphi(2).le.Zbox(2)))) &
    call find_RZ(sim%fields%node_list,sim%fields%element_list,RZphi(1),RZphi(2),&
    RZphi(1),RZphi(2),i_elm,st(1),st(2),ifail)
    call rngs(thread_id)%next(phi_array); 
    phi_array = phi_interval(1)+(phi_interval(2)-phi_interval(1))*phi_array;
    do jj=1,n_phi
      select type (p=>sim%groups(1)%particles((ii-1)*n_phi+jj))
        type is (particle_gc_relativistic)
        if(i_elm.le.0) then
          p%x = 0d0; p%st = 0d0; p%weight = 0d0; p%i_elm = 0; p%i_life = 0; p%t_birth = 0.;
          p%p = 0d0; p%q = int(0,kind=1);
        else
          p%x = [RZphi(1),RZphi(2),phi_array(jj)]; p%st = st; p%weight = 1d0;
          p%i_elm = i_elm; p%i_life = 0; p%t_birth = 0.;
          call sim%fields%calc_EBpsiU(time,i_elm,st,RZphi(3),E_field,B_field,psi,U)
          p%p = SPEED_OF_LIGHT*mass*[soft_orbit_ppar(ii),(5d-1*SPEED_OF_LIGHT*&
          (soft_orbit_pperp(ii)**2))/(norm2(B_field))]; p%q = int(charge,kind=1);
        endif
      end select
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  !> clean-up
  if(allocated(rngs)) deallocate(rngs)
end subroutine convert_soft_orbits_in_jorek_relativistic_gcs

!> scatter 2D array between MPI processes, the scattering is performed
!> along the second index
!> inputs:
!>   my_id:        (integer) MPI task rank
!>   n_cpus:       (integer) size of the MPI communicator
!>   n1:           (integer) first index size
!>   n2_loc:       (integer) local array second index size
!>   n2_glob:      (integer) global array second index size
!>   global_array: (n1,n2_glob) array to be scattered
!> outputs:
!>   local_array: (n1,n2_loc) array receiving the scattered global array
subroutine scatter_2D_arrays(my_id,n_cpus,n1,n2_loc,n2_glob,global_array,local_array)
  use mpi
  implicit none
  !> inputs:
  integer,intent(in) :: my_id,n_cpus,n1,n2_loc,n2_glob
  !> inputs-outputs:
  real*8,dimension(:,:),allocatable,intent(inout) :: global_array
  !> outputs:
  real*8,dimension(n1,n2_loc),intent(out) :: local_array
  !> variables:
  logical :: did_allocate
  integer :: ii,doublesize,subarraytype,resizedsubarraytype,errorcode,ierr
  integer(kind=MPI_Address_kind) :: startresized,extent
  integer,dimension(n_cpus)      :: disps,counts
  !> check consistency
  if(n2_loc.le.0) then
    write(*,*) 'Error: invalid size of the receiving array during scattering: abort!'
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  endif
  !> create a vector for each subblock and scatter them
  did_allocate = .false.
  if(my_id.eq.0) then
    counts = 1; disps = [(n2_loc*ii,ii=0,n_cpus-1)]; startresized = 0;
    call MPI_Type_size(MPI_REAL8,doublesize,ierr); extent = n1*doublesize;
    call MPI_Type_create_subarray(2,[n1,n2_glob],[n1,n2_loc],[0,0],MPI_ORDER_FORTRAN,MPI_REAL8,subarraytype,ierr)
    call MPI_Type_create_resized(subarraytype,startresized,extent,resizedsubarraytype,ierr)
    call MPI_Type_commit(resizedsubarraytype,ierr)
  endif
  local_array = 0d0; 
  if(my_id.ne.0) then
    if(.not.allocated(global_array)) allocate(global_array(1,1)); did_allocate = .true.;
  endif
  call MPI_scatterv(global_array,counts,disps,resizedsubarraytype,local_array,n1*n2_loc,MPI_REAL8,&
  0,MPI_COMM_WORLD,ierr)
  if(my_id.eq.0) call MPI_Type_free(resizedsubarraytype,ierr)
  if(did_allocate) deallocate(global_array)
end subroutine scatter_2D_arrays

!> read SOFT orbit file
!> inputs:
!>   soft_orbit_filename_in: (character)(*) name of the soft orbit file
!>   n_vec:                  (integer) size of the position vector
!> outputs:
!>   n_soft_points: (integer) total number of valid soft orbits
!>   x:             (real8)(3,n_soft_particles) soft positions in xyz coordinates
!>   ppar:          (real8)(n_soft_particles) soft parallel momentum
!>   pperp:         (real8)(n_soft_particles) soft perpendicular momentum
!>   Jdtdrho:       (real8)(n_soft_particles) jacobian*dpoloidal*dminor_radius
subroutine read_soft_orbit_file(soft_orbit_filename_in,n_vec,n_soft_points,x,&
ppar,pperp,Jdtdrho)
  use constants, only: ATOMIC_MASS_UNIT
  use hdf5
  use hdf5_io_module, only: HDF5_open,HDF5_close
  use hdf5_io_module, only: HDF5_allocatable_array2D_reading
  implicit none
  !> inputs:
  character(len=*),intent(in) :: soft_orbit_filename_in
  integer,intent(in)          :: n_vec
  !> outputs:
  integer,intent(out) :: n_soft_points
  real*8,dimension(:),allocatable,intent(out)   :: ppar,pperp,Jdtdrho
  real*8,dimension(:,:),allocatable,intent(out) :: x
  !> variables:
  integer(HID_T) :: file_id
  integer        :: ii,n_orbits,n_times,n_active_orbits,ierr
  integer,dimension(:),allocatable   :: valid_orbit_id
  real*8,dimension(:,:),allocatable  :: x_loc,ppar_loc,pperp_loc,Jdtdrho_loc
  !> open the soft orbit hdf5 file
  call HDF5_open(trim(soft_orbit_filename_in)//".h5",file_id,ierr)
  !> read sofit particles
  call HDF5_allocatable_array2D_reading(file_id,ppar_loc,"/ppar")       !< parallel momentum
  call HDF5_allocatable_array2D_reading(file_id,pperp_loc,"/pperp")     !< perpendicular momentum
  call HDF5_allocatable_array2D_reading(file_id,Jdtdrho_loc,"/Jdtdrho") !< jacobia*dpoloidal*dminorradius
  call HDF5_allocatable_array2D_reading(file_id,x_loc,"/x")             !< position in xyz coordinates
  !> close the soft orbit hdf5 file
  call HDF5_close(file_id)
  !> remove zero orbits
  n_orbits = size(ppar_loc,2); n_times  = size(ppar_loc,1);
  n_soft_points = n_orbits*n_times; n_active_orbits = 0;
  !> reshape arrays
  x_loc       = reshape(x_loc,[n_vec,n_soft_points])
  ppar_loc    = reshape(ppar_loc,[n_soft_points,1])
  pperp_loc   = reshape(pperp_loc,[n_soft_points,1])
  Jdtdrho_loc = reshape(Jdtdrho_loc,[n_soft_points,1])
  !> find id of active orbits
  allocate(valid_orbit_id(n_soft_points)); valid_orbit_id  = 0;
  do ii=1,n_soft_points
    if((x_loc(1,ii).eq.0d0).and.(x_loc(2,ii).eq.0d0).and.(x_loc(3,ii).eq.0d0)) cycle
    if((ppar_loc(ii,1).eq.0d0).and.(pperp_loc(ii,1).eq.0d0)) cycle
    n_active_orbits = n_active_orbits + 1
    valid_orbit_id(n_active_orbits) = ii;
  enddo
  !> copy only active orbits
  allocate(x(n_vec,n_active_orbits)); x     = 0d0;
  allocate(ppar(n_active_orbits));    ppar  = 0d0; 
  allocate(pperp(n_active_orbits));   pperp = 0d0;
  allocate(Jdtdrho(n_active_orbits)); Jdtdrho = 0d0;
  x       = x_loc(:,valid_orbit_id(1:n_active_orbits))
  ppar    = ppar_loc(valid_orbit_id(1:n_active_orbits),1)
  pperp   = pperp_loc(valid_orbit_id(1:n_active_orbits),1)
  Jdtdrho = Jdtdrho_loc(valid_orbit_id(1:n_active_orbits),1)
  n_soft_points = n_active_orbits
  if(allocated(valid_orbit_id))  deallocate(valid_orbit_id)
  if(allocated(ppar_loc))        deallocate(ppar_loc)
  if(allocated(pperp_loc))       deallocate(pperp_loc)
  if(allocated(Jdtdrho_loc))     deallocate(Jdtdrho_loc)
  if(allocated(x_loc))           deallocate(x_loc)
end subroutine read_soft_orbit_file

!> allocate soft pdf type and initialize to 0
!> inputs:
!>   soft_pdf: (type_soft_pdf) soft_pdf to be initialized
!> outputs:
!>   soft_pdf: (type_soft_pdf) initialized soft_pdf
subroutine init_soft_pdf(n_pxi,soft_pdf)
  implicit none
  !> inputs-outputs:
  type(type_soft_pdf),intent(inout) :: soft_pdf
  integer,dimension(2),intent(in)   :: n_pxi
  !> allocate soft_pdf and initialize to 0
  if(allocated(soft_pdf%p)) deallocate(soft_pdf%p); 
  allocate(soft_pdf%p(n_pxi(1))); soft_pdf%p=0d0;
  if(allocated(soft_pdf%xi)) deallocate(soft_pdf%xi); 
  allocate(soft_pdf%xi(n_pxi(2))); soft_pdf%xi=0d0;
  if(allocated(soft_pdf%pdf)) deallocate(soft_pdf%pdf); 
  allocate(soft_pdf%pdf(n_pxi(1),n_pxi(2))); soft_pdf%pdf=0d0;
end subroutine init_soft_pdf

!> read and distribute the pdf used in SOFT and its mesh (r,xi,p)
!> inputs:
!>   soft_pdf_filename_in: (charcater)(*) soft pdf filename
!> outputs:
!>   n_r_mesh: (integer) size of the pdf minor radius mesh
!>   r_mesh:   (real8)(n_r_mesh) pdf minor radius mesh
!>   soft_pdf: (type_soft_pdf)(n_r_mesh) the pdf used by soft with the xi,p meshes
subroutine read_soft_pdf_file(soft_pdf_filename_in,n_r_mesh,r_mesh,soft_pdf)
  use hdf5
  use hdf5_io_module, only: HDF5_open,HDF5_close
  use hdf5_io_module, only: HDF5_allocatable_array1D_reading
  use hdf5_io_module, only: HDF5_allocatable_array2D_reading
  implicit none
  !> inputs:
  character(len=*),intent(in) :: soft_pdf_filename_in
  !> outputs:
  integer,intent(out) :: n_r_mesh
  real*8,dimension(:),allocatable,intent(out)              :: r_mesh
  type(type_soft_pdf),dimension(:),allocatable,intent(out) :: soft_pdf
  !> variables:
  character(len=10) :: format_char
  character(len=:),allocatable  :: group_name
  integer(HID_T)    :: file_id
  integer           :: ii,ierr,r_id,n_r_id,group_name_len
  !> open hdf5 file
  call HDF5_open(trim(soft_pdf_filename_in)//".h5",file_id,ierr)
  !> read mesh datasets
  call HDF5_allocatable_array1D_reading(file_id,r_mesh,"r")
  !> allocate soft pdf data strucutre
  n_r_mesh = size(r_mesh); allocate(soft_pdf(n_r_mesh))
  !> read soft pdf
  do ii=1,n_r_mesh
    !> define the group name compatible with soft nomenclature
    r_id = ii-1
    n_r_id = int(log10(real(r_id)))+1
    if(r_id.eq.0) n_r_id = 1
    write(format_char,'(A,I1,A)') "(A,I",n_r_id,")"
    group_name_len = 2+n_r_id
    allocate(character(len=group_name_len)::group_name)
    write(group_name,trim(format_char)) "/r",r_id
    !> read the pdf data 
    call HDF5_allocatable_array1D_reading(file_id,soft_pdf(ii)%p,trim(group_name)//"/p")
    call HDF5_allocatable_array1D_reading(file_id,soft_pdf(ii)%xi,trim(group_name)//"/xi")
    call HDF5_allocatable_array2D_reading(file_id,soft_pdf(ii)%pdf,trim(group_name)//"/f")
    deallocate(group_name)
  enddo
  !> close hdf5 file
  call HDF5_close(file_id)
end subroutine read_soft_pdf_file

!> broadcast the soft pdf list and meshes
!>   n_r:           (integer) number of soft minor radii
!>   r_mesh:        (real8)(n_r) soft pdf minor radii
!>   soft_pdf_list: (type_pdf_list)(n_r) soft pdf for various minor radii
!> inputs:
!>   n_r:           (integer) number of soft minor radii
!>   r_mesh:        (real8)(n_r) soft pdf minor radii
!>   soft_pdf_list: (type_pdf_list)(n_r) soft pdf for various minor radii
!> outputs:
subroutine broadcast_soft_pdf_list(my_id,n_r,r_mesh,soft_pdf_list)
  use mpi
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: n_r
  real*8,dimension(:),allocatable,intent(inout) :: r_mesh
  type(type_soft_pdf),dimension(:),allocatable,intent(inout) :: soft_pdf_list
  !> inputs:
  integer,intent(in) :: my_id
  !> variables:
  integer :: ii
  integer,dimension(2)               :: n_pxi_mesh_global,n_pxi_sum
  real*8,dimension(:),allocatable    :: p_mesh_global,xi_mesh_global
  real*8,dimension(:,:),allocatable  :: pdf_global
  integer,dimension(:,:),allocatable :: dims_pxi
  !> broadcast the minor radius size and allocate p and xi mesh sizes
  call MPI_Bcast(n_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  allocate(dims_pxi(2,n_r)); dims_pxi = 0;
  !> extract and broadcast the pdf sizes
  if(my_id.eq.0) then
    do ii=1,n_r
      dims_pxi(:,ii) = [size(soft_pdf_list(ii)%p),size(soft_pdf_list(ii)%xi)]
    enddo
  endif
  call MPI_Bcast(dims_pxi,n_r*2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  !> allocate global pdf, r and xi arrays
  n_pxi_mesh_global = sum(dims_pxi,dim=2);
  allocate(p_mesh_global(n_pxi_mesh_global(1))); p_mesh_global    = 0d0;
  allocate(xi_mesh_global(n_pxi_mesh_global(2))); xi_mesh_global = 0d0;
  allocate(pdf_global(n_pxi_mesh_global(1),n_pxi_mesh_global(2))); pdf_global = 0d0;
  !> allocate soft_pdf_list for all ranks except master
  if(my_id.ne.0) then
    if(allocated(soft_pdf_list)) deallocate(soft_pdf_list); allocate(soft_pdf_list(n_r));
    do ii=1,n_r
      call init_soft_pdf(dims_pxi(:,ii),soft_pdf_list(ii))
    enddo
  endif
  !> construct the global arrays to be sent
  if(my_id.eq.0) then
    n_pxi_sum = 0;
    do ii=1,n_r
      p_mesh_global(n_pxi_sum(1)+1:n_pxi_sum(1)+dims_pxi(1,ii))  = soft_pdf_list(ii)%p 
      xi_mesh_global(n_pxi_sum(2)+1:n_pxi_sum(2)+dims_pxi(2,ii)) = soft_pdf_list(ii)%xi
      pdf_global(n_pxi_sum(1)+1:n_pxi_sum(1)+dims_pxi(1,ii),&
      n_pxi_sum(2)+1:n_pxi_sum(2)+dims_pxi(2,ii)) = soft_pdf_list(ii)%pdf
      n_pxi_sum = n_pxi_sum + dims_pxi(:,ii);
    enddo
  endif
  !> broadcast global arrays
  call MPI_Bcast(p_mesh_global,n_pxi_mesh_global(1),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(xi_mesh_global,n_pxi_mesh_global(2),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(pdf_global,n_pxi_mesh_global(1)*n_pxi_mesh_global(2),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  !> store the soft pdf into structure
  if(my_id.ne.0) then
    n_pxi_sum = 0;
    do ii=1,n_r
      soft_pdf_list(ii)%p   = p_mesh_global(n_pxi_sum(1)+1:n_pxi_sum(1)+dims_pxi(1,ii))
      soft_pdf_list(ii)%xi  = xi_mesh_global(n_pxi_sum(2)+1:n_pxi_sum(2)+dims_pxi(2,ii))
      soft_pdf_list(ii)%pdf = pdf_global(n_pxi_sum(1)+1:n_pxi_sum(1)+dims_pxi(1,ii),&
      n_pxi_sum(2)+1:n_pxi_sum(2)+dims_pxi(2,ii))
      n_pxi_sum = n_pxi_sum + dims_pxi(:,ii);
    enddo
  endif
  !> clean up
  if(allocated(dims_pxi))       deallocate(dims_pxi)
  if(allocated(p_mesh_global))  deallocate(p_mesh_global)
  if(allocated(xi_mesh_global)) deallocate(xi_mesh_global)
  if(allocated(pdf_global))     deallocate(pdf_global)
end subroutine broadcast_soft_pdf_list

!> read and broadcast the soft minor radii
!> inputs:
!>   my_id:       (integer) mpi task number
!>   filename_in: (character) name of the soft magnetic field file
!> outputs:
!>   R_mesh:  (real8)(nR) soft major radius mesh
!>   Z_mesh:  (real8)(nZ) soft vertical position mesh
!>   r_minor: (real8)(nZ,nR) soft minor radii
subroutine read_and_broadcast_soft_minor_radii(my_id,filename_in,R_mesh,Z_mesh,r_minor)
  use hdf5
  use hdf5_io_module, only: HDF5_open,HDF5_close
  use hdf5_io_module, only: HDF5_array1D_reading
  use hdf5_io_module, only: HDF5_allocatable_array1D_reading
  use hdf5_io_module, only: HDF5_allocatable_array2D_reading
  implicit none
  !> inputs:
  character(len=*),intent(in) :: filename_in
  integer,intent(in)          :: my_id
  !> outputs:
  real*8,dimension(:),allocatable,intent(out)   :: R_mesh,Z_mesh
  real*8,dimension(:,:),allocatable,intent(out) :: r_minor
  !> variables:
  integer(HID_T) :: file_id 
  integer        :: n_separatrix,ierr
  integer,dimension(2) :: id_R,id_Z
  integer,dimension(4) :: dims
  real*8               :: poloidal_flux_axis,poloidal_flux_bnd
  real*8,dimension(2)  :: RZ_axis,Rnodes,Znodes
  real*8,dimension(4)  :: poloidal_flux_nodes
  real*8,dimension(:,:),allocatable :: separatrix
  if(my_id.eq.0) then
    !> open hdf5 file
    call HDF5_open(trim(filename_in)//".h5",file_id,ierr)
    call HDF5_array1D_reading(file_id,RZ_axis,"/maxis")
    call HDF5_allocatable_array1D_reading(file_id,R_mesh,"/r")
    call HDF5_allocatable_array1D_reading(file_id,Z_mesh,"/z")
    call HDF5_allocatable_array2D_reading(file_id,r_minor,"/Psi")
    call HDF5_allocatable_array2D_reading(file_id,separatrix,"/separatrix")
    !> close hdf5 file
    call HDF5_close(file_id)
    dims(1:2) = shape(r_minor); dims(3:4) = shape(separatrix); 
  endif
  call MPI_Bcast(dims,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if(.not.allocated(Z_mesh))     allocate(Z_mesh(dims(1)))
  if(.not.allocated(R_mesh))     allocate(R_mesh(dims(2)))
  if(.not.allocated(r_minor))    allocate(r_minor(dims(1),dims(2)))
  if(.not.allocated(separatrix)) allocate(separatrix(dims(3),dims(4)))
  call MPI_Bcast(RZ_axis,2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Z_mesh,dims(1),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(R_mesh,dims(2),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(r_minor,dims(1)*dims(2),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(separatrix,dims(3)*dims(4),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  !> find poloidal flux at the magnetic axis
  call find_point_segment(RZ_axis(1),dims(2),R_mesh,id_R,Rnodes)
  call find_point_segment(RZ_axis(2),dims(1),Z_mesh,id_Z,Znodes)
  poloidal_flux_nodes = [r_minor(id_Z(1),id_R(1)),r_minor(id_Z(1),id_R(2)),\
                         r_minor(id_Z(2),id_R(2)),r_minor(id_Z(2),id_R(1))]
  poloidal_flux_axis = bilinear_interp_RZ(RZ_axis(1),RZ_axis(2),\
                       Rnodes,Znodes,poloidal_flux_nodes)
  !> find the poloidal flux at the plasma boundary
  call find_point_segment(separatrix(1,1),dims(2),R_mesh,id_R,Rnodes)
  call find_point_segment(separatrix(2,1),dims(1),Z_mesh,id_Z,Znodes)
  poloidal_flux_nodes = [r_minor(id_Z(1),id_R(1)),r_minor(id_Z(1),id_R(2)),\
                         r_minor(id_Z(2),id_R(2)),r_minor(id_Z(2),id_R(1))]
  poloidal_flux_bnd = bilinear_interp_RZ(separatrix(1,1),separatrix(2,1),\
                      Rnodes,Znodes,poloidal_flux_nodes)
  !> compute the magnetic minor radii
  r_minor = sqrt((r_minor-poloidal_flux_axis)/(poloidal_flux_bnd-poloidal_flux_axis))
  !> cleanup
  if(allocated(separatrix)) deallocate(separatrix)
end subroutine read_and_broadcast_soft_minor_radii

!> identify mesh element containing a target point
!> inputs:
!>   x:      (real8) target point
!>   nx:     (integer) number of mesh nodes
!>   x_mesh: (real8)(nx) 1D mesh
!> outputs:
!>   id_x:   (integer)(2) indexes of the nodes containing the target point
!>   xnodes: (real8)(2) min,max of the value
subroutine find_point_segment(x,nx,x_mesh,id_x,xnodes)
  implicit none
  !> inputs:
  integer,intent(in) :: nx
  real*8,intent(in)  :: x 
  real*8,dimension(nx),intent(in) :: x_mesh
  !> outputs:
  integer,dimension(2),intent(out) :: id_x
  real*8,dimension(2),intent(out)  :: xnodes
  !> variables:
  !> identify the nearest node
  id_x = [minloc(abs(x-x_mesh)),-1];
  if(id_x(1).le.1) then 
    id_x = [1,2];
  elseif(id_x(1).ge.nx) then
    id_x = [nx-1,nx]
  elseif(x.ge.x_mesh(id_x(1))) then
    id_x = [id_x(1),id_x(1)+1]
  elseif(x.lt.x_mesh(id_x(1))) then
    id_x = [id_x(1)-1,id_x(1)]
  endif
  xnodes = x_mesh(id_x(1):id_x(2))
end subroutine find_point_segment

!> perform a bilinear interpolation in the RZ coordinates
!> 4 - 3
!> |   |
!> 1 - 2
!> inputs:
!>   R:      (real8) target major radius
!>   Z:      (real8) target vertical position
!>   Rnodes: (real8)(2) min,max major radius value
!>   Znodes: (real8)(2) min,max vertical position value
!>   fnodes: (real8)(4) f-function value at the 4 nodes
!> outputs:
!>   f:      (real8) interpolated value
function bilinear_interp_RZ(R,Z,Rnodes,Znodes,fnodes) result(f)
  implicit none
  !> inputs:
  real*8,intent(in) :: R,Z
  real*8,dimension(2),intent(in) :: Rnodes,Znodes
  real*8,dimension(4),intent(in) :: fnodes
  !> outputs:
  real*8 :: f
  !> compute bilinear interpolation
  f = (fnodes(1)*(Rnodes(2)-R)*(Znodes(2)-Z) + &
      fnodes(2)*(R-Rnodes(1))*(Znodes(2)-Z) + &
      fnodes(3)*(R-Rnodes(1))*(Z-Znodes(1)) + &
      fnodes(4)*(Rnodes(2)-R)*(Z-Znodes(1)))/&
      ((Rnodes(2)-Rnodes(1))*(Znodes(2)-Znodes(1)))
end function bilinear_interp_RZ

!> write jorek particles in HDF5
!> inputs:
!>   my_id:    (integer)
!>   filename: (character) hdf5 filename
!>   n_cpus:   (integer) total number of mpi tasks
!>   n_vec:    (integer) size of the x position vector
!>   sim:      (particle_sim) jorek particle simulation
!> outputs:
subroutine write_particles_in_hdf5(my_id,filename,n_cpus,n_vec,sim)
  use constants,      only: SPEED_OF_LIGHT
  use mpi
  use hdf5
  use hdf5_io_module, only: HDF5_open_or_create,HDF5_close
  use hdf5_io_module, only: HDF5_array2D_saving,HDF5_array1D_saving
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  implicit none
  !> inputs:
  type(particle_sim),intent(in) :: sim
  character(len=*),intent(in)   :: filename
  integer,intent(in)            :: my_id,n_cpus,n_vec
  !> variables:
  integer(HID_T) :: file_id
  integer        :: ii,n_particles,ierr
  real*8         :: U,psi
  real*8,dimension(n_vec) :: Bvec,Evec
  real*8,dimension(:,:),allocatable :: x_pos,x_pos_glob
  real*8,dimension(:),allocatable   :: ppar,pperp,ppar_glob,pperp_glob
  !> initialisation
  n_particles = size(sim%groups(1)%particles);
  !> compute position in cartesian coordinates, ppar and pperp
  allocate(x_pos(n_vec,n_particles)); x_pos = 0d0;
  allocate(ppar(n_particles));        ppar  = 0d0;
  allocate(pperp(n_particles));       pperp = 0d0;
  allocate(x_pos_glob(n_vec,n_cpus*n_particles)); x_pos_glob = 0d0;
  allocate(ppar_glob(n_cpus*n_particles));        ppar_glob  = 0d0;
  allocate(pperp_glob(n_cpus*n_particles));       pperp_glob = 0d0; 
  !$omp parallel do default(shared) firstprivate(n_particles) &
  !$omp private(ii,Evec,Bvec,psi,U)
  do ii=1,n_particles
    if(sim%groups(1)%particles(ii)%i_elm.le.0) cycle
    x_pos(:,ii) = cylindrical_to_cartesian(sim%groups(1)%particles(ii)%x)
    call sim%fields%calc_EBpsiU(sim%time,sim%groups(1)%particles(ii)%i_elm,\
    sim%groups(1)%particles(ii)%st,sim%groups(1)%particles(ii)%x(3),Evec,Bvec,psi,U)
    select type (part=>sim%groups(1)%particles(ii))
    type is (particle_gc_relativistic)
      ppar(ii) = part%p(1)/(sim%groups(1)%mass*SPEED_OF_LIGHT)
      pperp(ii) = sqrt(2d0*sim%groups(1)%mass*norm2(Bvec)*part%p(2))/\
      (sim%groups(1)%mass*SPEED_OF_LIGHT)
    end select
  enddo
  !$omp end parallel do
  !> gather data from all mpi tasks 
  call MPI_gather(x_pos,n_vec*n_particles,MPI_REAL8,x_pos_glob,n_vec*n_particles,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_gather(ppar,n_particles,MPI_REAL8,ppar_glob,n_particles,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_gather(pperp,n_particles,MPI_REAL8,pperp_glob,n_particles,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  !> write data in hdf5
  if(my_id.eq.0) then
    call HDF5_open_or_create(trim(filename//'.h5'),H5P_DEFAULT_F,file_id,ierr,H5F_ACC_TRUNC_F)
    call HDF5_array2D_saving(file_id,x_pos_glob,n_vec,n_cpus*n_particles,'x')
    call HDF5_array1D_saving(file_id,ppar_glob,n_cpus*n_particles,'ppar')
    call HDF5_array1D_saving(file_id,pperp_glob,n_cpus*n_particles,'pperp')
    call HDF5_close(file_id)
  endif
  !> cleanup
  if(allocated(x_pos)) deallocate(x_pos); if(allocated(x_pos_glob)) deallocate(x_pos_glob);
  if(allocated(ppar))  deallocate(ppar);  if(allocated(ppar_glob))  deallocate(ppar_glob);
  if(allocated(pperp)) deallocate(pperp); if(allocated(pperp_glob)) deallocate(pperp_glob);
end subroutine write_particles_in_hdf5

!> -----------------------------------------------------------------------------------------

end program generate_particle_restart_from_SOFT
 
