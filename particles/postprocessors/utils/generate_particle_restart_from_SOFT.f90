!> generate_particle_restart_from_SOFT generates a 
!> JOREK particle restart file from a orbits file
!> produced by the SOFT code
program generate_particle_restart_from_SOFT
use constants,       only: PI,SPEED_OF_LIGHT,ATOMIC_MASS_UNIT,EL_CHG
use mod_mpi_tools,   only: init_mpi_threads,finalize_mpi_threads
use mod_particle_io, only: write_simulation_hdf5
use particle_tracer

implicit none

!> Variables -------------------------------------------------------------------------------
type(pcg32_rng)      :: rng_type
type(event)          :: field_reader
logical              :: do_write_particles_in_hdf5
integer              :: my_id,n_cpus,ierr,n_vec,n_groups,n_phi
integer,dimension(2) :: dims
real*8               :: time,mass,charge
real*8,dimension(2)  :: phi_interval
real*8,dimension(:),allocatable   :: soft_orbit_ppar_local,soft_orbit_pperp_local
real*8,dimension(:),allocatable   :: soft_orbit_ppar,soft_orbit_pperp
real*8,dimension(:,:),allocatable :: soft_orbit_x,soft_orbit_x_local
character(len=17)    :: fields_filename
character(len=125)   :: particle_filename
character(len=193)   :: soft_orbit_filename
character(len=25)    :: filename_jorek_hdf5
!> Variable definitions --------------------------------------------------------------------
do_write_particles_in_hdf5 = .true.       !< writeh particle in hdf5 if true
n_groups            = 1                   !< number of jorek particle groups
n_vec               = 3                   !< component of a vector
n_phi               = 1                   !< number of sampled toroidal positions to be sampled for each particle
time                = 0d0                 !< simulation time
mass                = 5.48579909065d-4    !< electron mass in AMU 
charge              = -1d0                !< electron charge
phi_interval        = [0d0,PI]            !< toroidal angle interval in which particles are sampled
fields_filename     = 'jorek_equilibrium' !< jorek restart filename
particle_filename   = 'part_restart_soft_orbits.h5' !< particle restart filename 
soft_orbit_filename = 'orbit_test_jorek_JET_pulse95135_t48dot54_parabolic_qprofile_q95_6dot8_press0_res1r5dot88en1m5_res2r4dot705en1m4_Ip612en1MA_Ekin20MeV_np10_theta1_2dot85_itheta2_pi_nitheta100_norbits100_a96_wall' !< soft orbit filename
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
write(*,*) "Reading SOFT orbit file ..."
if(my_id.eq.0) then
  call read_soft_orbit_file(soft_orbit_filename,n_vec,dims(1),soft_orbit_x,soft_orbit_ppar,soft_orbit_pperp)
  !> scatter the global arrays to each mpi process)
  dims(2) = dims(1)/n_cpus;
  if(dims(2)*n_cpus.lt.dims(1)) dims(2) = dims(2)+1
endif
call MPI_Bcast(dims,size(dims),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
allocate(soft_orbit_x_local(n_vec,dims(2))); soft_orbit_x_local     = 0d0;
allocate(soft_orbit_ppar_local(dims(2)));    soft_orbit_ppar_local  = 0d0;
allocate(soft_orbit_pperp_local(dims(2)));   soft_orbit_pperp_local = 0d0;
call scatter_2D_arrays(my_id,n_cpus,n_vec,dims(2),dims(1),soft_orbit_x,soft_orbit_x_local)
call MPI_Scatter(soft_orbit_ppar,dims(2),MPI_REAL8,soft_orbit_ppar_local,dims(2),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_Scatter(soft_orbit_pperp,dims(2),MPI_REAL8,soft_orbit_pperp_local,dims(2),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if(allocated(soft_orbit_x))     deallocate(soft_orbit_x)
if(allocated(soft_orbit_ppar))  deallocate(soft_orbit_ppar)
if(allocated(soft_orbit_pperp)) deallocate(soft_orbit_pperp)
write(*,*) "Reading SOFT orbit file: completed!"
!> Generate JOREK relativistic gc from SOFT orbits -----------------------------------------
write(*,*) 'Converting soft orbit to JOREK relativistic gc ...'
call convert_soft_orbits_in_jorek_relativistic_gcs(sim,rng_type,my_id,n_cpus,time,\
mass,charge,n_vec,dims(2),n_phi,phi_interval,soft_orbit_x_local,\
soft_orbit_ppar_local,soft_orbit_pperp_local)
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
if(allocated(soft_orbit_x))           deallocate(soft_orbit_x)
if(allocated(soft_orbit_ppar))        deallocate(soft_orbit_ppar)
if(allocated(soft_orbit_pperp))       deallocate(soft_orbit_pperp)
if(allocated(soft_orbit_x_local))     deallocate(soft_orbit_x_local)
if(allocated(soft_orbit_ppar_local))  deallocate(soft_orbit_ppar_local)
if(allocated(soft_orbit_pperp_local)) deallocate(soft_orbit_pperp_local)
call finalize_mpi_threads(ierr)
write(*,*) 'Program terminated: good bye!'

contains

!> generate jorek particle relativistic gc from soft orbits 
!> inputs:
!>   sim:              (particle_sim) particle simulation type to be initialised
!>   rng_base:         (type_rng) type of the random number generator
!>   my_id:            (integer) MPI task rank
!>   n_cpus:           (integer) number of MPI tasks
!>   time:             (real8) simulation time
!>   mass:             (real8) particle mass in AMU
!>   charge:           (real8) particle charge (RE: -1)
!>   n_vec:            (integer) size of the position vector
!>   n_points:         (integer) number of soft valid points
!>   n_phi:            (integer) number of toroidal samples per particles
!>   phi_interval:     (real8)(2) toroidal angle interval for sampling
!>   soft_orbit_x:     (real8)(n_vec,n_points) soft positions in xyz
!>   soft_orbit_ppar:  (real8)(n_points) soft parallel momentum
!>   soft_orbit_pperp: (real8)(n_points) soft perpendicular momentum
!> outpus:
!>   sim: (particle_sim) initialised particle simulation
subroutine convert_soft_orbits_in_jorek_relativistic_gcs(sim,rng_base,my_id,&
n_cpus,time,mass,charge,n_vec,n_points,n_phi,phi_interval,&
soft_orbit_x,soft_orbit_ppar,soft_orbit_pperp)
  use mpi
  use constants,                 only: SPEED_OF_LIGHT
  use mod_coordinate_transforms, only: cartesian_to_cylindrical
  use mod_random_seed,           only: random_seed
  use mod_rng
  !$ use omp_lib
  implicit none
  !> parameters:
  integer,parameter :: intkind=8
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
  real*8,dimension(n_vec,n_points),intent(in) :: soft_orbit_x

  !> variables:
  class(type_rng),dimension(:),allocatable :: rngs
  integer   :: ii,jj,i_elm,ifail,n_threads,thread_id
  integer(kind=intkind)   :: id,one
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
  one = int(1,kind=intkind)
  call domain_bounding_box(sim%fields%node_list,sim%fields%element_list,&
  Rbox(1),Rbox(2),Zbox(1),Zbox(2))
  sim%time = time; sim%groups(1)%mass = mass; 
  id = int(n_points,kind=intkind)*int(n_phi,kind=intkind);
  allocate(particle_gc_relativistic::sim%groups(1)%particles(id))
  !> loop on the soft orbits
  !$omp parallel default(shared) firstprivate(n_points,n_vec,time,charge,n_phi,mass,&
  !$omp Rbox,Zbox,phi_interval,one) private(ii,jj,RZphi,i_elm,st,B_field,E_field,U,&
  !$omp psi,ifail,phi_array,thread_id,id)
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
      id = (int(ii,kind=intkind)-one)*int(n_phi,kind=intkind)+int(jj,kind=intkind)
      select type (p=>sim%groups(1)%particles(id))
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
!>   soft_orbit_filename: (character)(*) name of the soft orbit file
!>   n_vec:               (integer) size of the position vector
!> outputs:
!>   n_soft_points: (integer) total number of valid soft orbits
!>   x:             (real8)(3,n_soft_particles) soft positions in xyz coordinates
!>   ppar:          (real8)(n_soft_particles) soft parallel momentum
!>   pperp:         (real8)(n_soft_particles) soft perpendicular momentum
subroutine read_soft_orbit_file(soft_orbit_filename,n_vec,n_soft_points,x,ppar,pperp)
  use constants, only: ATOMIC_MASS_UNIT
  use hdf5
  use hdf5_io_module, only: HDF5_open,HDF5_close
  use hdf5_io_module, only: HDF5_allocatable_array2D_reading
  implicit none
  !> inputs:
  character(len=*),intent(in) :: soft_orbit_filename
  integer,intent(in)          :: n_vec
  !> outputs:
  integer,intent(out) :: n_soft_points
  real*8,dimension(:),allocatable,intent(out)   :: ppar,pperp
  real*8,dimension(:,:),allocatable,intent(out) :: x
  !> variables:
  integer(HID_T)          :: file_id
  integer                 :: ii,jj,n_orbits,n_times,ierr
  integer,dimension(:),allocatable   :: n_active_orbits
  integer,dimension(:,:),allocatable :: valid_orbit_id
  real*8,dimension(:,:),allocatable  :: x_loc,ppar_loc,pperp_loc
  !> open the soft orbit hdf5 file
  call HDF5_open(trim(soft_orbit_filename//".h5"),file_id,ierr)
  !> read sofit particles
  call HDF5_allocatable_array2D_reading(file_id,ppar_loc,"/ppar")   !< parallel momentum
  call HDF5_allocatable_array2D_reading(file_id,pperp_loc,"/pperp") !< perpendicular momentum
  call HDF5_allocatable_array2D_reading(file_id,x_loc,"/x")         !< position in xyz coordinates
  !> close the soft orbit hdf5 file
  call HDF5_close(file_id)
  !> remove zero orbits
  n_orbits = size(ppar_loc,2); n_times  = size(ppar_loc,1); 
  !> find id of active orbits
  allocate(n_active_orbits(n_orbits));        n_active_orbits = 0;
  allocate(valid_orbit_id(n_times,n_orbits)); valid_orbit_id  = 0;
  do jj=1,n_orbits
    do ii=1,n_times
      if((x_loc((ii-1)*n_vec+1,jj).eq.0d0).and.(x_loc((ii-1)*n_vec+2,jj).eq.0d0).and.&
      (x_loc(ii*n_vec,jj).eq.0d0)) cycle
        n_active_orbits(jj) = n_active_orbits(jj) + 1
        valid_orbit_id(n_active_orbits(jj),jj) = ii;
    enddo
  enddo
  !> copy only active orbits
  n_soft_points = sum(n_active_orbits)
  allocate(x(n_vec,n_soft_points)); x     = 0d0;
  allocate(ppar(n_soft_points));    ppar  = 0d0; 
  allocate(pperp(n_soft_points));   pperp = 0d0;
  n_soft_points  = 0
  do jj=1,n_orbits   
    do ii = 1,n_active_orbits(jj)
      x(:,n_soft_points+ii)    = x_loc((valid_orbit_id(ii,jj)-1)*n_vec+1:valid_orbit_id(ii,jj)*n_vec,jj)
      ppar(n_soft_points+ii)   = ppar_loc(valid_orbit_id(ii,jj),jj)
      pperp(n_soft_points+ii)  = pperp_loc(valid_orbit_id(ii,jj),jj)
    enddo
    n_soft_points = n_soft_points + n_active_orbits(jj)
  enddo
  !> clean-up
  if(allocated(n_active_orbits)) deallocate(n_active_orbits)
  if(allocated(valid_orbit_id))  deallocate(valid_orbit_id)
  if(allocated(ppar_loc))        deallocate(ppar_loc)
  if(allocated(pperp_loc))       deallocate(pperp_loc)
  if(allocated(x_loc))           deallocate(x_loc)
end subroutine read_soft_orbit_file

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
 
