!> generate_particle_restart_from_SOFT generates a 
!> JOREK particle restart file from a orbits file
!> produced by the SOFT code
program generate_particle_restart_from_SOFT
use constants,       only: SPEED_OF_LIGHT,ATOMIC_MASS_UNIT,EL_CHG
use mod_mpi_tools,   only: init_mpi_threads,finalize_mpi_threads
use mod_particle_io, only: write_simulation_hdf5
use particle_tracer

implicit none

!> Variables -------------------------------------------------------------------------------
type(event)          :: field_reader
integer              :: my_id,n_cpus,ierr,n_vec,n_groups
integer,dimension(4) :: dims
real*8               :: time,mass,charge
real*8,dimension(:,:),allocatable :: soft_orbit_ppar,soft_orbit_pperp,soft_orbit_x
real*8,dimension(:,:),allocatable :: soft_orbit_ppar_local,soft_orbit_pperp_local,soft_orbit_x_local
character(len=17)    :: fields_filename
character(len=125)   :: particle_filename
character(len=193)   :: soft_orbit_filename
!> Variable definitions --------------------------------------------------------------------
n_groups            = 1                   !< number of jorek particle groups
n_vec               = 3                   !< component of a vector
time                = 0d0                 !< simulation time
mass                = 5.48579909065d-4    !< electron mass in AMU 
charge              = -1d0                !< electron charge
fields_filename     = 'jorek_equilibrium' !< jorek restart filename
particle_filename   = 'part_restart_soft_orbits.h5' !< particle restart filename 
soft_orbit_filename = 'orbit_test_jorek_JET_pulse95135_t48dot54_parabolic_qprofile_q95_6dot8_press0_res1r5dot88en1m5_res2r4dot705en1m4_Ip612en1MA_Ekin20MeV_np10_theta1_2dot85_itheta2_pi_nitheta100_norbits100_a96_wall' !< soft orbit filename
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
  call read_soft_orbit_file(soft_orbit_filename,soft_orbit_x,soft_orbit_ppar,soft_orbit_pperp)
  !> scatter the global arrays to each mpi process)
  dims(1) = size(soft_orbit_x,1); dims(2) = size(soft_orbit_ppar,1);
  dims(3) = size(soft_orbit_x,2); dims(4) = dims(3)/n_cpus;
endif
call MPI_Bcast(dims,size(dims),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
allocate(soft_orbit_x_local(dims(1),dims(4))); allocate(soft_orbit_ppar_local(dims(2),dims(4))); 
allocate(soft_orbit_pperp_local(dims(2),dims(4)));
call  scatter_2D_arrays(my_id,n_cpus,dims(1),dims(4),dims(3),soft_orbit_x,soft_orbit_x_local)
call  scatter_2D_arrays(my_id,n_cpus,dims(2),dims(4),dims(3),soft_orbit_ppar,soft_orbit_ppar_local)
call  scatter_2D_arrays(my_id,n_cpus,dims(2),dims(4),dims(3),soft_orbit_pperp,soft_orbit_pperp_local)
if(allocated(soft_orbit_x))     deallocate(soft_orbit_x)
if(allocated(soft_orbit_ppar))  deallocate(soft_orbit_ppar)
if(allocated(soft_orbit_pperp)) deallocate(soft_orbit_pperp)
write(*,*) "Reading SOFT orbit file: completed!"
!> Generate JOREK relativistic gc from SOFT orbits -----------------------------------------
write(*,*) 'Converting soft orbit to JOREK relativistic gc ...'
call convert_soft_orbits_in_jorek_relativistic_gcs(sim,time,mass,charge,n_vec,dims(2),&
dims(4),soft_orbit_x_local,soft_orbit_ppar_local,soft_orbit_pperp_local)
write(*,*) 'Converting soft orbit to JOREK relativistic gc: completed!'
write(*,*) 'Writing JOREK relativistic gc in ',particle_filename,' ...'
call write_simulation_hdf5(sim,trim(particle_filename))
write(*,*) 'Writing JOREK relativistic gc in ',particle_filename,' completed!'
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
!>   time:             (real8) simulation time
!>   mass:             (real8) particle mass in AMU
!>   charge:           (real8) particle charge (RE: -1)
!>   n_vec:            (integer) size of the position vector
!>   n_orbits:         (integer) number of soft orbits
!>   n_t_step:         (integer) number of time steps per soft orbit
!>   soft_orbit_x:     (real8)(n_vec*n_t_step,n_orbits) soft positions in xyz
!>   soft_orbit_ppar:  (real8)(n_t_step,n_orbits) soft parallel momentum
!>   soft_orbit_pperp: (real8)(n_t_step,n_orbits) soft perpendicular momentum
!> outpus:
!>   sim: (particle_sim) initialised particle simulation
subroutine convert_soft_orbits_in_jorek_relativistic_gcs(sim,time,mass,charge,&
n_vec,n_orbits,n_t_steps,soft_orbit_x,soft_orbit_ppar,soft_orbit_pperp)
  use constants,                 only: ATOMIC_MASS_UNIT
  use mod_coordinate_transforms, only: cartesian_to_cylindrical
  implicit none
  !> inputs-outputs:
  type(particle_sim),intent(inout) :: sim
  !> inputs:
  integer,intent(in)               :: n_vec,n_orbits,n_t_steps
  real*8,intent(in)                :: time,mass,charge
  real*8,dimension(n_vec*n_t_steps,n_orbits),intent(in) :: soft_orbit_x
  real*8,dimension(n_t_steps,n_orbits),intent(in)       :: soft_orbit_ppar
  real*8,dimension(n_t_steps,n_orbits),intent(in)       :: soft_orbit_pperp
  !> variables:
  integer :: ii,jj,i_elm,ifail
  real*8              :: U,psi
  real*8,dimension(2) :: st,Rbox,Zbox
  real*8,dimension(3) :: RZphi,B_field,E_field
  !> initialise and allocate particle simulation array
  call domain_bounding_box(sim%fields%node_list,sim%fields%element_list,&
  Rbox(1),Rbox(2),Zbox(1),Zbox(2))
  sim%time = time; sim%groups(1)%mass = mass;
  allocate(particle_gc_relativistic::sim%groups(1)%particles(n_orbits*n_t_steps))
  !> loop on the soft orbits
  !$omp parallel do default(shared) firstprivate(n_orbits,n_t_steps,n_vec,time,charge,&
  !$omp mass,Rbox,Zbox) private(ii,jj,RZphi,i_elm,st,B_field,E_field,U,psi,ifail) &
  !$omp collapse(2)
  do ii=1,n_orbits
    !> loop on the soft time steps
    do jj=1,n_t_steps
      !> transform the soft orbit coordinates in jorek global/local coordinates
      RZphi = cartesian_to_cylindrical(soft_orbit_x((jj-1)*n_vec+1:jj*n_vec,ii)); i_elm = -1;
      if(((RZphi(1).ge.Rbox(1)).and.(RZphi(1).le.Rbox(2))).and.((RZphi(2).ge.Zbox(1)).and.(RZphi(2).le.Zbox(2)))) &
      call find_RZ(sim%fields%node_list,sim%fields%element_list,RZphi(1),RZphi(2),&
      RZphi(1),RZphi(2),i_elm,st(1),st(2),ifail)
      select type (p=>sim%groups(1)%particles((ii-1)*n_t_steps+jj))
        type is (particle_gc_relativistic)
        if(i_elm.le.0) then
          p%x = 0d0; p%st = 0d0; p%weight = 0d0; p%i_elm = 0; p%i_life = 0; p%t_birth = 0.;
          p%p = 0d0; p%q = int(0,kind=1);
        else
          p%x = RZphi; p%st = st; p%weight = 1d0; p%i_elm = i_elm; p%i_life = 0; p%t_birth = 0.;
          call sim%fields%calc_EBpsiU(time,i_elm,st,RZphi(3),E_field,B_field,psi,U)
          p%p = ATOMIC_MASS_UNIT*[soft_orbit_ppar(jj,ii),&
          (5d-1*ATOMIC_MASS_UNIT*(soft_orbit_pperp(jj,ii)**2))/(mass*norm2(B_field))]; 
          p%q = int(charge,kind=1);
        endif
      end select
    enddo
  enddo
  !$omp end parallel do
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
!> outputs:
!>   x:       (real8)(,3*n_orbits) soft positions in xyz coordinates
!>   ppar:    (real8)(,n_orbits) soft parallel momentum
!>   pperp:   (real8)(,n_orbits) soft perpendicular momentum
subroutine read_soft_orbit_file(soft_orbit_filename,x,ppar,pperp)
  use constants, only: ATOMIC_MASS_UNIT
  use hdf5
  use hdf5_io_module, only: HDF5_open,HDF5_close
  use hdf5_io_module, only: HDF5_allocatable_array2D_reading
  implicit none
  !> inputs:
  character(len=*),intent(in) :: soft_orbit_filename
  !> outputs:
  real*8,dimension(:,:),allocatable,intent(out) :: x,ppar,pperp
  !> variables:
  integer(HID_T)          :: file_id
  integer                 :: ierr
  !> open the soft orbit hdf5 file
  call HDF5_open(trim(soft_orbit_filename//".h5"),file_id,ierr)
  !> read sofit particles
  call HDF5_allocatable_array2D_reading(file_id,ppar,"/ppar")   !< parallel momentum
  call HDF5_allocatable_array2D_reading(file_id,pperp,"/pperp") !< perpendicular momentum
  call HDF5_allocatable_array2D_reading(file_id,x,"/x")         !< position in xyz coordinates
  !> close the soft orbit hdf5 file
  call HDF5_close(file_id)
end subroutine read_soft_orbit_file

!> -----------------------------------------------------------------------------------------

end program generate_particle_restart_from_SOFT
 
