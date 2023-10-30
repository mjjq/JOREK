module mod_particle_projection_spec_mpi_test
use data_structure
use fruit
implicit none
private
public :: run_fruit_particle_projection_spec_mpi
!> Variables --------------------------------------
!> Set to true to wrtie restart files with projected density
logical,parameter :: write_projection_output=.false.
logical,parameter :: impose_dirichlet=.false.
integer,parameter :: message_len=100
integer,parameter :: filename_len=100
integer,parameter :: master_rank=0
real*8,parameter  :: R_particle_in=2.d0
real*8,parameter  :: Z_particle_in=1.d0
real*8,parameter  :: test_time=0.d0
integer,dimension(4),parameter :: n_particles=(/1000,10000,100000,1000000/)
integer,dimension(1),parameter :: nx=(/10/)
integer,dimension(1),parameter :: ny=(/10/)
type(type_node_list),pointer    :: test_nodes
type(type_element_list),pointer :: test_elements
integer :: rank_loc,n_tasks_loc,ifail_loc
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_particle_projection_spec_mpi(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  write(*,'(/A)') "  ... setting-up: particle projection spec mpi"
  call setup(rank,n_tasks,ifail)
  write(*,'(/A)') "  ... running: particle projection spec mpi"
  call run_test_case(test_particle_projection_square_10_10_pcg,&
  'test_particle_projection_square_10_10_pcg')
  write(*,'(/A)') "  ... tearing-down: particle projection spec mpi"
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_particle_projection_spec_mpi

!> Set-up and tear-down ---------------------------
subroutine setup(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc=rank; n_tasks_loc=n_tasks; ifail_loc=ifail;
  allocate(test_nodes); allocate(test_elements);
end subroutine setup

subroutine teardown(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc = -1; n_tasks_loc = 0; ifail=ifail_loc;
  deallocate(test_nodes); deallocate(test_elements);
end subroutine teardown

!> Tests ------------------------------------------
!> Project 10^3-10^5 particles generated with pcg ont square grid
subroutine test_particle_projection_square_10_10_pcg
  use constants,                         only: TWOPI
  use mod_pcg32_rng,                     only: pcg32_rng
  use mod_project_particles,             only: proj_one
  use mod_projection_helpers_test_tools, only: f_1,default_square_grid
  implicit none
  real*8,parameter              :: expect_mean=1.d0
  real*8,parameter              :: expect_rms=0.d0
  real*8,parameter              :: volume=TWOPI
  real*8,dimension(3),parameter :: tol_mean=[3.d-8,3.d-8,3.d-8]
  real*8,dimension(3),parameter :: tol_rms=[2.3d1/real(n_particles(1),kind=8),&
                                   2.3d1/real(n_particles(2),kind=8),&
                                   2.3d1/real(n_particles(3),kind=8)]
  character(len=message_len)    :: message
  character(len=filename_len)   :: filename
  write(message,'(A,I0,A,I0,A,I0,A)') 'Error particle projection square nx: ',&
  nx(1),' ny: ',ny(1),' pcg rank: ',rank_loc,':'
  write(filename,'(A,I0,A,I0,A,I0,A)') '_test_projection_square_rank',&
  rank_loc,'_nx',nx(1),'_ny',ny(1),'_pcg'
  call default_square_grid(rank_loc,n_tasks_loc,nx(1),nx(1),&
  test_nodes,test_elements,ifail_loc)
  call project_n(rank_loc,master_rank,test_nodes,test_elements,proj_one,f_1,&
  n_particles(1:3),pcg32_rng(),volume,expect_mean,expect_rms,tol_mean,&
  tol_rms,trim(adjustl(message)),trim(adjustl(filename)),ifail_loc,&
  apply_dirichlet_in=impose_dirichlet,write_particle_in=write_projection_output)
end subroutine test_particle_projection_square_10_10_pcg


!> Tool -------------------------------------------
!> Helper function to project n particles onto a grid in nod_list,
!> element list with optional smoothing. It creates a projection
!> type behind the scenes and uses that. We also need to create
!> a particle-sim here.
subroutine project_n(rank,master,node_list,element_list,proj_f_proj,&
f_proj,n,rng,volume,mean_expect,rms_expect,mean_tol,rms_tol,message,&
fname,ifail,smoothing_in,n_tor_local_in,i_tor_local_in,&
apply_dirichlet_in,write_particle_in,n_fields_write_in)
  use data_structure
  use mod_rng,                           only: type_rng
  use mod_initialise_particles,          only: initialise_particles
  use mod_particle_sim,                  only: particle_sim
  use mod_project_particles,             only: projection,new_projection,proj_f
  use mod_project_particles,             only: write_particle_distribution_to_h5
  use mod_particle_types,                only: particle_kinetic
  use mod_projection_helpers_test_tools, only: elements_mean_rms
  implicit none
  type(type_node_list),intent(inout)    :: node_list
  type(type_element_list),intent(inout) :: element_list
  class(type_rng),intent(in)            :: rng
  integer,intent(inout)                 :: ifail
  integer,intent(in)                    :: rank,master
  integer,dimension(:),intent(in)       :: n
  real*8,intent(in)                     :: volume,mean_expect,rms_expect
  real*8,dimension(:),intent(in)        :: rms_tol,mean_tol
  character(len=*),intent(in)           :: message,fname
  integer,optional,intent(in)           :: n_tor_local_in,i_tor_local_in
  integer,optional,intent(in)           :: n_fields_write_in
  real*8,intent(in),optional            :: smoothing_in
  logical,intent(in),optional           :: apply_dirichlet_in,write_particle_in
  real*8,external                       :: proj_f_proj,f_proj
  type(projection)                      :: project
  type(particle_sim)                    :: sim
  integer                               :: ii,jj,kk,ielm_out
  integer                               :: i_tor_local,n_tor_local,n_fields_write
  real*8                                :: smoothing 
  real*8                                :: R_out,Z_out,s_out,t_out
  real*8                                :: mean,rms_error
  character*8                           :: number_particles,tol_s
  character*8                           :: smooth_string,group_string
  logical                               :: apply_dirichlet,write_particle

  !> initialisations
  smoothing = 0.d0; smooth_string = '';
  if(present(smoothing_in)) then
    smoothing = smoothing_in; write(smooth_string,'(g8.1)') smoothing_in;
  endif
  write_particle =.false.; 
  if(present(write_particle_in)) write_particle = write_particle_in;
  apply_dirichlet = .true.
  if(present(apply_dirichlet_in)) apply_dirichlet = apply_dirichlet_in
  n_fields_write = 1; if(present(n_fields_write_in)) n_fields_write = n_fields_write_in;
  n_tor_local = 1;    if(present(n_tor_local_in)) n_tor_local = n_tor_local_in;
  i_tor_local = 1;    if(present(i_tor_local_in)) i_tor_local = i_tor_local_in;
  project = new_projection(node_list,element_list,filter_n0=smoothing,&
  f=[proj_f(proj_f_proj,group=1)],do_dirichlet=apply_dirichlet)
  project%n_tor_local = n_tor_local; project%i_tor_local = i_tor_local;
  allocate(sim%groups(1));
  !> fill the particle structure
  do kk=1,size(n)
    write(number_particles,'(I8)') n(kk)
    write(group_string,'(I8)') kk
    allocate(particle_kinetic::sim%groups(1)%particles(n(kk)))
    !> to prevent omp trouble (!?)
    call find_RZ(node_list,element_list,R_particle_in,Z_particle_in,&
    R_out,Z_out,ielm_out,s_out,t_out,ifail)
    call initialise_particles(sim%groups(1)%particles,node_list,element_list,rng)
    do ii=1,n(kk)
      sim%groups(1)%particles(ii)%weight = volume/real(n(kk),kind=8)
    enddo
    call project%do(sim) !< project particles, results in node_list
    deallocate(sim%groups(1)%particles) !< cleanup
    !> perform checks on the mean and rms
    call elements_mean_rms(project%node_list,project%element_list,&
    f_proj,mean,rms_error)
    write(tol_s,'(g8.1)') mean_tol(kk)
    call assert_equals(mean_expect,mean,mean_tol(kk),message//&
    ' mean value not matched for n='//trim(adjustl(number_particles))//&
    ', tol: '//trim(tol_s)//' smoothing: '//trim(smooth_string))
    write(tol_s,'(g8.1)') rms_tol(kk)
    call assert_equals(rms_expect,rms_error,rms_tol(kk),message//&
    ' RMS value not matched for n='//trim(adjustl(number_particles))//&
    ', tol: '//trim(tol_s)//' smoothing: '//trim(smooth_string))
    if(write_particle) then
      call write_particle_distribution_to_h5(project%node_list,project%element_list,&
      filename='part_'//trim(adjustl(number_particles))//'_group_'//&
      trim(adjustl(group_string))//trim(fname)//'.h5',n_fields=n_fields_write,time=test_time)
    endif
  enddo
  call project%close_mumps()
end subroutine project_n
!> ------------------------------------------------
end module mod_particle_projection_spec_mpi_test
