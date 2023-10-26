!> This module contains some testcases for projections, ensuring
!> that the projection matrix, RHS and MUMPS work for these cases.
!>
!> It contains tests of projecting zero, x, xy, x^4 onto square,
!> flux-alignes or circular grid
module mod_projection_spec_mpi_test
use fruit
use fruit_mpi
implicit none
private
public :: run_fruit_projection_spec_mpi
!> Variables --------------------------------------
!> Set to true to wrtie restart files with the projected density
logical,parameter :: write_proj_output=.true.
#ifdef PROJECTION_EXTRATEST
  logical,parameter :: EXTRATEST=.true.  !< Set to .true to do flux-alignes projectio test
#else
  logical,parameter :: EXTRATEST=.false. !< Set to .true to do flux-alignes projectio test
#endif
integer,parameter :: message_len=100
integer,parameter :: filename_len=100
integer,parameter :: n_fields_1=1
integer,parameter :: nx=10 !< x-dimension size of square mesh
integer,parameter :: ny=10 !< y-dimension size of square mesh
real*8,parameter  :: time_sol=0.d0;
real*8,dimension(1),parameter :: mean_sol=(/0.d0/)
real*8,dimension(1),parameter :: rms_sol=(/0.d0/)
integer           :: rank_loc,n_tasks_loc,ifail_loc

contains
!> Fruit basket -----------------------------------
subroutine run_fruit_projection_spec_mpi(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  write(*,'(/A)') "  ... setting-up: projection spec"
  call setup(rank,n_tasks,ifail)
  write(*,'(/A)') "  ... running: projection spec"
  call run_test_case(test_project_0_square,'test_project_0_square')
  write(*,'(/A)') "  ... tearing-down: projection spec"
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_projection_spec_mpi

!> Set-ups tear-downs -----------------------------
subroutine setup(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc=rank; n_tasks_loc=n_tasks; ifail_loc=ifail;
end subroutine setup

subroutine teardown(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc=-1; n_tasks_loc=0; ifail=ifail_loc;
end subroutine teardown

!> Tests ------------------------------------------
!> Project zero onto squre grid
subroutine test_project_0_square()
  use nodes_elements
  use mod_projection_helpers_test_tools, only: f_0,default_square_grid
  implicit none
  character(len=message_len),parameter :: message='Error project 0 square test'
  character(len=filename_len) :: filename
  write(filename,'(A,I0,A,I0,A,I0)') 'rank_',rank_loc,'_0_square_',nx,'_',ny
  call default_square_grid(rank_loc,n_tasks_loc,nx,ny,node_list,element_list,ifail_loc)
  call project_f_with_assert_and_write(node_list,element_list,f_0,mean_sol(1),&
  rms_sol(1),trim(filename),trim(message))
end subroutine test_project_0_square 

!> Tools ------------------------------------------
!> Project a function f onto grid in node_list and
!> element_list and test for mean and RMS value.
!> Optionally, write to file visual inspection.
!> inputs:
!>   node_list:    (type_node_list) JOREK node list
!>   element_list: (type_element_list) JOREK element list
!>   f:            (function,real8) function to project
!>   mean:         (real8) expected mean value
!>   rms:          (real8) expected root mean squared
!>   pfilename:    (character) particle file name
!>   message_root: (character) root of the error message to be printed
!>   mean_tol_in   (real8)(optional) tolerance for the mean check
!>   rms_tol_in:   (real8)(optional) tolerance for the rms check
!>   filter_in:    (real8)(optional) smoothing filter value to be used projection
!>   hyper_filter: (real8)(optional) smoothing hyper filter value to be use in projection
!> outputs:
!>   node_list:    (type_node_list) JOREK node list
!>   element_list: (type_element_list) JOREK element list
subroutine project_f_with_assert_and_write(node_list,element_list,f,mean,&
RMS,pfilename,message_root,mean_tol_in,rms_tol_in,filter_in,hyper_filter_in)
  use data_structure,                    only: type_node_list
  use data_structure,                    only: type_element_list
  use mod_project_particles,             only: write_particle_distribution_to_h5
  use mod_projection_helpers_test_tools, only: project_f,elements_mean_rms
  implicit none
  type(type_node_list),intent(inout)    :: node_list
  type(type_element_list),intent(inout) :: element_list
  real*8,external                       :: f
  real*8,intent(in)                     :: mean,rms
  character(len=*),intent(in)           :: pfilename,message_root
  real*8,intent(in),optional            :: mean_tol_in,rms_tol_in
  real*8,intent(in),optional            :: filter_in,hyper_filter_in
  real*8 :: mean_test,rms_test,rms_tol,mean_tol,filter,hyper_filter
  character(len=2*message_len) :: message
  !> initialisation
  filter=0.d0;       if(present(filter_in))       filter=filter_in;
  hyper_filter=0.d0; if(present(hyper_filter_in)) hyper_filter=hyper_filter_in;
  mean_tol=1.d-12;   if(present(mean_tol_in))     mean_tol=mean_tol_in;
  rms_tol=1.d-12;    if(present(rms_tol_in))      rms_tol=rms_tol_in;
  !> project the function and compute mean and rms
  call project_f(node_list,element_list,f,filter,hyper_filter)
  call elements_mean_rms(node_list,element_list,f,mean_test,rms_test)
  !> checks
  write(message,'(A,A)') trim(message_root),': unexpected mean value!'
  call assert_equals(mean,mean_test,mean_tol,trim(message))
  write(message,'(A,A)') trim(message_root),': unexpected RMS value!'
  call assert_equals(RMS,rms_test,rms_tol,trim(message))
  !> write projection in file if requried
  if(write_proj_output) then
    call write_particle_distribution_to_h5(node_list,element_list,filename=pfilename//'.h5',&
    n_fields=n_fields_1,time=time_sol)
  endif 
end subroutine project_f_with_assert_and_write
!> ------------------------------------------------
end module mod_projection_spec_mpi_test
