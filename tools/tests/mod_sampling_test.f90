!> mod_sampling_test contains all procedures used for
!> for testing the procedures contained in mod_samples
module mod_sampling_test
use fruit
use constants, only: PI,TWOPI
implicit none

private
public :: run_fruit_sampling

!> Variables -----------------------------------------
!> variables for testing the cone sampling methods
integer,parameter :: n_x=3
integer,parameter :: n_cos_half_angle=6
integer,parameter :: n_directions=4
integer,parameter :: n_origins=5
integer,parameter :: n_rays=21
real*8,parameter  :: tol_real8=5.d-15
real*8,dimension(2),parameter :: colat_int=(/-PI,PI/)
real*8,dimension(2),parameter :: azimuth_int=(/0.d0,TWOPI/)
real*8,dimension(3),parameter :: origin_min=(/-1.d0,2.d0,-7.d0/)
real*8,dimension(3),parameter :: origin_max=(/2.d0,5.d0,2.5d0/)
real*8,dimension(2),parameter :: half_angle_int=(/TWOPI/2.3d2,TWOPI/5.3d0/)
real*8,dimension(2),parameter :: length_int=(/1.d-1,3.5d1/)
real*8,dimension(n_cos_half_angle)        :: cos_half_angle_sol
real*8,dimension(n_x,n_directions)        :: directions
real*8,dimension(n_x,n_origins)           :: origins
!> Interfaces ----------------------------------------

contains
!> Fruit basket --------------------------------------
!> fruit basket performing all set-ups, tests and
!> clean-up procedured
subroutine run_fruit_sampling()
  implicit none
  write(*,*) "  ... setting-up: sampling tests"
  call setup()
  write(*,*) "  ... running: sampling tests"
  call test_sample_uniform_standard_cone
  call test_sample_uniform_direction_cone
  call test_sample_uniform_direction_length_cone
  write(*,*) "  ... tearing-down: sampling tests"
end subroutine run_fruit_sampling

!> Set-up and tear-down ------------------------------
!> set up tests parameters
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  integer :: ii
  real*8,dimension(n_directions) :: colatitude,azimuth,lengths
  !> generate set of half solid angle cosines
  call gnu_rng_interval(n_cos_half_angle,half_angle_int,cos_half_angle_sol)
  cos_half_angle_sol = cos(cos_half_angle_sol)
  !> generate set of origins
  do ii=1,n_origins
    call gnu_rng_interval(n_x,origin_min,origin_max,origins(:,ii))
  enddo
  !> generate set of direction vectors
  call gnu_rng_interval(n_directions,length_int,lengths)
  call gnu_rng_interval(n_directions,colat_int,colatitude)
  call gnu_rng_interval(n_directions,azimuth_int,azimuth)
  directions(1,:) = sin(colatitude)*cos(azimuth)
  directions(2,:) = sin(colatitude)*sin(azimuth)
  directions(3,:) = cos(colatitude)
  do ii=1,n_directions
    directions(:,ii) = lengths(ii)*directions(:,ii)
  enddo
end subroutine setup

!> Tests ---------------------------------------------
!> test the generation of random rays within a cone given
!> an origin and a direction vector and a length
subroutine test_sample_uniform_direction_length_cone()
  use mod_sampling, only: sample_uniform_cone
  implicit none
  integer :: ii,jj,kk,pp
  real*8,dimension(3) :: u
  real*8,dimension(3) :: ray
  real*8,dimension(n_rays) :: cos_half_angle,ray_length
  !> computes rays with different origins, directions
  !> and half angles
  do pp=1,n_origins
    do kk=1,n_directions
      do jj=1,n_cos_half_angle
        do ii=1,n_rays
          call random_number(u)
          ray = sample_uniform_cone(cos_half_angle_sol(jj),u,&
          directions(:,kk),origins(:,pp),length_int)
          ray = ray - origins(:,pp)
          ray_length(ii) = norm2(ray)
          cos_half_angle(ii) = dot_product(ray/ray_length(ii),directions(:,kk)/norm2(directions(:,kk)))
        enddo
        call assert_true(all((cos_half_angle.ge.cos_half_angle_sol(jj)).and.(cos_half_angle.le.1.d0)),&
        "Error uniform sampling direction length cone: rays cosinus half angle out-of-bound!")
        call assert_true(all((ray_length.ge.length_int(1)).and.(ray_length.le.length_int(2))),&
        "Error uniform sampling direction length cone: rays length out-of-bound!")
      enddo
    enddo
  enddo
end subroutine test_sample_uniform_direction_length_cone

!> test the generation of random rays within a cone given
!> an origin and a direction vector
subroutine test_sample_uniform_direction_cone()
  use mod_sampling, only: sample_uniform_cone
  implicit none
  integer :: ii,jj,kk,pp
  real*8,dimension(n_rays) :: ones=1.d0
  real*8,dimension(2) :: u
  real*8,dimension(3) :: ray
  real*8,dimension(n_rays) :: cos_half_angle,ray_length
  !> computes rays with different origins, directions
  !> and half angles
  do pp=1,n_origins
    do kk=1,n_directions
      do jj=1,n_cos_half_angle
        do ii=1,n_rays
          call random_number(u)
          ray = sample_uniform_cone(cos_half_angle_sol(jj),u,&
          directions(:,kk),origins(:,pp))
          ray = ray - origins(:,pp)
          ray_length(ii) = norm2(ray)
          cos_half_angle(ii) = dot_product(ray,directions(:,kk)/norm2(directions(:,kk)))
        enddo
        call assert_true(all((cos_half_angle.ge.cos_half_angle_sol(jj)).and.(cos_half_angle.le.1.d0)),&
        "Error uniform sampling direction cone: rays cosinus half angle out-of-bound!")
        call assert_equals(ray_length,ones,n_rays,tol_real8,&
        "Error uniform sampling direction cone: rays length not unitary!")
      enddo
    enddo
  enddo
end subroutine test_sample_uniform_direction_cone

!> test generation of random rays in the stantard cone
subroutine test_sample_uniform_standard_cone()
  use mod_sampling, only: sample_uniform_cone
  implicit none
  real*8,dimension(n_rays) :: ones=1.d0
  integer :: ii,jj
  real*8,dimension(3),parameter :: z_dir=(/0.d0,0.d0,1.d0/)
  real*8,dimension(2) :: u
  real*8,dimension(3) :: ray
  real*8,dimension(n_rays) :: cos_half_angle,ray_length
  !> compute rays and half angles
  do jj=1,n_cos_half_angle
    do ii=1,n_rays
      call random_number(u)
      ray = sample_uniform_cone(cos_half_angle_sol(jj),u)
      ray_length(ii) = norm2(ray)
      cos_half_angle(ii) = dot_product(ray,z_dir)
    enddo
    !> check correctness
    call assert_true(all((cos_half_angle.ge.cos_half_angle_sol(jj)).and.(cos_half_angle.le.1.d0)),&
    "Error uniform sampling standard cone: rays cosinus half angle out-of-bound!")
    call assert_equals(ray_length,ones,n_rays,tol_real8,&
    "Error uniform sampling standard cone: rays length not unitary!")
  enddo
end subroutine test_sample_uniform_standard_cone
!> Tools ---------------------------------------------
!>----------------------------------------------------
end module mod_sampling_test
