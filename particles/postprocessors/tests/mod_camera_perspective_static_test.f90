!> mod_camera_perspective_static_test contains all 
!> variables and procedures used for unit testing 
!> the camera perspective static model
module mod_camera_perspective_static_test
use fruit
use mod_assert_equals_tools,   only: assert_equals_allocatable_arrays
use constants,                 only: PI,TWOPI
use mod_pinhole_lens,          only: pinhole_lens
use mod_spectra_deterministic, only: spectrum_integrator_2nd
use mod_camera_perspective_static, only: camera_perspective_static
implicit none

private
public :: run_fruit_camera_perspective_static

!> Variable and data types -------------------------
integer,parameter :: n_int_param=3
integer,parameter :: n_real_param=8
integer,parameter :: n_points_per_pixel=11
integer,parameter :: n_planes=11
integer,parameter :: n_properties=1
integer,parameter :: n_st_sol=2
integer,parameter :: n_stq_sol=3
integer,parameter :: n_x_sol=3
integer,parameter :: n_times_sol=1
integer,parameter :: n_plane_vertices=3
integer,parameter :: n_pixels_x=512
integer,parameter :: n_pixels_y=512
integer,parameter :: n_points_on_lens_sol=2345
integer,parameter :: n_lines_per_spectrum=13
integer,parameter :: n_spectra=2
integer,parameter :: n_rays_sol=123
integer,parameter :: test_plane_pupil_id=1
real*8,parameter  :: tol_real8=5.d-15 
real*8,parameter  :: tol_real8_2=5.d-9
real*8,parameter  :: tol_real8_rel=3.d-6
real*8,parameter  :: plane_distance=2.5d1
real*8,parameter  :: accept_threshold=5.d-1
real*8,dimension(2),parameter :: costheta_interval=(/-1.d0,1.d0/)
real*8,dimension(2),parameter :: phi_interval=(/0.d0,TWOPI/)
real*8,dimension(2),parameter :: ray_q_interval=(/3.1d-2,9.5d-1/)
real*8,dimension(2),parameter :: half_angle_lowbnd=(/PI/1.d1,PI/4.d0/)
real*8,dimension(2),parameter :: half_angle_uppbnd=(/PI/1.9d0,PI/3.d0/)
real*8,dimension(3),parameter :: center_pos_lowbnd=(/-2.d-1,4.d1,-7.d0/)
real*8,dimension(3),parameter :: center_pos_uppbnd=(/3.4d1,3.d2,5.d0/)
real*8,dimension(3),parameter :: test_points_lowbnd=(/-7.4d-1,2.9d1,-5.d0/)
real*8,dimension(3),parameter :: test_points_uppbnd=(/1.4d1,4.5d2,9.d0/)
real*8,dimension(2),parameter :: st_false_lowbnd=(/-5.4d1,-9.d0/)
real*8,dimension(2),parameter :: st_false_uppbnd=(/1.4d1,4.5d2/)
real*8,dimension(n_spectra),parameter :: min_wlen=(/3.d0-6,2.5d-7/)
real*8,dimension(n_spectra),parameter :: max_wlen=(/3.5d0-6,4.2d-7/)
logical,dimension(n_rays_sol)          ::accept_ray_sol
integer,dimension(:,:,:,:),allocatable :: pixel_ids
integer,dimension(:,:,:,:),allocatable :: pixel_ids_sol
real*8,dimension(n_st_sol)             :: pixel_size_sol
real*8,dimension(2,n_planes)           :: half_angle_sol
real*8,dimension(n_x_sol,n_planes)     :: image_plane_coords
real*8,dimension(n_x_sol,n_planes)     :: pupil_positions
real*8,dimension(n_x_sol,n_planes)     :: test_points
real*8,dimension(n_x_sol,n_rays_sol)   :: test_ray_vertices
real*8,dimension(n_stq_sol,n_rays_sol) :: test_ray_stq_sol
real*8,dimension(:,:,:),allocatable    :: points_on_lens
real*8,dimension(:,:,:),allocatable    :: pdf_points_on_lens
real*8,dimension(:,:,:,:),allocatable  :: st_point_on_pixels
real*8,dimension(:,:,:,:),allocatable  :: st_point_on_pixels_sol
type(camera_perspective_static) :: camera_sol
type(pinhole_lens)              :: pinhole_sol
type(spectrum_integrator_2nd)   :: spectrum_sol

!> Interfaces --------------------------------------
interface compute_cos_angle_two_vectors
  module procedure compute_cos_angle_two_vectors_origin_points
end interface compute_cos_angle_two_vectors

contains
!> Fruit basket ------------------------------------
!> fruit basket executes all set-up, test and 
!> tear-down procedures
subroutine run_fruit_camera_perspective_static
  implicit none
  write(*,'(/A)') "  ... setting-up: camera perspective static tests"
  call setup
  write(*,'(/A)') "  ... running: camera perspective static tests"
  call test_points_on_lens_pdf_pinhole
  call test_de_allocation_camera_perspective_static
  call test_image_plane_pixel_size_definitions
  call test_computation_pixel_ids_st_plane_point
  call test_init_camera_perspective_static_pinhole
  call test_cosine_view_angle_static
  call test_material_funct_perspective_static
  call test_find_ray_image_plane_intersection 
  write(*,'(/A)') "  ... tearing-up: camera perspective static tests"
  call teardown
end subroutine run_fruit_camera_perspective_static

!> Set-up and tear-down ----------------------------
!> set-up unit test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  integer :: ii
  real*8,dimension(3) :: center
  !> compute the solution pixel size
  pixel_size_sol = (/1.d0,1.d0/)/real((/n_pixels_x,n_pixels_y/),kind=8)
  !> initialise the pinhole lens
  call gnu_rng_interval(3,center_pos_lowbnd,center_pos_uppbnd,center)
  call pinhole_sol%init_pinhole(n_x_sol,center)
  !> initialise camera spectrum
  spectrum_sol = spectrum_integrator_2nd(n_lines_per_spectrum,&
  n_spectra,min_wlen,max_wlen)
  !> initialise image plane variables
  call generate_image_planes_variables()
  !> initialise camera
  camera_sol%n_plane_points = n_plane_vertices
  !> initialise the test points
  do ii=1,n_planes
    call gnu_rng_interval(n_x_sol,test_points_lowbnd,&
    test_points_uppbnd,test_points(:,ii))
  enddo
  !> generate the ray variables for one image plane
  call generate_ray_variables_from_origin_plane(half_angle_sol(:,test_plane_pupil_id),&
  pupil_positions(:,test_plane_pupil_id),image_plane_coords(:,test_plane_pupil_id),center)
end subroutine setup

!> tearing-down unit test features
subroutine teardown()
  implicit none
  pixel_size_sol = 0.d0
  call pinhole_sol%deallocate_lens; call spectrum_sol%deallocate_spectrum;
  if(allocated(points_on_lens))     deallocate(points_on_lens)
  if(allocated(pdf_points_on_lens)) deallocate(pdf_points_on_lens)
  if(allocated(pixel_ids))          deallocate(pixel_ids)
  if(allocated(pixel_ids_sol))      deallocate(pixel_ids_sol)
  if(allocated(st_point_on_pixels)) deallocate(st_point_on_pixels)
  if(allocated(st_point_on_pixels)) deallocate(st_point_on_pixels_sol)
end subroutine teardown

!> Tests -------------------------------------------
!> test the initialisation camera perspective static
subroutine test_init_camera_perspective_static_pinhole()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  use mod_geometry, only: define_plane_from_half_angles
  implicit none
  !> variables
  integer,dimension(n_int_param)           :: int_param
  real*8,dimension(n_x_sol)                :: direction_sol,direction_test
  real*8,dimension(n_real_param)           :: real_param
  real*8,dimension(n_x_sol,n_plane_vertices) :: image_plane_sol
  !> initialisation, only one image plane is considered
  int_param = (/n_points_on_lens_sol,n_pixels_x,n_pixels_y/)
  real_param(1:2) = half_angle_sol(:,test_plane_pupil_id)
  real_param(3:5) = image_plane_coords(:,test_plane_pupil_id)
  real_param(6:8) = pupil_positions(:,test_plane_pupil_id)
  allocate(points_on_lens(n_x_sol,1,1))
  allocate(pdf_points_on_lens(1,1,1))
  call pinhole_sol%sampling(1,points_on_lens)
  call pinhole_sol%pdf(1,points_on_lens(:,:,1),pdf_points_on_lens(:,1,1))
  !> initialise the camera perspective static
  call camera_sol%init_camera(pinhole_sol,spectrum_sol,&
  n_int_param,n_real_param,int_param,real_param)
  !> perform tests
  call assert_equals(camera_sol%n_vertices,1,&
  "Error init camera perspective static pinhole: n vertices is not 1")
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%n_active_vertices,&
  "Error init camera perspective static pinhole: n active vertices")
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%times,&
  "Error init camera perspective static pinhole: times")
  call assert_equals_allocatable_arrays(n_x_sol,1,&
  n_times_sol,camera_sol%x,"Error init camera perspective static pinhole: x")
  call assert_equals_allocatable_arrays(n_properties,1,n_times_sol,&
  camera_sol%properties,"Error init camera perspective static pinhole: properties")
  call assert_equals(camera_sol%n_pixels_spectra,(/n_spectra,n_pixels_x,n_pixels_y/),&
  3,"Error init camera perspective static pinhole: n pixels / spectra")
  call assert_equals_allocatable_arrays(n_x_sol,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%x,points_on_lens,tol_real8,&
  ":Error init camera perspective static pinhole: points on lens")
  call assert_equals_allocatable_arrays(n_properties,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%properties,pdf_points_on_lens,tol_real8,&
  "Error init camera perspective static pinhole: pdf points on lens") 
  call assert_equals_allocatable_arrays(n_x_sol,camera_sol%image_plane_direction,&
  "Error init camera perspective static pinhole: image plane direction")
  if(allocated(camera_sol%image_plane_direction)) then
    direction_test = (/1.d0,acos(camera_sol%image_plane_direction(3)),&
    atan2(camera_sol%image_plane_direction(2),camera_sol%image_plane_direction(1))/)
    if(direction_test(3).lt.0.d0) direction_test(3) = TWOPI + direction_test(3)
    direction_sol = (/1.d0,image_plane_coords(2,1),image_plane_coords(3,1)/)
    call assert_equals(direction_test,direction_sol,n_x_sol,tol_real8,&
    "Error init camera perspective static pinhole: image plane mismatch!")
  endif
  call assert_equals_allocatable_arrays(n_x_sol,n_plane_vertices,&
  camera_sol%image_plane,"Error init camera perspective static pinhole: image plane")
  if(allocated(camera_sol%image_plane)) then
    call define_plane_from_half_angles(half_angle_sol(:,1),&
    image_plane_coords(:,1),pupil_positions(:,1),image_plane_sol)
    call assert_equals(camera_sol%image_plane,image_plane_sol,n_x_sol,&
    n_plane_vertices,tol_real8,&
    "Error init camera perspective static pinhole: image plane mismatch!")
  endif
  call assert_equals(camera_sol%pixel_size,pixel_size_sol,2,tol_real8,&
  "Error init camera perspective static pinhole: pixel size mismatch!")
  !> clean up
  deallocate(points_on_lens); deallocate(pdf_points_on_lens);
  call camera_sol%deallocate_camera_perspective_static
end subroutine test_init_camera_perspective_static_pinhole

!> test allocation and deallocation of camera_perspective_static
!> attributs (only)
subroutine test_de_allocation_camera_perspective_static()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  implicit none
  !> allocate camera perspective static
  call camera_sol%allocate_camera_perspective_static(spectrum_sol,3,&
  (/n_points_on_lens_sol,n_pixels_x,n_pixels_y/))
  !> test allocation
  call assert_equals(camera_sol%n_vertices,n_points_on_lens_sol,&
  "Error allocate camera perspective static: n vertices mismatch")
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%n_active_vertices,&
  "Error allocate camera perspective static: n active vertices")
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%times,&
  "Error allocate camera perspective static: times")
  call assert_equals_allocatable_arrays(n_x_sol,n_points_on_lens_sol,&
  n_times_sol,camera_sol%x,"Error allocate camera perspective static: x")
  call assert_equals_allocatable_arrays(n_properties,n_points_on_lens_sol,n_times_sol,&
  camera_sol%properties,"Error allocate camera perspective static: properties")
  call assert_equals(camera_sol%n_pixels_spectra,(/n_spectra,n_pixels_x,n_pixels_y/),&
  3,"Error allocate camera perspective static: n pixels / spectra")
  call assert_equals_allocatable_arrays(n_x_sol,camera_sol%image_plane_direction,&
  "Error allocate_camera perspective static: image plane direction")
  call assert_equals_allocatable_arrays(n_x_sol,n_plane_vertices,&
  camera_sol%image_plane,"Error allocate_camera perspective static: image plane") 
  !> deallocate camera perspective static
  call camera_sol%deallocate_camera_perspective_static
  call assert_equals(camera_sol%n_vertices,0,&
  "Error deallocate camera perspective static: n vertices not 0!")
  call assert_equals(camera_sol%n_pixels_spectra,(/0,0,0/),3,&
  "Error deallocate camera perspective static: n pixels spectra not 0!")
  call assert_false(allocated(camera_sol%n_active_vertices),&
  "Error deallocate camera perspective static: n active vertices not deallocated!")
  call assert_false(allocated(camera_sol%times),&
  "Error deallocate camera perspective static: times not deallocated!")
  call assert_false(allocated(camera_sol%x),&
  "Error deallocate camera perspective static: x not deallocated!")
  call assert_false(allocated(camera_sol%properties),&
  "Error deallocate camera perspective static: properties not deallocated!")
  call assert_false(allocated(camera_sol%image_plane_direction),&
  "Error deallocate camera perspective static: image plane direction not deallocated!")
  call assert_false(allocated(camera_sol%image_plane),&
  "Error deallocate camera perspective static: image plane not deallocated!")
end subroutine test_de_allocation_camera_perspective_static

!> test generation of points and pdf on lens using
!> a pinhole lens
subroutine test_points_on_lens_pdf_pinhole()
  implicit none
  !> allocate and initialise variables
  camera_sol%n_property_vertex = n_properties
  call camera_sol%allocate_vertices(1,1)
  allocate(points_on_lens(n_x_sol,1,1))
  allocate(pdf_points_on_lens(1,1,1))
  call pinhole_sol%sampling(1,points_on_lens)
  call pinhole_sol%pdf(1,points_on_lens(:,:,1),pdf_points_on_lens(:,1,1))
  !> test generation without input number of points
  call camera_sol%generate_points_on_lens_pdf(pinhole_sol)
  call assert_equals(camera_sol%n_vertices,1,&
  "Error camera perspective static generate points on pinhole: n vertices not 1!")
  call assert_equals_allocatable_arrays(n_x_sol,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%x,points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: x")
  call assert_equals_allocatable_arrays(n_properties,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%properties,pdf_points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: pdf points on lens")
  !> test generation with input number points
  call camera_sol%generate_points_on_lens_pdf(pinhole_sol,n_points_on_lens_sol)
  call assert_equals(camera_sol%n_vertices,1,&
  "Error camera perspective static generate points on pinhole: n vertices not 1!")
  call assert_equals_allocatable_arrays(n_x_sol,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%x,points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: points on lens")
  call assert_equals_allocatable_arrays(n_properties,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%properties,pdf_points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: pdf points on lens") 
  !> deallocate variable
  call camera_sol%deallocate_vertices
  deallocate(points_on_lens); deallocate(pdf_points_on_lens);
end subroutine test_points_on_lens_pdf_pinhole

!> test the generation of image planes and the calculation of the pixel size
subroutine test_image_plane_pixel_size_definitions()
  use constants,    only: TWOPI
  use mod_geometry, only: define_plane_from_half_angles
  implicit none
  !> variables
  integer                                    :: ii
  real*8,dimension(n_x_sol)                  :: direction_sol,direction_test
  real*8,dimension(n_x_sol,n_plane_vertices) :: image_plane_sol
  real*8,dimension(n_x_sol,n_plane_vertices) :: image_plane_std_sol
  real*8,dimension(n_real_param)             :: real_param
  !> initialisation
  camera_sol%n_property_vertex = n_properties;
  call camera_sol%allocate_camera_perspective_static(spectrum_sol,3,&
  (/n_points_on_lens_sol,n_pixels_x,n_pixels_y/))
  !> test the definition of the image plane and pixel size
  do ii=1,n_planes
    !> store plane value in parameters
    real_param(1:2) = half_angle_sol(:,ii)
    real_param(3:5) = image_plane_coords(:,ii)
    real_param(6:8) = pupil_positions(:,ii)
    !> generate camera and solution planes
    call camera_sol%define_image_plane_pixel_size(n_real_param,real_param)
    call define_plane_from_half_angles(half_angle_sol(:,ii),&
    image_plane_coords(:,ii),pupil_positions(:,ii),image_plane_sol)
    call define_plane_from_half_angles(half_angle_sol(:,ii),&
    image_plane_coords(:,ii),image_plane_std_sol)
    direction_test = (/1.d0,acos(camera_sol%image_plane_direction(3)),&
    atan2(camera_sol%image_plane_direction(2),camera_sol%image_plane_direction(1))/)
    if(direction_test(3).lt.0.d0) direction_test(3) = TWOPI + direction_test(3)
    direction_sol = (/1.d0,real_param(4),real_param(5)/)
    !> test plane and pixel size
    call assert_equals(direction_test,direction_sol,n_x_sol,tol_real8,&
    "Error camera perspective static define image plane: image plane direction mismatch!")
    call assert_equals(camera_sol%image_plane,image_plane_sol,n_x_sol,&
    n_plane_vertices,tol_real8,&
    "Error camera perspective static define image plane: image plane mismatch!")
    call assert_equals(camera_sol%pixel_size,pixel_size_sol,2,tol_real8,&
    "Error camera perspective static define image plane: pixel size mismatch!")
  enddo
  !> deallocate camera perspective static
  call camera_sol%deallocate_camera_perspective_static
end subroutine test_image_plane_pixel_size_definitions

!> test the calculation of the calculation of the position of a point on
!> a plane in the pixel coordinate system
subroutine test_computation_pixel_ids_st_plane_point()
  use mod_assert_equals_tools, only: assert_equals_extended
  use mod_assert_equals_tools, only: assert_equals_rel_error
  implicit none
  !> variables
  integer :: ii,jj,kk
  real*8,dimension(n_st_sol) :: st_plane
  !> initialisations
  camera_sol%pixel_size = pixel_size_sol
  allocate(st_point_on_pixels_sol(n_st_sol,n_points_per_pixel,n_pixels_x,n_pixels_y))
  allocate(st_point_on_pixels(n_st_sol,n_points_per_pixel,n_pixels_x,n_pixels_y))
  allocate(pixel_ids(n_st_sol,n_points_per_pixel,n_pixels_x,n_pixels_y))
  allocate(pixel_ids_sol(n_st_sol,n_points_per_pixel,n_pixels_x,n_pixels_y))
  !> compute local and global coordinates of the points
  do ii=1,n_pixels_y
    do jj=1,n_pixels_x
      do kk=1,n_points_per_pixel
        pixel_ids_sol(:,kk,jj,ii) = (/jj,ii/)
        call random_number(st_point_on_pixels_sol(:,kk,jj,ii))
        st_plane = (st_point_on_pixels_sol(:,kk,jj,ii) + &
        real(pixel_ids_sol(:,kk,jj,ii)-1,kind=8))*pixel_size_sol
        call camera_sol%plane_to_pixel_local_coord(st_plane,&
        pixel_ids(:,kk,jj,ii),st_point_on_pixels(:,kk,jj,ii))
      enddo
    enddo
  enddo
  !> test solutions
  call assert_equals_extended(n_st_sol,n_points_per_pixel,n_pixels_x,&
  n_pixels_y,pixel_ids,pixel_ids_sol,&
  "Error computation position in pixel coordinates: pixel ids mismatch!")
  call assert_equals_rel_error(n_st_sol,n_points_per_pixel,n_pixels_x,&
  n_pixels_y,st_point_on_pixels,st_point_on_pixels_sol,tol_real8_rel,&
  "Error computation position in pixel coordinates: pixel coords. mismatch!")
  !> cleanup
  camera_sol%pixel_size = 0.d0
  deallocate(st_point_on_pixels); deallocate(pixel_ids_sol);
  deallocate(st_point_on_pixels_sol); deallocate(pixel_ids);
end subroutine test_computation_pixel_ids_st_plane_point

!> test the calculation of the cosinus between the image plane direction and a ray
subroutine test_cosine_view_angle_static()
  use mod_geometry,     only: define_vertex_spherical_coord
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  !> variables
  type(pinhole_lens)             :: pinhole
  integer :: ii
  integer,dimension(n_int_param) :: int_param
  real*8,dimension(n_x_sol)      :: vertex_1
  real*8,dimension(n_real_param) :: real_param
  real*8,dimension(n_planes)     :: cos_view_angle_sol
  real*8,dimension(n_planes)     :: cos_view_angle_test
  !> initialisation
  camera_sol%n_property_vertex = n_properties;
  int_param = (/n_points_on_lens_sol,n_pixels_x,n_pixels_y/)
  call camera_sol%init_camera(pinhole_sol,spectrum_sol,&
  n_int_param,n_real_param,int_param,real_param)
  !> test the calculation of the view angle cosinus
  do ii=1,n_planes
    !> store plane value in parameters
    real_param(1:2) = half_angle_sol(:,ii)
    real_param(3:5) = image_plane_coords(:,ii)
    real_param(6:8) = pupil_positions(:,ii)
    !> compute solution
    call define_vertex_spherical_coord(real_param(3:5),real_param(6:8),vertex_1)
    call compute_cos_angle_two_vectors(real_param(6:8),vertex_1,test_points(:,ii),cos_view_angle_sol(ii))
    !> compute test value
    call pinhole%init_pinhole(n_x_sol,real_param(6:8))
    call camera_sol%generate_points_on_lens_pdf(pinhole)
    call camera_sol%define_image_plane_pixel_size(n_real_param,real_param)
    call camera_sol%cos_view_angle_static(test_points(:,ii),1,cos_view_angle_test(ii))
  enddo
  !> check solutions
  call assert_equals(cos_view_angle_test,cos_view_angle_sol,n_planes,tol_real8,&
  "Error computation cosine view angle perspective static: cosine view angles mismatch!")
  !> deallocate camera perspective static
  call pinhole%deallocate_lens
  call camera_sol%deallocate_camera_perspective_static
end subroutine test_cosine_view_angle_static 

!> test the calculation of the physical material function for the perspective static camera
subroutine test_material_funct_perspective_static()
  use mod_geometry,     only: define_vertex_spherical_coord
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  !> variables
  type(pinhole_lens)             :: pinhole
  integer :: ii
  integer,dimension(n_int_param) :: int_param
  real*8,dimension(n_x_sol)      :: vertex_1
  real*8,dimension(n_real_param) :: real_param
  real*8,dimension(n_planes)     :: material_sol
  real*8,dimension(n_planes)     :: material_test
  !> initialisation
  camera_sol%n_property_vertex = n_properties;
  int_param = (/n_points_on_lens_sol,n_pixels_x,n_pixels_y/)
  call camera_sol%init_camera(pinhole_sol,spectrum_sol,&
  n_int_param,n_real_param,int_param,real_param)
  !> test the calculation of the view angle cosinus
  do ii=1,n_planes
    !> store plane value in parameters
    real_param(1:2) = half_angle_sol(:,ii)
    real_param(3:5) = image_plane_coords(:,ii)
    real_param(6:8) = pupil_positions(:,ii)
    !> compute solution
    call define_vertex_spherical_coord(real_param(3:5),real_param(6:8),vertex_1)
    call compute_cos_angle_two_vectors(real_param(6:8),vertex_1,test_points(:,ii),material_sol(ii))
    !> compute test value
    call pinhole%init_pinhole(n_x_sol,real_param(6:8))
    call camera_sol%generate_points_on_lens_pdf(pinhole)
    call camera_sol%define_image_plane_pixel_size(n_real_param,real_param)
    call camera_sol%physical_material_funct(test_points(:,ii),1,material_test(ii),1)
  enddo
  !> check solutions
  call assert_equals(material_test,material_sol,n_planes,tol_real8,&
  "Error computation physical material funct perspective static: cosine view angles mismatch!")
  !> deallocate camera perspective static
  call pinhole%deallocate_lens
  call camera_sol%deallocate_camera_perspective_static
end subroutine test_material_funct_perspective_static

!> test the method for finding image plane - ray intersections
!> for simplicity, only one image plane is tested
subroutine test_find_ray_image_plane_intersection()
  implicit none
  !> variables
  integer :: ii
  integer,dimension(n_int_param)         :: int_param
  logical,dimension(n_rays_sol)          :: test_intersection
  real*8,dimension(n_real_param)         :: real_param
  real*8,dimension(n_stq_sol,n_rays_sol) :: test_ray_stq
  !> initialisation
  camera_sol%n_property_vertex = n_properties;
  int_param = (/n_points_on_lens_sol,n_pixels_x,n_pixels_y/)
  real_param(1:2) = half_angle_sol(:,test_plane_pupil_id)
  real_param(3:5) = image_plane_coords(:,test_plane_pupil_id)
  real_param(6:8) = pupil_positions(:,test_plane_pupil_id)
  call camera_sol%init_camera(pinhole_sol,spectrum_sol,&
  n_int_param,n_real_param,int_param,real_param)
  !> compute rays and intersections
  do ii=1,n_rays_sol
    call camera_sol%find_ray_image_plane_intersection(&
    test_ray_vertices(:,ii),test_plane_pupil_id,&
    test_intersection(ii),test_ray_stq(:,ii))
  enddo
  !> check solutions
  call assert_equals(test_intersection,accept_ray_sol,n_rays_sol,&
  "Error find ray image plane intersection: intersections mismatch!")
  call assert_equals(test_ray_stq,test_ray_stq_sol,n_stq_sol,n_rays_sol,&
  tol_real8_2,"Error find ray image plane intersection: local coordinates mismatch!")
  !> cleanup
  call camera_sol%deallocate_camera_perspective_static
end subroutine test_find_ray_image_plane_intersection 

!> Tools -------------------------------------------
!> sample the image plane variables
subroutine generate_image_planes_variables()
  use mod_gnu_rng,  only: gnu_rng_interval
  use mod_sampling, only: sample_uniform_sphere
  implicit none
  !> variables
  integer             :: ii
  real*8,dimension(3) :: rand
  !> generation routine
  do ii=1,n_planes
    call random_number(rand)
    image_plane_coords(:,ii) = sample_uniform_sphere(&
    plane_distance,costheta_interval,phi_interval,rand)
    call gnu_rng_interval(2,half_angle_lowbnd,&
    half_angle_uppbnd,half_angle_sol(:,ii))
    call gnu_rng_interval(n_x_sol,center_pos_lowbnd,&
    center_pos_uppbnd,pupil_positions(:,ii))
  enddo
end subroutine generate_image_planes_variables

!> sample the vertex of a ray given an origin and a plane
subroutine generate_ray_variables_from_origin_plane(half_width,&
origin,plane_coords,x_lens)
  use mod_geometry, only: define_plane_from_half_angles
  use mod_geometry, only: compute_global_cart_coord_plane_points
  use mod_gnu_rng,  only: gnu_rng_interval
  implicit none
  !> inputs
  real*8,dimension(2),intent(in)         :: half_width
  real*8,dimension(n_x_sol),intent(in)   :: origin,x_lens
  real*8,dimension(n_x_sol),intent(in)   :: plane_coords
  !> variables
  integer :: ii
  real*8 ::  rand
  real*8,dimension(n_st_sol)  :: st_value
  real*8,dimension(n_x_sol)   :: plane_pos
  real*8,dimension(n_x_sol,3) :: plane
  !> intiialisation
  call define_plane_from_half_angles(half_width,plane_coords,origin,plane)
  test_ray_stq_sol(1:2,:) = 5.d-1
  !> loop on the number of rays
  do ii=1,n_rays_sol
    !> check if a ray with or without intersection must be generated
    call random_number(rand)
    if(rand.gt.accept_threshold) then
      call random_number(test_ray_stq_sol(1:2,ii))
      accept_ray_sol(ii) = .true.
    else
      do while((all(test_ray_stq_sol(1:2,ii).ge.0.d0).and.all(test_ray_stq_sol(1:2,ii).le.1.d0)))
        call gnu_rng_interval(2,st_false_lowbnd,st_false_uppbnd,test_ray_stq_sol(1:2,ii))
      enddo
      accept_ray_sol(ii) = .false.
    endif
    !> sample a ray lenght
    call gnu_rng_interval(ray_q_interval,test_ray_stq_sol(3,ii))
    !> compute the position on the plane
    plane_pos = compute_global_cart_coord_plane_points(plane,test_ray_stq_sol(1:2,ii))
    !> compute and store the ray vertex and its local coordinate
    test_ray_vertices(:,ii) = x_lens+(plane_pos-x_lens)/test_ray_stq_sol(3,ii)
  enddo
end subroutine generate_ray_variables_from_origin_plane

!> compute angle between two vectors with same origin
subroutine compute_cos_angle_two_vectors_origin_points(&
origin,vertex_1,vertex_2,cos_angle)
  implicit none
  !> inputs:
  real*8,dimension(n_x_sol) :: origin,vertex_1,vertex_2
  !> outputs:
  real*8 :: cos_angle
  !> comput angle
  cos_angle = dot_product(vertex_2-origin,vertex_1-origin)
  cos_angle = cos_angle/(norm2(vertex_2-origin)*norm2(vertex_1-origin))
end subroutine compute_cos_angle_two_vectors_origin_points

!>--------------------------------------------------
end module mod_camera_perspective_static_test

