!------------------------------------------------------------------
module mod_wall_collision
use hdf5_io_module
use mod_coordinate_transforms
implicit none
contains

function cross(a,b)
  implicit none
  real(kind=8), dimension(3)                     :: cross
  real(kind=8), dimension(3), intent(in)         :: a,b
  
  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
end function cross

function cylindrical_to_cartesian_real8(cyl) result(xyz)
  real*8, intent(in)           :: cyl(3) !< The vector components in RZPhi coordinates
  real*8                       :: xyz(3) !< The vector components in xyz coordinates

  xyz(1) = cyl(1)*cos(-cyl(3))
  xyz(2) = cyl(1)*sin(-cyl(3))
  xyz(3) = cyl(2)
end function cylindrical_to_cartesian_real8

function  plucker(p,q)
  implicit none
  real(kind=8),dimension(3)                      :: p, q
  real(kind=8),dimension(6)                      :: plucker
  p = cylindrical_to_cartesian_real8(p)
  q = cylindrical_to_cartesian_real8(q)
  plucker(1) = p(1) * q(2) - p(2) * q(1)
  plucker(2) = p(1) * q(3) - p(3) * q(1)
  plucker(3) = p(1) - q(1)
  plucker(4) = p(2) * q(3) - p(3) * q(2)
  plucker(5) = p(3) - q(3)
  plucker(6) = p(2) - q(2)
end function plucker

function plane_check(p,q,v0,v1,v2)
  implicit none
  real(kind=8), dimension(3)                     :: p, q, v0, v1,v2
  real(kind=8), dimension(3)                     :: normal
  logical                                        :: plane_check
  real(kind=8)                                   :: dot
  p =  cylindrical_to_cartesian_real8(p)
  q =  cylindrical_to_cartesian_real8(q)
  v0 =  cylindrical_to_cartesian_real8(v0)
  v1 = cylindrical_to_cartesian_real8(v1)
  v2 = cylindrical_to_cartesian_real8(v2)
  normal = cross(v1-v0,v2-v0)
  plane_check = .false.
  dot = dot_product(q-p,normal)
  if(dot == 0.d0) then
    plane_check = .true.
  endif
end function plane_check

!subroutine wall_collision_init()
!  implicit none
!  integer(HID_T)    :: file, file_space, mem_space, dset, plist ! handles
!  integer(HID_T)    :: data_type
!  integer(HID_T)    :: group_id
!  integer(HID_T)    :: time_set_id
!  integer(HSIZE_T)  :: i_here
!  integer           :: n_here
!  integer           :: storage_type, max_corder
!  character(len=12) :: group_name
!  character(len=particle_type_name_length) :: particle_type_name
!  integer           :: i, j, n, hdferr, n_alive
!  integer, allocatable :: n_alive_all(:)
!  logical           :: exists

!  type(c_ptr) :: p_ptr
!  integer*8, dimension(1:2) :: tmp, maxdims
!  real*8, dimension(:,:), allocatable :: real8_2D
!  integer*4, dimension(:), allocatable :: int4_1D
!  real*4, dimension(:), allocatable :: real4_1D
!  real*8, dimension(:), allocatable :: real8_1D
  
!  call HDF5_open(iterwall_offset_20cm.h5)
  
!end subroutine wall_collision_init

subroutine read_array(wall)
  integer(HID_T), intent(in)        :: wall
  real*8        , dimension(:), allocatable      :: wall_array
  character(LEN=*), intent(in)      :: dsetname
  integer                           :: n
  n = size(wall%groups(1),1)
  allocate(wall_array(n))
  call HDF5_array1D_reading(wall,wall_array,group_name//"arr")
end subroutine read_array

subroutine side_operator(p,q,r,s,side)
  implicit none
  real(kind=8), dimension(3), intent(in)         :: p,q,r,s
  real(kind=8), dimension(6)                     :: a,b
  real(kind=8), intent(out)                      :: side
  a = plucker(p,q)
  b = plucker(r,s)
  side = a(1)*b(5) + a(2)*b(6) + a(3)*b(4) + a(4)*b(3) + a(5)*b(1) + a(6)*b(2)
end subroutine side_operator
  
subroutine intersection_lls(p,q,r,s,t,s0,s1,s2,intersect,check_sign)
  implicit none
  real(kind=8), dimension(3)        :: p, q, r, s, t
  real(kind=8)                      :: s0, s1, s2, check_sign
  real(kind=8), dimension(6)        :: L, L1, L2
  logical                           :: intersect
  
  
  t = cross(q-p,s-r)
  t = cylindrical_to_cartesian_real8(t)
  call side_operator(p,q,r,s,s0)
  if(s0 == 0.d0) then
    call side_operator(p,q,t,r,s1)
    call side_operator(p,q,t,s,s2)
    check_sign = sign(1.d0,s1*s2)

    if(check_sign < 0.d0) then
      intersect = .false.
      write(*,*) 'no intersection'
    else
      intersect = .true.
      write(*,*) 'passes through'
    endif
  endif


  
  if(s1*s2 == 0.d0 .and. s1/=s2) then
    write(*,*) 'line passes through one end point of the segment'
  elseif(s1*s2 == 0.d0 .and. s1==s2) then
    write(*,*) 'line contains the segment'
  end if
end subroutine intersection_lls
  
subroutine triangle_line_intersection_nc(p,q,v0,v1,v2,s1,s2,s3,sign_s1,sign_s2,sign_s3,proper,edge,vertex)    
  implicit none
  real(kind=8), dimension(3)        :: p, q, v0, v1, v2
  real(kind=8)                      :: s1, s2, s3, sign_s1, sign_s2, sign_s3
  logical                           :: proper, edge, vertex
  call side_operator(p,q,v0,v1,s1)
  call side_operator(p,q,v1,v2,s2)
  call side_operator(p,q,v2,v0,s3)
  
  proper = .false.
  edge   = .false.
  vertex = .false.

  sign_s1 = sign(1.d0,s1)
  sign_s2 = sign(1.d0,s2)
  sign_s3 = sign(1.d0,s3)

  if(sign_s1 == sign_s2 .and. sign_s1 == sign_s3) then
    proper = .true.
  endif
  if((s1*s2==0.d0 .and. s3==0.d0) .or. (s1*s3==0.d0 .and. s2==0.d0)) then
    vertex = .true.
  endif
  if(s1==0.d0 .or. s2*s3==0.d0) then
    edge = .true.
  endif  
end subroutine triangle_line_intersection_nc    

subroutine triangle_line_segment(p,q,v0,v1,v2,collision)
  implicit none
  real(kind=8), dimension(3)         :: p, q, t, v0, v1, v2
  real(kind=8)                       :: s1, s2, s3, s4, s5, s11, s12, s13, s21, s22, s23, s31, s32, s33, sign_s1, sign_s2, sign_s3, sign_s4, sign_s5, sign_s11, sign_s12, sign_s13, sign_s21, sign_s22, sign_s23, sign_s31, sign_s32, sign_s33
  logical                            :: collision, check, pr, ve, ed, pr1, ve1, ed1, pr2, ed2, ve2, pr3, ve3, ed3
  
  !calculate side operators
  call side_operator(p,q,v0,v1,s1)
  call side_operator(p,q,v1,v2,s2)
  call side_operator(p,q,v2,v0,s3)

  !check whether line and triangle are coplanar or not
  check = plane_check(p,q,v0,v1,v2)
  !Construct the side operators for side(L4,L3) (s4) and side(L4,L2) (s5)
  if(check == .false.) then

    call triangle_line_intersection_nc(p,q,v0,v1,v2,s1,s2,s3,sign_s1,sign_s2,sign_s3,pr,ve,ed)
    if(ve == .true.) then !intersection at vertex implies the common vertex is selected for the new triangles
      if(s1/=0.d0) then
        call side_operator(v0,v1,q,v2,s4)
        call side_operator(v0,v1,p,v2,s5)
      elseif(s2/=0.d0) then
        call side_operator(v1,v2,q,v0,s4)
        call side_operator(v1,v2,p,v0,s5)
      elseif(s3/=0.d0) then
        call side_operator(v2,v0,q,v1,s4)
        call side_operator(v2,v0,p,v1,s5)
      endif
    elseif(ve == .false. .and. (pr == .true. .or. ed == .true.)) then
      call side_operator(v0,v1,q,v2,s4)
      call side_operator(v0,v1,p,v2,s5)
    endif

    if(sign_s4 == sign_s5 .or. s4*s5==0.d0) then
      collision = .true.
    else
      collision = .false.
    endif
    
  else
    t = cross(v1-v0,v2-v0)   !define a point outside the plane
    call triangle_line_intersection_nc(p,q,v0,v1,t,s11,s12,s13,sign_s11,sign_s12,sign_s13,pr1,ve1,ed1)
    call triangle_line_intersection_nc(p,q,v1,v2,t,s21,s22,s23,sign_s21,sign_s22,sign_s23,pr2,ve2,ed2)
    call triangle_line_intersection_nc(p,q,v2,v0,t,s31,s32,s33,sign_s31,sign_s32,sign_s33,pr3,ve3,ed3)
    if(sign_s11 /= sign_s12 .and. sign_s21 /= sign_s22 .and. sign_s31 /= sign_s32) then
      collision = .false.
    else
      collision = .true.
    endif
  endif
end subroutine triangle_line_segment

end module mod_wall_collision
