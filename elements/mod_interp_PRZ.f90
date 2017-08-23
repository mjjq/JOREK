module mod_interp_PRZ
  implicit none
  private
  public :: interp_RZ, interp_PRZ, interp_PRZ_delta
contains
!> This subroutine interpolates space a specific position within one element at a given position (s,t)
pure subroutine interp_RZ(node_list, element_list, i_elm, s, t, R, R_s, R_t, Z, Z_s, Z_t)

use data_structure
use mod_basisfunctions
implicit none

! --- Routine parameters
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
real*8,                   intent(in)  :: s, t
real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t

! --- Local variables
real*8  :: H(4,4), H_s(4,4), H_t(4,4)
integer :: kv, iv
real*8  :: xR(n_order+1,n_vertex_max), xZ(n_order+1,n_vertex_max)
real*8  :: sizes(n_order+1)

call basisfunctions6(s,t,H,H_s,H_t)

R = 0.d0; R_s = 0.d0; R_t = 0.d0;
Z = 0.d0; Z_s = 0.d0; Z_t = 0.d0;

! Preload values and premultiply with sizes(:,kv)
do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  sizes(:) = element_list%element(i_elm)%size(kv,:)
  xR(:,kv) = node_list%node(iv)%x(:,1) * sizes(:)
  xZ(:,kv) = node_list%node(iv)%x(:,2) * sizes(:)
end do

R   = sum(xR*H)
R_s = sum(xR*H_s)
R_t = sum(xR*H_t)
Z   = sum(xZ*H)
Z_s = sum(xZ*H_s)
Z_t = sum(xZ*H_t)
end subroutine interp_RZ


!> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
pure subroutine interp_PRZ(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)

use data_structure
use mod_parameters, only: n_period, n_tor
use mod_basisfunctions
implicit none

! --- Routine parameters
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: n_v, i_v(n_v)
real*8,                   intent(in)  :: s, t, phi
real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v)
real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
real*8,                   intent(out) :: P_phi(n_v)

! --- Local variables
real*8  :: H(4,4), H_s(4,4), H_t(4,4), HZ(n_tor), dHZ(n_tor)
integer :: kv, iv, kf, i, i_tor
real*8  :: e1r, e1i, enr, eni, enrtmp
real*8  :: values(n_tor,n_order+1,n_v,n_vertex_max)
real*8  :: xR(n_order+1,n_vertex_max), xZ(n_order+1,n_vertex_max)
real*8  :: sizes(n_order+1), v, vp

! 7% exec time
call basisfunctions6(s,t,H,H_s,H_t)

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_phi = 0.d0

! 7% exec time
call sincosperiod_moivre(phi, HZ, dHZ)

! 30% exec time
! Preload values and premultiply with sizes(:,kv)
do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  sizes(:) = element_list%element(i_elm)%size(kv,:)
  do i = 1, n_v
    do kf=1,n_order+1
      values(1:n_tor,kf,i,kv) = node_list%node(iv)%values(1:n_tor,kf,i_v(i)) * sizes(kf)
    end do
  end do
  xR(:,kv) = node_list%node(iv)%x(:,1) * sizes(:)
  xZ(:,kv) = node_list%node(iv)%x(:,2) * sizes(:)
end do

! together 7%
R   = sum(xR*H)
R_s = sum(xR*H_s)
R_t = sum(xR*H_t)
Z   = sum(xZ*H)
Z_s = sum(xZ*H_s)
Z_t = sum(xZ*H_t)

! 40% exec time
do kv = 1, n_vertex_max
  do i = 1, n_v
    do kf = 1, n_order+1
      v = dot_product(values(1:n_tor,kf,i,kv),HZ(1:n_tor))
      P(i)     = P(i)     + v * H(kf, kv)
      P_s(i)   = P_s(i)   + v * H_s(kf, kv)
      P_t(i)   = P_t(i)   + v * H_t(kf, kv)
      vp = dot_product(values(1:n_tor,kf,i,kv),dHZ(1:n_tor))
      P_phi(i) = P_phi(i) + vp * H(kf, kv)
    enddo
  enddo
enddo
end subroutine interp_PRZ

! Apply De Moivre formula to calculate the series of sines.
! Assumes that mode is of the form [0 1 1 2 2 3 3 4 4] ([0 4 4 8 8 12 12])
! This is roughly 3-4 times faster in my tests than just calculating the sines
! and cosines (even when that is vectorized). Perhaps that changes for n_tor >> 10
! I tested n_tor = 17.
pure subroutine sincosperiod_moivre(phi,HZ,dHZ)
  use mod_parameters, only: n_tor, n_period
  integer, parameter :: n_mode = (n_tor-1)/2 ! number of modes excluding 0
  real*8, intent(in) :: phi
  real*8, intent(out) :: HZ(n_tor), dHZ(n_tor)

  integer :: i
  HZ(1) = 1.d0
  dHZ(1) = 0.d0
  if (n_mode .gt. 0) then
    HZ(2) = cos(n_period*phi)
    HZ(3) = sin(n_period*phi)
    dHZ(2) = HZ(3)*(-n_period)
    dHZ(3) = HZ(2)*(n_period)

    do i=2,n_mode
      call moivre(HZ(2),HZ(3), &
                  HZ(2),HZ(3), &
                  HZ(2*i),HZ(2*i+1))
      dHZ(2*i)   = HZ(2*i+1)*(-n_period*i)
      dHZ(2*i+1) = HZ(2*i)*(n_period*i)
    end do
  end if
end subroutine sincosperiod_moivre

pure subroutine moivre(ar,ai,br,bi,or,oi)
  real*8, intent(in) :: ar, ai !< real and imag part of e^(i x)
  real*8, intent(in) :: br, bi !< real and imag part of e^(i y)
  real*8, intent(out) :: or, oi !< real and imag part of e^(i (x+y))
  or = ar*br - ai*bi
  oi = ai*br + ar*bi
end subroutine moivre

!> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
pure subroutine interp_PRZ_delta(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)

use data_structure
use mod_parameters, only: n_period
use mod_basisfunctions
implicit none

! --- Routine parameters
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: n_v, i_v(n_v)
real*8,                   intent(in)  :: s, t, phi
real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v)
real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
real*8,                   intent(out) :: P_phi(n_v)

integer, parameter :: n_mode = (n_tor-1)/2 ! number of modes excluding 0

! --- Local variables
real*8  :: H(4,4), H_s(4,4), H_t(4,4), HZ(n_tor), dHZ(n_tor)
integer :: kv, iv, kf, i, i_tor
real*8  :: e1r, e1i, enr, eni, enrtmp
real*8  :: deltas(n_tor,n_order+1,n_v,n_vertex_max)
real*8  :: xR(n_order+1,n_vertex_max), xZ(n_order+1,n_vertex_max)
real*8  :: sizes(n_order+1), v, vp

call basisfunctions6(s,t,H,H_s,H_t)

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_phi = 0.d0

call sincosperiod_moivre(phi, HZ, dHZ)

! Preload deltas and premultiply with sizes(:,kv)
do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  sizes(:) = element_list%element(i_elm)%size(kv,:)
  do i = 1, n_v
    do kf=1,n_order+1
      deltas(1:n_tor,kf,i,kv) = node_list%node(iv)%deltas(1:n_tor,kf,i_v(i)) * sizes(kf)
    end do
  end do
  xR(:,kv) = node_list%node(iv)%x(:,1) * sizes(:)
  xZ(:,kv) = node_list%node(iv)%x(:,2) * sizes(:)
end do

R   = sum(xR*H)
R_s = sum(xR*H_s)
R_t = sum(xR*H_t)
Z   = sum(xZ*H)
Z_s = sum(xZ*H_s)
Z_t = sum(xZ*H_t)

do kv = 1,n_vertex_max
  do i = 1, n_v
    do kf = 1, n_order+1
      v = sum(deltas(1:n_tor,kf,i,kv)*HZ(1:n_tor))
      P(i)     = P(i)     + v * H(kf, kv)
      P_s(i)   = P_s(i)   + v * H_s(kf, kv)
      P_t(i)   = P_t(i)   + v * H_t(kf, kv)
      vp = sum(deltas(1:n_tor,kf,i,kv)*dHZ(1:n_tor))
      P_phi(i) = P_phi(i) + vp * H(kf, kv)
    enddo
  enddo
enddo
end subroutine interp_PRZ_delta
end module mod_interp_PRZ
