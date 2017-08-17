module mod_interp_PRZ
contains
!> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
pure subroutine interp_PRZ(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)

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
real*8  :: values(n_tor,n_order+1,n_v,n_vertex_max)
real*8  :: xR(n_order+1,n_vertex_max), xZ(n_order+1,n_vertex_max)
real*8  :: sizes(n_order+1), v, vp

call basisfunctions6(s,t,H,H_s,H_t)

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_phi = 0.d0

! Apply De Moivre formula to calculate the series of sines.
! Assumes that mode is of the form [0 1 1 2 2 3 3 4 4] ([0 4 4 8 8 12 12])
! This is roughly 3-4 times faster in my tests than just calculating the sines
! and cosines (even when that is vectorized). Perhaps that changes for n_tor >> 10
! I tested n_tor = 17.
e1r = cos(n_period*phi)
e1i = sin(n_period*phi)
enr = 1.d0
eni = 0.d0
HZ(1) = 1.d0
dHZ(1) = 0.d0
do i=1,n_mode
  enrtmp = enr
  enr = enr*e1r - eni*e1i
  eni = eni*e1r + enrtmp*e1i
  HZ(2*i)   = enr
  HZ(2*i+1) = eni
  dHZ(2*i) = eni*(-n_period*i)
  dHZ(2*i+1) = enr*(n_period*i)
end do

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

R   = sum(xR*H)
R_s = sum(xR*H_s)
R_t = sum(xR*H_t)
Z   = sum(xZ*H)
Z_s = sum(xZ*H_s)
Z_t = sum(xZ*H_t)

do kv = 1,n_vertex_max
  do i = 1, n_v
    do kf = 1, n_order+1
      v = sum(values(1:n_tor,kf,i,kv)*HZ(1:n_tor))
      P(i)     = P(i)     + v * H(kf, kv)
      P_s(i)   = P_s(i)   + v * H_s(kf, kv)
      P_t(i)   = P_t(i)   + v * H_t(kf, kv)
      vp = sum(values(1:n_tor,kf,i,kv)*dHZ(1:n_tor))
      P_phi(i) = P_phi(i) + vp * H(kf, kv)
    enddo
  enddo
enddo
end subroutine interp_PRZ

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

! Apply De Moivre formula to calculate the series of sines.
! Assumes that mode is of the form [0 1 1 2 2 3 3 4 4] ([0 4 4 8 8 12 12])
! This is roughly 3-4 times faster in my tests than just calculating the sines
! and cosines (even when that is vectorized). Perhaps that changes for n_tor >> 10
! I tested n_tor = 17.
e1r = cos(n_period*phi)
e1i = sin(n_period*phi)
enr = 1.d0
eni = 0.d0
HZ(1) = 1.d0
dHZ(1) = 0.d0
do i=1,n_mode
  enrtmp = enr
  enr = enr*e1r - eni*e1i
  eni = eni*e1r + enrtmp*e1i
  HZ(2*i)   = enr
  HZ(2*i+1) = eni
  dHZ(2*i) = eni*(-n_period*i)
  dHZ(2*i+1) = enr*(n_period*i)
end do

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
