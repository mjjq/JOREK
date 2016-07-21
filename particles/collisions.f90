module collisions

implicit none

contains

subroutine random_velocity(m_b, density_b, T_b, Grad_Tb, f_w, f_tht, f_phi)

use constants
use phys_module, only : central_density, central_mass

! density_b is background number density (not mass density)
! T_b is in JOREK units: kT = T_b / (mu_zero * density_b)

implicit none

real*8 :: m_b         ! mass in proton masss units of the background plasma ions
real*8 :: density_b   ! in units of central_density (i.e JOREK units)
real*8 :: T_b         ! temperature (in JOREK units)
real*8 :: grad_Tb     ! temperature gradient (in JOREK units)
real*8 :: f_w, f_tht, f_phi

real*8 :: R5(5), w_2(3), w(3), rotation(3,3)
real*8 :: wx, wy, wz, f_w2, A_kt, alfa, cos_tht2, tht_2, phi_2

!--- call init_random_seed()

call random_number(R5)

!wx = sqrt(-2*Log(R4(1))) * cos(TWOPI*R4(2))     ! box-mueller transform
!wy = sqrt(-2*Log(R4(1))) * sin(TWOPI*R4(2))
!wz = sqrt(-2*Log(R4(3))) * cos(TWOPI*R4(4))
!f_w2 = sqrt(wx*wx+wy*wy+wz*wz)

zn_b = density_b * 1d20      ! [1/m^3]
kb_T = T_b /(MU_ZERO*central_density*1d20)   ! [J]

V_thermal = sqrt(kb_T / (m_b * PROTON_MASS)  ! [m/s]

f_w = V_thermal * sqrt( -2*Log(R5(1)) -2 * Log(R5(2))* cos(TWOPI*R5(3))**2 ) ! checked

A_kT = 37.5 * PI**1.5 * EPS_ZERO**2 / 17.d0 * kb_T / zn_b / EL_CHG**4

alfa = A_kt * (1.d0  - f_w**2 / (5.d0 * v_thermal**2) ) * (f_w/V_thermal) * grad_Tb

write(*,'(A,e14.6)') ' Ion mass (background) : ',m_b
write(*,'(A,e14.6)') ' Density               : ',zn_b
write(*,'(A,e14.6)') ' Temperature           : ',kb_T
write(*,'(A,e14.6)') ' Thermal velocity      : ',V_thermal

write(*,*) A_kt,alfa

if (alfa .gt. 0.d0) then
  if alfa .le. 1.d0) then
    cos_tht2 = ( sqrt(4.d0*alfa*R5(4) + (1-alfa)**2) - 1.d0 ) / alfa
  else
    cos_tht2 = - 1.d0 +  2.d0 * sqrt(R5(4))
  endif
elseif (alfa .lt. 0.d0) then
  if alfa .ge. -1.d0) then
    cos_tht2 = ( sqrt(-4.d0*alfa*R5(4) + (1+alfa)**2) - 1.d0 ) / alfa
  else
    cos_tht2 = 1.d0 - 2.d0 * sqrt(R5(4))
  endif
else
  cos_tht2 = 2*R5(4) - 1.d0
endif

f_tht = acos(cos_tht2)  ! in coordinate system II (fig.7) of Homma 2013
f_phi = TWOPI * R5(5)

return
end subroutine

end module

subroutine collision(R, Z, phi, v_in,  v_out)

implicit none

real*8, intent(in) :: R, Z, phi ! position of the of the particle
real*8, intent(in) :: v_in(3)   ! (R,Z,phi) components of the particle velocity

real*8, intent(out) :: v_out(3)  (R,Z,phi) components of the particle velocity after a collision

real*8  :: P(5), P_s(5), P_t(5), R, R_s, R_t, Z, Z_s, Z_t, P_phi(5)
real*8  :: B(3), V_fluid(3)
integer :: i_var(5)


! find density, temperature, temperature gradient and velocity of the background fluid at this position

i_var = (/ 1, 2, 5, 6, 7 /)

!call interp_PRZ(node_list,element_list,i_elm,i_var,2,R,Z,P,P_s,P_t,R,R_s,R_t,Z,Z_s,Z_t)

call interp_PRZ(node_list,element_list,i_elm,i_var,2,s_elm,t_elm,phi,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)

psi_R    = (  psi_s * Z_t - psi_t * Z_s ) / st_jac
psi_Z    = (- psi_s * R_t + psi_t * R_s ) / st_jac
U_R      = (  u_s   * Z_t - u_t   * Z_s ) / st_jac
U_Z      = (- u_s   * R_t + u_t   * R_s ) / st_jac
U_phi    = 0.d0 ! TODO

B        =  (/ psi_Z, - psi_R, F0 /) / R

V_fluid  = (/ - R * U_Z, +R * U_R, 0.d0 /) + v_par * B



return
end

program test_random
use collisions
implicit none
real*8 :: mb, Tb, grad_Tb, f_w, f_tht, f_phi, average
real*8, allocatable :: x(:), f(:), occurence(:)
integer :: n_box, n_points, n_scale, i

call random_seed()

n_box = 1000
n_points = 10
n_scale = 5

allocate(x(n_box),occurence(n_box),f(n_points))

do i=1,n_points

  call random_distribution(mb, Tb, Grad_Tb, f_w, f_tht, f_phi)

  f(i) = f_w

  if (f_w .lt. n_scale) then
     occurence(int(n_box*f_w/n_scale)+1) =  occurence(int(n_box*f_w/n_scale)+1) + 1
  endif

enddo

do i=1,n_box
  x(i) = real(i-1)*n_scale
enddo

occurence = occurence / real(n_points)

average = sum(f)/real(n_points)

write(*,*) ' average  : ',average
write(*,*) ' variance : ',sum((f-average)**2)/real(n_points)

call begplt('random.ps')
call hplot6(11,1,x,occurence,n_box,'distribution')
call finplt

end
