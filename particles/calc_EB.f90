!> Calculates the electric and magnetic fields at a specific position
subroutine calc_EB(i_elm,st,phi,E,B,psi,U)
use parameters
use nodes_elements
use constants
use phys_module, only : F0, central_mass, central_density

implicit none

! Routine parameters
integer, intent(in) :: i_elm
real*8, intent(in)  :: st(2), phi

real*8, intent(out) :: E(3), B(3), psi, U

! Internal parameters
integer :: i_var(2)
real*8                    :: P(2), P_s(2), P_t(2), P_phi(2) ! Placeholder for evaluating variables and derivatives locally
real*8                    :: psi_R, psi_Z, U_R, U_Z, U_phi
real*8                    :: R, R_s, R_t, Z, Z_s, Z_t, inv_st_jac, R_inv

! Select psi and U
i_var = (/1,2/)

! Interpolate the fields to get psi and U at the current position
call interp_PRZ(node_list,element_list,i_elm,i_var,2,st(1),st(2),phi,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)

psi = P(1)
U   = P(2)
R_inv = 1.d0/R

inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
! these are the local derivatives to s and t of the flux and potential
!psi_s = P_s(1); psi_t = P_t(1)
!u_s   = P_s(2); u_t   = P_t(2)

! Calculate the derivatives to R and Z
psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
U_R      = (  P_s(2) * Z_t - P_t(2) * Z_s ) * inv_st_jac
U_Z      = (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac
! And to phi
!psi_phi  = P_phi(1) ! Not used now, assumed 0!
U_phi    = P_phi(2)

B     = (/ + psi_Z, - psi_R, F0 /) * R_inv

! The local electric field, obtained from E=-Grad (u F0)
E     = (/ - F0 * U_R, - F0 * U_Z, - F0 * U_phi * R_inv /)
end subroutine calc_EB
