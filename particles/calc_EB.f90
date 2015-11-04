!> Calculates the electric and magnetic fields at a specific position
subroutine calc_EB(i_elm,st,phi,E,B,psi,U)
use parameters
use nodes_elements
use constants
use phys_module, only : F0, central_mass, central_density

! Routine parameters
integer, intent(in) :: i_elm
real*8, intent(in)  :: st(2), phi

real*8, intent(out) :: E(3), B(3)

! Internal parameters
integer :: i_var(2)
real*8                    :: P(2), P_s(2), P_t(2) ! Placeholder for evaluating variables and derivatives locally
real*8                    :: qom, B02, psi_R, psi_Z, U_R, U_Z, U_phi, U
real*8                    :: psi, psi_s, psi_t, u_s, u_t, psi_prev, U_prev
real*8                    :: R, R_s, R_t, Z, Z_s, Z_t, st_jac

! Select psi and U
i_var = (/1,2/)

! Interpolate the fields to get psi and U at the current position
call interp_PRZ(node_list,element_list,i_elm,i_var,2,st(1),st(2),phi,P,P_s,P_t,R,R_s,R_t,Z,Z_s,Z_t)

psi = P(1)
U   = P(2)

st_jac = R_s * Z_t - R_t * Z_s
! these are the local derivatives to s and t of the flux and potential
psi_s = P_s(1); psi_t = P_t(1)
u_s   = P_s(2); u_t   = P_t(2)

! Calculate the derivatives to R and Z
psi_R    = (  psi_s * Z_t - psi_t * Z_s ) / st_jac
psi_Z    = (- psi_s * R_t + psi_t * R_s ) / st_jac
U_R      = (  u_s   * Z_t - u_t   * Z_s ) / st_jac
U_Z      = (- u_s   * R_t + u_t   * R_s ) / st_jac
! And assume for now no electric field in the phi-direction
U_phi    = 0.d0

B     = (/ + psi_Z, - psi_R, F0 /) / R

! The local electric field, obtained from E=-Grad (u F0)
E     = (/ - F0 * U_R, - F0 * U_Z, - F0 * U_phi / R /)
end subroutine calc_EB
