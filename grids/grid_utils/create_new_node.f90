subroutine create_new_node(node_list, element_list, newnode_list, index, i2, j2, nwpts)
!------------------------------------------------------------------------------------------
! subroutine defines a new node with the given index
!------------------------------------------------------------------------------------------

use tr_module 
use data_structure
use grid_xpoint_data
use mod_interp

implicit none

! --- Routine parameters
type (type_node_list)       , intent(inout) :: node_list
type (type_node_list)       , intent(inout) :: newnode_list
type (type_element_list)    , intent(inout) :: element_list
type (type_new_points)      , intent(in)    :: nwpts
integer,                      intent(in)    :: i2, j2, index

! --- local variables
integer             :: m
real*8              :: R_cub1d(4), Z_cub1d(4), dR_dt, dZ_dt, RZ_jac, PSI_R, PSI_Z, tmp1, tmp2
real*8              :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8              :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8              :: PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss


  m = nwpts%k_cross(i2,j2)
  
  R_cub1d = (/ nwpts%R_polar(m,1,j2), 3.d0/2.d0 *(nwpts%R_polar(m,2,j2)-nwpts%R_polar(m,1,j2)), &
  	       nwpts%R_polar(m,4,j2), 3.d0/2.d0 *(nwpts%R_polar(m,4,j2)-nwpts%R_polar(m,3,j2))  /)
  Z_cub1d = (/ nwpts%Z_polar(m,1,j2), 3.d0/2.d0 *(nwpts%Z_polar(m,2,j2)-nwpts%Z_polar(m,1,j2)), &
  	       nwpts%Z_polar(m,4,j2), 3.d0/2.d0 *(nwpts%Z_polar(m,4,j2)-nwpts%Z_polar(m,3,j2)) /)
  
  call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),nwpts%t_tht(i2,j2),tmp1, dR_dt)
  call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),nwpts%t_tht(i2,j2),tmp2, dZ_dt)

  call interp_RZ(node_list,element_list,nwpts%ielm_flux(i2,j2),nwpts%s_flux(i2,j2),nwpts%t_flux(i2,j2), &
  		 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
  		 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  call interp(node_list,element_list,nwpts%ielm_flux(i2,j2),1,1,nwpts%s_flux(i2,j2),nwpts%t_flux(i2,j2),&
  		 PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

  RZ_jac  = DRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr
  
  PSI_R = (   dPSg1_dr * dZZg1_ds - dPSg1_ds * dZZg1_dr ) / RZ_jac
  PSI_Z = ( - dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr ) / RZ_jac
  
  newnode_list%node(index)%x(1,:) = (/ nwpts%RR_new(i2,j2), nwpts%ZZ_new(i2,j2) /)
  newnode_list%node(index)%x(2,:) = (/ dR_dt, dZ_dt /)   / sqrt(dR_dt**2 + dZ_dt**2)
  newnode_list%node(index)%x(3,:) = (/ -PSI_Z, +PSI_R /) / sqrt(PSI_R**2 + PSI_Z**2)
  newnode_list%node(index)%x(4,:) = 0.d0
  newnode_list%node(index)%boundary = 0

return
end subroutine create_new_node
