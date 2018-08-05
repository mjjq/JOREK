subroutine flux_surface_add_point(node_list,element_list,surface_list,s,i_elm,iv,ifound,r_psi,s_psi,dpsi_dr,dpsi_ds)
use mod_interp, only: interp
use data_structure

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list

integer :: i_elm, iv, ifound, node1, node2, node3, node4
real*8  :: s, r_psi(*), s_psi(*), dpsi_dr(*), dpsi_ds(*), psi_test, P_st, P_ss, P_tt

if (iv .eq. 1) then
  r_psi(ifound) = s             ; s_psi(ifound) = 0.d0
elseif (iv .eq. 2) then
  r_psi(ifound) = 1.d0          ; s_psi(ifound) = s
elseif (iv .eq. 3) then
  r_psi(ifound) = 1.d0 - s      ; s_psi(ifound) = 1.d0
elseif (iv .eq. 4) then
  r_psi(ifound) = 0.d0          ; s_psi(ifound) = 1.d0 - s
endif

!write(*,*) ' add point : ',r_psi(ifound),s_psi(ifound)

!node1 = elm_list(1,i);  node2 = elm_list(2,i);  node3 = elm_list(3,i);  node4 = elm_list(4,i)
!call INTERP2(psi(:,node1),psi(:,node2),psi(:,node3),psi(:,node4), &
!             r_psi(ifound),s_psi(ifound),psi_test,dpsi_dr(ifound),dpsi_ds(ifound))

call interp(node_list,element_list,i_elm,1,1,r_psi(ifound),s_psi(ifound),psi_test,dpsi_dr(ifound),dpsi_ds(ifound),P_st,P_ss,P_tt)

!write(*,*) ' add point test : ',psi_test

return
end
