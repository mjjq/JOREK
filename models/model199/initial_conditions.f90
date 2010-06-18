subroutine initial_conditions(my_id,node_list,element_list,xpoint2)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list

integer    :: my_id, i, in, mm, i_elm_axis, i_elm_xpoint
real*8     :: amplitude, psi, psi_axis
real*8     :: zn, dn_dpsi, dn_dpsi2, dn_dz, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi2_dz, dn_dpsi_dz2
real*8     :: zT, dT_dpsi, dT_dpsi2, dT_dz, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi2_dz, dT_dpsi_dz2
real*8     :: zFFprime,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi_dz, dFFprime_dz2, dFFprime_dpsi2
real*8     :: R_axis, Z_axis, s_axis, t_axis, R, Z, BigR
real*8     :: psi_n, psi_bnd,psi_xpoint,R_xpoint,Z_xpoint,s_xpoint,t_xpoint
logical    :: xpoint2

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions  (199)      *'
  write(*,*) '***************************************'
endif

if (my_id .eq. 0) then
   i_elm_axis = 1 ! XL : uninitilised value... So let's say it'll be 1, to look in the first case of element array in interp.f90...
   call find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis)
  if (xpoint2) call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint)


  psi_bnd = 0.d0
  if (xpoint2) psi_bnd = psi_xpoint

  do i=1,node_list%n_nodes

    psi = node_list%node(i)%values(1,1,1)
    R   = node_list%node(i)%x(1,1)
    Z   = node_list%node(i)%x(1,2)

    call density(    xpoint2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,             &
                                                               dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    call temperature(xpoint2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,             &
                                                               dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)

    call FFprime(    xpoint2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zFFprime,dFFprime_dpsi,dFFprime_dz, &
                                                               dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz)

    node_list%node(i)%values(1,1,5) = zn
    node_list%node(i)%values(1,2,5) = dn_dpsi  * node_list%node(i)%values(1,2,1) + dn_dz * node_list%node(i)%x(2,2)
    node_list%node(i)%values(1,3,5) = dn_dpsi  * node_list%node(i)%values(1,3,1) + dn_dz * node_list%node(i)%x(3,2)
    node_list%node(i)%values(1,4,5) = dn_dpsi  * node_list%node(i)%values(1,4,1) + dn_dz * node_list%node(i)%x(4,2) &
                                    + dn_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dn_dz2   * node_list%node(i)%x(2,2)        * node_list%node(i)%x(3,2)

    node_list%node(i)%values(1,1,6) = zT
    node_list%node(i)%values(1,2,6) = dT_dpsi  * node_list%node(i)%values(1,2,1) + dT_dz * node_list%node(i)%x(2,2)
    node_list%node(i)%values(1,3,6) = dT_dpsi  * node_list%node(i)%values(1,3,1) + dT_dz * node_list%node(i)%x(3,2)
    node_list%node(i)%values(1,4,6) = dT_dpsi  * node_list%node(i)%values(1,4,1) + dT_dz * node_list%node(i)%x(4,2) &
                                    + dT_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dT_dz2   * node_list%node(i)%x(2,2)        * node_list%node(i)%x(3,2)

  enddo

endif

!---------------------------- initialise perturbations

amplitude = 1.d-12
mm = 2

do in=2,n_tor

  if (my_id .eq. 0) then

    do i=1,node_list%n_nodes

      psi = node_list%node(i)%values(1,1,1)
      Z   = node_list%node(i)%x(1,2)

      psi_bnd = 0.d0
      if (xpoint2) psi_bnd = psi_xpoint
      psi_n = (psi - psi_axis)/(psi_bnd - psi_axis)

      node_list%node(i)%values(in,1,4) = amplitude * psi_n * (1.d0 -psi_n)
      node_list%node(i)%values(in,2,4) = amplitude * (1. - 2.d0 * psi_n)/(psi_bnd - psi_axis) * node_list%node(i)%values(1,2,1)
      node_list%node(i)%values(in,3,4) = amplitude * (1. - 2.d0 * psi_n)/(psi_bnd - psi_axis) * node_list%node(i)%values(1,3,1)
      node_list%node(i)%values(in,4,4) = amplitude * (1. - 2.d0 * psi_n)/(psi_bnd - psi_axis) * node_list%node(i)%values(1,4,1)


      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. (Z .lt. Z_xpoint))) then
        node_list%node(i)%values(in,1:4,4) = 0.d0
      endif

      node_list%node(i)%deltas = 0.d0

    enddo

  endif

  call poisson(my_id,1,node_list,element_list,4,2,in,xpoint2)   !----------- for Poisson (toroidal) use 1

! call poisson(my_id,2,node_list,element_list,4,2,in,xpoint2) !----------- for cylinder use 2

enddo


!call add_pellet(node_list,element_list,50.,0.06,0.02,R_axis-0.96,Z_axis)

return
end
