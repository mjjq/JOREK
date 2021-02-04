subroutine initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint2, xcase2)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use mod_poiss
use equil_info
use mod_interp, only: interp

implicit none

type (type_node_list)        :: node_list
type (type_element_list)     :: element_list
type (type_surface_list)     :: surface_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

integer    :: my_id, i, in, mm, i_elm_axis, i_elm_xpoint(2), ifail, i_elm, xcase2
real*8     :: amplitude, psi, psi_n
real*8     :: zn
real*8     ::    dn_dpsi, dn_dz                                           ! 1st order derivatives
real*8     ::    dn_dpsi2, dn_dz2, dn_dpsi_dz                             ! 2nd order derivatives
real*8     ::    dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3              ! 2rd order derivatives
real*8     ::    dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz, dn_dz4 ! 4th order derivatives
real*8     :: zT
real*8     ::    dT_dpsi,  dT_dz                                          ! 1st order derivatives
real*8     ::    dT_dpsi2, dT_dz2, dT_dpsi_dz                             ! 2nd order derivatives
real*8     ::    dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,  dT_dz3              ! 2rd order derivatives
real*8     ::    dT_dpsi4, dT_dpsi_dz3, dT_dpsi2_dz2, dT_dpsi3_dz, dT_dz4 ! 4th order derivatives
real*8     :: zFFprime
real*8     ::    dFF_dpsi, dFF_dz                                                ! 1st order derivatives
real*8     ::    dFF_dpsi2, dFF_dz2, dFF_dpsi_dz                                 ! 2nd order derivatives
real*8     ::    dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3                 ! 2rd order derivatives
real*8     ::    dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4   ! 4th order derivatives
real*8     :: R, Z, BigR
real*8     :: R_out, Z_out, s_out, t_out, R_lim, Z_lim, s_lim, t_lim, psi_lim
real*8     :: p_s, p_t, p_ss, p_st, p_tt
logical    :: xpoint2

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions  (199)      *'
  write(*,*) '***************************************'
endif

if (my_id .eq. 0) then

  do i=1,node_list%n_nodes

    psi = node_list%node(i)%values(1,1,1)
    R   = node_list%node(i)%x(1,1,1)
    Z   = node_list%node(i)%x(1,1,2)

    call density(    xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zn,             &
                 dn_dpsi, dn_dz, &                                        ! 1st order derivatives
                 dn_dpsi2, dn_dz2, dn_dpsi_dz, &                          ! 2nd order derivatives
                 dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz,  dn_dz3, &           ! 2rd order derivatives
                 dn_dpsi4, dn_dpsi_dz3, dn_dpsi2_dz2, dn_dpsi3_dz, dn_dz4)! 4th order derivatives

    call temperature(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zT,             &
                     dT_dpsi,  dT_dz, &                                       ! 1st order derivatives
                     dT_dpsi2, dT_dz2, dT_dpsi_dz, &                          ! 2nd order derivatives
                     dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz,  dT_dz3, &           ! 2rd order derivatives
                     dT_dpsi4, dT_dpsi_dz3, dT_dpsi2_dz2, dT_dpsi3_dz, dT_dz4)! 4th order derivatives

    call FFprime(    xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zFFprime, &
                 dFF_dpsi, dFF_dz, &                                             ! 1st order derivatives
                 dFF_dpsi2, dFF_dz2, dFF_dpsi_dz, &                              ! 2nd order derivatives
                 dFF_dpsi3, dFF_dpsi_dz2, dFF_dpsi2_dz,  dFF_dz3, &              ! 2rd order derivatives
                 dFF_dpsi4, dFF_dpsi_dz3, dFF_dpsi2_dz2, dFF_dpsi3_dz, dFF_dz4, &! 4th order derivatives
                 .true.)

    node_list%node(i)%values(1,1,5) = zn
    node_list%node(i)%values(1,2,5) = dn_dpsi  * node_list%node(i)%values(1,2,1) + dn_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,5) = dn_dpsi  * node_list%node(i)%values(1,3,1) + dn_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,5) = dn_dpsi  * node_list%node(i)%values(1,4,1) + dn_dz * node_list%node(i)%x(1,4,2) &
                                    + dn_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dn_dz2   * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2) &
                                    + dn_dpsi_dz * ( &
                                      + node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,2,2) &
                                      + node_list%node(i)%values(1,2,1) * node_list%node(i)%x(1,3,2) &
                                      )

    node_list%node(i)%values(1,1,6) = zT
    node_list%node(i)%values(1,2,6) = dT_dpsi  * node_list%node(i)%values(1,2,1) + dT_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,6) = dT_dpsi  * node_list%node(i)%values(1,3,1) + dT_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,6) = dT_dpsi  * node_list%node(i)%values(1,4,1) + dT_dz * node_list%node(i)%x(1,4,2) &
                                    + dT_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dT_dz2   * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2) &
                                    + dT_dpsi_dz * ( &
                                      + node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,2,2) &
                                      + node_list%node(i)%values(1,2,1) * node_list%node(i)%x(1,3,2) &
                                      )

  enddo

endif

!---------------------------- initialise perturbations

amplitude = 1.d-10
mm = 2

do in=2,n_tor

  if (my_id .eq. 0) then

    do i=1,node_list%n_nodes

      node_list%node(i)%values(in,:,:) = 0.d0

      psi = node_list%node(i)%values(1,1,1)
      Z   = node_list%node(i)%x(1,1,2)

      psi_n = (psi - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)

      node_list%node(i)%values(in,1,4) = amplitude * psi_n * (1.d0 -psi_n)
      node_list%node(i)%values(in,2,4) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,2,1)
      node_list%node(i)%values(in,3,4) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,3,1)
      node_list%node(i)%values(in,4,4) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,4,1)

      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .lt. ES%Z_xpoint(1)) .and. (xcase2 .ne. 2)) ) ) then
        node_list%node(i)%values(in,1:4,4) = 0.d0
      endif
      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .gt. ES%Z_xpoint(2)) .and. (xcase2 .ne. 1)) ) ) then
        node_list%node(i)%values(in,1:4,4) = 0.d0
      endif

      node_list%node(i)%deltas = 0.d0

    enddo

  endif

  call Poisson(my_id,1,node_list,element_list,bnd_node_list,bnd_elm_list, &
               4,2,in, ES%psi_axis,ES%psi_bnd,xpoint2, xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)

enddo


!call add_pellet(node_list,element_list,50.,0.06,0.02,ES%R_axis-0.96,ES%Z_axis)

return
end
