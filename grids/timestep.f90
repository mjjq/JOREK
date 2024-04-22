subroutine timestep(my_id, node_list, element_list, local_elms, n_local_elms, cfl, tstep_cfl)
  use data_structure, only: type_node_list, type_element_list
  use phys_module, only: gamma
  use mod_parameters
  use mod_math_operators
  use basis_at_gaussian
  
  implicit none
  integer                                    :: my_id
  type(type_node_list),        intent(in)    :: node_list
  type(type_element_list),     intent(in)    :: element_list
  real*8,                      intent(in)    :: cfl
  real*8,                      intent(out)   :: tstep_cfl

  integer                                    :: local_elms(*)
  integer                                    :: n_local_elms

  ! --- local variables
  integer :: i, j, k
  integer :: i_tor, i_plane, ielm, iv, inode
  real*8  :: he_min
  real*8, dimension(n_vertex_max) :: xjac, lambda
  real*8, dimension(n_vertex_max) :: R,R_s,R_t,Z,Z_s,Z_t
  real*8, dimension(n_vertex_max) :: Temperature, Te, Ti, density, BR, BZ, Bphi, VR, VZ, Vphi, BB2, VV2
  real*8, dimension(n_vertex_max) :: AR, AZ, A3, AR_Z, AR_p, AZ_R, AZ_p, A3_R, A3_Z, F_prof
  real*8, dimension(n_vertex_max) :: AR_s, AR_t, AZ_s, AZ_t, A3_s, A3_t

  tstep_cfl = 1.d10
#ifdef fullmhd
!   do ielm = 1, element_list%n_elements
  do i = 1, n_local_elms
    ielm = local_elms(i)
    do i_plane = 1, n_plane
      Te = 0.d0; Ti = 0.d0; Temperature = 0.d0
      BR = 0.d0; BZ = 0.d0; Bphi = 0.d0
      VR = 0.d0; VZ = 0.d0; Vphi = 0.d0
      AR = 0.d0; AZ = 0.d0; A3 = 0.d0
      AR_Z = 0.d0; AR_p = 0.d0; AZ_R = 0.d0; AZ_p = 0.d0; A3_R = 0.d0; A3_Z = 0.d0
      F_prof = 0.d0
      do iv = 1, n_vertex_max
        inode = element_list%element(ielm)%vertex(iv)

        R(iv)    = node_list%node(inode)%x(1, 1, 1)
        R_s(iv)  = node_list%node(inode)%x(1, 2, 1)
        R_t(iv)  = node_list%node(inode)%x(1, 3, 1)
        Z(iv)    = node_list%node(inode)%x(1, 1, 2)
        Z_s(iv)  = node_list%node(inode)%x(1, 2, 2)
        Z_t(iv)  = node_list%node(inode)%x(1, 3, 2)
            
        xjac(iv) = op_xjac(R_s(iv), R_t(iv), Z_s(iv), Z_t(iv))
        if(dabs(xjac(iv)) .le. 1.d-10) xjac = 1.d10
            
        F_prof(iv)  = node_list%node(inode)%Fprof_eq(1)
        do i_tor = 1, n_tor
          density(iv) = density(iv) + node_list%node(inode)%values(i_tor, 1, var_rho) * HZ(i_tor,i_plane)
          if (with_TiTe) then
            Te(iv)      = Te(iv)      + node_list%node(inode)%values(i_tor, 1, var_Te)  * HZ(i_tor,i_plane)
            Ti(iv)      = Ti(iv)      + node_list%node(inode)%values(i_tor, 1, var_Ti)  * HZ(i_tor,i_plane)
          else
            temperature(iv)      = temperature(iv)      + node_list%node(inode)%values(i_tor, 1, var_T)  * HZ(i_tor,i_plane)
          endif
          VR(iv)      = VR(iv)      + node_list%node(inode)%values(i_tor, 1, var_uR)  * HZ(i_tor,i_plane)
          VZ(iv)      = VZ(iv)      + node_list%node(inode)%values(i_tor, 1, var_uZ)  * HZ(i_tor,i_plane)
          Vphi(iv)    = Vphi(iv)    + node_list%node(inode)%values(i_tor, 1, var_up)  * HZ(i_tor,i_plane)
                
          AR(iv)      = AR(iv)      + node_list%node(inode)%values(i_tor, 1, var_AR)  * HZ(i_tor,i_plane)
          AZ(iv)      = AZ(iv)      + node_list%node(inode)%values(i_tor, 1, var_AZ)  * HZ(i_tor,i_plane)
          A3(iv)      = A3(iv)      + node_list%node(inode)%values(i_tor, 1, var_A3)  * HZ(i_tor,i_plane)

          AR_s(iv)    = node_list%node(inode)%values(i_tor, 2, var_AR)
          AR_t(iv)    = node_list%node(inode)%values(i_tor, 3, var_AR)
          AZ_s(iv)    = node_list%node(inode)%values(i_tor, 2, var_AZ)
          AZ_t(iv)    = node_list%node(inode)%values(i_tor, 3, var_AZ)
          A3_s(iv)    = node_list%node(inode)%values(i_tor, 2, var_A3)
          A3_t(iv)    = node_list%node(inode)%values(i_tor, 3, var_A3)

          AR_Z(iv)    = AR_Z(iv)    + op_grad_Z(AR_s(iv), AR_t(iv), R_s(iv), R_t(iv), Z_s(iv), Z_t(iv), xjac(iv)) * HZ(i_tor,i_plane)
          AZ_R(iv)    = AZ_R(iv)    + op_grad_R(AZ_s(iv), AZ_t(iv), R_s(iv), R_t(iv), Z_s(iv), Z_t(iv), xjac(iv)) * HZ(i_tor,i_plane)
          A3_Z(iv)    = A3_Z(iv)    + op_grad_Z(A3_s(iv), A3_t(iv), R_s(iv), R_t(iv), Z_s(iv), Z_t(iv), xjac(iv)) * HZ(i_tor,i_plane)
          A3_R(iv)    = A3_R(iv)    + op_grad_R(A3_s(iv), A3_t(iv), R_s(iv), R_t(iv), Z_s(iv), Z_t(iv), xjac(iv)) * HZ(i_tor,i_plane)

          AR_p(iv)    = AR_p(iv)    + AR(iv) * HZ_p(i_tor,i_plane)
          AZ_p(iv)    = AZ_p(iv)    + AZ(iv) * HZ_p(i_tor,i_plane)

        enddo ! i_tor

        if (with_TiTe) then
          Temperature(iv)  = Te(iv) + Ti(iv)
        endif
        BR(iv)             = ( A3_Z(iv) - AZ_p(iv) ) / R(iv)
        BZ(iv)             = ( AR_p(iv) - A3_R(iv) ) / R(iv)
        Bphi(iv)           = ( AZ_R(iv) - AR_Z(iv) )    + F_prof(iv) / R(iv)

        BB2(iv)            = BR(iv) * BR(iv) + BZ(iv) * BZ(iv) + Bphi(iv) * Bphi(iv)
        VV2(iv)            = VR(iv) * VR(iv) + VZ(iv) * VZ(iv) + Vphi(iv) * Vphi(iv)
            
        lambda(iv)         = dsqrt( gamma * Temperature(iv) +  BB2(iv) / density(iv) ) + dsqrt(VV2(iv))
      enddo ! iv

      he_min     = dsqrt(min(  ((R(1) - R(3))**2 + (Z(1) - Z(3))**2), ((R(2) - R(4))**2 + (Z(2) - Z(4))**2)  ))
      tstep_cfl = min(tstep_cfl, cfl * he_min / maxval(lambda))

    enddo ! i_plane

  enddo ! i or ielm
  ! write(*,*) "my_id = ", my_id, ", tstep_cfl = ", tstep_cfl
#else
  write(*,*) "Warning, adaptive_time has not been implemented in reduced MHD models."
  stop
#endif

end subroutine timestep