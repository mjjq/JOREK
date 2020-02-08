module mod_boundary_conditions
contains
  !*******************************************************************************
  !* Subroutine: boundary_condition                                              *
  !*******************************************************************************
  !*                                                                             *
  !* Add boundary condition on the matrix.                                       *
  !*                                                                             *
  !* Parameters:                                                                 *
  !*   my_id        - Identifier of the node in MPI_COMM_WORLD                   *
  !*   node_list    - List of nodes                                              *
  !*   element_list - List of all elements                                       *
  !*   local_elms   - List of local elements                                     *
  !*   n_local_elms - Number of local elements                                   *
  !*   index_min    - Minimal index of local elements                            *
  !*   index_max    - Maximal index of local elements                            *
  !*   xpoint2      -                                                            *
  !*   xcase2       -                                                            *
  !*   psi_axis     -                                                            *
  !*   psi_bnd      -                                                            *
  !*   Z_xpoint     -                                                            *
  !*   gmres        - boolean indicating if we are using GMRES method            *
  !*   solve_only   - Indicate if we want to perform only solve                  *
  !*                                                                             *
  !*                                                                             *
  !* Note: At present only Dirichlet boundary conditions have been implemented   *
  !*                                                                             *
  !*******************************************************************************
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,    & 
       &                          n_local_elms, index_min, index_max, rhs_loc, xpoint2,   &
       &                          xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint, gmres, solve_only )

    use data_structure
    use global_distributed_matrix
    use phys_module, only: F0, GAMMA, n_pol, n_tht, Mach1_openBC, fix_axis_nodes
    use vacuum, only: is_freebound
    use mpi_mod
    use mod_locate_irn_jcn

    implicit none

    ! Subroutine parameters
    INTEGER                  :: my_id
    INTEGER                  :: local_elms(*)
    INTEGER                  :: n_local_elms
    INTEGER                  :: index_min
    INTEGER                  :: index_max
    INTEGER                  :: xcase2
    TYPE (type_node_list)    :: node_list
    TYPE (type_element_list) :: element_list
    TYPE (type_bnd_node_list):: bnd_node_list
    logical                  :: xpoint2
    REAL*8                   :: R_axis
    REAL*8                   :: Z_axis
    REAL*8                   :: psi_axis
    REAL*8                   :: psi_bnd
    REAL*8                   :: R_xpoint(2)
    REAL*8                   :: Z_xpoint(2)
    REAL*8                   :: psi_xpoint(2)
    logical                  :: gmres
    logical                  :: solve_only
    real*8                   :: rhs_loc(*)

    ! Internal parameters
    real*8  :: zbig
    integer :: i, in, iv, inode, k
    integer :: index_large_i, index_node, index_node2, ielm
    integer :: ijA_position, ijA_position2, ilarge2
    integer :: ilarge_v(n_var), ilarge_vs(n_var)

    integer :: loop_nbr, loop, cnt, cnt_prod
    integer :: ierr
    logical :: is_local, only_count

    real*8  :: R, R_s, R_t, R_mid, R_cnt
    real*8  :: Z, Z_s, Z_t, Z_mid, Z_cnt
    real*8  :: xjac

    real*8  :: T0, T0_corr, T0_s   
    real*8  :: uR0, uR0_s  
    real*8  :: uZ0, uZ0_s  
    real*8  :: up0, up0_s  

    real*8  :: AR0, AR0_s, AR0_t, AR0_p, AR0_R, AR0_Z
    real*8  :: AZ0, AZ0_s, AZ0_t, AZ0_p, AZ0_R, AZ0_Z
    real*8  :: A30, A30_s, A30_t, A30_p, A30_R, A30_Z
    real*8  :: Fprofile, BR0, BZ0, Bp0, BB2, B_dot_n

    real*8  :: normal(2), normal_direction(2)
    real*8  :: grad_s(2), grad_t(2)

    real*8  :: Cs,   Cs_T,   Cs_s,   Cs_s_T,   Cs_s_Ts
    real*8  :: beta, beta_T, beta_s, beta_s_T, beta_s_Ts

    real*8  :: Mach1, Mach1_U, Mach1_T, Mach1_s, Mach1_s_Us, Mach1_s_T, Mach1_s_Ts

    zbig = 1.d12

    do i=1, n_local_elms
      ielm = local_elms(i)

      ! --- Needed for Mach1 BCs
      R_cnt = 0.d0 ; Z_cnt = 0.d0
      do iv=1, n_vertex_max
        inode = element_list%element(ielm)%vertex(iv)
        R_cnt = R_cnt + node_list%node(inode)%x(1,1) / 4.d0     ! center point within element (approx.)
        Z_cnt = Z_cnt + node_list%node(inode)%x(1,2) / 4.d0
      enddo

      do iv=1, n_vertex_max
        inode = element_list%element(ielm)%vertex(iv)

        ! --- Needed for Mach1 BCs
        R_mid = node_list%node(inode)%x(1,1)
        Z_mid = node_list%node(inode)%x(1,2)
        normal_direction = (/R_mid - R_cnt, Z_mid - Z_cnt /) / norm2((/R_mid - R_cnt, Z_mid - Z_cnt /))

        ! A crude way of imposing partial regularity at the grid axis
        !---------------------------------------------------------------------------------------------
        do in=1, n_tor
          do k=1, n_var

            ! Restrain the coefficients of the 3rd basis functions on axis from changing
            if ( ( inode <= n_tht .or. ( n_tht < 1 .and. inode <= n_pol ) ) .and. (fix_axis_nodes) ) then

              index_node = node_list%node(inode)%index(3)
              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                index_large_i = n_tor * n_var * (index_node - 1)
                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                A_glob(ilarge2)    = zbig

              end if
            endif

          enddo
        enddo
        !---------------------------------------------------------------------------------------------


        if (node_list%node(inode)%boundary .ne. 0) then

          do in=1, n_tor
            do k=1, n_var

              !------------------------------------ the open field lines (in case of x-point grid)
              if ((node_list%node(inode)%boundary == 1) .or. (node_list%node(inode)%boundary == 3)) then

                if ((k .eq. var_AR) .or. (k .eq. var_AZ) .or. (k .eq. var_A3)) then

                  index_node = node_list%node(inode)%index(1)
                  if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                    call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                    index_large_i = n_tor * n_var * (index_node - 1)
                    ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                    irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                    jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                    A_glob(ilarge2)   = zbig

                  endif

                  index_node = node_list%node(inode)%index(2)
                  if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                    call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                    index_large_i = n_tor * n_var * (index_node - 1)
                    ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                    irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                    jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                    A_glob(ilarge2)    = zbig

                  endif

                endif

                        ! --- Mach-1 BCs
                if (.not. Mach1_openBC) then 
                  if ( (k == var_uR) .or. (k == var_uZ) .or. (k == var_up) ) then 
                    
                    index_node  = node_list%node(inode)%index(1)             ! position of value
                    index_node2 = node_list%node(inode)%index(2)             ! position of first deriative

                    R         = node_list%node(inode)%x(1,1)
                    R_s       = node_list%node(inode)%x(2,1)
                    R_t       = node_list%node(inode)%x(3,1)
                    Z         = node_list%node(inode)%x(1,2)
                    Z_s       = node_list%node(inode)%x(2,2)
                    Z_t       = node_list%node(inode)%x(3,2)
                    xjac      = R_s*Z_t - R_t*Z_s

                    T0        = node_list%node(inode)%values(1,1,var_T)
                    T0_corr   = max(T0, 1.d-12)
                    T0_s      = node_list%node(inode)%values(1,2,var_T)
                    uR0       = node_list%node(inode)%values(1,1,var_uR)
                    uR0_s     = node_list%node(inode)%values(1,2,var_uR)
                    uZ0       = node_list%node(inode)%values(1,1,var_uZ)
                    uZ0_s     = node_list%node(inode)%values(1,2,var_uZ)
                    up0       = node_list%node(inode)%values(1,1,var_up)
                    up0_s     = node_list%node(inode)%values(1,2,var_up)
                    
                    AR0       = node_list%node(inode)%values(1,1,var_AR)
                    AR0_s     = node_list%node(inode)%values(1,2,var_AR)
                    AR0_t     = node_list%node(inode)%values(1,3,var_AR)
                    AR0_p     = 0.d0
                    AZ0       = node_list%node(inode)%values(1,1,var_AZ)
                    AZ0_s     = node_list%node(inode)%values(1,2,var_AZ)
                    AZ0_t     = node_list%node(inode)%values(1,3,var_AZ)
                    AZ0_p     = 0.d0
                    A30       = node_list%node(inode)%values(1,1,var_A3)
                    A30_s     = node_list%node(inode)%values(1,2,var_A3)
                    A30_t     = node_list%node(inode)%values(1,3,var_A3)
                    A30_p     = 0.d0
                    
                    AR0_R = (   Z_t * AR0_s  - Z_s * AR0_t ) / xjac
                    AR0_Z = ( - R_t * AR0_s  + R_s * AR0_t ) / xjac
                    AZ0_R = (   Z_t * AZ0_s  - Z_s * AZ0_t ) / xjac
                    AZ0_Z = ( - R_t * AZ0_s  + R_s * AZ0_t ) / xjac
                    A30_R = (   Z_t * A30_s  - Z_s * A30_t ) / xjac
                    A30_Z = ( - R_t * A30_s  + R_s * A30_t ) / xjac
                    
                    Fprofile = node_list%node(inode)%Fprof_eq(1)
                    
                    BR0 = ( A30_Z - AZ0_p )/ R
                    BZ0 = ( AR0_p - A30_R )/ R
                    Bp0 = ( AZ0_R - AR0_Z )* R  +  Fprofile
                    BB2 = (BR0*BR0 + BZ0*BZ0 + Bp0*Bp0 / R**2)
                    
                    grad_s = (/   Z_t, - R_t /) / xjac
                    grad_t = (/ - Z_s,   R_s /) / xjac
                    normal = dot_product(grad_t,normal_direction) * grad_t
                    normal = normal / norm2(normal)

                    B_dot_n = (BR0 * normal(1) + BZ0 * normal(2))
                    Cs = (gamma * T0_corr)**0.5
                    beta = Cs * B_dot_n / BB2
                    
                    Cs_T = 0.5 * gamma * (gamma * T0_corr)**(-0.5)
                    beta_T   = Cs_T * B_dot_n / BB2

                    Cs_s = 0.5 * gamma * T0_s * (gamma * T0_corr)**(-0.5)
                    beta_s = Cs_s * B_dot_n / BB2
                    
                    Cs_s_T  = - 0.25 * gamma**2 * T0_s * (gamma * T0_corr)**(-1.5)
                    Cs_s_Ts = 0.5 * gamma * (gamma * T0_corr)**(-0.5)
                    beta_s_T  = Cs_s_T  * B_dot_n / BB2
                    beta_s_Ts = Cs_s_Ts * B_dot_n / BB2

                    if (k == var_uR) Mach1      = uR0 - beta      * BR0
                    if (k == var_uZ) Mach1      = uZ0 - beta      * BZ0
                    if (k == var_up) Mach1      = up0 - beta      * Bp0
                                     Mach1_U    = 1.0
                    if (k == var_uR) Mach1_T    =     - beta_T    * BR0
                    if (k == var_uZ) Mach1_T    =     - beta_T    * BZ0
                    if (k == var_up) Mach1_T    =     - beta_T    * Bp0
                    
                    if (k == var_uR) Mach1_s    = uR0_s - beta_s    * BR0
                    if (k == var_uZ) Mach1_s    = uZ0_s - beta_s    * BZ0
                    if (k == var_up) Mach1_s    = up0_s - beta_s    * Bp0
                    if (k == var_uR) Mach1_s_Us = 1.0
                    if (k == var_uR) Mach1_s_T  =     - beta_s_T  * BR0
                    if (k == var_uZ) Mach1_s_T  =     - beta_s_T  * BZ0
                    if (k == var_up) Mach1_s_T  =     - beta_s_T  * Bp0
                    if (k == var_uR) Mach1_s_Ts =     - beta_s_Ts * BR0
                    if (k == var_uZ) Mach1_s_Ts =     - beta_s_Ts * BZ0
                    if (k == var_up) Mach1_s_Ts =     - beta_s_Ts * Bp0

                    if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                      call locate_irn_jcn(index_node,index_node, index_min,index_max,ijA_position)
                      call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position2)

                      index_large_i = n_tor * n_var * (index_node - 1)
                      
                      ilarge_v (k)     = ijA_position  - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k    -1)*n_tor + in
                      ilarge_v (var_T) = ijA_position  - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (var_T-1)*n_tor + in

                      irn_glob(ilarge_v(k))      = n_tor * n_var * (index_node -1) + (k    -1)*n_tor + in
                      jcn_glob(ilarge_v(k))      = n_tor * n_var * (index_node -1) + (k    -1)*n_tor + in
                      A_glob  (ilarge_v(k))      = Zbig * Mach1_U

                      irn_glob(ilarge_v(var_T))  = n_tor * n_var * (index_node -1) + (k    -1)*n_tor + in
                      jcn_glob(ilarge_v(var_T))  = n_tor * n_var * (index_node -1) + (var_T-1)*n_tor + in
                      A_glob  (ilarge_v(var_T))  = Zbig * Mach1_T

                      RHS_loc(n_tor*n_var * (index_node-1) + (k-1)*n_tor + in) = 0.d0
                      if (in .eq. 1) then
                         RHS_loc(n_tor*n_var * (index_node-1) + (k-1)*n_tor + in) = - Zbig * Mach1
                      endif

                    endif

                    if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then

                      call locate_irn_jcn(index_node2,index_node, index_min,index_max,ijA_position)
                      call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2)

                      index_large_i = n_tor * n_var * (index_node2 - 1)

                      ilarge_vs(k)     = ijA_position2 - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k    -1)*n_tor + in
                      ilarge_v (var_T) = ijA_position  - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (var_T-1)*n_tor + in
                      ilarge_vs(var_T) = ijA_position2 - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (var_T-1)*n_tor + in

                      irn_glob(ilarge_vs(k))     = n_tor * n_var * (index_node2-1) + (k    -1)*n_tor + in
                      jcn_glob(ilarge_vs(k))     = n_tor * n_var * (index_node2-1) + (k    -1)*n_tor + in
                      A_glob  (ilarge_vs(k))     = Zbig * Mach1_s_Us

                      irn_glob(ilarge_v(var_T))  = n_tor * n_var * (index_node2-1) + (k    -1)*n_tor + in
                      jcn_glob(ilarge_v(var_T))  = n_tor * n_var * (index_node -1) + (var_T-1)*n_tor + in
                      A_glob  (ilarge_v(var_T))  = Zbig * Mach1_s_T

                      irn_glob(ilarge_vs(var_T)) = n_tor * n_var * (index_node2-1) + (k    -1)*n_tor + in
                      jcn_glob(ilarge_vs(var_T)) = n_tor * n_var * (index_node2-1) + (var_T-1)*n_tor + in
                      A_glob  (ilarge_vs(var_T)) = Zbig * Mach1_s_Ts

                      Rhs_loc(n_tor * n_var * (index_node2-1) + (k-1)*n_tor + in) = 0.d0
                      if (in .eq. 1) then
                         Rhs_loc(n_tor * n_var * (index_node2-1) + (k-1)*n_tor + in) = - Zbig * Mach1_s
                      endif

                    endif
                  
                  endif
                endif

              endif

              !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
              if ((node_list%node(inode)%boundary == 2) .or. (node_list%node(inode)%boundary == 3)) then

                if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                  index_node = node_list%node(inode)%index(1)
                  if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                    call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                    index_large_i = n_tor * n_var * (index_node - 1)
                    ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                    irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                    jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                    A_glob(ilarge2)   = zbig

                  endif

                  index_node = node_list%node(inode)%index(3)
                  if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                    call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)
                    index_large_i = n_tor * n_var * (index_node - 1)
                    ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                    irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                    jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                    A_glob(ilarge2)    = zbig
                  end if

                endif
              endif

            enddo

          enddo
        endif
      enddo
    enddo

    return
  end subroutine boundary_conditions

end module mod_boundary_conditions
