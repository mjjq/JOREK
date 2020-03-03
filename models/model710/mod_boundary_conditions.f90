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
    use phys_module, only: F0, GAMMA, n_pol, n_tht, linear_run
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
    real*8  :: zbig, zbig_backup
    integer :: i, in, iv, inode, k
    integer :: index_large_i, index_node, ielm
    integer :: ijA_position, ilarge2

    integer :: loop_nbr, loop, cnt, cnt_prod
    integer :: ierr
    logical :: is_local, only_count

    zbig = 1.d12
    zbig_backup = zbig
       do i=1, n_local_elms

          ielm = local_elms(i)

          do iv=1, n_vertex_max

             inode = element_list%element(ielm)%vertex(iv)


             ! A crude way of imposing partial regularity at the grid axis
             ! May be improved upon in the future ( see forthcoming paper on full MHD JOREK )
             !---------------------------------------------------------------------------------------------
             do in=1, n_tor
               if (linear_run  .and.  in .eq. 1 ) then
                 zbig = 1.d15
               else
                 zbig = zbig_backup
               endif
                 do k=1, n_var
                  
                   ! Restrain the coefficients of the 3rd basis functions on axis from changing
                   if ( ( inode <= n_tht .or. ( n_tht < 1 .and. inode <= n_pol ) ) .and. 1==1 ) then
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
                      !
                      ! open field line conditions not yet implemented
                      if ((node_list%node(inode)%boundary == 1) .or. (node_list%node(inode)%boundary == 3)) then

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
