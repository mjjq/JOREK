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
!*******************************************************************************
module mod_boundary_conditions
implicit none
contains
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,          &
                                  n_local_elms, index_min, index_max, rhs_loc, xpoint2, xcase2,       & 
                                  R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,  &  
                                  gmres, solve_only, ijA_index, ijA_size, irn_jcn, irn_glob, jcn_glob,& 
                                  A_glob, i_tor_min, i_tor_max )

    use data_structure
    use phys_module, only: F0, GAMMA
    use vacuum, only: is_freebound
    use mpi_mod
    use mod_locate_irn_jcn

    implicit none

    ! --- Routine parameters
    integer,                   intent(in)    :: my_id
    type (type_node_list),     intent(in)    :: node_list
    type (type_element_list),  intent(in)    :: element_list
    type (type_bnd_node_list), intent(in)    :: bnd_node_list
    integer,                   intent(in)    :: local_elms(*)
    integer,                   intent(in)    :: n_local_elms
    integer,                   intent(in)    :: index_min
    integer,                   intent(in)    :: index_max
    real*8,                    intent(inout) :: rhs_loc(*)
    logical,                   intent(in)    :: xpoint2
    integer,                   intent(in)    :: xcase2
    real*8,                    intent(in)    :: R_axis
    real*8,                    intent(in)    :: Z_axis
    real*8,                    intent(in)    :: psi_axis
    real*8,                    intent(in)    :: psi_bnd
    real*8,                    intent(in)    :: R_xpoint(2)
    real*8,                    intent(in)    :: Z_xpoint(2)
    real*8,                    intent(in)    :: psi_xpoint(2)
    logical,                   intent(in)    :: gmres
    logical,                   intent(in)    :: solve_only
    integer                                  :: irn_glob(:), jcn_glob(:) 
    real*8                                   :: A_glob(:) 
    integer,                   intent(in)    :: i_tor_min, i_tor_max
    integer,          intent(in), pointer    :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:)

    ! Internal parameters
    real*8  :: zbig
    integer :: i, in, iv, inode, k
    integer :: index_large_i, index_node, ielm
    integer :: ijA_position, ilarge2

    integer :: loop_nbr, loop, cnt, cnt_prod
    integer :: ierr
    logical :: is_local, only_count

    zbig = 1.d12

       do i=1, n_local_elms

          ielm = local_elms(i)

          do iv=1, n_vertex_max

             inode = element_list%element(ielm)%vertex(iv)

             if (node_list%node(inode)%boundary .ne. 0) then

                do in=i_tor_min, i_tor_max !1, n_tor

                   do k=1, n_var

                      !------------------------------------ the open field lines (in case of x-point grid)
                      if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then

                         if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
                              (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 6) ) then
 
                          if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                            index_node = node_list%node(inode)%index(1)
                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               index_large_i = (i_tor_max - i_tor_min + 1) * n_var * (index_node - 1)

                               ilarge2 = ijA_position - 1 + ((k-1)*(i_tor_max - i_tor_min + 1) + in-i_tor_min) * n_var*(i_tor_max - i_tor_min + 1)  & 
                                                          +  (k-1)*(i_tor_max - i_tor_min + 1) + in - i_tor_min + 1

                               irn_glob(ilarge2) = (i_tor_max - i_tor_min + 1) * n_var * (index_node-1) + (k-1)*(i_tor_max - i_tor_min + 1) + in - i_tor_min + 1
                               jcn_glob(ilarge2) = (i_tor_max - i_tor_min + 1) * n_var * (index_node-1) + (k-1)*(i_tor_max - i_tor_min + 1) + in - i_tor_min + 1
                               A_glob(ilarge2)   = zbig

                            endif

                            index_node = node_list%node(inode)%index(2)

                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               index_large_i = (i_tor_max - i_tor_min + 1) * n_var * (index_node - 1)

                               ilarge2 = ijA_position - 1 + ((k-1)*(i_tor_max - i_tor_min + 1) + in-i_tor_min) * n_var*(i_tor_max - i_tor_min + 1)   & 
                                                          +  (k-1)*(i_tor_max - i_tor_min + 1) + in - i_tor_min + 1

                               irn_glob(ilarge2) = (i_tor_max - i_tor_min + 1) * n_var * (index_node-1) + (k-1)*(i_tor_max - i_tor_min + 1) + in - i_tor_min + 1
                               jcn_glob(ilarge2) = (i_tor_max - i_tor_min + 1) * n_var * (index_node-1) + (k-1)*(i_tor_max - i_tor_min + 1) + in - i_tor_min + 1
                               A_glob(ilarge2)    = zbig

                            endif

                          endif
                        endif
                      endif

                      !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                      if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

                         if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                            index_node = node_list%node(inode)%index(1)

                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               index_large_i = (i_tor_max - i_tor_min + 1) * n_var * (index_node - 1)

                               ilarge2 = ijA_position - 1 + ((k-1)*(i_tor_max - i_tor_min + 1) + in-i_tor_min) * n_var*(i_tor_max - i_tor_min + 1)   & 
                                                          +  (k-1)*(i_tor_max - i_tor_min + 1) + in - i_tor_min + 1

                               irn_glob(ilarge2) = (i_tor_max - i_tor_min + 1) * n_var * (index_node-1) + (k-1)*(i_tor_max - i_tor_min + 1) + in - i_tor_min + 1
                               jcn_glob(ilarge2) = (i_tor_max - i_tor_min + 1) * n_var * (index_node-1) + (k-1)*(i_tor_max - i_tor_min + 1) + in - i_tor_min + 1
                               A_glob(ilarge2)   = zbig

                            endif

                            index_node = node_list%node(inode)%index(3)

                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               index_large_i = (i_tor_max - i_tor_min + 1) * n_var * (index_node - 1)

                               ilarge2 = ijA_position - 1 + ((k-1)*(i_tor_max - i_tor_min + 1) + in-i_tor_min) * n_var*(i_tor_max - i_tor_min + 1)   & 
                                                          +  (k-1)*(i_tor_max - i_tor_min + 1) + in - i_tor_min + 1

                               irn_glob(ilarge2) = (i_tor_max - i_tor_min + 1) * n_var * (index_node-1) + (k-1)*(i_tor_max - i_tor_min + 1) + in - i_tor_min + 1
                               jcn_glob(ilarge2) = (i_tor_max - i_tor_min + 1) * n_var * (index_node-1) + (k-1)*(i_tor_max - i_tor_min + 1) + in - i_tor_min + 1
                               A_glob(ilarge2)   = zbig
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
