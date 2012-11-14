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
!*   index_min    - Minimal index of local elements (not with murge assembly)  *
!*   index_max    - Maximal index of local elements (not with murge assembly)  *
!*   xpoint2      -                                                            *
!*   xcase2       -                                                            *
!*   psi_axis     -                                                            *
!*   psi_bnd      -                                                            *
!*   Z_xpoint     -                                                            *
!*   gmres        - boolean indicating if we are using GMRES method            *
!*   solve_only   - Indicate if we want to perform only solve                  *
!*                                                                             *
!* Authors:                                                                    *
!*   Xavier Lacoste - xavier.lacoste@inria.fr                                  *
!*                                                                             *
!*******************************************************************************
subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,    & 
     &                          n_local_elms, index_min, index_max, rhs_loc, xpoint2,   &
     &                          xcase2, psi_axis, psi_bnd, Z_xpoint, psi_xpoint, gmres, solve_only )

  use data_structure
  use global_distributed_matrix
  use phys_module, only: F0, GAMMA, freeboundary
  USE murge_module
  use mpi_mod

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
  REAL*8                   :: psi_axis
  REAL*8                   :: psi_bnd
  REAL*8                   :: Z_xpoint(2)
  REAL*8                   :: psi_xpoint(2)
  logical                  :: gmres
  logical                  :: solve_only
  real*8                   :: rhs_loc(*)

  ! Internal parameters
  real*8  :: zbig,  T0, Vpar0, bigR, dT0_ds, dVpar0_ds, dBigR_ds
  real*8  :: R_s, R_t, Z_s, Z_t, ps0_s, ps0_t, ps0_x, ps0_y, direction, xjac
  real*8  :: Btot
  real*8  :: grad_psi, u0_s, u0_t, u0_x, u0_y
  integer :: i, in, iv, inode, k
  integer :: index_large_i, index_node, index_node2, ielm
  integer :: ijA_position,ijA_position2, ilarge2, kv, kT, ku, ilarge_vv, ilarge_vT, ilarge_vus
  integer :: ilarge_vsvs, ilarge_vsTs, ilarge_vsT
  integer :: loop_nbr, loop, cnt, cnt_prod
  integer :: first_tor, last_tor, murge_ntor, ierr
  logical :: is_local, only_count
  
  zbig = 1.d10
  if (use_murge .and. use_murge_element) then
     ! when we use murge assembly we first count entries then we had them.
     loop_nbr   = 2
     cnt        = 0
     cnt_prod   = 0
     only_count = .true.
  else
     ! No need to do 2 loops when we build irn_glob, jcn_glob, A_glob.
     loop_nbr   = 1
     only_count = .false.
  end if

  do loop = 1, loop_nbr
     if (loop == 2)  then
        only_count = .false.
        write (*,*) my_id, ":: Murge Boundary Assembly phase :: ", cnt, " entries"
        if (.not. solve_only) then
           CALL MURGE_ASSEMBLYBEGIN(murge_id, cnt,                          &
                &                   MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW, &
                &                   MURGE_ASSEMBLY_FOOL, murge_sym, ierr)
        end if
        if (gmres) then
           CALL MURGE_ASSEMBLYBEGIN(murge_id_prod, cnt_prod,                &
                &                   MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW, &
                &                   MURGE_ASSEMBLY_FOOL, murge_sym, ierr)
        end if
     end if

     do i=1, n_local_elms

        ielm = local_elms(i)

        do iv=1, n_vertex_max

           inode = element_list%element(ielm)%vertex(iv)

           if (node_list%node(inode)%boundary .ne. 0) then

              do in=1, n_tor

                 do k=1, n_var

                    !------------------------------------ the open field lines (in case of x-point grid)
                    if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then

                       if ( (k .eq. 1) .or. (k .eq. 2) ) then

                          index_node = node_list%node(inode)%index(1)
                          if (use_murge .and. use_murge_element) then
                             call vertex_is_local(index_node, is_local)
                             if (is_local) then
                                call murge_add_one_entry( index_node, k, in, &
                                     &                    index_node, k, in, &
                                     &                    zbig, solve_only,  &
                                     &                    gmres,             &
                                     &                    cnt, cnt_prod,     &
                                     &                    only_count)

                             end if
                          else
                             if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                                index_large_i = n_tor * n_var * (index_node - 1)

                                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                A_glob(ilarge2)   = zbig

                             endif
                          end if

                          index_node = node_list%node(inode)%index(2)
                          if (use_murge .and. use_murge_element) then
                             call vertex_is_local(index_node, is_local)
                             if (is_local) then
                                call murge_add_one_entry( index_node, k, in,  &
                                     &                    index_node, k, in,  &
                                     &                    zbig, solve_only,   &
                                     &                    gmres,              &
                                     &                    cnt, cnt_prod,      &
                                     &                    only_count)
                             end if
                          else
                             if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                                index_large_i = n_tor * n_var * (index_node - 1)

                                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                A_glob(ilarge2)    = zbig

                             endif
                          end if
                       endif

                    end if

                    !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)

                    if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

                       if ( (k .eq. 1) .or. (k .eq. 2) ) then

                          index_node = node_list%node(inode)%index(1)
                          if (use_murge .and. use_murge_element) then
                             call vertex_is_local(index_node, is_local)
                             if (is_local) then
                                call murge_add_one_entry( index_node, k, in, &
                                     &                    index_node, k, in, &
                                     &                    zbig, solve_only,  &
                                     &                    gmres,             &
                                     &                    cnt, cnt_prod,     &
                                     &                    only_count)
                             end if
                          else
                             if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                                index_large_i = n_tor * n_var * (index_node - 1)

                                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                A_glob(ilarge2)   = zbig

                             endif
                          end if
                          index_node = node_list%node(inode)%index(3)

                          if (use_murge .and. use_murge_element) then
                             call vertex_is_local(index_node, is_local)
                             if (is_local) then
                                call murge_add_one_entry( index_node, k, in, &
                                     &                    index_node, k, in, &
                                     &                    zbig, solve_only,  &
                                     &                    gmres,             &
                                     &                    cnt, cnt_prod,     &
                                     &                    only_count)
                             end if
                          else
                             if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position)

                                index_large_i = n_tor * n_var * (index_node - 1)

                                ilarge2 = ijA_position - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k-1)*n_tor + in

                                irn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                jcn_glob(ilarge2) =  n_tor * n_var * (index_node-1) + (k-1)*n_tor + in
                                A_glob(ilarge2)    = zbig

                             endif
                          end if
                       endif

                    endif

                 enddo

              enddo
           endif
        enddo
     enddo
     if (loop == 2) then
        if (.not. solve_only) then
           CALL MURGE_ASSEMBLYEND(murge_id, ierr)
        end if
        if (gmres) then
           CALL MURGE_ASSEMBLYEND(murge_id_prod, ierr)
        end if
     end if
  end do
  return
end subroutine boundary_conditions
