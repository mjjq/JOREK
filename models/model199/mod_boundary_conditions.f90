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
module mod_boundary_conditions
implicit none
contains
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,         &
       n_local_elms, index_min, index_max, rhs_loc, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, &
       R_xpoint, Z_xpoint, psi_xpoint, gmres, solve_only )

    use data_structure
    use global_distributed_matrix
    use phys_module, only: F0, GAMMA
    use vacuum, only: is_freebound
    USE murge_module, ONLY : MURGE_ASSEMBLYBEGIN => MURGE_ASSEMBLYBEGIN_WRAPPER,     &
         use_murge, use_murge_element, murge_id, murge_global_n, MURGE_ASSEMBLY_OVW, &
         MURGE_ASSEMBLY_FOOL, murge_sym, murge_id_prod, murge_global_n_prod,         &
         MURGE_SUCCESS, murge_add_one_entry
    use murge_module, only : MURGE_ASSEMBLYEND, vertex_is_local
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

    ! Internal parameters
    real*8  :: zbig
    integer :: i, in, iv, inode, k
    integer :: index_large_i, index_node, ielm
    integer :: ijA_position, ilarge2

    integer :: loop_nbr, loop, cnt, cnt_prod
    integer :: ierr
    logical :: is_local, only_count

    zbig = 1.d12
    if (use_murge .and. use_murge_element) then
       ! when we use murge assembly we first count entries then we had them.
       loop_nbr = 2
       only_count = .true.
       cnt      = 0
       cnt_prod = 0
    else
       ! No need to do 2 loops when we build irn_glob, jcn_glob, A_glob.
       loop_nbr = 1
       only_count = .false.
    end if


    do loop = 1, loop_nbr
       if (loop == 2)  then
          only_count = .false.
#ifdef USE_MURGE
          if (.not. solve_only) then
             write (*,*) my_id, ":: Murge Boundary Assembly phase :: ", cnt, " entries"
             CALL MURGE_ASSEMBLYBEGIN(murge_id, murge_global_n, cnt, &
                  MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW,            &
                  MURGE_ASSEMBLY_FOOL, murge_sym, ierr)
          end if
          if (gmres) then
             write (*,*) &
                  my_id, ":: Murge Boundary Assembly phase :: ", &
                  cnt_prod, " product entries"
             CALL MURGE_ASSEMBLYBEGIN(murge_id_prod, murge_global_n_prod, cnt_prod, &
                  MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW,                           &
                  MURGE_ASSEMBLY_FOOL, murge_sym, ierr)
          end if
          cnt      = 0
          cnt_prod = 0
#else
          print *, "Binary built without murge"
          call abort()
#endif
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

                         if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
                              (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 6) ) then
 
                          if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                            index_node = node_list%node(inode)%index(1)
                            if (use_murge .and. use_murge_element) then
                               call vertex_is_local(index_node, is_local)
                               if (is_local) then
                                  call murge_add_one_entry( & 
                                       & index_node, k, in, &
                                       & index_node, k, in, &
                                       & zbig,              &
                                       & solve_only, gmres, &
                                       & cnt, cnt_prod, only_count)
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
                                  call murge_add_one_entry( & 
                                       & index_node, k, in, &
                                       & index_node, k, in, &
                                       & zbig,              &
                                       & solve_only, gmres, &
                                       & cnt, cnt_prod, only_count)
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
                      endif

                      !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                      if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

                         if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                            index_node = node_list%node(inode)%index(1)

                            if (use_murge .and. use_murge_element) then
                               call vertex_is_local(index_node, is_local)
                               if (is_local) then
                                  call murge_add_one_entry( & 
                                       & index_node, k, in, &
                                       & index_node, k, in, &
                                       & zbig,              &
                                       & solve_only, gmres, &
                                       & cnt, cnt_prod, only_count)
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
                                  call murge_add_one_entry( & 
                                       & index_node, k, in, &
                                       & index_node, k, in, &
                                       & zbig,              &
                                       & solve_only, gmres, &
                                       & cnt, cnt_prod, only_count)
                               end if
                            else
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

                      endif

                   enddo

                enddo
             endif
          enddo
       enddo
       if (loop == 2) then
#ifdef USE_MURGE
          if (.not. solve_only) then
             CALL MURGE_ASSEMBLYEND(murge_id, ierr)
             IF (ierr /= MURGE_SUCCESS) THEN
                WRITE (*,*)  "ERROR in MURGE_ASSEMBLYEND "
                STOP
             END IF
          end if
          if (gmres) then
             CALL MURGE_ASSEMBLYEND(murge_id_prod, ierr)
             IF (ierr /= MURGE_SUCCESS) THEN
                WRITE (*,*)  "ERROR in MURGE_ASSEMBLYEND "
                STOP
             END IF
          end if
#else
          print *, "Binary built without murge"
          call abort()
#endif
       end if
    end do

    return
  end subroutine boundary_conditions
end module mod_boundary_conditions
