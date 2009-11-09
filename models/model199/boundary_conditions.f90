!*
!* Important: If there are modifications on boundary condition, 
!*            do not forget to modify boundary_conditions_murge.f90.
!*
subroutine boundary_conditions(my_id,node_list,element_list,local_elms,n_local_elms,index_min,index_max, &
                               xpoint2,psi_axis,psi_bnd,Z_xpoint)			    
!---------------------------------------------------------------
! add the boundary condition to the global matrix
!---------------------------------------------------------------
use data_structure
use global_distributed_matrix
use phys_module

implicit none
include 'mpif.h'

! Subroutine parameters
integer                  :: my_id
type (type_node_list)    :: node_list
type (type_element_list) :: element_list
integer                  :: local_elms(*)
integer                  :: n_local_elms
integer                  :: index_min
integer                  :: index_max
logical                  :: xpoint2
real*8                   :: psi_axis
real*8                   :: psi_bnd
real*8                   :: Z_xpoint

! Internal parameters
type (type_element)      :: element
type (type_node_list)    :: nodes
real*8  :: zbig
integer :: i, in, ife, iv, inode, j, k, l
integer :: index_i, index_large_i, index_node, ielm
integer :: ijA_position, ilarge2


zbig = 1.d12

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

          endif

!------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
          if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

            if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
                (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 6) ) then

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

              endif

            endif

          endif

        enddo

      enddo
    endif
  enddo
enddo

return
end
