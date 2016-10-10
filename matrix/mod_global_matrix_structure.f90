module mod_global_matrix_structure
contains
subroutine global_matrix_structure(my_id,node_List,element_list,boundary_list,freeboundary,local_elms,n_local_elms,index_min,index_max)
  !***********************************************************************
  !* subroutine determines the position of the indices in the global     *
  !* matrix                                                              *
  !***********************************************************************
  use tr_module
  use data_structure
  use global_distributed_matrix
  use mod_ch_node_struct

  implicit none

  type (type_node_list)        :: node_list
  type (type_element_list)     :: element_list
  type (type_bnd_element_list) :: boundary_list
  type (type_surface_list)     :: flux_list
  type (type_element)          :: element
  type (type_node)             :: nodes(n_vertex_max)

  integer :: local_elms(*), index_min, index_max, my_id, n_local_elms
  integer :: i, ibnd, jbnd, idir, jdir, iv, ik, jv, jk, ielm, inode1, inode2, index1, index2, index1_local, index2_local
  integer :: j_larger, j, ibase, n_max
  integer :: inode,i_father
  integer, dimension(n_vertex_max) ::  node_out
  logical :: freeboundary

  if ( my_id == 0 ) then
    write(*,*) '**********************************'
    write(*,*) '* global_matrix_structure        *'
    write(*,*) '**********************************'
    if (freeboundary) write(*,*) ' FREEBOUNDARY is ON'
  end if
  
  n_matrix_block_size = n_tor * n_var

  ndof_glob = -1
  do inode1=1,node_list%n_nodes
     ndof_glob = max(ndof_glob,maxval(node_list%node(inode1)%index))
  enddo
  ndof_glob = ndof_glob * n_tor*n_var

  n_max = 8192

  call tr_allocate(ijA_size,1,index_max-index_min+1,"ijA_size",CAT_DMATRIX)
  call tr_allocate(irn_jcn,1,index_max-index_min+1,1,n_max,"irn_jcn",CAT_DMATRIX)

  ijA_size    = 0
  irn_jcn = 0

  do i=1,n_local_elms                 ! loop over the local elements

     ielm = local_elms(i)
     element = element_list%element(ielm)
     i_father= element_list%element(ielm)%father
     do iv = 1, n_vertex_max

        inode     = element%vertex(iv)

        nodes(iv) = node_list%node(inode)

     enddo


     !if( i_father.ne.0) then
     call Ch_node_struct(ielm, element,nodes,node_out) !  Processing  "constrained nodes"
     !else
     !  do j=1,4
     !  node_out(j)=element%vertex(j)
     ! enddo
     !endif

     if (element%n_sons .eq. 0) then

        do iv = 1, n_vertex_max                                           ! loop over the vertices

           inode1 =node_out(iv)! element_list%element(ielm)%vertex(iv)

           do ik = 1, n_degrees                                            ! loop over degrees of freedom

              index1       = node_list%node(inode1)%index(ik)
              index1_local = index1 - index_min + 1

              if ((index1 .ge. index_min) .and. (index1 .le. index_max)) then     ! keep contribution only if index belongs to local range of indices
                 ! distribution is by nodes (not by element)
                 do jv = 1,n_vertex_max

                    inode2 =node_out(jv)! element_list%element(ielm)%vertex(jv)

                    do jk = 1, n_degrees

                       index2       = node_list%node(inode2)%index(jk)
                       index2_local = index2 - index_min + 1

                       if (ijA_size(index1_local) .eq. 0) then                      ! if row index1_local still empty, fill with first value (index2)

                          ijA_size(index1_local) = 1
                          irn_jcn(index1_local,1) = index2

                       elseif (index2 .gt. irn_jcn(index1_local,ijA_size(index1_local))) then   ! if index2 larger than all previous indices, add at end of list

                          irn_jcn(index1_local,ijA_size(index1_local)+1) = index2
                          ijA_size(index1_local) = ijA_size(index1_local) + 1

                       else                                                         ! index2 falls somewhere in (between) existing values

                          do j = 1, ijA_size(index1_local)                           ! find the first index larger than index2

                             if (index2 .le. irn_jcn(index1_local,j) ) then

                                j_larger = j                                           ! j_larger is the position of the index larger than (or equal to)  index2
                                exit

                             endif

                          enddo

                          if (index2 .ne. irn_jcn(index1_local,j_larger) ) then      ! if index2 <> index(j_larger) add index2 and shift the following values

                             do j=ijA_size(index1_local), j_larger, -1                ! shift the higher indeices

                                irn_jcn(index1_local,j+1) = irn_jcn(index1_local,j)

                             enddo

                             irn_jcn(index1_local,j_larger) = index2                  ! fill the freed position with index2
                             ijA_size(index1_local) = ijA_size(index1_local) + 1      ! add one to the total number of contributions of this row (index1_local)

                             if (ijA_size(index1_local) .gt. n_max) then
                                write(*,*) ' FATAL error : irn_jcn too small ',ijA_size(index1_local)
                             endif

                          endif ! (index2 .ne. irn_jcn(index1_local,j_larger) )

                       endif !(ijA_size(index1_local) .eq. 0)

                    enddo !jk = 1, n_degrees
                 enddo !jv = 1,n_vertex_max
                 !endif	


              endif           ! check if node is local

           enddo             ! loop over degrees of freedom
        enddo               ! loop over vertices
     endif ! nsons eq 0
  enddo                 ! loop over local elements


  if (freeboundary) then      ! add contributions from all boundary nodes

  do ibnd = 1,boundary_list%n_bnd_elements                                 ! loop over the boundary elements

        do iv = 1, 2                                                           ! loop over the vertices

          inode1 = boundary_list%bnd_element(ibnd)%vertex(iv)

           do ik = 1, 2                                                         ! loop over degrees of freedom

              idir = boundary_list%bnd_element(ibnd)%direction(iv,ik)

              index1       = node_list%node(inode1)%index(idir)
              index1_local = index1 - index_min + 1


              if ((index1 .ge. index_min) .and. (index1 .le. index_max)) then     ! keep contribution only if index belongs to local range of indices

              do jbnd = 1,boundary_list%n_bnd_elements                          ! loop over all boundary elements

                    do jv = 1, 2                                                    ! loop over the nodes of the bounsdary element

                       inode2 = boundary_list%bnd_element(jbnd)%vertex(jv)

                       do jk = 1, 2

                          jdir = boundary_list%bnd_element(jbnd)%direction(jv,jk)

                          index2       = node_list%node(inode2)%index(jdir)
                          index2_local = index2 - index_min + 1

                          if (ijA_size(index1_local) .eq. 0) then

                             ijA_size(index1_local) = 1
                             irn_jcn(index1_local,1) = index2

                          elseif (index2 .gt. irn_jcn(index1_local,ijA_size(index1_local))) then

                             irn_jcn(index1_local,ijA_size(index1_local)+1) = index2
                             ijA_size(index1_local) = ijA_size(index1_local) + 1

                          else

                             do j = 1, ijA_size(index1_local)

                                if (index2 .le. irn_jcn(index1_local,j) ) then

                                   j_larger = j
                                   exit

                                endif

                             enddo

                             if (index2 .ne. irn_jcn(index1_local,j_larger) ) then

                                do j=ijA_size(index1_local), j_larger, -1

                                   irn_jcn(index1_local,j+1) = irn_jcn(index1_local,j)

                                enddo

                                irn_jcn(index1_local,j_larger) = index2
                                ijA_size(index1_local) = ijA_size(index1_local) + 1

                                if (ijA_size(index1_local) .gt. n_max) then
                                   write(*,*) ' FATAL error : irn_jcn too small ',ijA_size(index1_local)
                                endif

                             endif

                          endif

                       enddo  ! end loop over degrees of freedom (jk)
                    enddo    ! end loop over vertices (jv)
                 enddo      ! end loop over boundary elements (jbnd)

              endif        ! endif check if local index
              !endif Ce endif est en plus, pourquoi ?
           enddo          ! end loop over degrees of freedom (ik)
        enddo            ! end loop over vertices (iv)
     enddo              ! end loop over boundary elements (ibnd)
  endif                ! check if free boundary on

  call tr_allocate(ijA_index,1,index_max-index_min+1,1,n_max,"ijA_index",CAT_DMATRIX)

  ibase = 0
  do i=1,index_max-index_min+1

     do j=1,ijA_size(i)

        ijA_index(i,j) = ibase + 1

        ibase = ibase + (n_tor*n_var)**2

     enddo

  enddo

  n_glob  = (index_max-index_min+1) * n_tor * n_var

  nz_glob = ijA_index(index_max-index_min+1,ijA_size(index_max-index_min+1)) + (n_tor*n_var)**2 - 1

  write(*,*) my_id,' size matrices : n, nz = ',n_glob, nz_glob

  return
end subroutine global_matrix_structure
end module mod_global_matrix_structure
