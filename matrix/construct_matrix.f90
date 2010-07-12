subroutine construct_matrix(my_id,local_elms,n_local_elms,index_min,index_max, &
                            xpoint2,psi_axis,psi_bnd,Z_xpoint)
!---------------------------------------------------------------
! collect the element matrices into one large sparse matrix
! in coordinate format
!---------------------------------------------------------------
use parameters
use data_structure
use global_distributed_matrix
use phys_module
use nodes_elements

implicit none
include 'mpif.h'
#include "r3_info.h"


type (type_element)       :: element
type (type_node)          :: nodes(n_vertex_max)
type (type_element)       :: element_father
type (type_node)          :: nodes_father(n_vertex_max)
real*8, allocatable :: rhs_loc(:)
integer :: my_id, local_elms(*), n_local_elms, index_min, index_max,index_min_loc, index_max_loc
real*8  :: ELM(n_tor*n_vertex_max*(n_order+1)*n_var,n_tor*n_vertex_max*(n_order+1)*n_var)
real*8  :: RHS(n_tor*n_vertex_max*(n_order+1)*n_var)
real*8  :: ELM2(n_tor*n_vertex_max*(n_order+1)*n_var,n_tor*n_vertex_max*(n_order+1)*n_var)
real*8  :: RHS2(n_tor*n_vertex_max*(n_order+1)*n_var)
real*8  :: zbig, psi_axis, psi_bnd, Z_xpoint
integer :: i_bnd, i, in, ife, iv, iv2, inode, inode1, inode2, knode, j, k, l, index_ij, index_kl
integer :: index_i, index_node, index_node1, index_node2, i_order, k_order, ic, ielm, ierr
integer :: ijA_position, ijA_position2, nz_AA2, n_AA2, kv, kT, ku
integer :: index_large_i, index_large_k, ilarge2, vertex(2), direction(2)
integer :: omp_nthreads, omp_tid
logical :: xpoint2
integer, dimension(n_vertex_max) ::  node_out
integer, dimension(node_list%n_nodes) :: active_node
integer :: i_father,INODE_FATHER
integer, external :: omp_get_num_threads, omp_get_thread_num

call r3_info_begin (r3_info_index_0, 'construct_matrix')   ! timing
if (my_id .eq. 0) then
  write(*,*) '****************************************'
  write(*,*) '*  construct matrix                    *'
  write(*,*) '****************************************'
! write(*,*) ' n_elements (local)       : ',my_id,n_local_elms
! write(*,*) ' index_min,index_max      : ',my_id,index_min,index_max
endif

i_bnd = 0

do i=1, n_local_elms

  ielm = local_elms(i)

  do iv=1,n_vertex_max

    inode = element_list%element(ielm)%vertex(iv)

    if (node_list%node(inode)%boundary .eq. 1) i_bnd = i_bnd + 1
    if (node_list%node(inode)%boundary .eq. 2) i_bnd = i_bnd + 1
    if (node_list%node(inode)%boundary .eq. 3) i_bnd = i_bnd + 2
    if (i == 1 .and. iv == 1) then
       index_min_loc = minval(node_list%node(iv)%index)
       index_max_loc = maxval(node_list%node(iv)%index)
    else
       index_min_loc = min(index_min_loc, minval(node_list%node(iv)%index))
       index_max_loc = max(index_max_loc, maxval(node_list%node(iv)%index))
    end if
  enddo

enddo

if (.not. allocated(A_glob))    allocate(A_glob(nz_glob))
if (.not. allocated(irn_glob))  allocate(irn_glob(nz_glob))
if (.not. allocated(jcn_glob))  allocate(jcn_glob(nz_glob))

if (allocated(rhs_glob))        deallocate(rhs_glob)

allocate(rhs_glob(ndof_glob))

if (.not. allocated(rhs_loc))   allocate(rhs_loc(ndof_glob))

irn_glob = 0
jcn_glob = 0
A_glob   = 0.d0
RHS_glob = 0.d0
RHS_loc  = 0.d0




!$omp parallel default(none) &
!$omp   shared(n_local_elms,irn_glob,jcn_glob,A_glob,RHS_loc,local_elms,element_list,node_list,   &
!$omp          index_min, index_max,xpoint2,psi_axis,psi_bnd,Z_xpoint, my_id)                     &
!$omp   private(ife,ielm,iv,inode,element,nodes,ELM,RHS,ELM2,RHS2,i,inode1,i_order,index_node1,   &
!$omp           index_large_i,j,index_ij,k,knode,k_order,index_node2,index_large_k,ijA_position,  &
!$omp           l,index_kl,ilarge2,iv2,vertex,direction,inode2,omp_nthreads,omp_tid,              &
!$omp           i_father,element_father, nodes_father, inode_father,node_out)

!omp_nthreads = omp_get_num_threads()
!omp_tid      = omp_get_thread_num()
!write(*,*) my_id,' number of threads, my_tid : ',omp_nthreads,omp_tid


!$omp do
do ife =1, n_local_elms

  ielm = local_elms(ife)


  
  element = element_list%element(ielm)

  i_father= element_list%element(ielm)%father
 

   if( i_father.ne.0) then
    element_father = element_list%element(i_father)
   endif

  do iv = 1, n_vertex_max

    if( i_father.ne.0) then
      inode_father=element_father%vertex(iv)
      nodes_father(iv) = node_list%node(inode_father)
    endif

    inode     = element%vertex(iv)

    nodes(iv) = node_list%node(inode)

  enddo

  if (n_tor .gt. 3) then

    call element_matrix_fft(element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)      ! use fft for toroidal integration

  else

    call element_matrix(element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)          ! use direct integration

    do iv = 1, n_vertex_max                                                                     ! boundary integrals

      iv2  = mod(iv, n_vertex_max) + 1

      inode1 = element%vertex(iv)
      inode2 = element%vertex(iv2)

!      if (     ((node_list%node(inode1)%boundary .eq. 1) .or.(node_list%node(inode1)%boundary .eq. 3)) &
!         .and. ((node_list%node(inode2)%boundary .eq. 1) .or.(node_list%node(inode2)%boundary .eq. 3)) ) then

!	nodes(1)  = node_list%node(inode1)
!	nodes(2)  = node_list%node(inode2)
!	vertex    = (/ iv, iv2 /)
!        direction = (/  1, 2   /)

!	write(*,*) iv,iv2,'boundary_matrix_open : ',inode1,inode2

!        call boundary_matrix_open(vertex, direction, element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)    ! for open field lines

!      endif

      if (     ((node_list%node(inode1)%boundary .eq. 2) .or.(node_list%node(inode1)%boundary .eq. 3)) &
         .and. ((node_list%node(inode2)%boundary .eq. 2) .or.(node_list%node(inode2)%boundary .eq. 3)) ) then

	nodes(1)  = node_list%node(inode1)
	nodes(2)  = node_list%node(inode2)
	vertex    = (/ iv, iv2 /)
        direction = (/ 1, 3    /)

!	write(*,*) iv,iv2,'boundary_matrix : ',inode1,inode2

!        call boundary_matrix(vertex, direction, element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)    ! for closed field lines

      endif

    enddo

  endif

!------------------------------------------------------- comparing two versions of element_matrix
!  if (ife .eq. n_local_elms/2) then
!    call element_matrix_fft(element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM2, RHS2)
!    call element_matrix(element,nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)
!    do i=1,n_tor*n_vertex_max*(n_order+1)*n_var
!      if (abs(RHS(i)-RHS2(i))/(abs(RHS(i))+abs(RHS2(i))+1.d0) .gt. 1.d-12) then
!        write(*,'(i3,A,i6,3e16.8)') my_id,' RHS : ',i,RHS(i),RHS2(i),RHS(i)-RHS2(i)
!      endif
!    enddo
!    do i=1,n_tor*n_vertex_max*(n_order+1)*n_var
!      do j=1,n_tor*n_vertex_max*(n_order+1)*n_var
!        if (abs(ELM(i,j)-ELM2(i,j))/(abs(ELM(i,j))+abs(ELM2(i,j))+1.d0) .gt. 1.d-10) then
!          write(*,'(i3,A,2i6,3e16.8)') my_id,' ELM : ',i,j,ELM(i,j),ELM2(i,j),ELM(i,j)-ELM2(i,j)
!        endif
!     enddo
!    enddo
!  endif

     
   call ch_nod_rhs_elm(ielm,element,nodes,element_father,nodes_father,ELM,RHS,node_out) 

!$omp critical  

if (element%n_sons .eq. 0) then
  do i=1,n_vertex_max

    inode1         =node_out(i)! element%vertex(i)

    do i_order = 1, n_order+1

      index_node1 = node_list%node(inode1)%index(i_order)

      index_large_i = n_tor * n_var * (index_node1 - 1)

      if ((index_node1 .ge. index_min) .and. (index_node1 .le. index_max)) then

        do j = 1, n_var * n_tor

          index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j   ! index in the ELM matrix
       
          !if (index_large_i + j .eq. 2) then
	  	!write(*,*)  "RHS(2) = ",ife,inode1,j,index_ij, rhs_loc(index_large_i+j), "=>", rhs_loc(index_large_i+j) + RHS(index_ij)
         !endif
          rhs_loc(index_large_i+j) = rhs_loc(index_large_i+j) + RHS(index_ij)

          do k=1,n_vertex_max

            knode         =node_out(k)! element%vertex(k)

            do k_order = 1, n_order+1

              index_node2 = node_list%node(knode)%index(k_order)

              index_large_k = n_tor * n_var * (index_node2 - 1)

              call locate_irn_jcn(index_node1,index_node2,index_min,index_max,ijA_position)

              do l = 1, n_var * n_tor

                index_kl = n_tor * n_var * (n_order+1) * (k-1) + n_tor * n_var * (k_order-1) + l   ! index in the ELM matrix
               
                ilarge2 = ijA_position - 1 + (j-1) * n_var*n_tor + l

                irn_glob(ilarge2) = index_large_i   + j
                jcn_glob(ilarge2) = index_large_k   + l
                A_glob(ilarge2)   = A_glob(ilarge2) + ELM(index_ij,index_kl)

              enddo

            enddo

          enddo
        enddo

      endif

    enddo
   
  enddo

endif

!$omp end critical

enddo

!$omp end do
!$omp end parallel

call MPI_Reduce(RHS_loc,RHS_glob,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

deallocate(RHS_loc)

!------------------vacuum cannot be called on an element by element basis. The vacuum couples all boundary points!

if (freeboundary) then
  call vacuum(my_id,node_list,element_list,boundary_list,index_min,index_max,xpoint2,psi_axis,psi_bnd,Z_xpoint)
endif

call boundary_conditions(my_id,node_list,element_list,local_elms,n_local_elms,index_min,index_max, &
                         xpoint2,psi_axis,psi_bnd,Z_xpoint)

call r3_info_end (r3_info_index_0)    ! timing
return
end
