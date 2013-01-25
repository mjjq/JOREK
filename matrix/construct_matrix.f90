!> Construct the main matrix from the contributions of the Bezier elements.
!!
!! The element contributions are determined by element_matrix(_fft). Additional
!! contributions from boundary conditions and the free boundary extension are
!! added by external routine calls.
subroutine construct_matrix(my_id, local_elms, n_local_elms, index_min, index_max, &
                            xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint)
  
  use tr_module 
  use parameters
  use data_structure
  use global_distributed_matrix
  use phys_module
  use pellet_module
  use nodes_elements
  use vacuum_response, only: vacuum_boundary_integral
  use avoid_neg_dens,  only: identify_dens_problems, output_problem_pos, diffusivity_factor
  use mod_ch_nod_rhs_elm
  use mod_boundary_matrix_open
  use mod_elt_matrix
  use mod_elt_matrix_fft
  use mpi_mod
  implicit none
  
#include "r3_info.h"

  ! --- Routine parameters
  integer, intent(in) :: my_id
  integer, intent(in) :: local_elms(*)
  integer, intent(in) :: n_local_elms
  integer, intent(in) :: index_min
  integer, intent(in) :: index_max
  integer, intent(in) :: xcase2
  real*8,  intent(in) :: minRad
  real*8,  intent(in) :: R_axis
  real*8,  intent(in) :: Z_axis
  real*8,  intent(in) :: psi_axis
  real*8,  intent(in) :: psi_bnd
  real*8,  intent(in) :: R_xpoint(2),Z_xpoint(2)
  real*8,  intent(in) :: psi_xpoint(2)
  logical, intent(in) :: xpoint2

  !--- Internal variables
  type (type_element) :: element
  type (type_node)    :: nodes(n_vertex_max)
  type (type_element) :: element_father
  type (type_node)    :: nodes_father(n_vertex_max)
  real*8, pointer     :: rhs_loc(:)
  real*8, dimension (:,:), pointer  :: ELM
  real*8, dimension (:,:), pointer  :: ELM2
  real*8, dimension (:)  , pointer  :: RHS
  real*8, dimension (:)  , pointer  :: RHS2
  integer :: i_bnd, i, ife, iv, iv2, inode, inode1, inode2, knode, j, k, l, index_ij, index_kl
  integer :: index_node1, index_node2, i_order, k_order, ielm, ierr
  integer :: ijA_position, index_min_loc, index_max_loc
  integer :: index_large_i, index_large_k, ilarge2, vertex(2), direction(2)
  integer :: omp_nthreads, omp_tid
  integer :: node_out(n_vertex_max)
  integer :: i_father,INODE_FATHER, ios
  integer, external :: omp_get_num_threads, omp_get_thread_num
  integer :: ilarge_vp, in

  call r3_info_begin (r3_info_index_0, 'construct_matrix')   ! timing

  if (my_id .eq. 0) then
     write(*,*) '****************************************'
     write(*,*) '*  construct matrix                    *'
     write(*,*) '****************************************'
     ! write(*,*) ' n_elements (local)       : ',my_id,n_local_elms
     ! write(*,*) ' index_min,index_max      : ',my_id,index_min,index_max
  endif
  call tr_print_memsize("DebConstM")
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

  if (.not. allocated(A_glob))    call tr_allocate(A_glob,1,nz_glob,"A_glob",CAT_DMATRIX)
  if (.not. allocated(irn_glob))  call tr_allocate(irn_glob,1,nz_glob,"irn_glob",CAT_DMATRIX)
  if (.not. allocated(jcn_glob))  call tr_allocate(jcn_glob,1,nz_glob,"jcn_glob",CAT_DMATRIX)

  if (allocated(rhs_glob))        call tr_deallocate(rhs_glob,"rhs_glob",CAT_DMATRIX)

  call tr_allocate(rhs_glob,1,ndof_glob,"rhs_glob",CAT_DMATRIX)
  call tr_allocatep(rhs_loc,1,ndof_glob,"rhs_loc",CAT_DMATRIX)

  irn_glob = 0
  jcn_glob = 0
  A_glob   = 0.d0
  RHS_glob = 0.d0
  RHS_loc  = 0.d0

#ifdef AVOID_NEG_DENS     
  call identify_dens_problems(my_id)
  if ( my_id == 0 ) call output_problem_pos()
#endif

  !$omp parallel default(none) &
  !$omp   shared(n_local_elms,irn_glob,jcn_glob,A_glob,RHS_loc,local_elms,element_list,node_list,   &
  !$omp          index_min, index_max,xpoint2,xcase2,minRad,R_axis,Z_axis,psi_axis,psi_bnd,Z_xpoint,&
  !$omp          R_xpoint,my_id,bc_natural_open,refinement,thread_struct)                           &
  !$omp   private(ife,ielm,iv,inode,element,nodes,ELM,RHS,ELM2,RHS2,i,inode1,i_order,index_node1,   &
  !$omp           index_large_i,j,index_ij,k,knode,k_order,index_node2,index_large_k,ijA_position,  &
  !$omp           l,index_kl,ilarge2,iv2,vertex,direction,inode2,omp_nthreads,omp_tid,              &
  !$omp           i_father,element_father, nodes_father, inode_father, node_out)

#ifdef _OPENMP
  omp_nthreads = omp_get_num_threads()
  omp_tid      = 1+omp_get_thread_num()
#else
  omp_nthreads = 1
  omp_tid      = 1
#endif
  
  ELM  => thread_struct(omp_tid)%ELM
  RHS  => thread_struct(omp_tid)%RHS
  if (jorek_model .ne. 400) then
    ELM2 => thread_struct(omp_tid)%ELM2
    RHS2 => thread_struct(omp_tid)%RHS2
  endif

  !$omp do 
  do ife =1, n_local_elms
     ielm = local_elms(ife)

     element = element_list%element(ielm)
     
     if (refinement) then
       i_father = element_list%element(ielm)%father

       if( i_father.ne.0) then
  
         element_father = element_list%element(i_father)

         do iv = 1, n_vertex_max
           inode_father=element_father%vertex(iv)
           nodes_father(iv) = node_list%node(inode_father)
         enddo
  
       endif
  
     else
          
       do iv = 1, n_vertex_max

        inode     = element%vertex(iv)
        nodes(iv) = node_list%node(inode)

       enddo

     endif

     if (n_tor .gt. 3) then
       call element_matrix_fft(element,nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, ELM, RHS, omp_tid)      ! use fft for toroidal integration
     else
       call element_matrix(element,nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, ELM, RHS, omp_tid)          ! use direct integration
     endif

     do iv = 1, n_vertex_max									 ! boundary integrals

       iv2  = mod(iv, n_vertex_max) + 1

       inode1 = element%vertex(iv)
       inode2 = element%vertex(iv2)

       if (bc_natural_open) then

         if   (  ((node_list%node(inode1)%boundary .eq. 1) .or.(node_list%node(inode1)%boundary .eq. 3)) &
           .and. ((node_list%node(inode2)%boundary .eq. 1) .or.(node_list%node(inode2)%boundary .eq. 3)) ) then

           nodes(1)  = node_list%node(inode1)
           nodes(2)  = node_list%node(inode2)
           
           vertex    = (/ iv, iv2 /)
           direction = (/  1, 2   /)

!           call boundary_matrix_open(vertex, direction, element,nodes, xpoint2, xcase2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)    ! for open field lines

           call boundary_matrix_open(vertex, direction, element,nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                     psi_bnd, R_xpoint, Z_xpoint, ELM, RHS)    ! for open field lines
          endif
        
       endif
         
     enddo


     !------------------------------------------------------- comparing two versions of element_matrix
      ! if (ife .eq. n_local_elms/2) then
      !   call element_matrix_fft(element,nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, ELM2, RHS2, omp_tid)
      !   call element_matrix(element,nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, ELM, RHS, omp_tid)
      !   do i=1,n_tor*n_vertex_max*(n_order+1)*n_var
      !     if (abs(RHS(i)-RHS2(i))/(abs(RHS(i))+abs(RHS2(i))+1.d0) .gt. 1.d-12) then
      !       write(*,'(i3,A,i6,3e16.8)') my_id,' RHS : ',i,RHS(i),RHS2(i),RHS(i)-RHS2(i)
      !     endif
      !   enddo
      !   do i=1,n_tor*n_vertex_max*(n_order+1)*n_var
      !     do j=1,n_tor*n_vertex_max*(n_order+1)*n_var
      !       if (abs(ELM(i,j)-ELM2(i,j))/(abs(ELM(i,j))+abs(ELM2(i,j))+1.d0) .gt. 1.d-10) then
      !         write(*,'(i3,A,2i6,3e16.8)') my_id,' ELM : ',i,j,ELM(i,j),ELM2(i,j),ELM(i,j)-ELM2(i,j)
      !       endif
      !    enddo
      !   enddo
      ! endif

     if (refinement) then   
       call ch_nod_rhs_elm(ielm,element,nodes,element_father,nodes_father,ELM,RHS,node_out) 
     else
       do i=1, n_vertex_max
         node_out(i) = element%vertex(i)   
       enddo 
     endif

     !$omp critical  

     if ((.not. refinement) .or.(refinement .and. (element%n_sons .eq. 0))) then
     
        do i=1,n_vertex_max

           inode1 = node_out(i)

           do i_order = 1, n_order+1

              index_node1 = node_list%node(inode1)%index(i_order)

              index_large_i = n_tor * n_var * (index_node1 - 1)

              if ((index_node1 .ge. index_min) .and. (index_node1 .le. index_max)) then

                 do j = 1, n_var * n_tor

                    index_ij = n_tor * n_var * (n_order+1) * (i-1) + n_tor * n_var * (i_order-1) + j   ! index in the ELM matrix

                    rhs_loc(index_large_i+j) = rhs_loc(index_large_i+j) + RHS(index_ij)

                    do k=1,n_vertex_max

                       knode = node_out(k)

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

  end do
  !$omp end do
  !$omp end parallel


  ! --- Add vacuum response (boundary integral) for free boundary computations
  if (freeboundary) then
    call global_matrix_structure_vacuum(node_list, bnd_node_list, index_min, index_max) !###TODO### move somewhere else
    call vacuum_boundary_integral(my_id, bnd_node_list, node_list, bnd_elm_list,                   &
      freeboundary_equil, resistive_wall, index_min, index_max, rhs_loc, tstep, index_now)
  end if


! previous WRONG formulation
!!$
!!$  ! --- Form a global rhs from the rhss of the individual mpi threads.
!!$  call MPI_Reduce(RHS_loc,RHS_glob,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!!$  call tr_deallocatep(RHS_loc,"RHS_loc",CAT_DMATRIX)
!!$
!!$  ! --- Apply boundary conditions.
!!$  call boundary_conditions(my_id, node_list, element_list, bnd_node_list, local_elms, n_local_elms, index_min,      &
!!$       index_max, xpoint2, xcase2, psi_axis, psi_bnd, Z_xpoint, psi_xpoint, .false., .false.)
!!$
!!$


!! For debugging purpose
!  call MPI_Reduce(RHS_loc,RHS_glob,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!  call tr_locvnorms("cm_Rhs",RHS_glob,ndof_glob)
!  if (my_id .eq. 0) then
!     write(fname,'(A,I6.6)')"rhs",index_now
!     call tr_vdump(fname,RHS_glob,ndof_glob)
!  end if


! --- Apply boundary conditions.
  call tr_vnorms("cm_A_bef_bc",A_glob,nz_glob)

  ! --- Apply boundary conditions.
  call boundary_conditions(my_id, node_list, element_list,  bnd_node_list,local_elms, n_local_elms,            &
                           index_min, index_max, rhs_loc, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd,  &
                           R_xpoint, Z_xpoint, psi_xpoint, .false., .false.)

  call tr_vnorms("cm_A_aft_bc",A_glob,nz_glob)
 

  ! --- Form a global rhs from the rhss of the individual mpi threads.

  call MPI_Reduce(RHS_loc,RHS_glob,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call tr_deallocatep(RHS_loc,"RHS_loc",CAT_DMATRIX)

!! For debugging purpose
!  if (my_id .eq. 0) then
!     write(fname,'(A,I6.6)')"rhsbc",index_now
!     call tr_vdump(fname,RHS_glob,ndof_glob)
!  end if

  call tr_locvnorms("cm_BCRhs",RHS_glob,ndof_glob)
  call tr_debug_writei("ndof_glob",ndof_glob)
  !write(string, '(A8,I2.2,A1)') "matrice_",my_id,"\0"
  !open(unit=9, file=string, STATUS='replace')
  !do k = 1, nz_glob
  !   if (A_glob(k) /= 0.0_8) &
  !        write(9, '(I8.8,1X,I8.8,1X,E20.12)'), jcn_glob(k), irn_glob(k), A_glob(k)
  !end do
  !close(unit=9)


  call r3_info_end(r3_info_index_0) !timing
  call tr_print_memsize("EndConstM")
  return
end subroutine construct_matrix
