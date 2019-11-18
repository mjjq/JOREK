module construct_harmonic_matrix_mod

use mod_parameters, only : n_var
implicit none
logical  :: difference_found, rhs_problem(n_var), elm_problem(n_var,n_var)

contains

 subroutine construct_harmonic_matrix(my_id, my_id_n, n_cpu, m_cpu, local_elms, n_local_elms, index_min, index_max,& 
                                      xpoint2, xcase2, R_axis, Z_axis,& 
                                      psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,& 
                                      i_tor_min, i_tor_max, n_glob_harm, nz_glob_harm)

  use construct_matrix_mod!, only: elementary_matrix_build_new 
  use mumps_module 
  use tr_module 
!  use mod_parameters
  use data_structure
  use harmonic_distributed_matrix, only: ijA_size_harm, ijA_index_harm, irn_jcn_harm
  use global_distributed_matrix, only: ndof_glob
  use phys_module
  use pellet_module
  use nodes_elements
  use vacuum, only: sr
  use vacuum_response, only: vacuum_boundary_integral
  use mod_ch_nod_rhs_elm
  use mod_boundary_matrix_open
  use mod_elt_matrix
  use mod_elt_matrix_fft
  use mpi_mod
  use mod_boundary_conditions, only : boundary_conditions
  use mod_locate_irn_jcn
  !$ use omp_lib
  implicit none
  
#include "r3_info.h"

  ! --- Routine parameters
  integer, intent(in) :: my_id
  integer, intent(in) :: my_id_n         
  integer, intent(in) :: n_cpu         
  integer, intent(in) :: m_cpu         
  integer, intent(in) :: local_elms(*)
  integer, intent(in) :: n_local_elms
  integer, intent(in) :: index_min
  integer, intent(in) :: index_max
  integer, intent(in) :: n_glob_harm  
  integer, intent(in) :: nz_glob_harm
!  integer (8), intent(in) :: nz_total
  integer, intent(in) :: i_tor_min
  integer, intent(in) :: i_tor_max
  integer, intent(in) :: xcase2
  real*8,  intent(in) :: R_axis
  real*8,  intent(in) :: Z_axis
  real*8,  intent(in) :: psi_axis
  real*8,  intent(in) :: psi_bnd
  real*8,  intent(in) :: R_xpoint(2)
  real*8,  intent(in) :: Z_xpoint(2)
  real*8,  intent(in) :: psi_xpoint(2)
  logical, intent(in) :: xpoint2

  !--- Internal variables
!  type (type_node_list)             :: node_list
!  type (type_element_list)          :: element_list
  type (type_element)               :: element
  type (type_node)                  :: nodes(n_vertex_max)
  type (type_element)               :: element_father
  type (type_node)                  :: nodes_father(n_vertex_max)
  real*8,                  pointer  :: rhs_loc(:)
  integer                           :: i_bnd, i, ife, iv, iv2, inode, inode1, inode2, knode, j, k, l, index_ij, index_kl
  integer                           :: index_node1, index_node2, i_order, k_order, ielm, ierr
  integer                           :: ijA_position, index_min_loc, index_max_loc
  integer                           :: index_large_i, index_large_k, ilarge2, vertex(2), direction(2)
  integer                           :: omp_nthreads, omp_tid
  integer                           :: node_out(n_vertex_max)
  integer                           :: i_father,INODE_FATHER, ios
  integer                           :: ilarge_vp, in, ivertex, iorder, ivar, itor, jvertex, jorder, jvar, jtor
  integer                           :: random_element, n_var_reduced, v1, v2, im, index_ij_model400_e, index_kl_model400_e
  real*8                            :: tmp_rhs, tmp_elm, tmp_elm_v2_8
  real*8                            :: a_22, a_23, a_32, a_33 !**psv 
  integer                           :: cnt                    !**psv 
!  integer                           :: i_tor_min                    !**psv 
!  integer                           :: i_tor_max                    !**psv 
!  logical                           :: difference_found, rhs_problem(n_var), elm_problem(n_var,n_var)
  CHARACTER(LEN=128)                :: fname
!  integer                           :: index_min
!  integer                           :: index_max
!  integer                           :: n_local_elms
 
!  i_tor_min = 1
!  i_tor_max = n_tor
    
 
  ! --- Timing call
  call r3_info_begin (r3_info_index_0, 'construct_harmonic_matrix')

  ! --- Printout
  if (my_id .eq. 0) then
    write(*,*) '****************************************'
    write(*,*) '*  construct harmonic matrix           *'
    write(*,*) '****************************************'
  endif
  
  ! --- Memory tracking
  call tr_print_memsize("DebConstM")
 
  
!  if(my_id .lt. m_cpu) then 
!    mumps_par%nz = nz_total/(n_tor*n_tor) 
!    mumps_par%n  = ndof_glob  / n_tor  
!  else 
!    mumps_par%nz = 4*nz_total/(n_tor*n_tor)
!    mumps_par%n  = 2*ndof_glob  / n_tor  
!  endif  

 
  
!  n_local_elms = element_list%n_elements
 
!  mumps_par%nz = (i_tor_max-i_tor_min+1) * (i_tor_max-i_tor_min+1) * nz_total / (n_tor*n_tor) 
!  mumps_par%n  = (i_tor_max-i_tor_min+1) * ndof_glob  / n_tor  

  mumps_par%nz = nz_glob_harm  
  mumps_par%n  = n_glob_harm
!  write(*,*) 'size matrices : n, nz =', my_id, mumps_par%n, mumps_par%nz
!  write(*,*) 'my_id, index_min, index_max =', my_id, index_min, index_max !mumps_par%n, mumps_par%nz
!write(*,'(i3,A,12i8)') my_id,' recv_counts : ',recv_counts
!write(*,'(i3,A,12i8)') my_id,' recv_disp   : ',recv_disp
!write(*,'(i3,A,12i8)') my_id,' sizes       : ',sizes

  if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,"dh_mumps_par%A",CAT_DMATRIX)
  if (associated(mumps_par%irn))  call tr_deallocatep(mumps_par%irn,"dh_mumps_par%irn",CAT_DMATRIX)
  if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"dh_mumps_par%jcn",CAT_DMATRIX)

  call tr_allocatep(mumps_par%A,1,mumps_par%nz,"dh_mumps_par%A",CAT_DMATRIX)
  call tr_allocatep(mumps_par%irn,1,mumps_par%nz,"dh_mumps_par%irn",CAT_DMATRIX)
  call tr_allocatep(mumps_par%jcn,1,mumps_par%nz,"dh_mumps_par%jcn",CAT_DMATRIX)

  if (associated(mumps_par%rhs))  call tr_deallocatep(mumps_par%rhs,"dh_mumps_par%rhs",CAT_DMATRIX)
  call tr_allocatep(mumps_par%rhs,1,mumps_par%n,"dh_mumps_par%rhs",CAT_DMATRIX)


  ! --- Initialise internal variables
!  index_min     = 1
!  index_max     = 2609
  mumps_par%A   = 0.d0
  mumps_par%rhs = 0.d0
  mumps_par%irn = 0
  mumps_par%jcn = 0
  difference_found = .false.
  rhs_problem(:)   = .false.
  elm_problem(:,:) = .false.
!!mumps_par%irn,mumps_par%jcn,mumps_par%A,    
  
  ! --- Declare shared and private variables for omp
  !$omp parallel default(none) &
  !$omp   shared(n_local_elms,local_elms,element_list,node_list,& 
  !$omp          index_min, index_max,xpoint2,xcase2,R_axis,Z_axis,psi_axis,psi_bnd,Z_xpoint,       &
  !$omp          R_xpoint,my_id,bc_natural_open,bc_natural_flux,refinement,thread_struct,n_tor_fft_thresh, &
  !$omp          difference_found,rhs_problem,elm_problem, i_tor_min, i_tor_max, mumps_par,      & 
  !$omp                     ijA_index_harm, ijA_size_harm, irn_jcn_harm)                                                 &
  !$omp   private(ife,ielm,iv,inode,element,nodes,i,inode1,i_order,index_node1,          &
  !$omp           index_large_i,j,index_ij,k,knode,k_order,index_node2,index_large_k,ijA_position,         &
  !$omp           l,index_kl,ilarge2,iv2,vertex,direction,inode2,omp_nthreads,omp_tid,                     &
  !$omp           i_father,element_father, nodes_father, inode_father, node_out, ivertex, iorder,          &
  !$omp           ivar, itor, jvertex, jorder, jvar, jtor, random_element, n_var_reduced, v1, v2, im,      &
  !$omp           index_ij_model400_e, index_kl_model400_e,  tmp_rhs, tmp_elm, tmp_elm_v2_8        )

!! --- omp id
#ifdef _OPENMP
  omp_nthreads = omp_get_num_threads()
  omp_tid      = 1+omp_get_thread_num()
#else
  omp_nthreads = 1
  omp_tid      = 1
#endif

!   write(*,*) my_id, index_min, index_max 

!   write(*,*) 'PSV', my_id, n_local_elms
  ! --- Loop over local elements
  !$omp do schedule(runtime)
  do ife =1, n_local_elms!element_list%n_elements
  !do ife =1, element_list%n_elements
      
    ielm = local_elms(ife)
    ! --- Get element
    element = element_list%element(ielm)
       
      do iv = 1, n_vertex_max
       inode   = element%vertex(iv)
       nodes(iv) = node_list%node(inode)
      enddo
     
!    write(*,*) my_id, i_tor_min, i_tor_max, 'PSV'
     
    call elementary_matrix_build(element, nodes, xpoint2, xcase2, R_axis, &
         &                       Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint,   &
         &                       omp_tid, ife, n_local_elms, node_list, i_tor_min, i_tor_max)
    
    !if(my_id .eq. 0) write(*,*) 'PSV', i_tor_min, i_tor_max
!    write(*,*) 'PSV', i_tor_min, i_tor_max
    
    ! --- Define element nodes (depends if it's refined)
      do i=1, n_vertex_max
      node_out(i) = element%vertex(i)   
      enddo 

    ! --- We only look at non-refined elements
    
      do i=1,n_vertex_max

        inode1 = node_out(i)

        do i_order = 1, n_order+1

          index_node1 = node_list%node(inode1)%index(i_order)

          index_large_i = (i_tor_max - i_tor_min + 1) * n_var * (index_node1 - 1)
         
          if ((index_node1 .ge. index_min) .and. (index_node1 .le. index_max)) then
            do j = 1, n_var * (i_tor_max - i_tor_min + 1)

              index_ij = (i_tor_max - i_tor_min + 1) * n_var * (n_order+1) * (i-1) + (i_tor_max - i_tor_min + 1) * n_var * (i_order-1) + j   ! index in the ELM matrix
              !$omp atomic
              mumps_par%rhs(index_large_i+j) = mumps_par%rhs(index_large_i+j) + thread_struct(omp_tid)%RHS(index_ij)
              !$omp end atomic
            end do

            do k=1,n_vertex_max

              knode = node_out(k)

              do k_order = 1, n_order+1

                index_node2 = node_list%node(knode)%index(k_order)

                index_large_k = (i_tor_max - i_tor_min + 1) * n_var * (index_node2 - 1)
          !      if(my_id.eq.0) write(*,*)  index_large_k, 'PSV'
                call locate_irn_jcn(index_node1,index_node2,index_min,index_max,ijA_position,&
                                    ijA_index_harm, ijA_size_harm, irn_jcn_harm)
!                if(my_id.eq.0) write(*,*)  my_id, index_min, index_max, ijA_position, 'PSVERMA'
                  
                thread_struct(omp_tid)%synch_buff(:) = 0.d0

                do j = 1, n_var * (i_tor_max - i_tor_min + 1)
                  index_ij = (i_tor_max - i_tor_min + 1) * n_var * (n_order+1) * (i-1) + (i_tor_max - i_tor_min + 1) * n_var * (i_order-1) + j   ! index in the ELM matrix

                  do l = 1, n_var * (i_tor_max - i_tor_min + 1)

                    index_kl = (i_tor_max - i_tor_min + 1) * n_var * (n_order+1) * (k-1) + (i_tor_max - i_tor_min + 1) * n_var * (k_order-1) + l   ! index in the ELM matrix

                    ilarge2 = ijA_position - 1 + (j-1) * n_var*(i_tor_max - i_tor_min + 1) + l

                    mumps_par%irn(ilarge2) = index_large_i + j
                    mumps_par%jcn(ilarge2) = index_large_k + l
                    
                    thread_struct(omp_tid)%synch_buff((j-1)*n_var*(i_tor_max - i_tor_min + 1)+l) = &
                      thread_struct(omp_tid)%synch_buff((j-1)*n_var*(i_tor_max - i_tor_min + 1)+l) + thread_struct(omp_tid)%ELM(index_ij,index_kl)

                  enddo
                  
                enddo ! n_order+1

                !$omp critical
                mumps_par%A(ijA_position : ijA_position + n_var*(i_tor_max - i_tor_min + 1)*n_var*(i_tor_max - i_tor_min + 1) - 1) = &
                  mumps_par%A(ijA_position : ijA_position + n_var*(i_tor_max - i_tor_min + 1)*n_var*(i_tor_max - i_tor_min + 1) - 1) +  &
                  thread_struct(omp_tid)%synch_buff(:) 
                !$omp end critical                 

              enddo ! n_vertex_max
            enddo ! n_var * (i_tor_max - i_tor_min + 1)

          endif ! index_min < index < index_max
        enddo ! n_order+1

      enddo ! n_vertex_max 
       
!      if(my_id.eq.1)  write(*,*) 'my_id=', my_id, mumps_par%A
          !write(*,*) 'my_id=', my_id, 'PSV' 
  end do
  !$omp end do
  !$omp end parallel

 !   if(my_id.eq.1) write(*,*)  'mumps_par%A =', mumps_par%A

  ! --- Memory tracking
  call tr_vnorms("cm_A_bef_bc",mumps_par%A,mumps_par%nz)

  
  ! --- Apply boundary conditions.
  call boundary_conditions(my_id, node_list, element_list,  bnd_node_list,local_elms, n_local_elms,            &
                           index_min, index_max,  mumps_par%rhs, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd,  &
                           R_xpoint, Z_xpoint, psi_xpoint, .false., .false.,                & 
       &   ijA_index_harm, ijA_size_harm, irn_jcn_harm,        & 
       &   mumps_par%irn, mumps_par%jcn, mumps_par%A,          & 
       &   i_tor_min, i_tor_max )

  ! --- Memory tracking
  call tr_vnorms("cm_A_aft_bc",mumps_par%A,mumps_par%nz)

end subroutine construct_harmonic_matrix





end module construct_harmonic_matrix_mod
