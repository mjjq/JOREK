module direct_construction_mod
  implicit none
contains  
subroutine direct_construction_harmonic(my_id, my_id_n, m_cpu, n_cpu, MPI_COMM_N,  MPI_COMM_MASTER, my_id_master, & 
                                        node_list, element_list, xpoint2, xcase2, freeboundary, direct_construction)
  use data_structure 
  use global_distributed_matrix
  use mod_global_matrix_structure
  use equil_info
  use construct_matrix_mod, only : construct_matrix 
  use mpi_mod
  implicit none

type (type_node_list)        :: node_list
type (type_element_list)     :: element_list
type (type_bnd_element_list) :: bnd_elm_list
type (type_surface_list)     :: flux_list
type (type_element)          :: element
type (type_node)             :: nodes(n_vertex_max)



integer, intent(in)          :: my_id, my_id_n, n_cpu, m_cpu,  MPI_COMM_N,  MPI_COMM_MASTER, my_id_master, xcase2
logical, intent(in)          :: direct_construction, xpoint2, freeboundary
integer, allocatable         :: index_min_harm(:), index_max_harm(:)
integer, allocatable         :: local_elms_harm(:)
integer                      :: n_local_elms_harm
integer                      :: ndof 
integer                      :: i_tor_min, i_tor_max 
integer                      :: i, ierr
 


     
      
      if (allocated(local_elms_harm)) call tr_deallocate(local_elms_harm,"local_elms_harm",CAT_DMATRIX) 
      if (allocated(index_min_harm))  call tr_deallocate(index_min_harm,"index_min_harm",CAT_DMATRIX) 
      if (allocated(index_max_harm))  call tr_deallocate(index_max_harm,"index_max_harm",CAT_DMATRIX) 
 
      call tr_allocate(local_elms_harm,1,element_list%n_elements,"local_elms_harm",CAT_FEM)
      call tr_allocate(index_min_harm,1,n_cpu,"index_min_harm",CAT_FEM)
      call tr_allocate(index_max_harm,1,n_cpu,"index_max_harm",CAT_FEM)

      if(my_id .lt. m_cpu)  then
       i_tor_min = 1
       i_tor_max = 1
      else
       i_tor_min = 2*(my_id - MOD(my_id, m_cpu))/m_cpu
       i_tor_max = i_tor_min + 1
      endif
 
      call distribute_nodes_elements(my_id,m_cpu,n_cpu,node_list,element_list, direct_construction, & 
                                     local_elms_harm, n_local_elms_harm, ndof, index_min_harm,index_max_harm)

      call global_matrix_structure(my_id,my_id_n,node_List,element_list,bnd_elm_list, freeboundary, &
                                   local_elms_harm,n_local_elms_harm,index_min_harm(my_id+1),          & 
                                   index_max_harm(my_id+1), ijA_index_harm, ijA_size_harm,             &
                                   irn_jcn_harm, irn_glob_harm, jcn_glob_harm, i_tor_min, i_tor_max,   & 
                                   n_glob_harm, nz_glob_harm, ndof_glob_harm, n_matrix_block_size_harm)
 
 
      ! --- Memory allocation
      !if (allocated(A_glob_harm))    call tr_deallocate(A_glob_harm,"A_glob_harm",CAT_DMATRIX) 
      !call tr_allocate(A_glob_harm,1,nz_glob_harm,"A_glob_harm",  CAT_DMATRIX)

      !if (allocated(irn_glob_harm))  call tr_deallocate(irn_glob_harm,"irn_glob_harm",CAT_DMATRIX)
      !call tr_allocate(irn_glob_harm,1,nz_glob_harm,"irn_glob_harm",  CAT_DMATRIX)
 
      !if (allocated(jcn_glob_harm))  call tr_deallocate(jcn_glob_harm,"jcn_glob_harm",CAT_DMATRIX)
      !call tr_allocate(jcn_glob_harm,1,nz_glob_harm,"jcn_glob_harm",  CAT_DMATRIX) 

      !if (allocated(rhs_glob_harm))  call tr_deallocate(rhs_glob_harm,"rhs_glob_harm",CAT_DMATRIX)
      !call tr_allocate (rhs_glob_harm,1,ndof_glob_harm,"rhs_glob_harm",CAT_DMATRIX)

      if (allocated(A_glob_harm)) deallocate(A_glob_harm) 
      allocate(A_glob_harm(1:nz_glob_harm))
      if (allocated(irn_glob_harm)) deallocate(irn_glob_harm) 
      allocate(irn_glob_harm(1:nz_glob_harm))
      if (allocated(jcn_glob_harm)) deallocate(jcn_glob_harm) 
      allocate(jcn_glob_harm(1:nz_glob_harm))
      if (allocated(rhs_glob_harm)) deallocate(rhs_glob_harm) 
      allocate(rhs_glob_harm(1:ndof_glob_harm))


      A_glob_harm     = 0.0d0 
      rhs_glob_harm   = 0.0d0 
      irn_glob_harm   = 0
      jcn_glob_harm   = 0

      call construct_matrix(my_id, MPI_COMM_N, my_id_n, MPI_COMM_MASTER, my_id_master,              &
      local_elms_harm, n_local_elms_harm, index_min_harm(my_id+1), index_max_harm(my_id+1), xpoint2,&
      xcase2, ES%R_axis, ES%Z_axis, ES%psi_axis, ES%psi_bnd, ES%R_xpoint, ES%Z_xpoint,              &
      ES%psi_xpoint, i_tor_min, i_tor_max, n_glob_harm, nz_glob_harm, ndof_glob_harm, A_glob_harm,  &
      rhs_glob_harm, irn_glob_harm, jcn_glob_harm, ijA_index_harm, ijA_size_harm, irn_jcn_harm,     &
      direct_construction)

      ! if(my_id .eq. 0) then
      !   do i = 1, ndof_glob_loc !mumps_par%n 
      !      print*, 'i, mumps_par%rhs:', i, rhs_glob_harm(i)
      !   enddo 
      ! endif
      ! call MPI_Barrier(MPI_COMM_WORLD, ierr)  



end subroutine direct_construction_harmonic
end module direct_construction_mod
