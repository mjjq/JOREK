module direct_construction_mod
  implicit none
contains  
subroutine direct_construction_harmonic(my_id, my_id_n, m_cpu, n_cpu, MPI_COMM_N,  MPI_COMM_MASTER, my_id_master, node_list, element_list, & 
                                        index_min_loc, index_max_loc, local_elms_loc, n_local_elms_loc, ijA_index_loc, ijA_size_loc,       & 
                                        irn_jcn_loc, irn_glob_loc, jcn_glob_loc, i_tor_min, i_tor_max, n_glob_loc, nz_glob_loc,            & 
                                        ndof_glob_loc, n_matrix_block_size_loc, direct_construction)
  use data_structure 
  use global_distributed_matrix
  use mod_global_matrix_structure
  use equil_info
  use construct_matrix_mod, only : construct_matrix
  implicit none

type (type_node_list)        :: node_list
type (type_element_list)     :: element_list
type (type_bnd_element_list) :: bnd_elm_list
type (type_surface_list)     :: flux_list
type (type_element)          :: element
type (type_node)             :: nodes(n_vertex_max)



integer, intent(in)          :: my_id, my_id_n, n_cpu, m_cpu,  MPI_COMM_N,  MPI_COMM_MASTER, my_id_master
logical, intent(in)          :: direct_construction
integer, intent(out)         :: index_min_loc(*), index_max_loc(*)
integer                      :: i_tor_min, i_tor_max 
integer                      :: local_elms_loc(*)
integer                      :: n_local_elms_loc
integer                      :: ndof 
integer,  allocatable        :: ijA_index_loc(:,:), ijA_size_loc(:), irn_jcn_loc(:,:)
integer,  allocatable        :: irn_glob_loc(:), jcn_glob_loc(:)
logical                      :: freeboundary
integer                      :: n_glob_loc, nz_glob_loc, ndof_glob_loc, n_matrix_block_size_loc
integer                      :: xcase2
logical                      :: xpoint2
 


      if(my_id .lt. m_cpu)  then
       i_tor_min = 1
       i_tor_max = 1
      else
       i_tor_min = 2*(my_id - MOD(my_id, m_cpu))/m_cpu
       i_tor_max = i_tor_min + 1
      endif
 
      call distribute_nodes_elements(my_id,m_cpu,n_cpu,node_list,element_list, direct_construction, & 
                                     local_elms_loc, n_local_elms_loc, ndof, index_min_loc,index_max_loc)

      call global_matrix_structure(my_id,my_id_n,node_List,element_list,bnd_elm_list, freeboundary, &
                                   local_elms_loc,n_local_elms_loc,index_min_loc(my_id+1),          & 
                                   index_max_loc(my_id+1), ijA_index_loc, ijA_size_loc,             &
                                   irn_jcn_loc, irn_glob_loc, jcn_glob_loc, i_tor_min, i_tor_max,   & 
                                   n_glob_loc, nz_glob_loc, ndof_glob_loc, n_matrix_block_size_loc)
 
      ! --- Memory allocation
      if (allocated(A_glob_harm))    call tr_deallocate(A_glob_harm,"A_glob_harm",CAT_DMATRIX) 
      call tr_allocate(A_glob_harm,1,nz_glob_harm,"A_glob_harm",  CAT_DMATRIX)

      if (allocated(irn_glob_harm))  call tr_deallocate(irn_glob_harm,"irn_glob_harm",CAT_DMATRIX)
      call tr_allocate(irn_glob_harm,1,nz_glob_harm,"irn_glob_harm",  CAT_DMATRIX)
 
      if (allocated(jcn_glob_harm))  call tr_deallocate(jcn_glob_harm,"jcn_glob_harm",CAT_DMATRIX)
      call tr_allocate(jcn_glob_harm,1,nz_glob_harm,"jcn_glob_harm",  CAT_DMATRIX) 

      if (allocated(rhs_glob_harm))  call tr_deallocate(rhs_glob_harm,"rhs_glob_harm",CAT_DMATRIX)
      call tr_allocate (rhs_glob_harm,1,ndof_glob_harm,"rhs_glob_harm",CAT_DMATRIX)


      A_glob_harm     = 0.0d0 
      rhs_glob_harm   = 0.0d0 
      irn_glob_harm   = 0
      jcn_glob_harm   = 0

      call construct_matrix(my_id, MPI_COMM_N, my_id_n, MPI_COMM_MASTER, my_id_master,          &
      local_elms_loc, n_local_elms_loc, index_min_loc(my_id+1), index_max_loc(my_id+1), xpoint2,&
      xcase2, ES%R_axis, ES%Z_axis, ES%psi_axis, ES%psi_bnd, ES%R_xpoint, ES%Z_xpoint,          &
      ES%psi_xpoint, i_tor_min, i_tor_max, n_glob_loc, nz_glob_loc, ndof_glob_loc, A_glob_harm, &
      rhs_glob_harm, irn_glob_harm, jcn_glob_harm, ijA_index_harm, ijA_size_harm, irn_jcn_harm, &
      direct_construction)


end subroutine direct_construction_harmonic
end module direct_construction_mod
