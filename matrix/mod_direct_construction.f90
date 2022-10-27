module mod_direct_construction
#ifdef DIRECT_CONSTRUCTION

  implicit none
  public update_pc_mat

contains

  subroutine update_pc_mat(pc, a_mat, sim)

    use mod_parameters, only : n_tor, n_var
    use mpi_mod
    use mod_integer_types
    use data_structure, only: type_SP_MATRIX, type_PRECOND, type_RHS, type_MHD_SIM    
    
    implicit none
    
    type(type_PRECOND)                 :: pc
    type(type_SP_MATRIX)               :: a_mat
    type(type_MHD_SIM)                 :: sim
    
    if (.not.pc%structured) call set_pc_structure(pc, a_mat, sim)
    
    !call construct_matrix(sim%my_id, pc%local_elms, pc%n_local_elms, xpoint, xcase, ES%R_axis, ES%Z_axis,&
    !                      ES%psi_axis, ES%psi_bnd, ES%R_xpoint, ES%Z_xpoint, ES%psi_xpoint, &
    !                      a_mat, rhs_vec, harmonic_matrix=.false.)
    
  end subroutine update_pc_mat
  
  subroutine set_pc_structure(pc, a_mat, sim)
    use mpi_mod
    use tr_module
    use mod_integer_types
    use mod_parameters, only : n_tor, n_var
    use data_structure, only: type_SP_MATRIX, type_PRECOND, type_MHD_SIM
    use mod_global_matrix_structure, only: global_matrix_structure
    use global_distributed_matrix, only: global_matrix_structure_vacuum
    
    implicit none
    
    type(type_PRECOND)                 :: pc
    type(type_MHD_SIM)                 :: sim
    type(type_SP_MATRIX)               :: a_mat
    integer                            :: i_tor_min, i_tor_max
    integer                            :: ierr
    
    if (pc%my_id.eq.0) write(*,*) "Analyzing preconditioner"
    
    pc%mat%ng = (pc%mode_set_n)*(a_mat%ng)/n_tor ! rank of local PC matrix
    pc%mat%nr = pc%mat%ng
    pc%mat%nc = pc%mat%ng
    pc%n_glob = a_mat%ng
    
    i_tor_min = pc%mode_set(1)
    i_tor_max = pc%mode_set(pc%mode_set_n)
    
    call tr_allocatep(pc%local_elms,1,sim%element_list%n_elements,"local_elms_harm",CAT_FEM)
    
    call distribute_nodes_elements(sim%my_id, pc%n_cpu_n, sim%n_cpu, sim%node_list, sim%element_list, .true., pc%local_elms, & 
                                   pc%n_local_elms, sim%restart, sim%freeboundary, pc%mat)
                                   
    call global_matrix_structure(sim%node_list, sim%element_list, sim%bnd_elm_list, sim%freeboundary, &
                                 pc%local_elms, pc%n_local_elms, pc%mat, i_tor_min=i_tor_min, i_tor_max=i_tor_max)
                                 
    if (sim%freeboundary .and. (sim%sr_n_tor /= 0)) then 
      call global_matrix_structure_vacuum(sim%node_list, sim%bnd_node_list, pc%mat, i_tor_min=i_tor_min, i_tor_max=i_tor_max) 
    endif                                 

    write(*,*) sim%my_id, i_tor_min, i_tor_max, pc%n_local_elms
                                   
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Finalize(ierr)
    
  end subroutine set_pc_structure
  
  !
  !!> Constructing harmonic matrix directly from the elementary matrix 
  !subroutine direct_construction_harmonic(my_id, my_id_n, m_cpu, n_cpu, MPI_COMM_N,  MPI_COMM_MASTER, my_id_master, & 
  !  node_list, element_list, bnd_elm_list, bnd_node_list, xpoint, xcase, restart, freeboundary, direct_construction, a_mat)
  !
  !use data_structure 
  !use global_distributed_matrix
  !use mod_global_matrix_structure
  !use equil_info
  !use construct_matrix_mod, only : construct_matrix 
  !use mpi_mod
  !use vacuum
  !use mod_integer_types
  !
  !implicit none
  !
  !! --- Routine parameters
  !type (type_node_list),        intent(in)    :: node_list
  !type(type_bnd_node_list),     intent(inout) :: bnd_node_list
  !type (type_element_list),     intent(in)    :: element_list
  !type (type_bnd_element_list), intent(in)    :: bnd_elm_list
  !type(type_SP_MATRIX)                        :: a_mat
  !
  !integer, intent(in) :: my_id, my_id_n, n_cpu, m_cpu, MPI_COMM_N, MPI_COMM_MASTER, my_id_master, xcase
  !logical, intent(in) :: direct_construction, xpoint, restart, freeboundary
  !
  !! --- Local variables
  !integer, dimension(:), pointer :: index_min_harm => null(), index_max_harm => null() !< division of work across processes  
  !
  !
  !integer,               allocatable :: local_elms_harm(:)
  !integer                            :: n_local_elms_harm
  !integer(kind=int_all)              :: ndof 
  !integer                            :: i_tor_min, i_tor_max 
  !integer                            :: i, ierr
  !type(type_RHS)                     :: rhs_vec
  !  
  !! --- Memory allocation 
  !if (allocated(local_elms_harm)) call tr_deallocate(local_elms_harm,"local_elms_harm",CAT_DMATRIX) 
  !if (associated(index_min_harm)) call tr_deallocatep(index_min_harm,"index_min_harm",CAT_DMATRIX) 
  !if (associated(index_max_harm)) call tr_deallocatep(index_max_harm,"index_max_harm",CAT_DMATRIX) 
  !
  !call tr_allocate(local_elms_harm,1,element_list%n_elements,"local_elms_harm",CAT_FEM)
  !call tr_allocatep(index_min_harm,1,n_cpu,"index_min_harm",CAT_FEM)
  !call tr_allocatep(index_max_harm,1,n_cpu,"index_max_harm",CAT_FEM)
  !
  !if(my_id .lt. m_cpu)  then
  !  i_tor_min = 1
  !  i_tor_max = 1
  !else
  !  i_tor_min = 2*(my_id - MOD(my_id, m_cpu))/m_cpu
  !  i_tor_max = i_tor_min + 1
  !endif
  !
  !a_mat%comm = MPI_COMM_N
  !
  !call distribute_nodes_elements(my_id, m_cpu, n_cpu, node_list, element_list, direct_construction, local_elms_harm, & 
  !       n_local_elms_harm, restart, freeboundary, a_mat)  
  !
  !if ( .not. matrix_structure_initialized ) then
  !    call global_matrix_structure(node_list, element_list, bnd_elm_list, freeboundary, &
  !       local_elms_harm, n_local_elms_harm, a_mat, i_tor_min=i_tor_min, i_tor_max=i_tor_max)
  !       
  !  call MPI_Barrier(a_mat%comm, ierr)
  !
  !  if ( freeboundary .and. ( sr%n_tor /= 0 ) ) then 
  !    call global_matrix_structure_vacuum(node_list, bnd_node_list, a_mat, i_tor_min=i_tor_min, i_tor_max=i_tor_max) 
  !  endif         
  !  
  !  matrix_structure_initialized = .true.
  !endif
  !
  !call construct_matrix(my_id, local_elms_harm, n_local_elms_harm, xpoint, xcase, ES%R_axis, ES%Z_axis,&
  !       ES%psi_axis, ES%psi_bnd, ES%R_xpoint, ES%Z_xpoint, ES%psi_xpoint, &
  !       a_mat, rhs_vec, harmonic_matrix=.true.)
  !
  !end subroutine direct_construction_harmonic
  !
#endif
end module mod_direct_construction
