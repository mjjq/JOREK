!> Demonstration of the diagnostic framework mod_position / mod_expression / mod_four_filter / mod_straight_field_line / ...
program test_rhs_diagno 
  
  use mod_parameters
  use data_structure
  use tr_module 
  use phys_module
  use mod_boundary
  use mod_new_diag
  use basis_at_gaussian
  use mod_import_restart
  use equil_info
!!  use construct_matrix_mod
  use mod_elt_matrix_fft
  !$ use omp_lib

  implicit none
  
  type(type_node_list),         pointer :: node_list
  type(type_element_list),      pointer :: element_list
  type (type_bnd_element_list), pointer :: bnd_elm_list
  type (type_bnd_node_list),    pointer :: bnd_node_list
  type (type_element)                   :: element
  type (type_node)                      :: nodes(n_vertex_max)

  type(t_pol_pos_list) :: pol_pos_list
  type(t_tor_pos_list) :: tor_pos_list
  type(t_four_filter)  :: filter
  type(t_expr_list)    :: expr_list
  integer :: my_id, ierr, k_tor, i, j, k, n(4), ife, iv, inode, index_node, index_total
  integer :: omp_nthreads, omp_tid, n_tor_local, i_order, index_large_i, index_ij
  integer :: index_RHS, k_var, i_tor, only_term(2), iterm, term_count
  real*8, allocatable :: result(:,:,:,:), res0d(:), res1d(:,:), res2d(:,:,:)
  real*8,     allocatable :: rhs(:)
  integer, parameter      :: max_terms=100
  character(len=64)       :: file_name, label 

  type(type_node_list) :: node_list_rhs

 
  ! --- Normal initialization
  allocate(node_list)
  allocate(element_list)
  allocate(bnd_elm_list)
  allocate(bnd_node_list)
  my_id = 0
  call initialise_parameters(my_id, "__NO_FILENAME__")
  call det_modes()
  call import_restart(node_list, element_list, 'jorek_restart',  rst_format, ierr, .true.)
  call initialise_basis()
  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
  
  ! --- Initialize the plasma equilibrium data structure
  call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
  call print_equil_state(.false.)
  
  ! --- Initialize the new_diag framework and print some information (.true.)
  call init_new_diag(.true.)
  
 
  call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread


  call tr_meminit(my_id, 1)
  call new_thread_buffers()

#ifdef _OPENMP
  omp_nthreads = omp_get_num_threads()
  omp_tid      = 1+omp_get_thread_num()
#else
  omp_nthreads = 1
  omp_tid      = 1
#endif

  write(*,*) 'OpenMP threads = ', omp_nthreads
  ! --- initialize properly!
  tstep = 1
  tstep_prev = 1
  time_evol_zeta =1


  index_total = -1
  do inode=1,node_list%n_nodes
    index_total = max(index_total,maxval(node_list%node(inode)%index))
  enddo
  node_list%n_dof = index_total * n_tor * n_var

  allocate(rhs(node_list%n_dof))

  node_list_rhs%n_nodes = node_list%n_nodes 
  do i=1, node_list%n_nodes 
    node_list_rhs%node(i)%x              = node_list%node(i)%x
    node_list_rhs%node(i)%boundary       = node_list%node(i)%boundary
    node_list_rhs%node(i)%boundary_index = node_list%node(i)%boundary_index
    node_list_rhs%node(i)%deltas         = 0.d0                            
  enddo


  do k_var=1, n_var

    term_count = 0

    if (k_var/=6) cycle  ! JUST to get T RHS

    do iterm=1, max_terms

      rhs = 0.0d0 

      ! --- Get RHS of ZK_perp of the T equation
      only_term = (/ k_var, iterm /)


      do ife = 1, element_list%n_elements 
        
        element = element_list%element(ife)
          
        do iv = 1, n_vertex_max
         inode     = element%vertex(iv)
         nodes(iv) = node_list%node(inode)
        enddo
    
    
        call element_matrix_fft(element,nodes, xpoint, xcase, ES%R_axis, ES%Z_axis, ES%psi_axis, ES%psi_bnd,   &
         ES%R_xpoint, ES%Z_xpoint, thread_struct(omp_tid)%ELM, thread_struct(omp_tid)%RHS, omp_tid,       &
         thread_struct(omp_tid)%ELM_p, thread_struct(omp_tid)%ELM_n, thread_struct(omp_tid)%ELM_k,  &
         thread_struct(omp_tid)%ELM_kn, thread_struct(omp_tid)%RHS_p, thread_struct(omp_tid)%RHS_k, &
         thread_struct(omp_tid)%eq_g, thread_struct(omp_tid)%eq_s, thread_struct(omp_tid)%eq_t,     &
         thread_struct(omp_tid)%eq_p, thread_struct(omp_tid)%eq_ss, thread_struct(omp_tid)%eq_st,   &
         thread_struct(omp_tid)%eq_tt, thread_struct(omp_tid)%delta_g,                              &
         thread_struct(omp_tid)%delta_s, thread_struct(omp_tid)%delta_t, 1, n_tor, only_term)
    
        do iv=1,n_vertex_max
    
          inode = element%vertex(iv)
    
          do i_order = 1, n_order+1
    
            index_node = node_list%node(inode)%index(i_order)
    
            index_large_i = n_tor * n_var * (index_node - 1)
    
            do j = 1, n_var * n_tor
    
              index_ij = n_tor * n_var * (n_order+1) * (iv-1) + n_tor * n_var * (i_order-1) + j   ! index in the ELM matrix
               
              !$omp atomic
              rhs(index_large_i+j) = rhs(index_large_i+j) + thread_struct(omp_tid)%RHS(index_ij) 
              !$omp end atomic
            enddo 
    
          enddo ! order
        enddo ! vertex
    
      enddo ! --- elements

      if (sum(rhs) == 0.d0) cycle
   
      ! --- Save rhs term into node_list_rhs 
      term_count = term_count + 1

      do inode=1,node_list%n_nodes
        do i_order=1, n_order+1
        
          index_node = node_list%node(inode)%index(i_order)
          
          do i_tor=1, n_tor
      
            index_RHS = n_tor*n_var*(index_node - 1) + n_tor*(k_var-1) + i_tor 
      
            node_list_rhs%node(inode)%values(i_tor, i_order, 1) = rhs(index_RHS)
            
          end do
        end do
      end do

      ! Export vtk for variable k with several plotted terms   
      ! --- Plot expressions in vtk 
      expr_list = exprs((/'R         ', 'Z         ',  'Psi       '/), 3, 2)

      call create_pol_pos(pol_pos_list, ierr, node_list_rhs, element_list, ES, grid=.true., nsub=4)
      call create_tor_pos(tor_pos_list, ierr, nphi=2)
  
      call eval_expr(ES, JOREK_UNITS, expr_list, pol_pos_list, tor_pos_list, result, ierr)
  
      call reduce_result_to_2d(ierr, result, res2d, i1=1)
      write (file_name,'(a, i2.2, a, i3.3, a)') 'RHS_', k_var, '_', iterm, '.vtk'
      write (label,'(a, i2.2, a, i3.3)') 'RHS_', k_var, '_', iterm
      expr_list%expr(3)%name=label
      call write_vtk_2d(ierr, expr_list, res2d, file_name, (/1,2/), close1=.true.)
 
    enddo    ! --- max terms



  enddo !--- variables


 
 
end program test_rhs_diagno
