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
  use mod_vtk
  use equil_info
  use mod_elt_matrix_fft
  use omp_lib
  use mpi_mod
  use mod_impurity, only: init_imp_adas

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
  integer,allocatable  :: ien (:,:)
  real*4,allocatable   :: xyz (:,:), scalars(:,:), scalars_o(:,:)
  character*12, allocatable :: scalar_names(:), scalar_names_o(:)

  integer :: my_id, ierr, k_tor, i, j, k, n(4), ife, iv, inode, index_node, index_total
  integer :: omp_nthreads, omp_tid, n_tor_local, i_order, index_large_i, index_ij
  integer :: index_RHS, k_var, i_tor, only_term(2), iterm, term_count, nsub, nnos, ielm, it_o 
  real*8, allocatable :: result(:,:,:,:), res2d(:,:,:)
  real*8,     allocatable :: rhs(:)
  integer, parameter      :: max_terms=20
  character(len=64)       :: file_name, label 
integer   :: required,provided,StatInfo
#ifdef USE_FFTW
  real*8     :: in_fft(1:n_plane)
  complex*16 :: out_fft(1:n_plane)
#endif

  ! --- Initialize FFTW
#ifdef USE_FFTW
  call dfftw_plan_dft_r2c_1d(fftw_plan,n_plane,in_fft,out_fft,FFTW_PATIENT)
#endif

  call MPI_Init_thread(required, provided, StatInfo)
  call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread

  ! --- Normal initialization
  allocate(node_list)
  allocate(element_list)
  allocate(bnd_elm_list)
  allocate(bnd_node_list)
  my_id = 0
  call initialise_parameters(my_id, "__NO_FILENAME__")
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
  ! --- Read ADAS data and generate coronal equilibrium if needed
  call init_imp_adas(my_id)
#endif
  call det_modes()
  call import_restart(node_list, element_list, 'jorek_restart',  rst_format, ierr, .true.)
  call initialise_basis()
  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

  
  ! --- Initialize the plasma equilibrium data structure
  call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
  call print_equil_state(.false.)
  
  ! --- Initialize the new_diag framework and print some information (.true.)
  call init_new_diag(.true.)
  
  call tr_meminit(my_id, 1)
  call new_thread_buffers()

  write(*,*) 'OpenMP threads = ', nbthreads

  ! --- Initialize time stepping parameters
  call update_time_evol_params()

  if ( index_start <= 1 ) then
    tstep_prev = tstep
  else
    tstep_prev = xtime(index_start) - xtime(index_start-1)
  end if


  ! --- Calculate DOFs
  index_total = -1
  do inode=1,node_list%n_nodes
    index_total = max(index_total,maxval(node_list%node(inode)%index))
  enddo
  node_list%n_dof = index_total * n_tor * n_var

  allocate(rhs(node_list%n_dof))

  ! --- Create grid points and save them for the vtk
  nsub    = 4
  nnos    = nsub*nsub*element_list%n_elements
  allocate(xyz(3,nnos))
  xyz     = 0


  call create_pol_pos(pol_pos_list, ierr, node_list, element_list, ES, grid=.true., nsub=nsub)
  call create_tor_pos(tor_pos_list, ierr, nphi=1)

  do i = 1, pol_pos_list%n_pos(1)
    do j = 1, pol_pos_list%n_pos(2)
      xyz(1:3,i) = (/ pol_pos_list%pos(i,j)%R, pol_pos_list%pos(i,j)%Z, 0.d0 /)
    enddo
  enddo
 
  ! --- Create vtk grid indices
  allocate(ien(4,(nsub-1)*(nsub-1)*element_list%n_elements))
  ielm  = 0
  inode = 0
  ien   = 0

  do i=1,element_list%n_elements

    do j=1,nsub
      do k=1,nsub
        inode       = inode +1
      enddo
    enddo

    do j=1,nsub-1
      do k=1,nsub-1
        ielm        = ielm  +1
        ien(1,ielm) = inode - nsub*nsub + nsub*(j-1) + k-1       ! indices for VTK
        ien(2,ielm) = inode - nsub*nsub + nsub*(j  ) + k-1
        ien(3,ielm) = inode - nsub*nsub + nsub*(j  ) + k
        ien(4,ielm) = inode - nsub*nsub + nsub*(j-1) + k
      enddo
    enddo
  enddo

  term_count = 0  ! Counts terms with non-zero RHS

  ! --- Go over variables 
  do k_var=1, n_var

    do iterm=1, max_terms

      rhs = 0.0d0 

      only_term = (/ k_var, iterm /)

      ! --- Declare shared and private variables for omp
      !$omp parallel default(none) &
      !$omp   shared(element_list,node_list, ES, only_term,        &
      !$omp          xpoint,xcase, rhs, my_id, thread_struct,n_tor_fft_thresh)                      &
      !$omp   private(ife,iv,inode,element,nodes,i,i_order,               &
      !$omp           index_large_i,j,index_ij,k, index_node,         &
      !$omp           omp_nthreads,omp_tid  )

! --- omp id
#ifdef _OPENMP
      omp_nthreads = omp_get_num_threads()
      omp_tid      = 1+omp_get_thread_num()
#else
      omp_nthreads = 1
      omp_tid      = 1
#endif

      !$omp do schedule(runtime)
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
      !$omp end do
      !$omp end parallel
 
      if (sum(abs(rhs)) < 1.d-30) cycle
   
      term_count = term_count + 1

      write(*,*) 'term counter  = ', term_count
      write(*,*) 'sum RHS       = ', sum(rhs)
      write(*,*) 'variable      = ', k_var
      Write(*,*) 'term number   = ', iterm 
      Write(*,*) ' ' 
      
      ! --- Replace node psi values in pol_pos_list by RHS
      do i = 1, pol_pos_list%n_pos(1)

        do i_order=1, n_order+1

          do iv=1, n_vertex_max
          
            index_node = pol_pos_list%pos(i,1)%nodes(iv)%index(i_order)
  
            do i_tor=1, n_tor
  
              index_RHS = n_tor*n_var*(index_node - 1) + n_tor*(k_var-1) + i_tor 
        
              pol_pos_list%pos(i,1)%nodes(iv)%values(i_tor, i_order, 1) = rhs(index_RHS)
  
            enddo 
          enddo
        enddo
      enddo

      
      it_o = term_count - 1
      if (term_count/=1) then 
        scalars_o(:, 1:it_o) = scalars(:,1:it_o)
        scalar_names_o(1:it_o) = scalar_names(1:it_o)
      endif

      if (allocated(scalars))      deallocate(scalars)
      if (allocated(scalar_names)) deallocate(scalar_names)
      allocate(  scalars(nnos,1:term_count),   scalar_names(term_count))

      ! --- Evaluate RHS term (Psi) in the grid points
      expr_list = exprs((/'Psi       '/), 1, 0)
      call eval_expr(ES, JOREK_UNITS, expr_list, pol_pos_list, tor_pos_list, result, ierr)
      call reduce_result_to_2d(ierr, result, res2d, i1=1)

      ! --- Save RHS values into scalar vtk vector 
      scalars(:,term_count)    = res2d(:,1,1)
      write (label,'(a, i2.2, a, i3.3)') 'RHS_', k_var, '_', iterm
      scalar_names(term_count) = label

      ! Recover old values (to append the array in fortran...) 
      if (term_count/=1) then
        scalars(:, 1:it_o) = scalars_o(:,1:it_o)
        scalar_names(1:it_o) = scalar_names_o(1:it_o)
      endif

      if (allocated(scalars_o))      deallocate(scalars_o)
      if (allocated(scalar_names_o)) deallocate(scalar_names_o)
      allocate(scalars_o(nnos,1:term_count), scalar_names_o(term_count))
      scalars_o(:, 1:term_count) = scalars(:,1:term_count)
      scalar_names_o(1:term_count) = scalar_names(1:term_count)

    enddo    ! --- max terms

  enddo  ! --- variables


  write (file_name,'(a, i5.5, a)') 'RHS.', index_start, '.vtk'
  call write_vtk(file_name,xyz,ien,9,scalar_names,scalars)
 
  deallocate(rhs, scalars, scalar_names, scalars_o, scalar_names_o, result, res2d)

#ifdef USE_FFTW
  call dfftw_destroy_plan(fftw_plan)
#endif


end program test_rhs_diagno
