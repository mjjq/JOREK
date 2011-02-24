subroutine resistive_wall_starwall(my_id, node_list, bnd_elm_list, bnd_node_list)
!-------------------------------------------------------------------
! Reads the STARWALL resistive wall vacuum response matrices from files
!-------------------------------------------------------------------

  use data_structure
  use vacuum_response
  use phys_module
  
  
  implicit none


  integer,                     intent(in) :: my_id            ! MPI thread number of current thread
  type(type_node_list),        intent(in) :: node_list        ! List of boundary nodes
  type(type_bnd_element_list), intent(in) :: bnd_elm_list    ! List of boundary elements
  type(type_bnd_node_list),    intent(in) :: bnd_node_list    ! List of boundary nodes

  integer :: dim(2)
  integer :: dim1, i, j, k
  real*8  :: TWOPI, a, b
  
  TWOPI = 8.d0 * atan(1.d0)
  
  write(*,*) '@@> resistive_wall_starwall'
  
  ! --- Read number of wall currents.
  write(*,*) 'Determining number of wall currents.'
  open(42, FILE='starwall_d_yy', status='old', action='read')
  read(42,*) n_wall_curr
  write(*,'(1x,a,i7)') 'n_wall_curr=', n_wall_curr
  close(42)
  
  
  ! --- Read 'EE' matrix.
  dim(:) = response_index(bnd_node_list%n_bnd_nodes, n_tor, 2)
  call read_response_matrix( matrix_ee, dim, 'starwall_m_ee' )
  matrix_ee(:,:) = matrix_ee(:,:) * TWOPI


  ! --- Read 'EY' matrix.
  dim(:) = (/ response_index(bnd_node_list%n_bnd_nodes, n_tor, 2), n_wall_curr /)
  call read_response_matrix( matrix_ey, dim, 'starwall_m_ey' )
  matrix_ey(:,:) = matrix_ey(:,:) * TWOPI
  
  
  ! --- Read 'YE' matrix.
  dim(:) = (/ n_wall_curr, response_index(bnd_node_list%n_bnd_nodes, n_tor, 2) /)
  call read_response_matrix( matrix_ye, dim, 'starwall_m_ye' )
  
  
  ! --- Read 'YY' matrix (diagonal matrix).
  call read_response_diagonal( diagonal_yy, n_wall_curr, 'starwall_d_yy' )
  
  
  ! --- Prepare the wall current array.
  if ( allocated(wall_curr) ) deallocate(wall_curr)
  allocate( wall_curr(n_wall_curr) )
  wall_curr = 0.d0
  
  
  ! --- Prepare additional matrices required for implicit wall current time-evolution.
  if ( .not. wall_curr_treatment == 'explicit' ) then
  
    write(*,*) 'Determining derived response matrices for implicit wall-current evolution.'
    
    dim1 = response_index(bnd_node_list%n_bnd_nodes, n_tor, 2)
    allocate( diagonal_r(n_wall_curr), matrix_s(n_wall_curr,dim1), &
              matrix_t(dim1, dim1), matrix_u(dim1 ,n_wall_curr)    )
    
    diagonal_r(:) = 1.d0 / ( 1.d0 + wall_thickness / ( wall_resistivity * tstep * diagonal_yy(:) ) )
    
    matrix_s = 0.d0
    do j = 1, dim1
      matrix_s(:,j) = matrix_ye(:,j) / ( 1.d0 + tstep * wall_resistivity * diagonal_yy(:) / wall_thickness )
    end do
    
    matrix_t(:,:) = matrix_ee(:,:) - matmul( matrix_ey(:,:), matrix_s(:,:) )
    
    matrix_u = 0.d0
    do i = 1, dim1
      matrix_u(i,:) = matrix_ey(i,:) * ( 1.d0 - diagonal_r(:) )
    end do
    
    allocate( matrix_v(dim1, dim1) )
    matrix_v = 0.d0
    do i = 1, dim1
      do j = 1, dim1
        do k = 1, n_wall_curr
          matrix_v(i,j) = matrix_v(i,j) + matrix_ey(i,k) * matrix_ye(k,j)
        end do
      end do
    end do    
    
  end if
  
  
  ! --- Check matrix consistency: Compare resistive wall with zero resistivity to ideal wall
  if ( wall_resistivity == 0.d0 ) then
    
    21 format(1x,a,1x,l,1x,2ES14.5)
    write(*,*) 'Checking matrix consistency.'
    call ideal_wall_starwall(my_id,node_list,bnd_elm_list,bnd_node_list)
    
    a = sum(abs( diagonal_r ))
    write(*,21) 'diagonal_r:', a == 0.d0, a
    
    a = sum(abs( matrix_s - matrix_ye ))
    b = sum(abs( matrix_ye ))
    write(*,21) 'matrix_s:  ', a < 1.d-5*b, a/b
    
    a = sum(abs( matrix_u - matrix_ey ))
    b = sum(abs( matrix_ey ))
    write(*,21) 'matrix_u:  ', a < 1.d-5*b, a/b
    
    a = sum(abs( matrix_t - vac_response ))
    b = sum(abs( vac_response ))
    write(*,21) 'matrix_t:  ', a < 1.d-5*b, a/b
    
    a = sum(abs( matrix_ee-matmul(matrix_ey,matrix_ye) - vac_response ))
    b = sum(abs( vac_response ))
    write(*,21) 'ee-ey*ye:  ', a < 1.d-5*b, a/b
    
  end if

  write(*,*) '@@< resistive_wall_starwall'
  
end subroutine resistive_wall_starwall

