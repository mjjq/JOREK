module high_resolution_wall
  real*8,  parameter    :: wall_accuracy = 1.d-3
  integer               :: n_wall_HR
  real*8,  allocatable  :: R_wall_HR(:), Z_wall_HR(:), psi_wall_HR(:)
  integer, allocatable  :: initial_wall_piece(:)

  contains  


  !> This routine increases the resolution of the limiter
  subroutine get_high_resolution_wall(node_list, element_list)

    use data_structure
    use phys_module
    use mod_interp, only: interp
    implicit none
    
    ! --- Routine variables
    type (type_node_list),    intent(in)        :: node_list
    type (type_element_list), intent(in)        :: element_list
    
    ! --- Internal variables
    integer                     :: i, j
    integer                     :: n_wall_tmp, n_tmp
    real*8,  allocatable        :: R_wall_tmp(:), Z_wall_tmp(:)
    integer, allocatable        :: initial_wall_piece_tmp(:)
    real*8                      :: total_length, length
    real*8                      :: R_out, Z_out
    real*8                      :: s_out, t_out
    real*8                      :: P_s,P_t,P_st,P_ss,P_tt
    integer                     :: i_elm_out,ier
    
    ! --- Initialise
    n_wall_HR = 0
    
    ! --- Get total length of wall
    total_length = 0.d0
    do i=1,n_limiter-1
      total_length = total_length + sqrt( (R_limiter(i+1) - R_limiter(i))**2.d0 + (Z_limiter(i+1) - Z_limiter(i))**2.d0 )
    enddo
    
    ! --- Deduce approximate number of new wall points and allocate temporary vectors (take times 2 to make sure we don't step off array...)
    n_wall_tmp = 2 * total_length / wall_accuracy
    allocate(R_wall_tmp(n_wall_tmp), Z_wall_tmp(n_wall_tmp), initial_wall_piece_tmp(n_wall_tmp))
  
    ! --- Now step on each piece and get required accuracy
    do i=1,n_limiter-1
      length = sqrt( (R_limiter(i+1) - R_limiter(i))**2.d0 + (Z_limiter(i+1) - Z_limiter(i))**2.d0 )
      n_tmp = length / wall_accuracy
      n_tmp = min(1, n_tmp)
      do j=1,n_tmp
        n_wall_HR = n_wall_HR + 1
        if (n_wall_HR .gt. n_wall_tmp) then
          write(*,*)'Warning! Stepping off array in get_high_resolution_wall, aborting...'
          return
        endif
        R_wall_tmp(n_wall_HR) = R_limiter(i) + (R_limiter(i+1) - R_limiter(i)) * real(j-1)/real(n_tmp)
        Z_wall_tmp(n_wall_HR) = Z_limiter(i) + (Z_limiter(i+1) - Z_limiter(i)) * real(j-1)/real(n_tmp)
        initial_wall_piece_tmp(n_wall_HR) = i
      enddo
    enddo
    
    ! --- Finally allocate HR_wall data and copy
    allocate(R_wall_HR(n_wall_HR), Z_wall_HR(n_wall_HR), psi_wall_HR(n_wall_HR), initial_wall_piece(n_wall_HR))
    do i=1,n_wall_HR
      R_wall_HR(i) = R_wall_tmp(i)
      Z_wall_HR(i) = Z_wall_tmp(i)
      initial_wall_piece(i) = initial_wall_piece_tmp(i)
      call find_RZ(node_list,element_list,R_wall_HR(i),Z_wall_HR(i),R_out,Z_out,i_elm_out,s_out,t_out,ier)
      if (ier .eq. 0) then
        call interp(node_list,element_list,i_elm_out,1,1,s_out,t_out,psi_wall_HR(i),P_s,P_t,P_st,P_ss,P_tt)
      else
        psi_wall_HR(i) = 1.d10
      endif
    enddo
        
    
    ! --- Deallocate temporary arrays and exit
    deallocate(R_wall_tmp, Z_wall_tmp, initial_wall_piece_tmp)
    return

  end subroutine get_high_resolution_wall



end module high_resolution_wall




