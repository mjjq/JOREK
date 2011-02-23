module global_distributed_matrix

  implicit none
  
  real*8,  allocatable, target  :: A_glob(:),   rhs_glob(:)            ! the distributed global matrix and rhs
  integer, allocatable, target  :: irn_glob(:), jcn_glob(:)            ! the row and column indices for coordinate format sparse matrix/ (or CSR)
  integer, allocatable :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:)    ! contains the structure of the sparse matrix (to fill in CSR format)
  real*8,  allocatable :: deltas(:)                                    ! solution from previous step
  real*8,  allocatable :: column_scaling(:)                            ! column scaling of the global matrix
  integer, allocatable :: local_index_start(:), local_index_end(:)     ! range of indices local to one MPI process 
  integer              :: ndof_glob, n_glob, nz_glob
  integer              :: n_matrix_block_size                          ! Size of a matrix block (n_var x n_tor)
  
  
  
  contains
  
  
  
    !> Determine the matrix row or column for given values of ::i_index, ::i_var, and ::i_tor.
  integer recursive function det_row_col(i_index, i_var, i_tor)
    
    use parameters, only: n_tor, n_var
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: i_index   !< node%index property for the node and node degree of freedom
    integer, intent(in) :: i_var     !< Variable number
    integer, intent(in) :: i_tor     !< Toroidal mode number
    
    det_row_col =  n_tor * n_var * (i_index-1) + n_tor * (i_var-1) + i_tor
    
  end function det_row_col
  
  
  
  
  
  
  !> Determine the position of a matrix entry given by its row and column positions (::i_row and
  !! ::j_col) in the sparse matrix structure. @todo put somewhere else
  integer recursive function det_sparse_pos(i_row, j_col, index_min)
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: i_row     !< Matrix row
    integer, intent(in) :: j_col     !< Matrix column
    integer, intent(in) :: index_min !< Smallest block index dealt with by current MPI proc
    
    ! --- Local variables
    integer :: i_block, j_block           ! Block indices
    integer :: i_row_block, j_col_block   ! Row and column in the block
    integer :: ij_sparse_block            ! Position of the block in the sparse matrix structure
    integer :: i_block_local              ! Block index at local MPI proc
    logical :: found_index
    integer :: i
    
    i_block       = (i_row-1) / n_matrix_block_size + 1
    i_block_local = i_block - index_min + 1
    i_row_block   = i_row - (i_block-1) * n_matrix_block_size
    
    j_block       = (j_col-1) / n_matrix_block_size + 1
    j_col_block   = j_col - (j_block-1) * n_matrix_block_size
    
    ! --- Determine the position of the block in the sparse matrix.
    found_index = .false.
    do i = 1, ijA_size(i_block_local)
      
      if ( irn_jcn(i_block_local,i) == j_block ) then ! Block index found?
        ij_sparse_block = ijA_index(i_block_local,i)
        det_sparse_pos  = ij_sparse_block-1 + n_matrix_block_size*(i_row_block-1) + j_col_block
        found_index     = .true.
        exit
      end if
      
    end do
    
    if ( .not. found_index ) det_sparse_pos = -999999999 !###
    
  end function det_sparse_pos
  
  
  
  
  
  
end module global_distributed_matrix
