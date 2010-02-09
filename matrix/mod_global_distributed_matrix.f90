module global_distributed_matrix
  real*8,  allocatable, target  :: A_glob(:),   rhs_glob(:)            ! the distributed global matrix and rhs
  integer, allocatable, target  :: irn_glob(:), jcn_glob(:)            ! the row and column indices for coordinate format sparse matrix/ (or CSR)
  integer, allocatable :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:)    ! contains the structure of the sparse matrix (to fill in CSR format)
  real*8,  allocatable :: deltas(:)                                    ! solution from previous step
  real*8,  allocatable :: column_scaling(:)                            ! column scaling of the global matrix
  integer, allocatable :: local_index_start(:), local_index_end(:)     ! range of indices local to one MPI process 
  integer              :: ndof_glob, n_glob, nz_glob
end module
