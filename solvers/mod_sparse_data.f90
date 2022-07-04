module mod_sparse_data

  use iso_c_binding
  use mpi
  use phys_module, only:    use_pastix, use_mumps, use_strumpack, iter_precon, gmres_max_iter
#ifdef USE_PASTIX
  use mod_pastix, only:     type_PASTIX_SOLVER
#endif
#ifdef USE_STRUMPACK
  use mod_strumpack, only:  type_STRUMPACK_SOLVER
#endif
  use data_structure, only: type_PRECOND

  type type_SP_SOLVER
#ifdef USE_PASTIX
    type(type_PASTIX_SOLVER)    :: ptss
#endif
#ifdef USE_STRUMPACK
    type(type_STRUMPACK_SOLVER) :: spss
#endif
    type(type_PRECOND)          :: pc
    integer                     :: iter_prev, iter_gmres
  end type type_SP_SOLVER
  
  integer, parameter :: MHD_EQUILI = 0
  integer, parameter :: MHD_DIRECT = 1
  integer, parameter :: MHD_PRECON = 2
   
  private
  public :: MHD_EQUILI, MHD_DIRECT, MHD_PRECON, type_SP_SOLVER

end module mod_sparse_data
