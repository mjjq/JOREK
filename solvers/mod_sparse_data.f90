module mod_sparse_data

  use iso_c_binding
  use mpi
  use phys_module, only:    use_pastix, use_mumps, use_strumpack, iter_precon, gmres_max_iter
#ifdef USE_MUMPS
  use mod_mumps, only:      type_MUMPS_SOLVER
#endif  
#if (defined USE_PASTIX) || (defined USE_PASTIX6)
  use mod_pastix, only:     type_PASTIX_SOLVER
#endif
#ifdef USE_STRUMPACK
  use mod_strumpack, only:  type_STRUMPACK_SOLVER
#endif
  use data_structure, only: type_PRECOND

  
  type type_SP_SOLVER
#ifdef USE_MUMPS
    type(type_MUMPS_SOLVER)     :: mmss
#endif  
#if (defined USE_PASTIX) || (defined USE_PASTIX6)
    type(type_PASTIX_SOLVER)    :: ptss
#endif
#ifdef USE_STRUMPACK
    type(type_STRUMPACK_SOLVER) :: spss
#endif
    type(type_PRECOND)          :: pc
    
    integer                     :: index_now                           !< current time step index (absolute)
    real(kind=8)                :: tstep                               !< current time step value
    real(kind=8)                :: tstep_prev                          !< previous time step
    integer                     :: istep                               !< curent time step index within jstep group
    
    integer                     :: iter_precon                         !< maximum number of iterations without pc update (input)
    integer                     :: max_steps_noUpdate                  !< maximum number of time steps without pc update (input)
    integer                     :: iter_max                            !< maximum allowed number of iterations (input)
    
    integer                     :: n_since_update = 0                  !< number of time steps since last pc update
    integer                     :: iter_prev = 0                       !< number of iterations in the previous step
    integer                     :: iter_gmres                          !< number of iterations in the current step
    real(kind=8)                :: iter_tol                            !< iterative convergence criteria
    logical                     :: solve_only = .false.
    logical                     :: step_success = .false.
    logical                     :: iterative = .false.
    logical                     :: equilibrium = .false.
    
  end type type_SP_SOLVER
  
  integer, parameter :: MHD_EQUILI = 0
  integer, parameter :: MHD_DIRECT = 1
  integer, parameter :: MHD_PRECON = 2
   
  private
  public :: MHD_EQUILI, MHD_DIRECT, MHD_PRECON, type_SP_SOLVER

end module mod_sparse_data
