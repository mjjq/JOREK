!> New pastix module to be used with core version
module mod_pastix
#ifdef USE_PASTIX6
  use mod_integer_types
  
  type type_PASTIX_SOLVER
    type(C_PTR) :: spm, pastix_data ! sparse solver (distributed)
    type(C_PTR) :: iparm, dparm
    
    integer                                      :: comm = 0
    real(kind=8), dimension(:), pointer          :: solution_scaling => Null()    !< matrix column scaling to be applied to solution vector
    logical                                      :: initialized = .false.
    logical                                      :: analyzed    = .false.
    logical                                      :: equilibrium = .false.    
    logical                                      :: scaled      = .false.
    logical                                      :: refine      = .false.    

    integer(kind=int_all), dimension(:), pointer :: loc2glob, glob2loc  ! mapping for column distribution
    
    integer(kind=int_all), dimension(:), pointer :: irn => Null()
    integer(kind=int_all), dimension(:), pointer :: jcn => Null()
    real(kind=8), dimension(:), pointer          :: val => Null()
    real(kind=8), dimension(:), pointer          :: rhs_val => Null()    

  end type type_PASTIX_SOLVER
  
  private
  public :: type_PASTIX_SOLVER, pastix_initialize, pastix_analyze, pastix_factorize, pastix_solve, pastix_finalize

  interface
    subroutine ptx() bind(C)
      use iso_c_binding
    end subroutine ptx

    subroutine ptx_init(pastix_data,spm,iparm,dparm,comm) bind(C)

      use iso_c_binding
      implicit none

      type(C_PTR), intent(out) :: pastix_data, spm
      type(C_PTR), intent(out) :: iparm, dparm
      integer, intent(in) :: comm

    end subroutine ptx_init
  end interface

  contains
!> Initialize PaStiX solver instance
  subroutine pastix_initialize(ptss)
    use, intrinsic :: iso_c_binding
    use mpi_mod
    implicit none
      
    type(type_PASTIX_SOLVER)          :: ptss        

    integer ierr

    call ptx_init(ptss%pastix_data,ptss%spm,ptss%iparm,ptss%dparm,ptss%comm)
    call MPI_Barrier(ptss%comm,ierr)

    return
  end subroutine pastix_initialize
    
  subroutine pastix_finalize(ptss)
    implicit none
    
    type(type_PASTIX_SOLVER)          :: ptss
    
    return
  end subroutine pastix_finalize

  subroutine pastix_analyze(ptss)
    implicit none
    
    type(type_PASTIX_SOLVER)          :: ptss
    return
  end subroutine pastix_analyze  

  subroutine pastix_factorize(ptss)
    implicit none
    
    type(type_PASTIX_SOLVER)          :: ptss

    return
  end subroutine pastix_factorize
  
  subroutine pastix_solve(ptss,rhs_vec)
    use data_structure, only: type_RHS
  
    implicit none
    
    type(type_PASTIX_SOLVER)          :: ptss
    type(type_RHS)                    :: rhs_vec    
   
    return
  end subroutine pastix_solve
    
    
    
    
    
    



#elif USE_PASTIX
#include "pastix_fortran.h"

  use mod_integer_types

  type type_PASTIX_SOLVER
    integer                                      :: comm = 0
    real(kind=8), dimension(:), pointer          :: solution_scaling => Null()    !< matrix column scaling to be applied to solution vector
    logical                                      :: initialized = .false.
    logical                                      :: analyzed    = .false.
    logical                                      :: equilibrium = .false.    
    logical                                      :: scaled      = .false.
    logical                                      :: refine      = .false.
    
    integer(kind=int_all), dimension(:), pointer :: irn => Null()
    integer(kind=int_all), dimension(:), pointer :: jcn => Null()
    real(kind=8), dimension(:), pointer          :: val => Null()
    real(kind=8), dimension(:), pointer          :: rhs_val => Null()

    integer(kind=int_all)                        :: iparm(IPARM_SIZE)
    real*8                                       :: dparm(DPARM_SIZE)
    
    integer(kind=int_all), dimension(:), pointer :: perm_vars => Null()
    integer(kind=int_all), dimension(:), pointer :: iperm_vars => Null()

    integer(kind=8)                              :: idata    = 0
    integer(kind=int_all)                        :: sym      = API_SYM_NO    
    integer(kind=int_all)                        :: iter     = 250
    integer(kind=int_all)                        :: ricar    = 0
    integer(kind=int_all)                        :: iluk     = 3
    integer(kind=int_all)                        :: amalg    = 5 
    real*8                                       :: eps      = 1.d-12
    real*8                                       :: pivot    = 1.d-64
    integer(kind=int_all)                        :: maxthrd  = 1024
    integer(kind=int_all)                        :: verb     = API_VERBOSE_NO
    integer(kind=int_all)                        :: facto    = API_FACT_LU
    integer(kind=int_all)                        :: rhs      = 0
    integer(kind=int_all)                        :: nblock   = 0
    
  end type type_PASTIX_SOLVER

  private
  public :: type_PASTIX_SOLVER, pastix_initialize, pastix_set_mat, pastix_analyze, &
            pastix_factorize, pastix_solve, pastix_finalize

  contains
  
  subroutine pastix_set_mat(ptss, a_mat)
  end subroutine pastix_set_mat
  
  subroutine pastix_initialize(ptss)
    use mpi_mod
    use mod_clock
  
    implicit none
    
    type(type_PASTIX_SOLVER)          :: ptss
    integer                           :: my_id, ierr
    type(clcktype)                    :: t_itstart, t0, t1, t2, t3
    real*8                            :: tsecond
    
    allocate(ptss%perm_vars(ptss%nblock))
    allocate(ptss%iperm_vars(ptss%nblock))
    ptss%perm_vars(1:ptss%nblock) = 0
    ptss%iperm_vars(1:ptss%nblock) = 0    
    
    call MPI_COMM_RANK(ptss%comm, my_id, ierr)
    
    call pastix_init_nthreads(ptss)
    
    if (my_id .eq. 0) write(*,'(i5,A,i5)') my_id,' PastiX n_threads : ', ptss%iparm(IPARM_THREAD_NBR)
    
    ptss%iparm(IPARM_MODIFY_PARAMETER)  = API_NO          ! insert default values
    ptss%iparm(IPARM_START_TASK)        = API_TASK_INIT   ! initializse
    ptss%iparm(IPARM_END_TASK)          = API_TASK_INIT  
  
    if (my_id .eq. 0) write(*,*) "Initializing PaStiX solver"
                          
    call pastix_fortran(ptss%idata, ptss%comm, ptss%nblock, ptss%jcn, ptss%irn, ptss%val, &
                        ptss%perm_vars, ptss%iperm_vars, ptss%rhs_val, int1, ptss%iparm, ptss%dparm)
  
    ptss%iparm(IPARM_VERBOSE)               = ptss%verb              
    ptss%iparm(IPARM_ITERMAX)               = ptss%iter
  
    ptss%iparm(IPARM_FACTORIZATION)         = ptss%facto
    ptss%iparm(IPARM_INCOMPLETE)            = ptss%ricar
    ptss%iparm(IPARM_LEVEL_OF_FILL)         = ptss%iluk
    ptss%dparm(DPARM_EPSILON_REFINEMENT)    = ptss%eps
    ptss%dparm(DPARM_EPSILON_MAGN_CTRL)     = ptss%pivot
    ptss%iparm(IPARM_RHS_MAKING)            = ptss%rhs
    ptss%iparm(IPARM_SYM)                   = ptss%sym
    ptss%iparm(IPARM_AMALGAMATION_LEVEL)    = ptss%amalg
    
#ifdef FUNNELED
    ptss%iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_FUNNELED
#else
    ptss%iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_MULTIPLE
#endif
  
    ptss%initialized = .true.      
    
    return
    
  end subroutine pastix_initialize
  
  subroutine pastix_finalize(ptss)
    implicit none
    
    type(type_PASTIX_SOLVER)          :: ptss
 
    
    ptss%iparm(IPARM_START_TASK) = API_TASK_CLEAN
    ptss%iparm(IPARM_END_TASK)   = API_TASK_CLEAN
    
    call pastix_fortran(ptss%idata, ptss%comm, ptss%nblock, ptss%jcn, ptss%irn, ptss%val, &
                        ptss%perm_vars, ptss%iperm_vars, ptss%rhs_val, int1, ptss%iparm, ptss%dparm)

    if (associated(ptss%perm_vars)) deallocate(ptss%perm_vars)
    if (associated(ptss%iperm_vars)) deallocate(ptss%iperm_vars)
    
    ptss%perm_vars  => Null()
    ptss%iperm_vars => Null()
    
    !if (associated(ptss%irn)) deallocate(ptss%irn)
    !if (associated(ptss%jcn)) deallocate(ptss%jcn)
    !if (associated(ptss%val)) deallocate(ptss%val)
    !if (associated(ptss%rhs_val)) deallocate(ptss%rhs_val)
    
    ptss%irn => Null()
    ptss%jcn => Null()
    ptss%val => Null()
    ptss%rhs_val => Null()
    
    ptss%idata = 0
    ptss%comm  = 0
    
    ptss%scaled      = .false.
    ptss%refine      = .false.
    ptss%equilibrium = .false.    
                        
    ptss%initialized = .false.
    ptss%analyzed    = .false.

    
    return

  end subroutine pastix_finalize

  
  subroutine pastix_analyze(ptss)
    use mpi_mod
    use mod_clock
  
    implicit none
    
    type(type_PASTIX_SOLVER)          :: ptss
    integer                           :: my_id, ierr
    type(clcktype)                    :: t_itstart, t0, t1, t2, t3
    real*8                            :: tsecond
    
    call MPI_COMM_RANK(ptss%comm, my_id, ierr)
    
    if (my_id .eq. 0) write(*,*) "PaStiX solver: analyzing matrix"
    
    ptss%iparm(IPARM_START_TASK) = API_TASK_ORDERING
    ptss%iparm(IPARM_END_TASK)   = API_TASK_ANALYSE
  
    call clck_time(t0)
  
    call pastix_fortran(ptss%idata, ptss%comm, ptss%nblock, ptss%jcn, ptss%irn, ptss%val, &
                        ptss%perm_vars, ptss%iperm_vars, ptss%rhs_val, int1, ptss%iparm, ptss%dparm)
                        
    call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time analysis :', tsecond
  
    ptss%analyzed = .true.
    
    return
    
  end subroutine pastix_analyze
  
  subroutine pastix_factorize(ptss)
    use mpi_mod
    use mod_clock
  
    implicit none
    
    type(type_PASTIX_SOLVER)          :: ptss
    integer                           :: my_id, ierr
    type(clcktype)                    :: t_itstart, t0, t1, t2, t3
    real*8                            :: tsecond
    
    call MPI_COMM_RANK(ptss%comm, my_id, ierr)

    call clck_time(t0)
    
    ptss%iparm(IPARM_START_TASK) = API_TASK_NUMFACT
    ptss%iparm(IPARM_END_TASK)   = API_TASK_NUMFACT
    
    call pastix_fortran(ptss%idata, ptss%comm, ptss%nblock, ptss%jcn, ptss%irn, ptss%val, &
                        ptss%perm_vars, ptss%iperm_vars, ptss%rhs_val, int1, ptss%iparm, ptss%dparm)
                        
    call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time factorize:', tsecond
    
    return
    
  end subroutine pastix_factorize
  
  subroutine pastix_solve(ptss,rhs_vec)
    use mpi_mod
    use data_structure, only: type_RHS
    use mod_clock
  
    implicit none
    
    type(type_PASTIX_SOLVER)          :: ptss
    type(type_RHS)                    :: rhs_vec    
    integer                           :: my_id, ierr
    type(clcktype)                    :: t_itstart, t0, t1, t2, t3
    real*8                            :: tsecond
    integer(kind=int_all)             :: i
    
    ! dummy arguments to be used with pastix interface
    integer(kind=int_all), pointer    :: irn(:), jcn(:)
    real(kind=8), pointer             :: val(:)    
    
    call MPI_COMM_RANK(ptss%comm, my_id, ierr)

    call clck_time(t0)
  
    ptss%iparm(IPARM_START_TASK) = API_TASK_SOLVE
    ptss%iparm(IPARM_END_TASK)   = API_TASK_SOLVE
    if (ptss%refine) ptss%iparm(IPARM_END_TASK) = API_TASK_REFINE
 
    call pastix_fortran(ptss%idata, ptss%comm, ptss%nblock, ptss%jcn, ptss%irn, ptss%val, &
                        ptss%perm_vars, ptss%iperm_vars, rhs_vec%val, int1, ptss%iparm, ptss%dparm)
                      
    if (ptss%scaled) then
      do i=1,rhs_vec%n
        rhs_vec%val(i) =  rhs_vec%val(i)/ptss%solution_scaling(i)
      enddo
    endif
  
    call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time solve:', tsecond
   
    return
    
  end subroutine pastix_solve
  

  subroutine pastix_init_nthreads(ptss)
    use mpi_mod
!$  use omp_lib    
    
    implicit none
    type(type_PASTIX_SOLVER) :: ptss
    integer :: nthrd = 1
    integer :: n_cpu, ierr
    

    call MPI_COMM_SIZE(ptss%comm, n_cpu, ierr)

!$omp parallel default(none) shared(nthrd)
!$omp master
!$      nthrd = omp_get_num_threads()    
!$omp end master
!$omp end parallel
    
    if (nthrd*n_cpu > ptss%maxthrd) then
      nthrd = max(ptss%maxthrd/n_cpu, 1)
    endif
    ptss%iparm(IPARM_THREAD_NBR) = nthrd
    
    return
  end subroutine pastix_init_nthreads

#endif
end module mod_pastix

