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
  public :: type_PASTIX_SOLVER, pastix_initialize, pastix_set_mat, pastix_analyze, pastix_factorize, pastix_solve, pastix_finalize

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

  subroutine pastix_set_mat(ptss, ad_mat, ac_mat)
    use mpi_mod
    use mod_coicsr
    use data_structure, only: type_SP_MATRIX

    implicit none

    type(type_PASTIX_SOLVER)          :: ptss
    type(type_SP_MATRIX)              :: ad_mat, ac_mat
    integer                           :: n_cpu, my_id, ierr, comm
    integer                           :: i, j
    integer(kind=int_all)             :: k, nnz
    integer*8                         :: check_data
    integer                           :: block_size2
    integer(kind=int_all)             :: nblock, nnz_block
    integer(kind=int_all), allocatable         :: sparskit_work(:)

    comm = ad_mat%comm

    call MPI_COMM_RANK(comm, my_id, ierr)
    call MPI_COMM_SIZE(comm, n_cpu, ierr)

    return

  end subroutine pastix_set_mat


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

!> Prepares matrix for PaStiX5 solver:
!  rescale columns; convert to csc format
  subroutine pastix_set_mat(ptss, ad_mat, ac_mat)
    use mpi_mod
    use mod_clock
    use mod_coicsr
    use data_structure, only: type_SP_MATRIX

    implicit none

    type(type_PASTIX_SOLVER)          :: ptss
    type(type_SP_MATRIX)              :: ad_mat, ac_mat
    type(clcktype)                    :: t_itstart, t0, t1, t2, t3
    real*8                            :: tsecond
    integer                           :: n_cpu, my_id, ierr, comm
    integer                           :: i, j
    integer(kind=int_all)             :: k, nnz
    integer*8                         :: check_data
    integer                           :: block_size2
    integer(kind=int_all)             :: nblock, nnz_block
    integer(kind=int_all), allocatable         :: sparskit_work(:)

    comm = ad_mat%comm

    call MPI_COMM_RANK(comm, my_id, ierr)
    call MPI_COMM_SIZE(comm, n_cpu, ierr)

    if (.not.ptss%equilibrium) then

      call scale_by_cols(ad_mat)
      if (associated(ptss%solution_scaling)) then
        deallocate(ptss%solution_scaling); ptss%solution_scaling => Null()
      endif
      allocate(ptss%solution_scaling(ad_mat%ng))
      do k = 1, ad_mat%ng
        ptss%solution_scaling(k) = ad_mat%column_scaling(k)
      enddo
      ptss%scaled = .true.

    endif

    if (n_cpu>1) then

      call MPI_Allreduce(ad_mat%nnz,ac_mat%nnz,1,MPI_INTEGER_ALL,MPI_SUM,comm,ierr)

      ac_mat%ng = ad_mat%ng
      ac_mat%block_size = ad_mat%block_size
      ac_mat%comm = ad_mat%comm

      allocate(ac_mat%irn(ac_mat%nnz))
      allocate(ac_mat%jcn(ac_mat%nnz))
      allocate(ac_mat%val(ac_mat%nnz))

      call clck_time(t0)

      call split_allgathersolve(n_cpu,my_id,ad_mat,ac_mat)

      call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
      if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time mpi_gather :', tsecond

    else
      ac_mat = ad_mat
    endif

    call clck_time(t0)

    block_size2 = ac_mat%block_size**2

    nblock   = ac_mat%ng/ac_mat%block_size
    nnz_block = ac_mat%nnz/block_size2

    ac_mat%nblock = nblock
    ac_mat%nzblock = nnz_block

    if (ac_mat%block_size > 1) then
      do i = 1,nnz_block
        ac_mat%irn(i) = (ac_mat%irn((i-1)*block_size2+1) - 1)/ac_mat%block_size + 1
        ac_mat%jcn(i) = (ac_mat%jcn((i-1)*block_size2+1) - 1)/ac_mat%block_size + 1
      enddo
    endif

    allocate(sparskit_work(nblock+1))

    call coicsr2(nblock,nnz_block,ac_mat%val,ac_mat%irn(1:nnz_block),ac_mat%jcn(1:nnz_block),ac_mat%block_size,sparskit_work)

    deallocate(sparskit_work)

    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time coicsr :', tsecond

    ! End of matrix preparation
    ptss%irn => ac_mat%irn
    ptss%jcn => ac_mat%jcn
    ptss%val => ac_mat%val
    ptss%nblock = ac_mat%nblock
    ptss%iparm(IPARM_DOF_NBR) = ac_mat%block_size

    if (ptss%equilibrium) then
    ! combine duplicated values
      nnz = ac_mat%jcn(ac_mat%ng + 1) - 1
      call pastix_fortran_checkmatrix(check_data, ac_mat%comm, &
       Int1, ptss%sym, Int1, ac_mat%ng, ac_mat%jcn, ac_mat%irn, ac_mat%val, -Int1, Int1)

      ac_mat%nnz = ac_mat%jcn(ac_mat%ng+1) - 1
      if (ac_mat%nnz < nnz) then
         call pastix_fortran_checkmatrix_end(check_data, Int1, ac_mat%irn, ac_mat%val, Int1)
      endif
    endif

    return

  end subroutine pastix_set_mat

!> Initializes PaStiX5 solver
  subroutine pastix_initialize(ptss)
    use mpi_mod
    use mod_clock

    implicit none

    type(type_PASTIX_SOLVER)          :: ptss
    integer                           :: my_id, ierr
    type(clcktype)                    :: t_itstart, t0, t1, t2, t3
    real*8                            :: tsecond


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

!> Performs matrix analyzis/reordering with PaStiX5 solver
  subroutine pastix_analyze(ptss)
    implicit none

    type(type_PASTIX_SOLVER)          :: ptss

    allocate(ptss%perm_vars(ptss%nblock))
    allocate(ptss%iperm_vars(ptss%nblock))
    ptss%perm_vars(1:ptss%nblock) = 0
    ptss%iperm_vars(1:ptss%nblock) = 0

    ptss%iparm(IPARM_START_TASK) = API_TASK_ORDERING
    ptss%iparm(IPARM_END_TASK)   = API_TASK_ANALYSE

    call pastix_fortran(ptss%idata, ptss%comm, ptss%nblock, ptss%jcn, ptss%irn, ptss%val, &
                        ptss%perm_vars, ptss%iperm_vars, ptss%rhs_val, int1, ptss%iparm, ptss%dparm)

    ptss%analyzed = .true.

    return

  end subroutine pastix_analyze

!> Performs matrix LU factorization with PaStiX5 solver
  subroutine pastix_factorize(ptss)

    implicit none

    type(type_PASTIX_SOLVER)          :: ptss

    ptss%iparm(IPARM_START_TASK) = API_TASK_NUMFACT
    ptss%iparm(IPARM_END_TASK)   = API_TASK_NUMFACT

    call pastix_fortran(ptss%idata, ptss%comm, ptss%nblock, ptss%jcn, ptss%irn, ptss%val, &
                        ptss%perm_vars, ptss%iperm_vars, ptss%rhs_val, int1, ptss%iparm, ptss%dparm)

    return

  end subroutine pastix_factorize

!> Finds the solution with PaStiX5 solver for given RHS
  subroutine pastix_solve(ptss,rhs_vec)
    use data_structure, only: type_RHS

    implicit none

    type(type_PASTIX_SOLVER)          :: ptss
    type(type_RHS)                    :: rhs_vec

    integer(kind=int_all)             :: i
    ! dummy arguments to be used with pastix interface
    integer(kind=int_all), pointer    :: irn(:), jcn(:)
    real(kind=8), pointer             :: val(:)

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

