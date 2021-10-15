module mod_pastix
#ifdef WITH_PASTIX62

  use iso_c_binding
  use mpi

  implicit none
  
  type(C_PTR) :: spm, pastix_data ! sparse solver (distributed)
  type(C_PTR) :: iparm, dparm

  integer(kind=C_INT), dimension(:), pointer :: loc2glob, glob2loc  ! mapping for column distribution
  logical :: spm_initialized, spm_analyzed

  private
  public :: pastix_init, pastix_set_mat, pastix_analyze, &
            pastix_factorize, pastix_solve, get_residual, &
            pastix_finalize, test_reallocate_array, spm_initialized, spm_analyzed

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

    subroutine ptx_set_mat(spm, indx, n, nnz, n_d, nnz_d, rptr, cptr, values, &
               loc2glob, glob2loc, comm, update, check) bind(C)

      use iso_c_binding
      implicit none

      type(C_PTR), intent(inout) :: spm
      integer(kind=C_INT) :: indx, comm
      integer(kind=C_INT) :: n, n_d, nnz, nnz_d
      
      !type(C_PTR), intent(inout) :: rptr, cptr, values, loc2glob

      integer(kind=C_INT), dimension(:), pointer :: rptr, cptr, loc2glob, glob2loc
      real(kind=C_DOUBLE), dimension(:), pointer :: values
      logical :: update, check

    end subroutine ptx_set_mat

    subroutine ptx_analyze(pastix_data, spm) bind(C)

      use iso_c_binding
      implicit none

      type(C_PTR), intent(inout) :: spm, pastix_data

    end subroutine ptx_analyze

    subroutine ptx_factorize(pastix_data, spm) bind(C)

      use iso_c_binding
      implicit none

      type(C_PTR), intent(inout) :: spm, pastix_data

    end subroutine ptx_factorize

    subroutine ptx_solve(pastix_data, spm, rhsc, refine) bind(C)

      use iso_c_binding
      implicit none

      type(C_PTR), intent(inout) :: spm, pastix_data
      type(C_PTR), intent(inout) :: rhsc
      logical :: refine

    end subroutine ptx_solve
    
    subroutine ptx_finalize(pastix_data, spm, iparm, dparm) bind(C)

      use iso_c_binding
      implicit none

      type(C_PTR), intent(inout) :: spm, pastix_data
      type(C_PTR), intent(inout) :: iparm, dparm

    end subroutine ptx_finalize 

    subroutine get_residual(n, nnz, irn, jcn, val, x, rhs, indexing) bind(C)

      use iso_c_binding
      implicit none

      integer(kind=C_INT), dimension(:), pointer, intent(in) :: irn, jcn
      real(kind=C_DOUBLE), dimension(:), pointer, intent(in) :: val, rhs, x
      integer(kind=C_INT) :: n, nnz, indexing

    end subroutine get_residual
    
    !subroutine free_double(val) bind(C)
    !  use iso_c_binding
    !  implicit none
    !  
    !  !real(kind=C_DOUBLE), pointer :: val
    !  type(C_PTR), intent(inout) :: val
    !end subroutine free_double


  end interface

  contains
    subroutine pastix_init(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none

        integer comm,ierr

        call ptx_init(pastix_data,spm,iparm,dparm,comm)
        call MPI_Barrier(comm,ierr)

        return
    end subroutine pastix_init

!    subroutine pastix_set_mat(indexing, n, nnz, n_d, nnz_d, rptr, cptr, values, loc2glob, comm, update) bind(C)
    subroutine pastix_set_mat(n, nnz, irn, jcn, val, block_size, comm, &
                update,distributed,equilibrium) bind(C)

      use iso_c_binding
      use mod_coicsr, only: coicsr
      use sorting_module, only: remove_duplicates, convert2csr
      implicit none

      integer(kind=C_INT), dimension(:), pointer :: irn, jcn
      real(kind=C_DOUBLE), dimension(:), pointer :: val
      
      logical,intent(in),optional :: update, distributed, equilibrium

      integer(kind=C_INT) :: indexing=1, block_size
      integer(kind=C_INT) :: n, n_d, nnz, nnz_d
      
      integer :: my_id_n, ncpu_n, ierr, comm
     
      integer(kind=C_INT), dimension(:), allocatable :: iwk

      logical :: upd, dflag, eql
      upd = .false.
      dflag = .false.
      eql = .false.
      
      if(present(update)) upd = update
      if(present(distributed)) dflag = distributed
      if(present(equilibrium)) eql = equilibrium      
      
      ! if not already distributed distribute matrix column-wise
      ! add n_cpu>1
      if (.not.dflag) then

        call distribute_matrix(indexing, block_size, n, nnz, n_d, nnz_d, jcn, irn, val, loc2glob, glob2loc, comm)
       
      endif

      ! Prepare matrix for using block-structure
      !block_size2 = block_size**2
      !n_block   = n_d/block_size
      !nnz_block = nnz_d/block_size2
      
      !if (block_size>1) then
      !  do i=1,nnz_block  
      !    irn(i) = (irn((i-1)*block_size2+1) - 1)/block_size + 1 
      !    jcn(i) = (jcn((i-1)*block_size2+1) - 1)/block_size + 1 
      !  enddo
      !endif

      !if (allocated(iwk)) deallocate(iwk)
      !allocate(iwk(n_block+1))
      !call coicsr2(n_block,nnz_block,val,irn(1:nnz_block),jcn(1:nnz_block),block_size,iwk)
      !deallocate(iwk)
     
      ! allocating local arrays to be used to store matrix in csc format
      !allocate(rptr(nnz_d),cptr(n_d+1),values(nnz_d))

      if (eql) then
#if (defined(USEMKL))
        call convert2csr(indexing, n_d, n, nnz_d, jcn, irn, val)
#else
        call remove_duplicates(n, nnz_d, jcn, irn, val)
        call convert2csr(indexing, n_d, n, nnz_d, jcn, irn, val)
#endif
      else
        if (allocated(iwk)) deallocate(iwk)
        allocate(iwk(n+1))
        call coicsr(n,nnz_d,1,val,irn,jcn,iwk)
        deallocate(iwk)
      endif

      call MPI_Barrier(comm,ierr)      
     
      call ptx_set_mat(spm, indexing, n, nnz, n_d, nnz_d, irn, jcn, val, loc2glob, glob2loc, comm, upd, eql)      
      
      return
    
    end subroutine pastix_set_mat

    subroutine pastix_analyze() bind(C)

      use iso_c_binding
      implicit none

      call ptx_analyze(pastix_data, spm)

    end subroutine pastix_analyze

    subroutine pastix_factorize() bind(C)

      use iso_c_binding
      implicit none

      call ptx_factorize(pastix_data, spm)

    end subroutine pastix_factorize

    subroutine pastix_solve(rhs, refine) bind(C)

      use iso_c_binding
      implicit none

      real(kind=C_DOUBLE), pointer, intent(inout) :: rhs(:)
      logical, intent(in), optional :: refine
      type(C_PTR) :: rhsc
      
      logical :: ref
      ref = .false.
      
      if(present(refine)) ref = refine

      rhsc = c_loc(rhs)

      call ptx_solve(pastix_data, spm, rhsc, ref)

    end subroutine pastix_solve
    
    subroutine pastix_finalize() bind(C)

      use iso_c_binding
      implicit none

      call ptx_finalize(pastix_data, spm,  iparm, dparm)
      deallocate(loc2glob, glob2loc)      

    end subroutine pastix_finalize
    
    subroutine distribute_matrix(indexing,block_size,n,nnz,n_d,nnz_d,irn,jcn,val,loc2glob,glob2loc,comm)
    ! distibute by the first array (irn or jcn)
      implicit none
    
      integer(kind=C_INT), dimension(:), pointer :: irn, jcn
      real(kind=C_DOUBLE), dimension(:), pointer :: val
      integer(kind=C_INT), dimension(:), pointer :: loc2glob, glob2loc
    
      
      integer, intent(in):: indexing, block_size
      integer, intent(in) :: n, nnz
      integer, intent(out) :: n_d, nnz_d
      
      integer(kind=C_INT), dimension(:), pointer :: irn_d, jcn_d
      real(kind=C_DOUBLE), dimension(:), pointer :: val_d
      integer :: ncpu_n, my_id_n, comm, ierr
      integer :: i, j, indx
      
      integer(kind=C_INT), allocatable :: dist(:), myelm(:)
      
      call MPI_Comm_rank(comm, my_id_n, ierr)
      call MPI_Comm_size(comm, ncpu_n, ierr)  
    
      if (allocated(dist)) deallocate(dist)
      allocate(dist(ncpu_n+1))
      
      ! number of rows/columns per cpu with the last one getting extra
      dist(2:ncpu_n+1) = block_size*((n/block_size)/ncpu_n)
      dist(ncpu_n+1) = dist(ncpu_n+1) + (n - sum(dist(2:ncpu_n+1)))
      
      dist(1) = indexing
      do i=2, ncpu_n+1
        dist(i)= dist(i) + dist(i-1)
      enddo      

      allocate(myelm(nnz))
      j = 1
      do i=1, nnz
        if ((irn(i)>= dist(my_id_n+1)).and.(irn(i)<=(dist(my_id_n+2)-1))) then
          myelm(j) = i
          j = j + 1
        endif
      enddo

      nnz_d = j - 1
      n_d = dist(my_id_n+2) - dist(my_id_n+1)
      
      allocate(irn_d(nnz_d),jcn_d(nnz_d),val_d(nnz_d),loc2glob(n_d))

      do i = 1, nnz_d
        irn_d(i) = irn(myelm(i)) - dist(my_id_n+1) + indexing
        jcn_d(i) = jcn(myelm(i))
        val_d(i) = val(myelm(i))
      enddo
      
      do i = 1, n_d
        loc2glob(i) = i - 1 + dist(my_id_n + 1);
      enddo
            
      deallocate(myelm,dist)
      deallocate(irn,jcn,val)

      irn => irn_d
      jcn => jcn_d
      val => val_d
      
      allocate(glob2loc(n))
      glob2loc = 0
      
      if (indexing == 0) then
        indx = 1
      else
        indx = 0
      endif
      
      do i = 1, n_d
        glob2loc(loc2glob(i) + indx)= - my_id_n - 1;
      enddo
      
      call MPI_Allreduce(MPI_IN_PLACE, glob2loc, n, MPI_INTEGER, MPI_SUM, comm, ierr)
      
      do i = 1, n_d
        glob2loc(loc2glob(i) + indx)= loc2glob(i) - loc2glob(1)
      enddo
      
      return
      
    end subroutine distribute_matrix
    
  subroutine test_reallocate_array(n, irn, jcn, val, dum) bind(C)
  
    implicit none
    
    integer(kind=C_INT), dimension(:), pointer :: irn, jcn, irn_d, jcn_d
    real(kind=C_DOUBLE), dimension(:), pointer :: val, val_d
    integer :: n, dum
    
    integer :: ncpu_n, my_id_n, comm, ierr
    
    allocate(irn_d(n),jcn_d(n)); allocate(val_d(n))
    irn_d = dum; jcn_d = dum; val_d = dum
    
    deallocate(irn); deallocate(jcn); deallocate(val)
    
    irn => irn_d
    jcn => jcn_d
    val => val_d
    
    return  
  end subroutine test_reallocate_array

#endif
end module mod_pastix

