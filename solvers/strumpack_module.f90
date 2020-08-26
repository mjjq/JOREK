module strumpack_module
#ifdef USE_STRUMPACK        
     
  use iso_c_binding
  use mpi

  implicit none
  type(C_PTR) :: spss ! STRUMPACK sparse solver
  integer(kind=C_INT), dimension(:), pointer :: dist  ! row distribution
  logical :: spss_initialized, spss_analyzed

  private
  public :: strumpack_init, strumpack_set_mat, strumpack_analyze, &
            strumpack_factorize, strumpack_solve, strumpack_finalize, &
            spss_initialized, spss_analyzed

  interface
    subroutine spk() bind(C)
      use iso_c_binding
    end subroutine spk
    
    subroutine spk_init(spss,comm) bind(C)
      use iso_c_binding
      use mpi            
      implicit none
      
      type(c_ptr), intent(inout) :: spss
      integer, intent(in) :: comm
    end subroutine spk_init
    
    subroutine spk_set_mat(n,dist,irn,jcn,val,spss,comm,upd) bind(C)
      use iso_c_binding
      use mpi            
      implicit none
      
      integer(kind=C_INT), dimension(:), pointer, intent(in) :: irn,jcn,dist
      real(kind=C_DOUBLE), dimension(:), pointer, intent(in) :: val
      integer(kind=C_INT), intent(in) :: n
      integer, intent(in) :: comm
      type(c_ptr), intent(inout) :: spss
      logical :: upd
    end subroutine spk_set_mat    

    subroutine spk_reord(spss,comm) bind(C)
      use iso_c_binding
      use mpi            
      implicit none
      
      type(c_ptr), intent(inout) :: spss
      integer, intent(in) :: comm
    end subroutine spk_reord    
    
    subroutine spk_fact(spss,comm) bind(C)
      use iso_c_binding
      use mpi            
      implicit none
      
      type(c_ptr), intent(inout) :: spss
      integer, intent(in) :: comm
    end subroutine spk_fact
    
    subroutine spk_solve(n,dist,rhsc,spss,comm) bind(C)
      use iso_c_binding
      use mpi            
      implicit none
      
      integer(kind=C_INT), intent(in) :: n
      integer(kind=C_INT), dimension(:), pointer, intent(in) :: dist
      type(c_ptr), intent(inout) :: spss, rhsc      
      integer, intent(in) :: comm
    end subroutine spk_solve    
    
    subroutine spk_finalize(spss,comm) bind(C)
      use iso_c_binding
      use mpi            
      implicit none
      
      type(c_ptr), intent(inout) :: spss
      integer, intent(in) :: comm
    end subroutine spk_finalize

    subroutine convert2csr(indx, n, m, nnz, irn, jcn, val) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=C_INT), dimension(:), pointer, intent(in) :: irn,jcn
      real(kind=C_DOUBLE), dimension(:), pointer, intent(in) :: val
      integer(kind=C_INT), intent(in) :: n, m, indx
      integer(kind=C_INT), intent(inout) :: nnz       
    end subroutine convert2csr

  end interface

  contains
    subroutine strumpack_init(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none
        
        integer comm,ierr
        
        call spk_init(spss,comm)
        call MPI_Barrier(comm,ierr)
        
        return
    end subroutine strumpack_init
    
    subroutine strumpack_set_mat(n,nnz,irn,jcn,val,comm,update,distributed) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        use sorting_module, only : remove_duplicates
        implicit none

        integer comm,ierr
        integer(kind=C_INT), dimension(:), pointer :: irn, jcn, irnl, jcnl
        real(kind=C_DOUBLE),  dimension(:), pointer :: val, vall
        integer(kind=C_INT), intent(in) :: n
        integer(kind=C_INT), intent(inout) :: nnz 
        logical,intent(in),optional :: update, distributed

        integer(kind=C_INT), dimension(:), pointer :: myelm
        logical :: upd=.false., dflag=.false.
        
        integer :: rank, ncpu, nnzloc, nloc, i, j, indx=1

        if(present(update)) upd = update
        if(present(distributed)) dflag = distributed

        call MPI_COMM_RANK(comm, rank, ierr)
        call MPI_COMM_SIZE(comm, ncpu, ierr)

        if ((minval(irn).lt.1).or.(maxval(irn).gt.n)) then
          write(*,*) rank, ": Error: inconsistent row indices", minval(irn), maxval(irn)
          !call exit(0)
        endif

        if ((.not. dflag).and.(ncpu>1)) then

          call distribute_rows(n,ncpu)

          allocate(myelm(nnz))
          j = 1
          do i=1, nnz
            if ((irn(i)>= dist(rank+1)).and.(irn(i)<= (dist(rank+2)-1))) then
              myelm(j) = i
              j = j + 1
            endif
          enddo
        
          nnzloc = j - 1 
          nloc = dist(rank+2) - dist(rank+1)

          allocate(irnl(nnzloc), jcnl(nnzloc), vall(nnzloc))

          do i = 1, nnzloc
            irnl(i) = irn(myelm(i)) - dist(rank+1) + 1       ! irn starts from 1
            jcnl(i) = jcn(myelm(i))                          ! jcn remains the same
            vall(i) = val(myelm(i))
          enddo
#ifndef USE_MKL
          call remove_duplicates(n,nnzloc,irnl,jcnl,vall)
#endif
          call convert2csr(indx,nloc,n,nnzloc,irnl,jcnl,vall)
          dist(:) = dist(:) - indx                           ! convert to c-indexing
          call spk_set_mat(nloc,dist,irnl,jcnl,vall,spss,comm,upd)
          deallocate(myelm,irnl,jcnl,vall)

        elseif (dflag.and.(ncpu>1)) then
           ! get row distribution from irn

          if (associated(dist)) dist=>null()
          allocate(dist(ncpu+1))
          dist(:) = 0
          dist(rank+1) = minval(irn)

          if (rank.eq.(ncpu-1)) dist(rank+2) = maxval(irn) + 1

          call MPI_Allreduce(MPI_IN_PLACE,dist,ncpu+1,MPI_INTEGER,MPI_SUM,comm,ierr)

          ! check for consistency
          ierr = 0
          if ((dist(1).ne.1)) ierr = 1
          do i = 2, ncpu+1
            if (.not.(dist(i)>dist(i-1))) ierr = 1
          enddo
         
          if (ierr.ne.0) then
            write(*,*) "Error in harmonic matrix distribution"
            call exit(2)
          endif

          nloc = dist(rank+2) - dist(rank+1)
          irn(:) = irn(:) - dist(rank+1) + 1 ! irn starts from 1
#ifndef USE_MKL
          call remove_duplicates(n,nnz,irn,jcn,val)
#endif
          call convert2csr(indx,nloc,n,nnz,irn,jcn,val)          
          dist(:) = dist(:) - indx
        
          call spk_set_mat(nloc,dist,irn,jcn,val,spss,comm,upd)

        else

          call distribute_rows(n,1)
          dist(:) = dist(:) - indx
#ifndef USE_MKL
          call remove_duplicates(n,nnz,irn,jcn,val)
#endif
          call convert2csr(indx,n,n,nnz,irn,jcn,val)
          call spk_set_mat(n,dist,irn,jcn,val,spss,comm,upd)

        endif

        call MPI_Barrier(comm,ierr)

        return  
    end subroutine strumpack_set_mat

    subroutine strumpack_analyze(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none
        
        integer comm,ierr
        
        call spk_reord(spss,comm)
        call MPI_Barrier(comm,ierr)
        
        return
    end subroutine strumpack_analyze

    subroutine strumpack_factorize(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none
        
        integer comm,ierr
        
        call spk_fact(spss,comm)
        call MPI_Barrier(comm,ierr)
        
        return
    end subroutine strumpack_factorize
    
    subroutine strumpack_solve(n,rhs,comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none

        real(kind=C_DOUBLE),  dimension(:), pointer :: rhs        
        integer(kind=C_INT), intent(in) :: n
        integer, intent(in) :: comm 
        
        integer :: rank, ncpu, ierr, nloc
        type(C_PTR) :: rhsc

        call MPI_COMM_RANK(comm, rank, ierr)
        call MPI_COMM_SIZE(comm, ncpu, ierr)

        rhsc = c_loc(rhs);        

        call spk_solve(n,dist,rhsc,spss,comm)
        call MPI_Barrier(comm,ierr)

        return  
    end subroutine strumpack_solve      
    
    subroutine strumpack_finalize(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none
        
        integer comm,ierr
        
        call spk_finalize(spss,comm)
        call MPI_Barrier(comm,ierr);
        
        return
    end subroutine strumpack_finalize

    subroutine distribute_rows(n,ncpu) bind(C)
      !> Distribute rows between members of MPI group

        use, intrinsic :: iso_c_binding
        use mpi
        implicit none

        integer, intent(in) :: n, ncpu
        integer(kind=C_INT), dimension(:), allocatable :: nl
        integer :: ierr, i

        if (associated(dist)) dist=>null()
        allocate(dist(ncpu+1))
        allocate(nl(ncpu))

        do i=1, ncpu
          nl(i) = n/ncpu
          if (i<mod(n,ncpu)+1) nl(i) = nl(i) + 1
        enddo

        dist(1) = 1
        do i=1, ncpu
          dist(i+1)= dist(i) + nl(i)
        enddo

        deallocate(nl)

        return
    end subroutine distribute_rows

#endif
end module strumpack_module

