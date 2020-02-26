module spk_module
#ifdef USE_STRUMPACK        
     
  use iso_c_binding
  use mpi

  implicit none
  type(C_PTR) :: spss ! STRUMPACK sparse solver
  logical :: spss_initialized, spss_analyzed

  private
  public :: f2spk, f2spk_finalize, getptr, fprintmem, spss_initialized, spss_analyzed

  interface
    subroutine spk(n,nnz,irn,jcn,val,rhsc,spss,comm,phase) bind(C)
      use iso_c_binding
      use mpi            
      implicit none
      
      integer(kind=C_INT), dimension(:), pointer, intent(in) :: irn,jcn
      real(kind=C_DOUBLE), dimension(:), pointer, intent(in) :: val
      integer(kind=C_INT), intent(in) :: n, phase
      integer(kind=C_INT), intent(inout) :: nnz      
      integer, intent(in) :: comm
      type(c_ptr), intent(inout) :: spss, rhsc
    end subroutine spk
    
    subroutine spk_finalize(spss,comm) bind(C)
      use iso_c_binding
      use mpi            
      implicit none
      
      type(c_ptr), intent(inout) :: spss
      integer, intent(in) :: comm
    end subroutine spk_finalize
    
    subroutine printptr(cptr) bind(C)
      use iso_c_binding
      implicit none
      type(C_PTR) cptr
    end subroutine printptr

    subroutine printmem(rank,cstr) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=C_INT) :: rank
      type(C_PTR) :: cstr
    end subroutine printmem    

  end interface  

  contains
    subroutine f2spk(n,nnz,irn,jcn,val,rhs,comm,phase) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none

        integer comm,ierr
        integer(kind=C_INT), intent(in) :: phase
        integer(kind=C_INT), dimension(:), pointer, intent(in) :: irn,jcn
        real(kind=C_DOUBLE),  dimension(:), pointer, intent(in) :: val
        real(kind=C_DOUBLE),  dimension(:), pointer, intent(inout) :: rhs        
        integer(kind=C_INT), intent(in) :: n
        integer(kind=C_INT), intent(inout) :: nnz           
        type(C_PTR) :: rhsc

        rhsc = c_loc(rhs);
        call spk(n,nnz,irn,jcn,val,rhsc,spss,comm,phase)
        call MPI_Barrier(comm,ierr)

        return          
    end subroutine f2spk
    
    subroutine f2spk_finalize(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none
        
        integer comm,ierr
        
        call spk_finalize(spss,comm)
        call MPI_Barrier(comm,ierr);
        
        return
    end subroutine f2spk_finalize    

    subroutine getptr() bind(C)
        use, intrinsic :: iso_c_binding

        call printptr(spss)
        return
    end subroutine getptr

    subroutine fprintmem(rank,msg)
        use iso_c_binding
        implicit none
        integer(kind=C_INT) :: rank
        character(len=24), target :: msg
        character(len=24), pointer :: fstr
        type(C_PTR) :: cstr
        integer :: i, mlen
        
        mlen = len(trim(msg))
        allocate(fstr)
        fstr=>msg(1:mlen)
        cstr = c_loc(fstr)
        call printmem(rank,cstr)
        deallocate(fstr)

        return
    end subroutine fprintmem

#endif
end module spk_module

