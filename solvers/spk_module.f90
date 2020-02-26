module spk_module
#ifdef USE_STRUMPACK        
     
  use iso_c_binding
  use mpi

  implicit none
  type(C_PTR) :: spss ! STRUMPACK sparse solver
  logical :: spss_initialized, spss_analyzed

  private
  public :: f2spk_init, f2spk_set_mat, f2spk_fact, f2spk_solve, f2spk_finalize
  public :: spss_initialized, spss_analyzed

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
    
    subroutine spk_set_mat(n,nnz,irn,jcn,val,spss,comm) bind(C)
      use iso_c_binding
      use mpi            
      implicit none
      
      integer(kind=C_INT), dimension(:), pointer, intent(in) :: irn,jcn
      real(kind=C_DOUBLE), dimension(:), pointer, intent(in) :: val
      integer(kind=C_INT), intent(in) :: n
      integer(kind=C_INT), intent(inout) :: nnz      
      integer, intent(in) :: comm
      type(c_ptr), intent(inout) :: spss
    end subroutine spk_set_mat
    
    subroutine spk_fact(spss,comm) bind(C)
      use iso_c_binding
      use mpi            
      implicit none
      
      type(c_ptr), intent(inout) :: spss
      integer, intent(in) :: comm
    end subroutine spk_fact
    
    subroutine spk_solve(n,rhsc,spss,comm) bind(C)
      use iso_c_binding
      use mpi            
      implicit none
      
      integer(kind=C_INT), intent(in) :: n
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

  end interface  

  contains
    subroutine f2spk_init(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none
        
        integer comm,ierr
        
        call spk_init(spss,comm)
        call MPI_Barrier(comm,ierr)
        
        return
    end subroutine f2spk_init
    
    subroutine f2spk_set_mat(n,nnz,irn,jcn,val,comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none

        integer comm,ierr
        integer(kind=C_INT), dimension(:), pointer :: irn,jcn
        real(kind=C_DOUBLE),  dimension(:), pointer :: val
        integer(kind=C_INT), intent(in) :: n
        integer(kind=C_INT), intent(inout) :: nnz           

        call spk_set_mat(n,nnz,irn,jcn,val,spss,comm)
        call MPI_Barrier(comm,ierr)

        return  
    end subroutine f2spk_set_mat
    
    subroutine f2spk_fact(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none
        
        integer comm,ierr
        
        call spk_fact(spss,comm)
        call MPI_Barrier(comm,ierr)
        
        return
    end subroutine f2spk_fact
    
    subroutine f2spk_solve(n,rhs,comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none

        integer comm,ierr
        real(kind=C_DOUBLE),  dimension(:), pointer :: rhs        
        integer(kind=C_INT), intent(in) :: n
        type(C_PTR) :: rhsc

        rhsc = c_loc(rhs);        

        call spk_solve(n,rhsc,spss,comm)
        call MPI_Barrier(comm,ierr)

        return  
    end subroutine f2spk_solve    
    
    subroutine f2spk_finalize(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none
        
        integer comm,ierr
        
        call spk_finalize(spss,comm)
        call MPI_Barrier(comm,ierr);
        
        return
    end subroutine f2spk_finalize    

#endif
end module spk_module

