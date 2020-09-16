module sorting_module
!> Contains subroutines to sort (in continuos ij index) 
!  and remove duplicates from the sparse matrix as needed for STRUMPACK solver

  use iso_c_binding
  implicit none
  private
  public remove_duplicates

#define INTSIZE 8
#define CINT c_int64_t

interface
  subroutine qsort(array,elem_count,elem_size,compare) bind(C,name="qsort")
  !> Interface to C-function qsort
    import
    type(c_ptr),value       :: array
    integer(c_size_t),value :: elem_count
    integer(c_size_t),value :: elem_size
    type(c_funptr),value    :: compare
  end subroutine qsort
end interface  

contains
  integer(c_int) function compar(a, b) bind(C)
    use iso_c_binding
    integer(CINT) a, b

    if ( a .lt. b ) compar = -1
    if ( a .eq. b ) compar = 0
    if ( a .gt. b ) compar = 1
  end function compar

  subroutine unique_sorted(list,n)
  !> Sort and remove duplicates from 1D list
  ! replace list with uniquelly sorted entries; return number of unique elements

    implicit none
    integer,intent(inout) :: n
    integer(kind=INTSIZE), allocatable, intent(inout), target :: list(:)
    logical,allocatable :: duplicates(:)
    integer(c_size_t) l,isize
    integer :: i, j, m

    l = n; isize = INTSIZE
    call qsort(c_loc(list(1)),l,isize,c_funloc(compar))

    ! removing duplicates
    allocate(duplicates(n)); duplicates=.false.
    duplicates(1:n)=list(1:n-1).eq.list(2:n)

    m = count(duplicates)
    if (m.gt.0) then
      j = 1
      do i=1, n
        if (.not.duplicates(i)) then
          list(j) = list(i)
          j = j + 1
        endif
      enddo
      n = j - 1
    endif
    return

  end subroutine unique_sorted

  recursive function find_index(list,low,high,x) result(idx)
  !> Find index of element x in the list
    integer, intent(in) :: low, high
    integer(kind=INTSIZE), intent(in) :: list(:), x
    integer :: mid
    integer(kind=INTSIZE) :: idx

    if (low.gt.high) then
      idx = 0
      write(*,*) "Error in find_index: element not found", x, low, high
      call exit(0)
    endif

    mid = (low + high)/2

    ! target value is found
    if (x .eq. list(mid)) then
      idx = mid

    ! discard all elements in the right search space
    ! including the mid element
    elseif (x .lt. list(mid)) then
      idx = find_index(list, low,  mid - 1, x)

    ! discard all elements in the left search space
    ! including the mid element
    else
      idx = find_index(list, mid + 1, high, x)
    endif

  end function find_index

  subroutine remove_duplicates(n,nnz,irn,jcn,val)
  !> Sort and remove duplicates from sparse matrix
    use, intrinsic :: iso_c_binding

    integer, intent(in) :: n
    integer, intent(inout) :: nnz
    integer(kind=C_INT), dimension(:), pointer :: irn, jcn
    real(kind=C_DOUBLE), dimension(:), pointer :: val

    real(kind=C_DOUBLE), allocatable :: val_new(:)
    ! long integer is required for 1d representation of coordinate index
    integer(kind=INTSIZE), allocatable :: ij(:), ij_new(:), new_ind(:)
    integer(kind=INTSIZE) :: dum, i1, i2, i3
    integer :: i, j, nnz0

    allocate(ij(nnz), new_ind(nnz))
    do i = 1, nnz
      i1 = int(irn(i)-1,kind=INTSIZE)
      i2 = int(n,kind=INTSIZE)
      i3 = int(jcn(i),kind=INTSIZE)
      ij(i) = i1*i2 + i3
    enddo
    
    nnz0 = nnz
    call unique_sorted(ij,nnz)

    if (nnz.ne.nnz0) write(*,*) "Number of nnz changed: nnz_old, nnz_new = ", nnz0, nnz
    
    ! find index of original element in the new (ordered) list
    do i = 1, nnz0
      i1 = int(irn(i)-1,kind=INTSIZE)
      i2 = int(n,kind=INTSIZE)
      i3 = int(jcn(i),kind=INTSIZE)
      dum = i1*i2 + i3
      new_ind(i) = find_index(ij,1,nnz,dum)
    enddo

    allocate(val_new(nnz)); val_new = 0.0
    do i =1, nnz0
      val_new(new_ind(i)) = val_new(new_ind(i)) + val(i)
    enddo

    do i = 1, nnz
      irn(i) = int((ij(i)-1)/n) + 1
      jcn(i) = mod(ij(i)-1,n) + 1
      val(i) = val_new(i)
    enddo

    deallocate(ij,new_ind,val_new)

  end subroutine remove_duplicates

end module sorting_module
