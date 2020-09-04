module sorting_module
  use iso_c_binding
  implicit none
  private
  public remove_duplicates

#define INTSIZE 8
#define CINT c_int64_t

interface
  subroutine qsort(array,elem_count,elem_size,compare) bind(C,name="qsort")
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

  function unique_sorted(list)
    implicit none
    integer :: n
    integer(kind=INTSIZE), intent(inout), target :: list(:)
    integer(kind=INTSIZE), allocatable :: unique_sorted(:)
    logical,allocatable :: duplicates(:)
    integer(c_size_t) l,isize

    n=size(list)
    
    l = n; isize = INTSIZE
    call qsort(c_loc(list(1)),l,isize,c_funloc(compar))

    ! removing duplicates
    allocate(duplicates(n))
    duplicates=.false.
    !removes duplicates and zeros
    duplicates(1:n)=list(1:n-1).eq.list(2:n)
    !duplicates(1:n-1)=((list(1:n-1).eq.list(2:n)).or.(list(1:n).eq.0))
    unique_sorted = pack(list,.not.duplicates)

  end function unique_sorted

  recursive function find_index(list,low,high,x) result(idx)
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
  ! sort and remove duplicates from sparse matrix
    use, intrinsic :: iso_c_binding

    integer, intent(in) :: n
    integer, intent(inout) :: nnz
    integer(kind=C_INT), dimension(:), pointer :: irn, jcn
    real(kind=C_DOUBLE), dimension(:), pointer :: val

    real(kind=C_DOUBLE), allocatable :: val_new(:)
    ! long integer is required for 1d representation of coordinate index
    integer(kind=INTSIZE), allocatable :: ij(:), ij_new(:), new_ind(:)
    integer(kind=INTSIZE) :: dum, i1, i2, i3
    integer :: i, j, nnz_new

    allocate(ij(nnz), new_ind(nnz))
    do i = 1, nnz
      i1 = int(irn(i)-1,kind=INTSIZE)
      i2 = int(n,kind=INTSIZE)
      i3 = int(jcn(i),kind=INTSIZE)
      ij(i) = i1*i2 + i3
    enddo

    ij_new=unique_sorted(ij)
    deallocate(ij)
    nnz_new = size(ij_new)
    if (nnz.ne.nnz_new) write(*,*) "Number of nnz changed: nnz_old, nnz_new = ", nnz, nnz_new
    
    ! find index of original element in the new (ordered) list
    do i = 1, nnz
      i1 = int(irn(i)-1,kind=INTSIZE)
      i2 = int(n,kind=INTSIZE)
      i3 = int(jcn(i),kind=INTSIZE)
      dum = i1*i2 + i3
      new_ind(i) = find_index(ij_new,1,nnz_new,dum)
    enddo

    allocate(val_new(nnz_new)); val_new = 0.0
    do i =1, nnz
      val_new(new_ind(i)) = val_new(new_ind(i)) + val(i)
    enddo

    do i = 1, nnz_new
      irn(i) = int((ij_new(i)-1)/n) + 1
      jcn(i) = mod(ij_new(i)-1,n) + 1
      val(i) = val_new(i)
    enddo

    deallocate(ij_new,new_ind,val_new)
    nnz = nnz_new

  end subroutine remove_duplicates

end module sorting_module
