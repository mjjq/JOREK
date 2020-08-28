module sorting_module
  implicit none
  private
  public remove_duplicates

contains
  recursive subroutine sort(temp, start, finish, list)
    implicit none
    integer(kind=8), intent(inout) :: list(:), temp(:)
    integer, intent(inout) :: start
    integer, intent(in) :: finish
    integer :: middle
    if (finish-start<2) then
      return
    else
      middle = (start + finish)/2
      call sort(list, start, middle, temp)
      call sort(list, middle, finish, temp)
      call merge(temp, start, middle, finish, list)
    endif
  end subroutine sort

  subroutine merge(list, start, middle, finish, temp)
    implicit none
    integer(kind=8),intent(in) :: list(:)
    integer(kind=8),intent(inout) :: temp(:)
    integer,intent(in) ::start, middle, finish
    integer :: i, i1, i2
    i1 = start
    i2 = middle
    do i=start, finish-1
      if (i1.lt.middle.and.(i2.ge.finish.or.list(i1).le.list(i2))) then
        temp(i)=list(i1)
        i1 = i1 + 1
      else
        temp(i)=list(i2)
        i2 = i2 + 1
      endif
    enddo
  end subroutine merge

  function unique_sorted(list)
    implicit none
    integer :: start, finish, n
    integer(kind=8), intent(inout) :: list(:)
    integer(kind=8), allocatable :: unique_sorted(:), work(:)
    logical,allocatable :: duplicates(:)
    ! sorting
    work=list
    start=1
    n=size(list)
    finish=n+1
    call sort(work,start,finish,list)
    ! removing duplicates
    allocate(duplicates(n))
    duplicates=.false.
    !removes duplicates and zeros
    duplicates(1:n-1)=list(1:n-1).eq.list(2:n)
    !duplicates(1:n-1)=((list(1:n-1).eq.list(2:n)).or.(list(1:n).eq.0))
    unique_sorted=pack(list,.not.duplicates)

  end function unique_sorted

  recursive function find_index(list,low,high,x) result(idx)
    integer, intent(in) :: low, high
    integer(kind=8), intent(in) :: list(:), x
    integer :: mid
    integer(kind=8) :: idx

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
    integer(kind=8), allocatable :: ij(:), ij_new(:), new_ind(:)
    integer(kind=8) :: dum, i1, i2, i3
    integer :: i, j, nnz_new

    allocate(ij(nnz), new_ind(nnz))
    do i = 1, nnz
      i1 = int(irn(i)-1,kind=8)
      i2 = int(n,kind=8)
      i3 = int(jcn(i),kind=8)
      ij(i) = i1*i2 + i3
    enddo

    ij_new=unique_sorted(ij)
    deallocate(ij)
    nnz_new = size(ij_new)
    if (nnz.ne.nnz_new) write(*,*) "Number of nnz changed: nnz_old, nnz_new = ", nnz, nnz_new
    
    ! find index of original element in the new (ordered) list
    do i = 1, nnz
      i1 = int(irn(i)-1,kind=8)
      i2 = int(n,kind=8)
      i3 = int(jcn(i),kind=8)
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
