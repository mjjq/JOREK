subroutine reduce(my_id,A,IRN,JCN,NZ,N,IHWB,NZNEW,first_row,last_row)
!--------------------------------------------------------------------------
! subroutine to change format of matrix from MUMPS to compressed column
! format (and a MUMPS with multiple entries added).
!
!   version : bricolage
!   status  : OK
!--------------------------------------------------------------------------
real*8    :: A(*)
integer :: IRN(*),JCN(*),IHWB(*),N,NZ
integer :: my_id
integer, allocatable :: ITOT(:),ISTART(:),IPOS(:),JCOL(:),IRNK(:)
integer, allocatable :: INEW(:),JNEW(:),  IJold(:)
integer :: irn_min, irn_max, first_row , last_row, nrows
integer :: i, irow, index, n_IRNK, index_new, nj, index_old, k, nznew

real*8,allocatable     :: Anew(:)

write(*,*) '***************************************'
write(*,*) '*   reduce from coordinate to CSC     *'
write(*,*) '***************************************'

allocate(JCOL(NZ))
allocate(IJold(NZ))

irn_min = minval(IRN(1:NZ))  ! distributed matrices do not necesarily start at 1
irn_max = maxval(IRN(1:NZ))

nrows = irn_max - irn_min + 1
allocate(itot(nrows),ipos(nrows),istart(nrows+1))

itot = 0
do i=1,NZ
  itot(irn(i)-irn_min+1)   = itot(irn(i)-irn_min+1) + 1     ! the total number of entries for row (i)
enddo

istart(1) = 1                  ! the starting index
do i=2,nrows+1
  istart(i) =  istart(i-1) + itot(i-1)
enddo

first_row = irn_min
last_row  = irn_max
!---------------- collect column indices for each row
ipos = 0
do i=1,nz
  irow        = irn(i) - irn_min + 1
  ipos(irow)  = ipos(irow)   + 1
  index       = istart(irow) + ipos(irow) - 1
  jcol(index) = jcn(i)
  IJold(index) = i
enddo
n_IRNK = maxval(itot)

deallocate(itot,ipos)

!-------------- combine multiple entries
allocate(IRNK(n_IRNK))
allocate(Anew(NZ))
allocate(Inew(NZ))
allocate(Jnew(NZ))

anew      = 0.
index_new = 0

ihwb(1) = 1
do i=1,nrows

  nj = istart(i+1)-istart(i)

  if (nj .gt. 0 ) then

    call I_mrgrnk(jcol(istart(i):istart(i+1)-1),irnk,nj)  ! find index of sorted column

    index_new = index_new + 1
    index_old = IJold(istart(i)+irnk(1)-1)

    anew(index_new) = A(index_old)
    inew(index_new) = i + first_row - 1
    jnew(index_new) = jcn(index_old)

    do k=2,nj

      index_old = IJold(istart(i)+irnk(k)-1)

      if (jcol(istart(i)+irnk(k)-1) .gt. jcol(istart(i)+irnk(k-1)-1)) then
        index_new = index_new + 1
        INEW(index_new) = i + first_row - 1
        JNEW(index_new) = JCOL(ISTART(i)+IRNK(k)-1)
      endif

      ANEW(index_new) = ANEW(index_new) + A(index_old)

    enddo

    ihwb(i+1) = index_new + 1

  endif
enddo

nznew = ihwb(nrows+1)-1

write(*,'(A,3i12)') ' N, NZ, NZNEW : ',N,NZ,NZNEW

do k=1,nznew
  A(k)   = ANEW(k)
  IRN(k) = INEW(k)
  JCN(k) = JNEW(k)
enddo

deallocate(JCOL)
deallocate(IJold)
deallocate(ISTART)
deallocate(IRNK)
deallocate(Anew)
deallocate(Inew)
deallocate(Jnew)

end
