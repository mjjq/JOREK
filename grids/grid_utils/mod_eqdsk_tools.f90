module mod_eqdsk_tools

character*9, parameter :: eqdsk_filename = 'eqdsk.dat' ! this could become an input parameter, but that's really irrelevant detail at this point...

contains


subroutine get_eqdsk_style(normal_eqdsk, normal_eqdsk_wall, ier)
  ! note: "normal" means the strict eqdsk format, where 5 floats are written for each line
  ! however, if the radial grid size nr is not a multiple of 5, then the last line of each profile
  ! will not be 5 floats. But this is not true for the psi-grid. So instead of blocks of 5-floats-lines
  ! plus one line with the remainder, the strict eqdsk just keeps going with 5-float lines. Note that the
  ! same applies for the boundary and the wall at the end of the file. This is what I call "normal"
  ! The "not normal" one is when the psi-grid is written as blocks of 5-floats-lines plus one with the remainder,
  ! and the boundary and the wall are written as lines of 2-floats (ie. R & Z)

  implicit none
  
  ! --- Routine variables
  logical, intent(inout) :: normal_eqdsk, normal_eqdsk_wall
  integer, intent(out)   :: ier
  
  ! --- eqdsk variables
  character        :: AA*52
  integer          :: nr, nz
  integer          :: i, n_skip, nbbs, n_wall, str_length
  character*256    :: skip
  
  ! --- Open eqdsk fileget_data_from_eqdsk(nR, nZ, R_grid, Z_grid, psi_grid, ier)
  open(unit=5,file=eqdsk_filename, ACTION = 'read', iostat=ier)
  if (ier .ne. 0) return
  
  read(5,'(A52,2i4)') AA,nr,nz
  
  normal_eqdsk = .true.
  if (mod(nr,5) .ne. 0) then
    ! --- The 4 lines of basic info, the 4 profiles F, P, dF2, dP and the first chunk of the psi_RZ grid
    n_skip = 4 + nr/5 + nr/5 + nr/5 + nr/5 + nr/5 + 4
    do i=1,n_skip
      read(5,*) skip
    enddo
    read(5,'(Q,A)') str_length,skip
    if (str_length .ne. 80) normal_eqdsk = .false.
  endif
  
  ! --- The last profile q and the rest of the psi_RZ grid
  if (.not. normal_eqdsk) then
    n_skip = nr/5 + (nr/5)*(nz-1)
  else
    n_skip = nr/5 + (nr*nz)/5 - 1
  endif
  if (mod(nr,5)    .ne. 0) n_skip = n_skip + 1
  if (mod(nr*nz,5) .ne. 0) n_skip = n_skip + 1
  if ( (.not. normal_eqdsk) .and. (mod(nr,5) .ne. 0) ) n_skip = n_skip - 1 + nz - 1
  
  do i=1,n_skip
    read(5,*) skip
  enddo
  
  read(5,*) nbbs,n_wall
  
  read(5,'(Q,A)') str_length,skip
  normal_eqdsk_wall = .true.
  if (str_length .ne. 80) normal_eqdsk_wall = .false.
  
  close(5)
  
  return

end subroutine get_eqdsk_style








subroutine get_eqdsk_dimensions(normal_eqdsk, nR, nZ, n_wall, ier)

  implicit none
  
  ! --- Routine variables
  logical, intent(in)    :: normal_eqdsk
  integer, intent(inout) :: nR, nZ, n_wall
  integer, intent(out)   :: ier
  
  ! --- eqdsk variables
  character        :: AA*52
  integer          :: i, n_skip, nbbs
  character*256    :: skip
  
  ! --- Open eqdsk fileget_data_from_eqdsk(nR, nZ, R_grid, Z_grid, psi_grid, ier)
  open(unit=5,file=eqdsk_filename, ACTION = 'read', iostat=ier)
  if (ier .ne. 0) return
  
  read(5,'(A52,2i4)') AA,nr,nz
  
  ! --- The 4 lines of basic info, the 5 profiles F, P, dF2, dP, q and the psi_RZ grid
  if (.not. normal_eqdsk) then
    n_skip = 4 + nr/5 + nr/5 + nr/5 + nr/5 + nr/5 + (nr/5)*nz
  else
    n_skip = 4 + nr/5 + nr/5 + nr/5 + nr/5 + nr/5 + (nr*nz)/5
  endif
  if (mod(nr,5)    .ne. 0) n_skip = n_skip + 5
  if (mod(nr*nz,5) .ne. 0) n_skip = n_skip + 1
  if ( (.not. normal_eqdsk) .and. (mod(nr,5) .ne. 0) ) n_skip = n_skip - 1 + nz
  
  do i=1,n_skip
    read(5,*) skip
  enddo
  
  read(5,*) nbbs,n_wall
  
  close(5)
  
  return

end subroutine get_eqdsk_dimensions







subroutine get_data_from_eqdsk(normal_eqdsk, normal_eqdsk_wall, nR, nZ, R_grid, Z_grid, psi_grid, n_wall, R_wall, Z_wall, ier)

  implicit none
  
  ! --- Routine variables
  logical, intent(in)    :: normal_eqdsk, normal_eqdsk_wall
  integer, intent(inout) :: nR, nZ, n_wall
  real*8,  intent(inout) :: R_grid(nR), Z_grid(nZ), psi_grid(nR,nZ)
  real*8,  intent(inout) :: R_wall(n_wall), Z_wall(n_wall)
  integer, intent(out)   :: ier
  
  ! --- eqdsk variables
  integer          :: i, j, k
  real,allocatable :: psi(:),p(:),f(:),q(:),rlim(:),zlim(:),rbnd(:), zbnd(:)
  real,allocatable :: dpsi(:),dp(:),df(:)
  real,allocatable :: dpr(:),df2(:),dg(:),work(:),psirz(:,:)
  real,allocatable :: xx(:),yy(:)
  real             :: xip,xdum1,xdum2,xdum3,xdum4,xdum5
  real             :: dummy(3), xdim,zdim,rzero,rgrid1,zmid,rmaxis,zmaxis,ssimag,ssibry,bcentr
  integer          :: mod_lines, n_lines
  integer          :: n_psi, nbbs, nlim
  character        :: AA*52
  
  ! --- Open eqdsk file
  open(unit=5,file=eqdsk_filename, ACTION = 'read', iostat=ier)
  if (ier .ne. 0) return
  
  read(5,'(A52,2i4)') AA,nr,nz
  read(5,'(5e16.9)') xdim,zdim,rzero,rgrid1,zmid
  read(5,'(5e16.9)') rmaxis,zmaxis,ssimag,ssibry,bcentr
  read(5,'(5e16.9)') xip,ssimag,xdum1,rmaxis,xdum2
  read(5,'(5e16.9)') zmaxis,xdum3,ssibry,xdum4,xdum5
  
  ! --- reading profiles
  n_psi=nr
  mod_lines = mod(n_psi,5)
  n_lines   = (n_psi-mod(n_psi,5))/5
  
  allocate(f(n_psi),p(n_psi),df2(n_psi),dpr(n_psi),psirz(nr,nz),q(n_psi))
  
  read(5,'(5e16.9)') (f(i),i=1,5*n_lines)
  if (mod_lines .ne. 0) read(5,'(4e16.9)') (f(i),i=5*n_lines+1,n_psi)
  
  read(5,'(5e16.9)') (p(i),i=1,5*n_lines)
  if (mod_lines .ne. 0) read(5,'(4e16.9)') (p(i),i=5*n_lines+1,n_psi)
  
  read(5,'(5e16.9)') (df2(i),i=1,5*n_lines)
  if (mod_lines .ne. 0) read(5,'(4e16.9)') (df2(i),i=5*n_lines+1,n_psi)
  
  read(5,'(5e16.9)') (dpr(i),i=1,5*n_lines)
  if (mod_lines .ne. 0) read(5,'(4e16.9)') (dpr(i),i=5*n_lines+1,n_psi)
  
  ! --- reading psi-map
  if (normal_eqdsk) then
    read(5,'(5e16.9)') ((psirz(i,j),i=1,nr),j=1,nz)
  else
    do k=1,nz
      do i=1,n_lines
        j = 5*(i-1)
        read(5,'(5e16.9)') psirz(j+1,k),psirz(j+2,k),psirz(j+3,k),psirz(j+4,k),psirz(j+5,k)
      enddo
      j = 5*n_lines
      read(5,'(4e16.9)') psirz(j+1,k),psirz(j+2,k),psirz(j+3,k),psirz(j+4,k)
    enddo
  endif
  
  write(*,*) 'reading q-profile'
    
  read(5,'(5e16.9)') (q(i),i=1,5*n_lines)
  if (mod_lines .ne. 0) read(5,'(4e16.9)') (q(i),i=5*n_lines+1,n_psi)
  
  write(*,*) 'reading limiter'
  
  read(5,*)  nbbs,nlim
  allocate(rbnd(nbbs),zbnd(nbbs))
  allocate(rlim(nlim),zlim(nlim))
  if (normal_eqdsk_wall) then
    read(5,'(5e16.9)') (rbnd(i),zbnd(i),i=1,nbbs)
    read(5,'(5e16.9)') (rlim(i),zlim(i),i=1,nlim)
  else
    read(5,'(2e16.9)') (rbnd(i),zbnd(i),i=1,nbbs)
    read(5,'(2e16.9)') (rlim(i),zlim(i),i=1,nlim)
  endif
  
  close(5)
  
  ! --- allocate grid data
  allocate(xx(nr),yy(nz))
  do i=1,nr
    xx(i) = rgrid1 + xdim*real(i-1)/real(nr-1)
  enddo     
  do i=1,nz
    yy(i) = zmid + zdim*(real(i-1)/real(nz-1)-0.5)
  enddo     
  
  ! --- Copy into arguments
  do i=1,nR
    R_grid(i) = xx(i)
  enddo
  do j=1,nZ
    Z_grid(j) = yy(j)
  enddo
  do i=1,nR
    do j=1,nZ
      psi_grid(i,j) = psirz(i,j)
    enddo
  enddo
  do i=1,nlim
    R_wall(i) = rlim(i)
    Z_wall(i) = zlim(i)
  enddo
  
  deallocate(f,p,df2,dpr,psirz,q)
  deallocate(rbnd,zbnd)
  deallocate(rlim,zlim)
  deallocate(xx,yy)
  
  return
end subroutine get_data_from_eqdsk











subroutine interpolate_psi_from_eqdsk_grid(nr_eqdsk, nz_eqdsk, xx, yy, psirz, R_find, Z_find, psi, psi_R, psi_Z)

  implicit none
  
  ! --- Input variables
  integer, intent(in)  :: nr_eqdsk, nz_eqdsk
  real*8,  intent(in)  :: xx(nr_eqdsk), yy(nz_eqdsk), psirz(nr_eqdsk, nz_eqdsk), R_find, Z_find
  real*8,  intent(out) :: psi, psi_R, psi_Z
  
  ! --- Local variables
  integer :: ii, jj, i_save, j_save, count
  real*8  :: psi1, psi2, psi3, psi4
  real*8  :: R_elm, Z_elm
  real*8  :: psi_left, psi_right, deriv_left, deriv_right
  
  if ( (R_find .gt. maxval(xx)) .or. (R_find .lt. minval(xx)) ) then
    write(*,*)'Warning, asking for point outside eqdsk grid, this should not happen...'
    psi = 0.d0 ; psi_R = 0.d0 ; psi_Z = 0.d0
    return
  endif
  
  if ( (Z_find .gt. maxval(yy)) .or. (Z_find .lt. minval(yy)) ) then
    write(*,*)'Warning, asking for point outside eqdsk grid, this should not happen...'
    psi = 0.d0 ; psi_R = 0.d0 ; psi_Z = 0.d0
    return
  endif
  
  
  i_save = 0
  j_save = 0
  do ii = 2,nr_eqdsk
    if ( (R_find .ge. xx(ii-1)) .and. (R_find .le. xx(ii)) ) then
      i_save = ii
      exit
    endif
  enddo
  do jj = 2,nz_eqdsk
    if ( (Z_find .ge. yy(jj-1)) .and. (Z_find .le. yy(jj)) ) then
      j_save = jj
      exit
    endif
  enddo
  ! --- Extrapolate psi value from the 4 eqdsk nodes around our point
  psi1 = psirz(i_save-1,j_save-1)
  psi2 = psirz(i_save  ,j_save-1)
  psi3 = psirz(i_save  ,j_save  )
  psi4 = psirz(i_save-1,j_save  )
  R_elm = (R_find-xx(i_save-1)) / (xx(i_save)-xx(i_save-1))
  Z_elm = (Z_find-yy(j_save-1)) / (yy(j_save)-yy(j_save-1))
  psi = psi1 + R_elm*(psi2-psi1) + Z_elm*(psi4-psi1) + R_elm*Z_elm*(psi3+psi1-psi2-psi4)
  
  ! --- The R derivative
  psi_R = 0.d0
  count = 0
  psi_left    = psi1 + Z_elm*(psi4-psi1)
  psi_right   = psi2 + Z_elm*(psi3-psi2)
  if (abs(R_find-xx(i_save-1)) .gt. 1.d-10) then
    count = count + 1
    deriv_left  = (psi   - psi_left)/(R_find-xx(i_save-1))
    psi_R = psi_R + deriv_left
  endif
  if (abs(xx(i_save)-R_find) .gt. 1.d-10) then
    count = count + 1
    deriv_right = (psi_right - psi )/(xx(i_save)-R_find)
    psi_R = psi_R + deriv_right
  endif
  if (count .gt. 0) psi_R = psi_R / real(count)
  
  ! --- The Z derivative
  psi_Z = 0.d0
  count = 0
  psi_left    = psi1 + R_elm*(psi2-psi1)
  psi_right   = psi4 + R_elm*(psi3-psi4)
  if (abs(Z_find-yy(j_save-1)) .gt. 1.d-10) then
    count = count + 1
    deriv_left  = (psi   - psi_left)/(Z_find-yy(j_save-1))
    psi_Z = psi_Z + deriv_left
  endif
  if (abs(yy(j_save)-Z_find) .gt. 1.d-10) then
    count = count + 1
    deriv_right = (psi_right - psi )/(yy(j_save)-Z_find)
    psi_Z = psi_Z + deriv_right
  endif
  if (count .gt. 0) psi_Z = psi_Z / real(count)

  ! --- min sign because of JOREK definition of psi
  psi   = - psi  
  psi_R = - psi_R
  psi_Z = - psi_Z

  return

end subroutine interpolate_psi_from_eqdsk_grid













subroutine get_wall_from_eqdsk(n_wall, R_wall, Z_wall, ier)

  use grid_xpoint_data, only: n_wall_max
  
  implicit none
  
  ! --- Routine variables
  integer, intent(inout) :: n_wall
  real*8,  intent(inout) :: R_wall(n_wall_max), Z_wall(n_wall_max)
  integer, intent(out)   :: ier
  
  ! --- eqdsk variables
  integer          :: i, j, k
  real,allocatable :: psi(:),p(:),f(:),q(:),rlim(:),zlim(:),rbnd(:), zbnd(:)
  real,allocatable :: dpsi(:),dp(:),df(:)
  real,allocatable :: dpr(:),df2(:),dg(:),work(:),psirz(:,:)
  real             :: xip,xdum1,xdum2,xdum3,xdum4,xdum5
  real             :: dummy(3), xdim,zdim,rzero,rgrid1,zmid,rmaxis,zmaxis,ssimag,ssibry,bcentr
  integer          :: mod_lines, n_lines
  integer          :: nr, nz, n_psi, nbbs, nlim
  integer          :: nr_eqdsk, nz_eqdsk
  character        :: AA*52
  
  ! --- Open eqdsk file
  !write(*,*) 'Opening file eqdsk.dat'
  open(unit=5,file=eqdsk_filename, ACTION = 'read', iostat=ier)
  if (ier .ne. 0) return
  
  read(5,'(A52,2i4)') AA,nr,nz
  read(5,'(5e16.9)') xdim,zdim,rzero,rgrid1,zmid
  read(5,'(5e16.9)') rmaxis,zmaxis,ssimag,ssibry,bcentr
  read(5,'(5e16.9)') xip,ssimag,xdum1,rmaxis,xdum2
  read(5,'(5e16.9)') zmaxis,xdum3,ssibry,xdum4,xdum5
  
  ! --- reading profiles
  n_psi=nr
  mod_lines = mod(n_psi,5)
  n_lines   = (n_psi-mod(n_psi,5))/5
  
  allocate(f(n_psi),p(n_psi),df2(n_psi),dpr(n_psi),psirz(nr,nz),q(n_psi))
  
  read(5,'(5e16.9)') (f(i),i=1,5*n_lines)
  if (mod_lines .ne. 0) read(5,'(4e16.9)') (f(i),i=5*n_lines+1,n_psi)
  
  read(5,'(5e16.9)') (p(i),i=1,5*n_lines)
  if (mod_lines .ne. 0) read(5,'(4e16.9)') (p(i),i=5*n_lines+1,n_psi)
  
  read(5,'(5e16.9)') (df2(i),i=1,5*n_lines)
  if (mod_lines .ne. 0) read(5,'(4e16.9)') (df2(i),i=5*n_lines+1,n_psi)
  
  read(5,'(5e16.9)') (dpr(i),i=1,5*n_lines)
  if (mod_lines .ne. 0) read(5,'(4e16.9)') (dpr(i),i=5*n_lines+1,n_psi)
  
  ! --- reading psi-map
  do k=1,nz
    do i=1,n_lines
      j = 5*(i-1)
      read(5,'(5e16.9)') psirz(j+1,k),psirz(j+2,k),psirz(j+3,k),psirz(j+4,k),psirz(j+5,k)
    enddo
    j = 5*n_lines
    read(5,'(4e16.9)') psirz(j+1,k),psirz(j+2,k),psirz(j+3,k),psirz(j+4,k)
  enddo
  
  ! --- reading q-profile
  read(5,'(5e16.9)') (q(i),i=1,5*n_lines)
  if (mod_lines .ne. 0) read(5,'(4e16.9)') (q(i),i=5*n_lines+1,n_psi)
  
  ! --- reading limiter
  read(5,*)  nbbs,nlim
  allocate(rbnd(nbbs),zbnd(nbbs))
  read(5,'(2e16.9)') (rbnd(i),zbnd(i),i=1,nbbs)
  allocate(rlim(nlim),zlim(nlim))
  read(5,'(2e16.9)') (rlim(i),zlim(i),i=1,nlim)
  
  close(5)
  
  n_wall = nlim
  R_wall(1:n_wall) = rlim(1:nlim)
  Z_wall(1:n_wall) = zlim(1:nlim)
  
  return
end subroutine get_wall_from_eqdsk




end module mod_eqdsk_tools


