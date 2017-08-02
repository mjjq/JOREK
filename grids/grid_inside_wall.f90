!> Create a square grid based on the Bezier/modified cubic Hermite representation
subroutine grid_inside_wall(n_R,n_Z,R_begin,R_end,Z_begin,Z_end,boundary,node_list,element_list)

  use mod_parameters
  use data_structure
  use mod_export_restart
  use phys_module, only: xshift

  implicit none
  
  ! --- Routine parameters
  integer,		   intent(in)	 :: n_R 	    !< Number of horizontal nodes (square grid)
  integer,		   intent(in)	 :: n_Z 	    !< Number of vertical nodes (square grid)
  real*8,		   intent(in)	 :: R_begin	   !< R-min (square grid)
  real*8,		   intent(in)	 :: R_end	   !< R-max (square grid)
  real*8,		   intent(in)	 :: Z_begin	   !< Z-min (square grid)
  real*8,		   intent(in)	 :: Z_end	   !< Z-max (square grid)
  logical,		   intent(in)	 :: boundary	   !< Fill boundary information?
  type(type_node_list),    intent(inout) :: node_list	   !< list of grid nodes
  type(type_element_list), intent(inout) :: element_list   !< list of finite elements
  
  ! --- Local variables
  integer		   :: m, n, i, j, k, ii, jj, iv, inode_0, inode_p, inode, n_elements, i_element, ip, iuv
  integer		   :: n_element_start, n_node_start, n_index_start
  real*8		   :: xx_0(n_dim),xx_p(n_dim),uv_0(n_dim),uv_p(n_dim)
  real*8, external	   :: dlength, ddot
  
  ! --- eqdsk variables
  real,allocatable :: psi(:),p(:),f(:),q(:),rlim(:),zlim(:),rbnd(:), zbnd(:)
  real,allocatable :: dpsi(:),dp(:),df(:)
  real,allocatable :: dpr(:),df2(:),dg(:),work(:),psirz(:,:)
  integer	   :: mod_lines, n_lines
  integer	   :: nr, nz, n_psi, nbbs, nlim
  integer	   :: nr_eqdsk, nz_eqdsk
  
  ! --- Dierckx variables
  real,allocatable :: xx(:),yy(:),zc(:), r_bnd(:), z_bnd(:), psi_bnd(:)
  real,allocatable :: tx(:),ty(:),c(:,:),wrk(:)
  integer,allocatable :: iwrk(:)
  real  	   :: dummy(3), xdim,zdim,rzero,rgrid1,zmid,rmaxis,zmaxis,ssimag,ssibry,bcentr
  real  	   :: xip,xdum1,xdum2,xdum3,xdum4,xdum5
  real  	   :: psi_sep, sig_sep, tanh1, zmu0, zn0, zmd
  real  	   :: xb ,xe, yb, ye, smth, fp, fout
  integer	   :: mx,my,kx,ky,nxest,nyest,lwrk,kwrk,ier,iopt,nx,ny, i1, j1
  character	   :: AA*52
  
  ! --- Grid variables
  integer :: nR_save, n_off, i_seg, nR_tmp, i_save, j_save
  integer :: iRp1, iZp1
  integer :: inside_boot, n_extra, inbetweenSegment
  real*8  :: Rmin, Rmax, r_min, r_max
  real*8  :: Zmin, Zmax, z_min, z_max, z_tmp, Z_first_jump
  real*8  :: width, width_prev, diff, psi_tmp
  real*8  :: accuracy
  real*8  :: psi1, psi2, psi3, psi4
  real*8  :: R_elm, Z_elm
  real*8  :: PSI_R, PSI_Z
  real*8  :: psi_right, psi_left
  real*8  :: deriv_right, deriv_left
  integer :: count
  integer, allocatable :: nR_grid(:), node_index(:,:)
  real*8,  allocatable :: R_grid(:,:), Z_grid(:,:), Zlines(:)
  
  
  write(*,*) '*************************************'
  write(*,*) '*       grid_inside_wall  	  *'
  write(*,*) '*************************************'
  write(*,*) ' creating a rectangular grid inside wall:  n_R, n_Z = ',nR,nZ
  

  ! -------------------------------------
  ! --- Read equilibrium from eqdsk file
  ! -------------------------------------
  
  ! --- Open eqdsk file
  write(*,*) 'Opening file eqdsk.dat'
  open(unit=5,file='eqdsk.dat', ACTION = 'read', iostat=ier)
  if (ier .ne. 0) then
    write(*,*) 'Error opening eqdsk file. You need a file named : eqdsk.dat'
    return
  endif
  
  read(5,'(A52,2i4)') AA,nr,nz
  
  write(*,*) AA
  write(*,'(A,2i5)') ' nr, nz : ',nr,nz
  
  read(5,'(5e16.9)') xdim,zdim,rzero,rgrid1,zmid
  read(5,'(5e16.9)') rmaxis,zmaxis,ssimag,ssibry,bcentr
  read(5,'(5e16.9)') xip,ssimag,xdum1,rmaxis,xdum2
  read(5,'(5e16.9)') zmaxis,xdum3,ssibry,xdum4,xdum5
  
  write(*,'(A,2f10.5,A)') ' xdim,  zdim : ',xdim,zdim,' m'
  write(*,'(A,2f10.5,A)') ' rzero, zmid : ',rzero,zmid, ' m'
  write(*,'(A,f10.5,A)')  ' xip 	: ',xip/1e6,' MA'
  write(*,'(A,f10.5,A)')  ' zmaxis	: ',zmaxis,' m'
  write(*,'(A,f10.5,A)')  ' Bvac	: ',bcentr,' T'
  
 	
  write(*,*) 'reading profiles'
    
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
  
  write(*,*) 'reading psi-map'
    
  do k=1,nz
    do i=1,n_lines
      j = 5*(i-1)
      read(5,'(5e16.9)') psirz(j+1,k),psirz(j+2,k),psirz(j+3,k),psirz(j+4,k),psirz(j+5,k)
    enddo
    j = 5*n_lines
    read(5,'(4e16.9)') psirz(j+1,k),psirz(j+2,k),psirz(j+3,k),psirz(j+4,k)
  enddo
  
  write(*,*) 'reading q-profile'
    
  read(5,'(5e16.9)') (q(i),i=1,5*n_lines)
  if (mod_lines .ne. 0) read(5,'(4e16.9)') (q(i),i=5*n_lines+1,n_psi)
  
  write(*,*) 'reading limiter'
  
  read(5,*)  nbbs,nlim
  allocate(rbnd(nbbs),zbnd(nbbs))
  read(5,'(2e16.9)') (rbnd(i),zbnd(i),i=1,nbbs)
  allocate(rlim(nlim),zlim(nlim))
  read(5,'(2e16.9)') (rlim(i),zlim(i),i=1,nlim)
  
  !write(*,*)'hello!',rbnd(1),zbnd(1)
  !write(*,*)'hello!',rbnd(nbbs),zbnd(nbbs)
  !write(*,*)'hello!',rlim(1),zlim(1)
  !write(*,*)'hello!',rlim(nlim),zlim(nlim)
  !stop
  close(5)
  write(*,*) 'done reading'
  
  ! ----------------------------------------------------------------------------------------------
  ! --- Spline psi-map (this doesn't work, I think it might be due to the shape of the eqdsk map)
  ! ----------------------------------------------------------------------------------------------
  
  write(*,*) 'splining psi-map'
  
  ! --- interpolate flux using Dierckx spline routine
  allocate(xx(nr),yy(nz),psi(n_psi))
  nr_eqdsk = nr
  nz_eqdsk = nz
  do i=1,n_psi
    psi(i) = real(i-1)/real(n_psi-1)
  enddo     
  do i=1,nr
    xx(i) = rgrid1 + xdim*real(i-1)/real(nr-1)
  enddo     
  do i=1,nz
    yy(i) = zmid + zdim*(real(i-1)/real(nz-1)-0.5)
  enddo     
  
  iopt= 0 
  mx = nr
  my = nz
  xb = xx(1)
  xe = xx(nr)
  yb = yy(1)
  ye = yy(nz)
  kx = 3
  ky = 3
  smth = 1.d-6
  nxest = nr-5
  nyest = nz-5
  lwrk  = 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+my*(ky+1)+my+nxest
  kwrk  = 3+mx+my+nxest+nyest
  
  allocate(tx(nxest),ty(nyest),c(nxest,nyest),wrk(lwrk),iwrk(kwrk))
  
  call regrid(iopt,mx,xx,my,yy,transpose(psirz),xb,xe,yb,ye,kx,ky,smth,nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
  
  write(*,*) 'Dierckx ier   : ',ier
  write(*,*) 'Dierckx fp    : ',fp
  write(*,*) 'Dierckx nx,ny : ',nx,ny
  
  lwrk = mx*(kx+1)+my*(ky+1)
  kwrk = mx+my
  deallocate(wrk,iwrk)
  allocate(wrk(lwrk),iwrk(kwrk))
  
  write(*,*)'done splining psi-map'
  
  
  ! --------------------------------------------------------------------------
  ! --- Build a grid inside the wall (not the case below has jumps for MAST-U)
  ! --------------------------------------------------------------------------
  
  write(*,*)'Building grid inside wall'
  
  ! --- Input parameters
  Z_first_jump = 1.0!meters
  nR = n_R
  nZ = n_Z - mod(n_Z,2) ! This needs to be even!!!
  accuracy = 1.d-5
  
  ! --- We cut the domain with horizontal lines
  Rmin = minval(rlim)-1.e-3 ; Rmax = maxval(rlim)+1.e-3 ! want them slightly outside wall
  Zmin = minval(zlim)+1.e-3 ; Zmax = maxval(zlim)-1.e-3 ! want them slightly inside wall
  write(*,*)'Rminmax:',Rmin, Rmax
  write(*,*)'Zminmax:',Zmin, Zmax
  
  ! --- nR is the resolution we want in the core, but in MAST-U there are jumps
  ! --- We start at midplane and go down (or up)
  allocate(nR_grid(nZ+1),node_index(4*nR,nZ+1))
  allocate(R_grid(4*nR,nZ+1),Z_grid(4*nR,nZ+1)) ! need some margin for nR due to jumps
  R_grid = 0.d0
  Z_grid = 0.d0
  
  ! --- First lower part
  nR_save = nR
  allocate(Zlines(nZ/2+1))
  do i=1,nZ/2+1
    Zlines(i) = real(i-1)/real(nZ/2) * Zmin
  enddo
  width_prev  = 0.0
  inside_boot = 0
  n_extra     = 0
  do i = 1,nZ/2+1
    !write(*,*)'Doing lower horizontal line ',nZ/2+1-i
    n_off = nZ/2
    call RintersectPolygon(nlim, rlim, zlim, Rmin, Rmax, Zlines(i), accuracy, r_min, r_max)
    width = r_max-r_min
    ! First jump
    if ( (width .lt. 0.85*width_prev) .and. (Zlines(i) .gt. -Z_first_jump) .and. (i .gt. 1) ) then
      diff = width/width_prev
      nR = int(real(nR) * diff) + 1
    endif
    ! Second jump
    if ( (width .gt. 1.2*width_prev) .and. (i .gt. 1) ) then
      diff = width/width_prev
      n_extra = int(real(nR) * diff) - nR
      inside_boot = 1
    endif
    nR_grid(n_off+i) = nR + n_extra
    if (inside_boot .eq. 1) then
      do j = 1,nR
  	R_grid(j,n_off+i) = r_min + (R_grid(j,n_off+i-1) - R_grid(1,n_off+i-1))
  	Z_grid(j,n_off+i) = Zlines(i)
      enddo
      do j = 1,n_extra
  	R_grid(nR+j,n_off+i) = R_grid(nR,n_off+i) + real(j)/real(n_extra) * (r_max - R_grid(nR,n_off+i))
  	Z_grid(nR+j,n_off+i) = Zlines(i)
      enddo
      i_seg = inbetweenSegment(nlim, rlim, zlim, r_max, Zlines(i))
      r_max = min(rlim(i_seg),rlim(i_seg+1))
      z_tmp = max(zlim(i_seg),zlim(i_seg+1))
      nR_grid(n_off+i-1) = nR + n_extra
      do j = 1,n_extra
  	R_grid(nR+j,n_off+i-1) = R_grid(nR,n_off+i-1) + real(j)/real(n_extra) * (r_max - R_grid(nR,n_off+i-1))
  	Z_grid(nR+j,n_off+i-1) = Z_grid(nR,n_off+i-1) + real(j)/real(n_extra) * (z_tmp - Z_grid(nR,n_off+i-1))
      enddo
      inside_boot = 2
    else
      do j = 1,nR_grid(n_off+i)
  	R_grid(j,n_off+i) = r_min + real(j-1)/real(nR_grid(n_off+i)-1) * (r_max - r_min)
  	Z_grid(j,n_off+i) = Zlines(i)
      enddo
    endif
    width_prev = width
  enddo
  write(*,*)'Finished lower part of grid'
  
  ! --- Then upper part
  nR = nR_save
  do i=1,nZ/2+1
    Zlines(i) = real(i-1)/real(nZ/2) * Zmax
  enddo
  width_prev  = 0.0
  inside_boot = 0
  n_extra     = 0
  do i = 1,nZ/2
    n_off = nZ/2+1
    !write(*,*)'Doint upper horizontal line ',nZ/2-i+1
    call RintersectPolygon(nlim, rlim, zlim, Rmin, Rmax, Zlines(i+1), accuracy, r_min, r_max)
    width = r_max-r_min
    ! First jump
    if ( (width .lt. 0.85*width_prev) .and. (Zlines(i+1) .lt. Z_first_jump) .and. (i .gt. 1) ) then
      diff = width/width_prev
      nR = int(real(nR) * diff) + 1
    endif
    ! Second jump
    if ( (width .gt. 1.2*width_prev) .and. (i .gt. 1) ) then
      diff = width/width_prev
      n_extra = int(real(nR) * diff) - nR
      inside_boot = 1
    endif
    nR_grid(n_off-i) = nR + n_extra
    if (inside_boot .eq. 1) then
      do j = 1,nR
  	R_grid(j,n_off-i) = r_min + (R_grid(j,n_off-i+1)-R_grid(1,n_off-i+1))
  	Z_grid(j,n_off-i) = Zlines(i+1)
      enddo
      do j = 1,n_extra
  	R_grid(nR+j,n_off-i) = R_grid(nR,n_off-i) + real(j)/real(n_extra) * (r_max - R_grid(nR,n_off-i))
  	Z_grid(nR+j,n_off-i) = Zlines(i+1)
      enddo
      i_seg = inbetweenSegment(nlim, rlim, zlim, r_max, Zlines(i+1))
      r_max = min(rlim(i_seg),rlim(i_seg+1))
      z_tmp = min(zlim(i_seg),zlim(i_seg+1))
      nR_grid(n_off-i+1) = nR + n_extra
      do j = 1,n_extra
  	R_grid(nR+j,n_off-i+1) = R_grid(nR,n_off-i+1) + real(j)/real(n_extra) * (r_max - R_grid(nR,n_off-i+1))
  	Z_grid(nR+j,n_off-i+1) = Z_grid(nR,n_off-i+1) + real(j)/real(n_extra) * (z_tmp - Z_grid(nR,n_off-i+1))
      enddo
      inside_boot = 2
    else
      do j = 1,nR_grid(n_off-i)
  	R_grid(j,n_off-i) = r_min + real(j-1)/real(nR_grid(n_off-i)-1) * (r_max - r_min)
  	Z_grid(j,n_off-i) = Zlines(i+1)
      enddo
    endif
    width_prev = width
  enddo
  write(*,*)'Finished upper part of grid'
  
  ! --- Print to file?
  if (.true.) then
    open(21,file='RZ_grid_inside_wall.txt')
    do i = 1,nZ
      !write(*,*)'Doint horizontal line ',nZ-i
      nR_tmp = min(nR_grid(i),nR_grid(i+1))
      do j = 1,nR_tmp-1
        write(21,'(e16.9,A,e16.9,A,e16.9,A,e16.9,A,e16.9,A,e16.9,A,e16.9,A,e16.9)') &
          R_grid(j,i),     ' ',Z_grid(j,i),    ' ', &
          R_grid(j,i+1),   ' ',Z_grid(j,i+1),  ' ', &
          R_grid(j+1,i+1), ' ',Z_grid(j+1,i+1),' ', &
          R_grid(j+1,i),   ' ',Z_grid(j+1,i)
      enddo
    enddo
    close(21)
  endif
    
  ! --- Convert into jorek grid
  node_list%n_nodes = 0
  do i = 1,nZ+1
    do j = 1,nR_grid(i)
      node_list%n_nodes = node_list%n_nodes + 1
      node_index(j,i) = node_list%n_nodes
      ! --- Position
      node_list%node(node_list%n_nodes)%x(1,1) = R_grid(j,i)
      node_list%node(node_list%n_nodes)%x(1,2) = Z_grid(j,i)
      
      ! --- vectors
      iRp1 = j+1
      if (iRp1 .gt. nR_grid(i)) iRp1 = j-1
      iZp1 = i+1
      if (iZp1 .gt. nZ+1) iZp1 = i-1
      if ( (iRp1 .gt. nR_grid(iZp1)) .and. (iZp1 .eq. i+1) ) then
        iZp1 = i-1
      else if ( (iRp1 .gt. nR_grid(iZp1)) .and. (iZp1 .eq. i-1) ) then
        write(*,*)'something wrong in defining vectors...'
      endif
      ! --- vector u
      width = sqrt((R_grid(iRp1,i) - R_grid(j,i))**2 + (Z_grid(iRp1,i) - Z_grid(j,i))**2)
      node_list%node(node_list%n_nodes)%x(2,1) = (R_grid(iRp1,i) - R_grid(j,i)) / width
      node_list%node(node_list%n_nodes)%x(2,2) = (Z_grid(iRp1,i) - Z_grid(j,i)) / width
      ! --- vector v
      width = sqrt((R_grid(j,iZp1) - R_grid(j,i))**2 + (Z_grid(j,iZp1) - Z_grid(j,i))**2)
      node_list%node(node_list%n_nodes)%x(3,1) = (R_grid(j,iZp1) - R_grid(j,i)) / width
      node_list%node(node_list%n_nodes)%x(3,2) = (Z_grid(j,iZp1) - Z_grid(j,i)) / width
      ! --- vector w
      node_list%node(node_list%n_nodes)%x(4,1) = 0.d0
      node_list%node(node_list%n_nodes)%x(4,2) = 0.d0
      
      ! --- boundary
      node_list%node(node_list%n_nodes)%boundary = 0
      
      ! --- matrix index
      do k=1, n_order+1
        node_list%node(node_list%n_nodes)%index(k) = (n_order+1)*(node_list%n_nodes-1) + k
      enddo
      
      ! --- psi values
      if (ier .eq. 0) then
        call bispev(tx,nx,ty,ny,c,kx,ky,R_grid(j,i),1,Z_grid(j,i),1,psi_tmp,wrk,lwrk,iwrk,kwrk,ier)
      else
        i_save = 0
        j_save = 0
        do ii = 2,nr_eqdsk
          if ( (R_grid(j,i) .ge. xx(ii-1)) .and. (R_grid(j,i) .le. xx(ii)) ) then
            i_save = ii
            exit
          endif
        enddo
        do jj = 2,nz_eqdsk
          if ( (Z_grid(j,i) .ge. yy(jj-1)) .and. (Z_grid(j,i) .le. yy(jj)) ) then
            j_save = jj
            exit
          endif
        enddo
        ! --- Extrapolate value from 4 nodes
        psi1 = psirz(i_save-1,j_save-1)
        psi2 = psirz(i_save  ,j_save-1)
        psi3 = psirz(i_save  ,j_save  )
        psi4 = psirz(i_save-1,j_save  )
        R_elm = (R_grid(j,i)-xx(i_save-1)) / (xx(i_save)-xx(i_save-1))
        Z_elm = (Z_grid(j,i)-yy(j_save-1)) / (yy(j_save)-yy(j_save-1))
        psi_tmp = psi1 + R_elm*(psi2-psi1) + Z_elm*(psi4-psi1) + R_elm*Z_elm*(psi3+psi1-psi2-psi4)
        ! --- The R derivative
        psi_R = 0.d0
        count = 0
        psi_left    = psi1 + Z_elm*(psi4-psi1)
        psi_right   = psi2 + Z_elm*(psi3-psi2)
        if (abs(R_grid(j,i)-xx(i_save-1)) .gt. 1.d-10) then
          count = count + 1
          deriv_left  = (psi_tmp   - psi_left)/(R_grid(j,i)-xx(i_save-1))
          psi_R = psi_R + deriv_left
        endif
        if (abs(xx(i_save)-R_grid(j,i)) .gt. 1.d-10) then
          count = count + 1
          deriv_right = (psi_right - psi_tmp )/(xx(i_save)-R_grid(j,i))
          psi_R = psi_R + deriv_right
        endif
        if (count .gt. 0) psi_R = - psi_R / real(count) ! min sign because of JOREK definition of psi
        ! --- The Z derivative
        psi_Z = 0.d0
        count = 0
        psi_left    = psi1 + R_elm*(psi2-psi1)
        psi_right   = psi4 + R_elm*(psi3-psi4)
        if (abs(Z_grid(j,i)-yy(j_save-1)) .gt. 1.d-10) then
          count = count + 1
          deriv_left  = (psi_tmp   - psi_left)/(Z_grid(j,i)-yy(j_save-1))
          psi_Z = psi_Z + deriv_left
        endif
        if (abs(yy(j_save)-Z_grid(j,i)) .gt. 1.d-10) then
          count = count + 1
          deriv_right = (psi_right - psi_tmp )/(yy(j_save)-Z_grid(j,i))
          psi_Z = psi_Z + deriv_right
        endif
        if (count .gt. 0) psi_Z = - psi_Z / real(count) ! min sign because of JOREK definition of psi
      endif
      node_list%node(node_list%n_nodes)%values(1,1,1) = -psi_tmp ! min sign because of JOREK definition of psi
      ! --- Derivatives?
      node_list%node(node_list%n_nodes)%values(1,2,1) = &
        psi_R * node_list%node(node_list%n_nodes)%x(2,1) + psi_Z * node_list%node(node_list%n_nodes)%x(2,2)
      node_list%node(node_list%n_nodes)%values(1,3,1) = &
        psi_R * node_list%node(node_list%n_nodes)%x(3,1) + psi_Z * node_list%node(node_list%n_nodes)%x(3,2)
      node_list%node(node_list%n_nodes)%values(1,4,1) = &
        psi_R * node_list%node(node_list%n_nodes)%x(4,1) + psi_Z * node_list%node(node_list%n_nodes)%x(4,2)
                         
      
    enddo
  enddo
  
  
  element_list%n_elements = 0
  do i=1,nZ
    nR_tmp = min(nR_grid(i),nR_grid(i+1))
    do j=1,nR_tmp-1
      element_list%n_elements = element_list%n_elements + 1
      element_list%element(element_list%n_elements)%vertex(1) = node_index(j  ,i  )
      element_list%element(element_list%n_elements)%vertex(2) = node_index(j+1,i  )
      element_list%element(element_list%n_elements)%vertex(3) = node_index(j+1,i+1)
      element_list%element(element_list%n_elements)%vertex(4) = node_index(j  ,i+1)
    enddo
  enddo

  
  do k=1, element_list%n_elements   ! fill in the size of the elements
    do iv=1,4			    ! loop over the vertices

      iuv = mod(iv+1,2)+1	    ! the direction vector corresponding to this edge (i)

      inode_0 = element_list%element(k)%vertex(iv)
      xx_0    = node_list%node(inode_0)%x(1,:)
      uv_0    = node_list%node(inode_0)%x(iuv+1,:)

      ip      = mod(iv,4)+1
      inode_p = element_list%element(k)%vertex(ip)
      xx_p    = node_list%node(inode_p)%x(1,:)
      uv_p    = node_list%node(inode_p)%x(iuv+1,:)

      element_list%element(k)%size(iv,1)     = 1.
      element_list%element(k)%size(iv,iuv+1) = sign(dlength(xx_p,xx_0),ddot(n_dim,xx_p - xx_0,1,uv_0,1)) /3.d0
      element_list%element(k)%size(ip,iuv+1) = sign(dlength(xx_p,xx_0),ddot(n_dim,xx_0 - xx_p,1,uv_p,1)) /3.d0

    enddo

    do iv=1,4
      element_list%element(k)%size(iv,4) = element_list%element(k)%size(iv,2) * element_list%element(k)%size(iv,3)
    enddo
  enddo
  
  call update_boundary_types(element_list,node_list, .false., .false.)
  
  ! --- MAST-U specific (due to stupid equilibria with nonsense p and ffp profiles...)
  ! --- Shift plasma away from solenoid
  if (abs(xshift) .gt. 0.d0) then
    write(*,*)'WARNING!!! shifting plasma, use xshift=0.0 if you do not want to!' 
    do inode=1, node_list%n_nodes
      if ( (node_list%node(inode)%boundary .gt. 0) .and. (node_list%node(inode)%x(1,1) .lt. 0.37) ) then
      !if ( (node_list%node(inode)%boundary .gt. 0) .and. (abs(node_list%node(inode)%x(1,2)) .le. 1.2) ) then
  	Zmin  = -0.7
        z_min = -0.4
  	Zmax  = +0.7
        z_max = +0.4
  	R_elm = node_list%node(inode)%x(1,1)
  	Z_elm = node_list%node(inode)%x(1,2)
  	diff  = xshift
  	if ( (Z_elm .le. z_min) .and. (Z_elm .ge. Zmin) ) then
  	  diff = (Z_elm-Zmin)/(z_min-Zmin)*diff
        elseif ( (Z_elm .ge. z_max) .and. (Z_elm .le. Zmax) ) then
  	  diff = (Z_elm-Zmax)/(z_max-Zmax)*diff
        endif
  	if ( (Z_elm .lt. Zmin) .or. (Z_elm .gt. Zmax) ) then
  	  diff = 0.d0
        endif
  	node_list%node(inode)%values(1,1,1) = node_list%node(inode)%values(1,1,1) + diff * node_list%node(inode)%values(1,1,1)
      endif
    enddo
  endif

  ! --- Update neighbours
  ! --- This may be long, but is necessary for the rest
  call update_neighbours_basic(element_list,node_list)
  
  deallocate(nR_grid, node_index)
  deallocate(R_grid, Z_grid) ! need some margin for nR due to jumps
  deallocate(Zlines)
  call export_restart(node_list, element_list, 'jorek_restart')
  !stop
  
  return
 
end subroutine grid_inside_wall
  






! --- Find the polygon segment a point belongs to
integer function inbetweenSegment(N, rpol, zpol, R, Z)
  
  implicit none
  
  ! --- Input variables
  integer, intent(in) :: N
  real*8,  intent(in) :: rpol(N), zpol(N), R, Z
  
  ! --- Local variables
  integer :: i
  real*8  :: rmin, rmax, zmin, zmax
  
  inbetweenSegment = 0
  do i = 1,N-1
    rmin = min(rpol(i),rpol(i+1))
    rmax = max(rpol(i),rpol(i+1))
    zmin = min(zpol(i),zpol(i+1))
    zmax = max(zpol(i),zpol(i+1))
    if ( (R .ge. rmin) .and. (R .le. rmax) .and. (Z .ge. zmin) .and. (Z .le. zmax) ) then
      inbetweenSegment = i
      exit
    endif
  enddo
  
  return
end function inbetweenSegment




! --- Find intersection between horizontal segment (Rmin,Rmax)@height=Z and polygon
subroutine RintersectPolygon(N, rpol, zpol, Rmin, Rmax, Z, accuracy, new_Rmin, new_Rmax)

  implicit none
  
  ! --- Input variables
  integer, intent(in)  :: N
  real*8,  intent(in)  :: rpol(N), zpol(N), Rmin, Rmax, Z, accuracy
  real*8,  intent(out) :: new_Rmin, new_Rmax
  
  ! --- Local variables
  integer, parameter :: n_iter = 100
  integer :: inside, i1, i2, i, insidePolygon
  real*8  :: R_iter(n_iter), r_tmp, z_tmp
  real*8  :: Rleft, Rmid, Rright
  real*8  :: diff
  
  ! --- Find entry and exit points of horizontal line (in/out of polygon)
  inside = 0
  do i = 1,n_iter
    R_iter(i) = Rmin + real(i-1)/real(n_iter-1) * (Rmax-Rmin)
  enddo
  z_tmp  = Z
  i1 = -1
  i2 = -1
  do i=1,n_iter
    r_tmp = R_iter(i)
    if (insidePolygon(N, rpol, zpol, r_tmp, z_tmp) .eq. 1) then
      if (inside .eq. 0) then
  	i1 = i
        inside = 1
      endif
    endif
    if (insidePolygon(N, rpol, zpol, r_tmp, z_tmp) .eq. 0) then
      if (inside .eq. 1) then
  	i2 = i
        exit
      endif
    endif
  enddo
  if (i1 .eq. -1) then
    write(*,*)'Warning: did not find entry'
    i1 = 1
  endif
  if (i2 .eq. -1) then
    write(*,*)'Warning: did not find exit'
    i2 = 1
  endif
  if (i1 .eq. 1) i1 = 2
  
  ! --- Converge to solution at entry point
  Rleft  = R_iter(i1-1)
  Rright = R_iter(i1)
  Rmid   = 0.5*(Rleft+Rright)
  do i = 1,n_iter
    diff = abs(Rmid - 0.5*(Rleft+Rright))
    Rmid = 0.5*(Rleft+Rright)
    if ( (diff .lt. accuracy) .and. (i .gt. 1) ) then
      new_Rmin = Rmid
      exit
    endif
    if (insidePolygon(N, rpol, zpol, Rmid, z_tmp) .eq. 1) then
      Rright = Rmid
    else
      Rleft  = Rmid
    endif
    if (i  .eq. n_iter) then
      write(*,*)'WARNING: RintersectPolygon did not converge in entry search!'
    endif
  enddo
  
  ! --- Converge to solution at exit point
  Rleft  = R_iter(i2-1)
  Rright = R_iter(i2)
  Rmid   = 0.5*(Rleft+Rright)
  do i = 1,n_iter
    diff = abs(Rmid - 0.5*(Rleft+Rright))
    Rmid = 0.5*(Rleft+Rright)
    if ( (diff .lt. accuracy) .and. (i .gt. 1) ) then
      new_Rmax = Rmid
      exit
    endif
    if (insidePolygon(N, rpol, zpol, Rmid, z_tmp) .eq. 1) then
      Rleft = Rmid
    else
      Rright  = Rmid
    endif
    if (i  == n_iter-1) then
      write(*,*)'WARNING: RintersectPolygon did not converge in exit search!'
    endif
  enddo
  
  return
end subroutine RintersectPolygon





! --- Check if point is inside polygon
integer function insidePolygon(N, rpol, zpol, r, z)

  implicit none
  
  ! --- Input variables
  integer, intent(in)  :: N
  real*8,  intent(in)  :: rpol(N), zpol(N), r, z
  
  ! --- Local variables
  integer :: counter, i
  real*8  :: r1, z1, r2, z2, xinters
  
  insidePolygon = 0
  
  counter = 0
  r1 = rpol(N)
  z1 = zpol(N)
  do i = 1,N
    r2 = rpol(i)
    z2 = zpol(i)
    if ( (z .ge. min(z1,z2)) .and. (z .le. max(z1,z2)) .and. (r .lt. max(r1,r2)) ) then
      if (z1 .ne. z2) then
    	xinters = (z-z1)*(r2-r1)/(z2-z1) + r1
    	if ((r .lt. xinters) .and. (z .ne. z2)) then
	  counter = counter + 1
        endif
      endif
    endif
    r1 = r2
    z1 = z2
  enddo

  ! --- Outside
  if (mod(counter,2) .eq. 0) then
    insidePolygon = 0
  ! --- Inside
  else
    insidePolygon = 1
  endif
  
  return

end function insidePolygon







