program convert

  use hdf5
  use hdf5_io_module
  
implicit none
integer, parameter :: ind_max =    1000000
integer, parameter :: n_nodes_max = 100000
integer, parameter :: nodes_per_elem = 8
logical, parameter :: ascii = .true.
real*8,  parameter :: PI = 3.141592653589793

real*8  :: xyznode(n_nodes_max,3), xyzav(3), maxdist
integer :: elemnode(n_nodes_max/2, nodes_per_elem)
integer :: tmp_elemnode(nodes_per_elem)
integer :: elem_component(n_nodes_max/2)
integer :: n_nodes, n_elems, ierr, i, j, ifile = 143, min_elem, max_elem
character(len=80) :: buffer
character(len=12) :: str1, str2
character(len=1), parameter :: lf = char(10)
character(len=20), parameter :: node_file='x.dat', elem_file='ix.dat'

logical :: include_struct(50)

! for current conversion
integer(HID_T)     :: file_id
integer :: n_wall_curr,l1,l2,ind_last, error
integer :: ij_glob(ind_max,2)=0, ind_glob(ind_max)=0, dof
real*8  :: xyz_e(3), phi
real*8, allocatable, dimension(:)     :: wall_curr, wall_curr_real, wall_curr_el
real*8, allocatable, dimension(:)     :: jx, jy, jz, jphi
real*8, allocatable, dimension(:,:)   :: S, gmat_sp
real*8, allocatable, dimension(:,:,:) :: gmat

integer, allocatable :: ndofel(:), ind_g(:,:), full_glob(:,:)
include_struct(:) = .true.
include_struct(1:50) = .true.
min_elem = 1
max_elem = 100000

! --- Read nodes
open(42, file=trim(node_file), status='old', action='read')
i = 0
do
  i = i + 1
  read(42,*, iostat=ierr) xyznode(i,1), xyznode(i,2), xyznode(i,3)
  if ( ierr /= 0 ) exit
  write(98,'(3ES18.9)') xyznode(i,1), xyznode(i,3), xyznode(i,2) !###
end do
close(42)
n_nodes = i-1
write(*,*) 'Read ', n_nodes, ' nodes.'

! --- Read elements
open(42, file=trim(elem_file), status='old', action='read')
i = 0
do
  i = i + 1
  read(42,*, iostat=ierr) elemnode(i,:), elem_component(i)
  if ( ierr /= 0 ) exit
  if ( nodes_per_elem == 8 ) then
    tmp_elemnode(:) = elemnode(i,:)
    elemnode(i,1) = tmp_elemnode(1)
    elemnode(i,2) = tmp_elemnode(5)
    elemnode(i,3) = tmp_elemnode(8)
    elemnode(i,4) = tmp_elemnode(4)
    elemnode(i,5) = tmp_elemnode(2)
    elemnode(i,6) = tmp_elemnode(6)
    elemnode(i,7) = tmp_elemnode(7)
    elemnode(i,8) = tmp_elemnode(3)
  else
    write(*,*) 'STOP, not implemented'
    stop
  end if
  
  write(99,'(99i10)') elemnode(i,:), elem_component(i) !###
end do
close(42)
n_elems = i-1
write(*,*) 'Read ', n_elems, ' elements.'

! --- Check elements
write(*,*) '  Min node index in elements: ', minval(elemnode(1:n_elems,:))
write(*,*) '  Max node index in elements: ', maxval(elemnode(1:n_elems,:))
write(*,*) '  Indices first element: ', elemnode(1,:)
write(*,*) '  Indices last  element: ', elemnode(n_elems,:)
write(*,*)


! --- export to VTK
if ( ascii ) then
   open(ifile, file='3dwall.vtk', form='formatted', status='replace')
else
  open(ifile, file='3dwall.vtk', form='binary', convert='BIG_ENDIAN', status='replace')
end if

if ( ascii ) then
  write(ifile,'(a)')'# vtk DataFile Version 3.0'
  write(ifile,'(a)')'vtk output'
  write(ifile,'(a)')'ASCII'
  write(ifile,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(str1,'(i12)') n_nodes
  write(ifile,'(a)')'POINTS '//str1//'  float'
else
  buffer = '# vtk DataFile Version 3.0'//lf    ; write(ifile) trim(buffer)
  buffer = 'vtk output'//lf                    ; write(ifile) trim(buffer)
  buffer = 'BINARY'//lf                        ; write(ifile) trim(buffer)
  buffer = 'DATASET UNSTRUCTURED_GRID'//lf     ; write(ifile) trim(buffer)
  
  write(str1,'(i12)') n_nodes
  buffer = 'POINTS '//str1//'  float'//lf      ; write(ifile) trim(buffer)
end if
if ( ascii ) then
  do j = 1, n_nodes
    write(ifile,'(99f12.5)') xyznode(j,:)
  end do
else
  write(ifile) ((real(xyznode(j,i),4),i=1,3),j=1,n_nodes)
end if

if ( ascii ) then
  write(str1(1:12),'(i12)') n_elems
  write(str2(1:12),'(i12)') (nodes_per_elem+1)*n_elems
  write(ifile,'(a)') 'CELLS '//str1//' '//str2
  do j = 1, n_elems
    write(ifile,'(99i10)') nodes_per_elem, elemnode(j,:)-1
  end do
else
  write(str1(1:12),'(i12)') n_elems
  write(str2(1:12),'(i12)') (nodes_per_elem+1)*n_elems
  buffer = lf//'CELLS '//str1//' '//str2//lf  ; write(ifile) trim(buffer)
  write(ifile) (int(nodes_per_elem,4), (int(elemnode(j,i),4),i=1,nodes_per_elem), j=1,n_elems)
end if

if ( ascii ) then
  write(str1(1:12),'(i12)') n_elems
  write(ifile,'(a)') 'CELL_TYPES'//str1
  if ( nodes_per_elem == 8 ) then
    write(ifile,'(4i10)') (12, i=1,n_elems)
  else
    write(ifile,'(4i10)') (13, i=1,n_elems)
  end if
else
  write(str1(1:12),'(i12)') n_elems
  buffer = lf//'CELL_TYPES'//str1//lf         ; write(ifile) trim(buffer)
  if ( nodes_per_elem == 8 ) then
    write(ifile) (int(12,4),i=1,n_elems)
  else
    write(ifile) (int(13,4),i=1,n_elems)
  end if
end if

! ================= Convert wall currents =================================
call HDF5_open('jorek_restart.h5',file_id,error)
call HDF5_integer_reading(file_id,n_wall_curr,"n_wall_curr")
allocate(wall_curr(n_wall_curr),S(n_wall_curr,n_wall_curr),wall_curr_real(n_wall_curr))
call HDF5_array1D_reading(file_id,wall_curr,"wall_curr")
call HDF5_close(file_id)
n_wall_curr=20320
!open(13, file='wall_curr.dat', action='read')
!i=0
!
!do
!  i = i + 1
!  read(13,*, iostat=ierr) wall_curr(i)
!  if ( ierr /= 0 ) exit
!end do
!close(13)
!write(*,*) 'wall_curr', sum(abs(wall_curr))
!#============== LOAD number of DOFS per element
open(13, file='ndofel.dat', action='read')
allocate(ndofel(n_elems))
do i=1,n_elems
  read(13,'(i)') ndofel(i)
end do
close(13)

!=== gives current at baricenter of each element, i_element , i_edge, jx, jy, jz
open(13, file='gmat.dat', action='read')
allocate(gmat_sp(ind_max,3),ind_g(ind_max,2))
ind_g=0; gmat_sp=0.d0
i=0
do
   i = i + 1
   read(13,*,iostat=ierr) ind_g(i,:), gmat_sp(i,:)
   if ( ierr /= 0 ) exit
   ind_last = i
end do
close(13)
! Transform into full matrix
l1=maxval(ind_g(:,1));l2=maxval(ind_g(:,2))
allocate(gmat(l1,l2,3))
gmat=0.d0
do i = 1,ind_last
   gmat(ind_g(i,1),ind_g(i,2),:) = gmat_sp(i,:)
end do
deallocate(gmat_sp)

! relationship between dof and element
open(13, file='iglobdof.dat', action='read')
i=0
do
   i = i + 1
   read(13,*,iostat=ierr) ij_glob(i,:), ind_glob(i)
   if ( ierr /= 0 ) exit
   ind_last=i
end do
close(13)
l1=maxval(ij_glob(:,1));l2=maxval(ij_glob(:,2))

allocate(full_glob(l1,l2))
full_glob=0.d0
do i =1,ind_last
   full_glob(ij_glob(i,1),ij_glob(i,2)) = ind_glob(i)
end do


open(13, file='S_mat.bin', status='old', form='unformatted', access='stream', action='read')
read(13) S
close(13)
do i = 1, n_wall_curr
   wall_curr_real(i) = sum(S(i,:)*wall_curr)
end do
deallocate(wall_curr)

open(13, file='wall_curr_test.dat')
do i=1,n_wall_curr
   write(13, '(ES16.8)') wall_curr_real(i)/4e-7/3.141592653589793
end do
close(13)
allocate(wall_curr_el(n_elems),jx(n_elems),jy(n_elems),jz(n_elems), jphi(n_elems))
wall_curr_el=0.d0

jx=0.d0;   jy=0.d0;    jz=0.d0
do i = 1, n_elems
   do j = 1, ndofel(i)
      dof = full_glob(i,j)
      if (dof > (n_wall_curr-32)) cycle
      jx(i) =jx(i) + gmat(i,j,1)*wall_curr_real(dof)/4e-7/PI
      jy(i) =jy(i) + gmat(i,j,2)*wall_curr_real(dof)/4e-7/PI
      jz(i) =jz(i) + gmat(i,j,3)*wall_curr_real(dof)/4e-7/PI
   end do
   do j = 1,nodes_per_elem
      xyz_e(:) = xyznode(elemnode(i,j),:)/nodes_per_elem
   end do
   phi = atan2(xyz_e(2),xyz_e(1))
   jphi(i) = jx(i) * sin(phi) - jy(i) * cos(phi)
   wall_curr_el(i) = (jx(i)**2+jy(i)**2+jz(i)**2)**.5
end do
if ( ascii ) then
   write(ifile,'(a,i)')'CELL_DATA ',n_elems
   write(ifile,'(a)')'Scalars j_w(Am^-3) float 1'
   write(ifile,'(a)')'LOOKUP_TABLE default'
   write(ifile,'(4es28.9)') (wall_curr_el(i), i=1,n_elems)
else
   write(*,*) 'not implemented'
endif
if ( ascii ) then
   write(ifile,'(a)')'Scalars I*e_phi float 1'
   write(ifile,'(a)')'LOOKUP_TABLE default'
   write(ifile,'(4es28.9)') (jphi(i), i=1,n_elems)
else
   write(*,*) 'not implemented'
endif

if ( ascii ) then
   write(ifile,'(a,i)')'VECTORS vectors float'
   write(ifile,'(4es28.9)') ((/jx(i),jy(i),jz(i)/), i=1,n_elems)
else
   write(*,*) 'not implemented'
endif
close(ifile)
stop



write(str1(1:12),'(i12)') n_elems
buffer = lf//'CELL_DATA '//str1            ; write(ifile) trim(buffer)

buffer = lf//'SCALARS '//'blank'//' float'//lf ; write(ifile) trim(buffer)
buffer = 'LOOKUP_TABLE default'//lf;                                write(ifile) trim(buffer)
end program convert


