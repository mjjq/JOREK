!> module to do the spline-linear interpolation in 2D and normal 2D linear
!> interpolation.
module mod_interp_splinear
implicit none

type Fspline
  real*8, allocatable :: Aspline(:,:),Bspline(:,:),Cspline(:,:),Dspline(:,:) !< Four arrays of spline functions used for linear-spline interpolation, the first dimension is for linear, the second for spline
  real*8, allocatable :: xspline(:) !< The spline direction
  real*8, allocatable :: ylinear(:) !< The linear direction
  integer             :: n_x, n_y   !< Size on x and y direction
end type Fspline

contains

! This is a standard rountine to allocate and initialize Fspline type
subroutine AllocFspline(f,nx,ny)
integer, intent(in)          :: nx, ny !< The dimension of the spline and linear direction
type(Fspline)                :: f      !< Fspline type to be initialized

  allocate(f%Aspline(ny,nx))
  allocate(f%Bspline(ny,nx))
  allocate(f%Cspline(ny,nx))
  allocate(f%Dspline(ny,nx))
  allocate(f%xspline(nx))
  allocate(f%ylinear(ny))
  f%n_x     = nx
  f%n_y     = ny
  f%Aspline = 0.
  f%Bspline = 0.
  f%Cspline = 0.
  f%Dspline = 0.
  f%xspline = 0.
  f%ylinear = 0.

end subroutine AllocFspline

!> Linear 2D interpolation on a rectangular grid
!> x2y1       xy1    x1y1
!>  *----------*------*
!>             |
!>             * xy
!>             |
!>             |
!>  *----------*------*
!> x2y2       xy2    x1y2
!>
!> Calculates the interpolation using two intermediate values
!> fxy1 and fxy2.
!> Equations used are:
!> \[fx1  = \frac{f_{11}-f_{21}}{x_1-x_2} (x-x_1) + f_{11}\]
!> \[fx2  = \frac{f_{12}-f_{22}}{x_1-x_2} (x-x_1) + f_{12}\]
!> \[fout = \frac{f_{x1}-f_{x2}}{y_1-y_2} (y-y_1) + f_{x1}\]
!> x1,2 and y1,2 are chosen in order of closeness
!> This algorithm can also be used for extrapolation
pure function L2Dinterp(tx,ty,f,x,y) result(fout)
real*8, intent(in), dimension(:)                 :: tx !< Grid points in x
real*8, intent(in), dimension(:)                 :: ty !< Grid points in y
real*8, intent(in), dimension(size(tx),size(ty)) :: f !< Function values at these points
real*8, intent(in)  :: x, y !< Points at which to interpolate
real*8              :: fout

integer :: ix1, iy1 !< Index of closest point
integer :: ix2, iy2 !< Index of other (usually next closest) point
real*8  :: fx1, fx2 ! Temporary variables
ix1 = minloc(abs(tx - x), dim=1)
if (x .ge. tx(ix1)) ix2 = ix1 + 1 ! find other index
if (x .lt. tx(ix1)) ix2 = ix1 - 1
if (ix2 .gt. size(tx)) ix2 = size(tx) - 1 ! if it does not exist, extrapolate
if (ix2 .lt. 1       ) ix2 = 2
iy1 = minloc(abs(ty - y), dim=1)
if (y .ge. ty(iy1)) iy2 = iy1 + 1
if (y .lt. ty(iy1)) iy2 = iy1 - 1
if (iy2 .gt. size(ty)) iy2 = size(ty) - 1
if (iy2 .lt. 1       ) iy2 = 2

fx1  = (f(ix1,iy1) - f(ix2,iy1))/(tx(ix1) - tx(ix2)) * (x - tx(ix1)) + f(ix1,iy1)
fx2  = (f(ix1,iy2) - f(ix2,iy2))/(tx(ix1) - tx(ix2)) * (x - tx(ix1)) + f(ix1,iy2)
fout = (fx1 - fx2) / (ty(iy1) - ty(iy2)) * (y - ty(iy1)) + fx1
end function L2Dinterp

!> Like [[L2Dinterp]], but interpolate a vector (3rd dimension) of size nz
pure function L2D2interp(tx,ty,nz,f,x,y) result(fout)
real*8, intent(in), dimension(:)                    :: tx !< Grid points in x
real*8, intent(in), dimension(:)                    :: ty !< Grid points in y
integer, intent(in)                                 :: nz !< number of scalars
real*8, intent(in), dimension(size(tx),size(ty),nz) :: f !< Function values at these points
real*8, intent(in)  :: x, y !< Points at which to interpolate
real*8, dimension(nz) :: fout

integer :: ix1, iy1 !< Index of closest point
integer :: ix2, iy2 !< Index of other (usually next closest) point
real*8  :: fx1(nz), fx2(nz) ! Temporary variables
ix1 = minloc(abs(tx - x), dim=1)
if (x .ge. tx(ix1)) ix2 = ix1 + 1 ! find other index
if (x .lt. tx(ix1)) ix2 = ix1 - 1
if (ix2 .gt. size(tx)) ix2 = size(tx) - 1 ! if it does not exist, extrapolate
if (ix2 .lt. 1       ) ix2 = 2
iy1 = minloc(abs(ty - y), dim=1)
if (y .ge. ty(iy1)) iy2 = iy1 + 1
if (y .lt. ty(iy1)) iy2 = iy1 - 1
if (iy2 .gt. size(ty)) iy2 = size(ty) - 1
if (iy2 .lt. 1       ) iy2 = 2

fx1  = (f(ix1,iy1,:) - f(ix2,iy1,:))/(tx(ix1) - tx(ix2)) * (x - tx(ix1)) + f(ix1,iy1,:)
fx2  = (f(ix1,iy2,:) - f(ix2,iy2,:))/(tx(ix1) - tx(ix2)) * (x - tx(ix1)) + f(ix1,iy2,:)
fout = (fx1 - fx2) / (ty(iy1) - ty(iy2)) * (y - ty(iy1)) + fx1
end function L2D2interp

!> Like [[L2D2interp]], but calculate the gradient in direction `dim`.
!> Dim must be 1 or 2.
pure function L2D2interp_grad(tx,ty,nz,f,x,y,dim) result(fout)
real*8, intent(in), dimension(:)                    :: tx !< Grid points in x
real*8, intent(in), dimension(:)                    :: ty !< Grid points in y
integer, intent(in)                                 :: nz !< number of scalars
real*8, intent(in), dimension(size(tx),size(ty),nz) :: f !< Function values at these points
real*8, intent(in)  :: x, y !< Points at which to interpolate
integer, intent(in) :: dim !< which dimension to interpolate in
real*8, dimension(nz) :: fout

integer :: ix1, iy1 !< Index of closest point
integer :: ix2, iy2 !< Index of other (usually next closest) point
real*8  :: fx1(nz), fx2(nz) ! Temporary variables
ix1 = minloc(abs(tx - x), dim=1)
if (x .ge. tx(ix1)) ix2 = ix1 + 1 ! find other index
if (x .lt. tx(ix1)) ix2 = ix1 - 1
if (ix2 .gt. size(tx)) ix2 = size(tx) - 1 ! if it does not exist, extrapolate
if (ix2 .lt. 1       ) ix2 = 2
iy1 = minloc(abs(ty - y), dim=1)
if (y .ge. ty(iy1)) iy2 = iy1 + 1
if (y .lt. ty(iy1)) iy2 = iy1 - 1
if (iy2 .gt. size(ty)) iy2 = size(ty) - 1
if (iy2 .lt. 1       ) iy2 = 2

if (dim .eq. 1) then
  fx1 = (f(ix1,iy1,:) - f(ix2,iy1,:))/(tx(ix1) - tx(ix2)) * (x - tx(ix1)) + f(ix1,iy1,:)
  fx2 = (f(ix1,iy2,:) - f(ix2,iy2,:))/(tx(ix1) - tx(ix2)) * (x - tx(ix1)) + f(ix1,iy2,:)
  fout = (fx1 - fx2) / (ty(iy1) - ty(iy2))
elseif (dim .eq. 2) then
  fx1 = (f(ix1,iy1,:) - f(ix1,iy2,:))/(ty(iy1) - ty(iy2)) * (y - ty(iy1)) + f(ix1,iy1,:)
  fx2 = (f(ix2,iy1,:) - f(ix2,iy2,:))/(ty(iy1) - ty(iy2)) * (y - ty(iy1)) + f(ix1,iy2,:)
  fout = (fx1 - fx2) / (tx(ix1) - tx(ix2))
else
  fout = -1d99
end if

end function L2D2interp_grad


!> Like [[L2Dinterp]], but do spline interpolation on one dimension and linear
!> interpolation on the other.
subroutine SL2Dinterp(f,x,y,fout,dfout_dx,d2fout_dx2,dfout_dy)
real*8, intent(in)           :: x, y !< Points at which to interpolate
type(Fspline), intent(in)    :: f  
real*8, intent(out),optional :: fout
real*8, intent(out),optional :: dfout_dx
real*8, intent(out),optional :: d2fout_dx2
real*8, intent(out),optional :: dfout_dy

integer :: ix1, iy1 !< Index of closest point
integer :: ix2, iy2 !< Index of other (usually next closest) point
real*8  :: fy1, fy2 ! Temporary variables on each y index

real*8, external :: spwert ! Evaluation of spline
real*8  :: ABLTG1(3), ABLTG2(3)        ! The evaluated first, second and third derivatives

! Find the neighboring index in linear direction
iy1 = minloc(abs(f%ylinear - y), dim=1) 
if (y .ge. f%ylinear(iy1)) iy2 = iy1 + 1
if (y .lt. f%ylinear(iy1)) iy2 = iy1 - 1
if (iy2 .gt. size(f%ylinear)) iy2 = size(f%ylinear) - 1 ! if it does not exist, extrapolate
if (iy2 .lt. 1              ) iy2 = 2

fy1  = spwert(f%n_x,x,f%Aspline(iy1,:),f%Bspline(iy1,:),f%Cspline(iy1,:),&
                      f%Dspline(iy1,:),f%xspline,ABLTG1)
fy2  = spwert(f%n_x,x,f%Aspline(iy2,:),f%Bspline(iy2,:),f%Cspline(iy2,:),&
                      f%Dspline(iy2,:),f%xspline,ABLTG2)

if (present(fout)) &
  fout = (fy1 - fy2) / (f%ylinear(iy1) - f%ylinear(iy2)) * (y - f%ylinear(iy1)) + fy1
if (present(dfout_dx)) &
  dfout_dx   = (ABLTG1(1) - ABLTG2(1)) / (f%ylinear(iy1) - f%ylinear(iy2)) * (y - f%ylinear(iy1)) + ABLTG1(1)
if (present(d2fout_dx2)) &
  d2fout_dx2 = (ABLTG1(2) - ABLTG2(2)) / (f%ylinear(iy1) - f%ylinear(iy2)) * (y - f%ylinear(iy1)) + ABLTG1(2)
if (present(dfout_dy)) &
  dfout_dy   = (fy1 - fy2) / (f%ylinear(iy1) - f%ylinear(iy2))
end subroutine SL2Dinterp

end module mod_interp_splinear
