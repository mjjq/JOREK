!> The mod_fast_camera_io module contains procedures
!> for writing and reading fast camera data 
!> in from HDF5 file
module mod_fast_camera_io
#ifdef USE_HDF5
implicit none

private 
public :: write_pixel_intensity_hdf5

contains

!> Procedures first order integrator-----------------------------
!> Write pixel intensity / filter structure in HDF5 file
!> inputs:
!> outputs:
subroutine write_pixel_intensity_hdf5(filename,n_spectra,&
n_values,n_pixel_x,n_pixel_y,n_times,pixel_filter_array,ierr)
  use hdf5
  use hdf5_io_module, only: HDF5_open_or_create,HDF5_close
  use hdf5_io_module, only: HDF5_array5D_saving
  implicit none
  !> Inputs:
  character(len=*),intent(in) :: filename
  integer,intent(in)          :: n_spectra,n_values,n_pixel_x
  integer,intent(in)          :: n_pixel_y,n_times
  real*8,dimension(n_spectra,n_values,n_pixel_x,n_pixel_y,n_times),intent(out) :: pixel_filter_array
  !> Outputs:
  integer,intent(out)         :: ierr
  !> Variables:
  integer(HID_T)              :: file_id
  !> open / store / close image
  call HDF5_open_or_create((trim(filename)//'.h5'),H5P_DEFAULT_F,file_id,ierr,H5F_ACC_TRUNC_F)
  call HDF5_array5D_saving(file_id,pixel_filter_array,n_spectra,&
  n_values,n_pixel_x,n_pixel_y,n_times,'pixel_filter_intensities')
  call HDF5_close(file_id)
end subroutine write_pixel_intensity_hdf5

!>---------------------------------------------------------------
#endif
end module mod_fast_camera_io

