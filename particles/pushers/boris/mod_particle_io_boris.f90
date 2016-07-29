!> Particle pusher module with the Boris scheme
module mod_particle_io_boris
  use mod_particle_boris
  use mod_particle_io

  !> HDF5 IO type for boris method particles.
  !> Usage example: 
  !>```fortran
  !> class(particle_hdf5_io), allocatable :: io
  !> allocate(particle_hdf5_io_boris::io); allocate(particle_boris::io%particle_type)
  !> call io%set_data_type
  !>```
  !=> <<<
  !> @note expect no output from this test
  type, extends(particle_hdf5_io) :: particle_hdf5_io_boris
  contains
    procedure, pass, public :: set_data_type => set_data_type_boris
  end type
contains

!> Create hdf5 data type for a passed-in particle of type particle_boris
subroutine set_data_type_boris(this)
use hdf5
implicit none
integer        :: hdferr
integer(HID_T) :: v_array

call this%parent%set_data_type

call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, (/3_HSIZE_T/), v_array, hdferr)

call h5tinsert_f(this%data_type, "v [m/s] at time t-1/2dt", &
     H5OFFSETOF(C_LOC(this%particle),C_LOC(this%particle%v)), v_array, hdferr)
end subroutine set_data_type_boris
end module mod_particle_io_boris
