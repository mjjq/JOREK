!> Module for filling hdf5 datatype with particle_base parameters
module mod_particle_base_hdf5
  use mod_particle_base
contains

!> Create hdf5 data type for a passed-in particle
subroutine h5_datatype_base(particles, datatype)
use iso_c_binding
use hdf5
implicit none

class(particle_base), intent(in), dimension(:), target :: particles
integer(HID_T), intent(out) :: datatype !< HDF5 type

integer                     :: hdferr
integer(HID_T)              :: st_array, x_array

integer, parameter :: n_dim = 2

! Create the compound datatype
call h5tcreate_f(H5T_COMPOUND_F, H5OFFSETOF(C_LOC(particles(1)), &
                                            C_LOC(particles(2))), &
                                 datatype, hdferr)

call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, (/int(n_dim,HSIZE_T)/), st_array, hdferr)
call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, (/3_HSIZE_T/), x_array, hdferr)


! Fill type
call h5tinsert_f(datatype, "x [m] at time t", &
     H5OFFSETOF(C_LOC(particle(1)),C_LOC(particle(1)%x)), x_array_mem, hdferr)

call h5tinsert_f(datatype, "mass [atomic mass units]", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%mass)), &
     H5T_NATIVE_REAL, hdferr)

call h5tinsert_f(datatype, "weight (number of particles)", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%weight)), &
     H5T_NATIVE_REAL, hdferr)


call h5tinsert_f(datatype, "st", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%st)), st_array_mem, hdferr)

call h5tinsert_f(datatype, "i_elm", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%i_elm)), &
     H5T_NATIVE_INTEGER, hdferr)


call h5tinsert_f(datatype, "q [electron charges]", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%q)), &
     h5kind_to_type(kind(p%particle(1)%q),H5_INTEGER_KIND), hdferr)

call h5tinsert_f(datatype, "label", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%label)), &
     h5kind_to_type(kind(p%particle(1)%label),H5_INTEGER_KIND), hdferr)

call h5tinsert_f(datatype, "lost", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%lost)), &
     h5kind_to_type(kind(p%particle(1)%lost),H5_INTEGER_KIND), hdferr)

end subroutine h5_dtypes
end module mod_particle_base_hdf5
