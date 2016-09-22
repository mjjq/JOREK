!> Particle input-output module, containing hdf5 data_type and writing routines
!> TODO: add metadata and/or use H5MD format (http://nongnu.org/h5md/h5md.html)
module mod_particle_io
use hdf5, only: HSIZE_T
use mod_constants, only: CARTESIAN, CYLINDRICAL
public :: write_simulation_hdf5
public :: read_simulation_hdf5
private

integer(HSIZE_T), parameter :: geometry_name_length = 20 !< length of the string used to write geometries (cartesian or cylindrical)
integer(HSIZE_T), parameter :: particle_type_name_length = 40 !< length of the string used to identify a specific type of particle
character(len=*), parameter :: particle_type_name_field_name = 'particle_type' !< Name of the field containing the particle_type_name
contains

!> Create hdf5 data type.
!> The code below ([[import_particles]] and [[export_particles]]) is slightly 
!> illegal according to the standard.
!> - interoperability between C and fortran is not supported for polymorphism
!> - we don't know if all of the particles will follow eachother in memory
!> From "15.2.3.6 C_LOC(X)": (see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=56305)
!> > Argument. X shall have either the POINTER or TARGET attribute. It shall not be a coindexed object. It shall either be a variable with interoperable type and kind type parameters, or be a scalar, nonpolymorphic variable with no length type parameters. If it is allocatable, it shall be allocated. If it is a pointer, it shall be associated. If it is an array, it shall be contiguous and have nonzero size. It shall not be a zero-length string.
!> but it seems to work in ifort and in gfortran with a workaround for C_LOC.
!> TODO:
!> 
!>* Check portability when using only a single datatype instead of a filetype and memtype
function get_hdf5_particle_data_type(particles) result(data_type)
use iso_c_binding
use hdf5
use mod_boris !< Gfortran workaround for C_LOC not allowing polymorphism
implicit none
integer(HID_T)              :: data_type

class(particle_base), dimension(:), intent(in) :: particles
class(particle_base), dimension(:), allocatable :: particles_2
integer                     :: hdferr
integer(HID_T)              :: st_array, x_array
integer(HSIZE_T), parameter :: st_dim(1) = (/2/) !< JOREK integration
integer(HSIZE_T), parameter :: x_dim(1) = (/3/) !< JOREK integration
integer(HSIZE_T), dimension(0:9) :: offsets

! Reallocate to a fixed-size list to allow for single-particle lists
if (size(particles,1) .eq. 0) then
  write(*,*) "ERROR: no particles given"
  call exit(1)
end if
allocate(particles_2(1:2), source=particles(1)) 

! ugly workaround, remove as soon as gfortran supports the 2012 interop TS
! for now, copy this code for each particle type expected
select type(p => particles_2)
type is (particle_boris)
  offsets(0) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(2))) ! full size
  offsets(1) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%x))
  offsets(2) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%m))
  offsets(3) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%weight))
  offsets(4) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%st))
  offsets(5) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%i_elm))
  offsets(7) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%label))
  offsets(6) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%q))
  offsets(8) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%lost))
class default
  write(*,*) "ERROR: unknown particle type for creating hdf5 type"
  call exit(1)
end select

! Reinitialize the library
call h5open_f(hdferr)

! Create the compound data_type
call h5tcreate_f(H5T_COMPOUND_F, offsets(0), &
                                 data_type, hdferr)

call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, st_dim, st_array, hdferr)
call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, x_dim, x_array, hdferr)

! Fill type
call h5tinsert_f(data_type, "x [m] at time t", &
     offsets(1), x_array, hdferr)
call h5tinsert_f(data_type, "m [atomic mass units]", &
     offsets(2), H5T_NATIVE_REAL, hdferr)
call h5tinsert_f(data_type, "weight (number of particles)", &
     offsets(3), H5T_NATIVE_REAL, hdferr)


call h5tinsert_f(data_type, "st", offsets(4), st_array, hdferr)
call h5tinsert_f(data_type, "i_elm", offsets(5), H5T_NATIVE_INTEGER, hdferr)


call h5tinsert_f(data_type, "q [electron charges]", offsets(6), &
     h5kind_to_type(kind(particles(1)%q),H5_INTEGER_KIND), hdferr)

call h5tinsert_f(data_type, "label", offsets(7), &
     h5kind_to_type(kind(particles(1)%label),H5_INTEGER_KIND), hdferr)

call h5tinsert_f(data_type, "lost", offsets(8), &
     h5kind_to_type(kind(particles(1)%lost),H5_INTEGER_KIND), hdferr)


! type-specific fields
select type(p => particles)
type is (particle_boris)
  offsets(9) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%v))
  call h5tinsert_f(data_type, "v [m/s] at time t-1/2 dt", &
     offsets(9), x_array, hdferr)
! add new particle types here
end select
end function get_hdf5_particle_data_type




!> Export all particles using HDF5 Parallel File IO
subroutine write_simulation_hdf5(sim, filename)
use mpi
use hdf5
use mod_particle_sim
use mod_boris !< Gfortran workaround for C_LOC not allowing polymorphism
implicit none

type(particle_sim)   , intent(in) :: sim
character*(*)        , intent(in) :: filename

integer :: my_id, n_cpu, ierr
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 writing
integer(HID_T)                :: file, file_space, mem_space, dset, plist ! handles
integer(HID_T)                :: group_id, attr_id, aspace_id, atype_id
integer(HID_T)                :: data_type
integer(HID_T)                :: time_space_id, time_set_id
integer(HID_T)                :: geometry_space_id, geometry_set_id
character(len=80)             :: dataset_name
character(len=particle_type_name_length) :: particle_type_name
integer                       :: i, hdferr
type(c_ptr) :: p_ptr

! Preparation
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)
allocate(particles_per_proc(0:n_cpu-1))

! Create file property list for parallel access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdferr)
call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)


! Create file
call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr, access_prp=plist)
if (hdferr .gt. 0) then
  write(*,*) "file open failed:", hdferr
  call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
endif
call h5pclose_f(plist, hdferr)


! Create group to write particle groups in
call h5gcreate_f(file, "/groups", group_id, hdferr)
call h5gclose_f(group_id, hdferr)


! Create attribute space to store particle type
call h5screate_simple_f(1, [1_HSIZE_T], aspace_id, hdferr)
! Create a character type of length particle_type_name_length
call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
call h5tset_size_f(atype_id, particle_type_name_length, hdferr)


! Write the time
call h5screate_simple_f(1, [1_HSIZE_T], time_space_id, hdferr)
call h5dcreate_f(file, '/time', H5T_NATIVE_DOUBLE, time_space_id, time_set_id, hdferr)
call h5dwrite_f(time_set_id, H5T_NATIVE_DOUBLE, sim%time, [1_HSIZE_T], hdferr)
call h5dclose_f(time_set_id, hdferr)
call h5sclose_f(time_space_id, hdferr)


! Write the geometry used
! Create a character type of length geometry_name_length
if (allocated(sim%fields)) then
  call h5screate_simple_f(1, [1_HSIZE_T], geometry_space_id, hdferr)
  call h5dcreate_f(file, '/geometry', H5T_NATIVE_INTEGER, geometry_space_id, geometry_set_id, hdferr)
  call h5dwrite_f(geometry_set_id, H5T_NATIVE_INTEGER, sim%fields%geometry, [1_HSIZE_T], hdferr)
  call h5sclose_f(geometry_space_id, hdferr)
  call h5dclose_f(geometry_set_id, hdferr)
end if


do i=1,size(sim%groups,1)
  ! Find the number of particles on each node
  call MPI_AllGather(size(sim%groups(i)%particles,1),1,MPI_INTEGER,&
      particles_per_proc,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

  ! Create dataspace for file and memory separately
  call h5screate_simple_f(1, (/int(sum(particles_per_proc),HSIZE_T)/), file_space, hdferr)
  call h5screate_simple_f(1, (/size(sim%groups(i)%particles,dim=1,kind=HSIZE_T)/), mem_space, hdferr)

  ! Create dataset for file
  write(dataset_name,"(A,i0.3)") "/groups/", i
  data_type = get_hdf5_particle_data_type(sim%groups(i)%particles)
  call h5dcreate_f(file, trim(dataset_name), data_type, file_space, dset, hdferr)
  call h5sclose_f(file_space, hdferr)

  ! Create an attribute for this set with the particle type
  call h5acreate_f(dset, particle_type_name_field_name, atype_id, aspace_id, attr_id, hdferr)
  select type (p => sim%groups(i)%particles)
  type is (particle_boris)
    particle_type_name = 'particle_boris'
  class default
    write(*,*) "error: missing type name declaration for write"
    call exit(1)
  end select
  call h5awrite_f(attr_id, atype_id, particle_type_name, [1_HSIZE_T], hdferr)
  call h5aclose_f(attr_id, hdferr)

  ! Select hyperslab in the file (offset only, no stride)
  call h5dget_space_f(dset, file_space, hdferr)
  call h5sselect_hyperslab_f(file_space, H5S_SELECT_SET_F, &
      start=(/int(sum(particles_per_proc(0:my_id-1)),HSIZE_T)/), &
      count=(/size(sim%groups(i)%particles,dim=1,kind=HSIZE_T)/), &
      hdferr=hdferr, stride=(/1_HSIZE_T/), block=(/1_HSIZE_T/))

  ! ugly workaround, remove as soon as gfortran supports the 2012 interop TS
  ! for now, copy this code for each particle type expected
  select type (p => sim%groups(i)%particles)
  type is (particle_boris)
    p_ptr = C_LOC(p(1))
  end select
  ! Write the dataset independently
  call h5dwrite_f(dset, data_type, p_ptr, &
       hdferr, file_space_id = file_space, mem_space_id = mem_space)
  call h5sclose_f(mem_space, hdferr)
  call h5dclose_f(dset, hdferr)
end do

! Close everything
call h5sclose_f(aspace_id, hdferr)
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)
end subroutine write_simulation_hdf5


!> Import all particles using MPI File this
!> Reads the number of particles from a file, determines a
!> particle distribution over all processors and read this many
!> particles per processor.
subroutine read_simulation_hdf5(sim, filename)
use mod_particle_sim
use mpi
use hdf5
use mod_boris !< Gfortran workaround for C_LOC not allowing polymorphism
implicit none

type(particle_sim) , intent(out) :: sim
character*(*)      , intent(in)  :: filename

integer                            :: my_id, n_cpu, ierr, rank, n_particles
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 reading
integer(HID_T)    :: file, file_space, mem_space, dset, plist ! handles
integer(HID_T)    :: data_type
integer(HID_T)    :: group_id
integer(HID_T)    :: attr_id, atype_id
integer           :: storage_type, max_corder
character(len=80) :: dataset_name
character(len=particle_type_name_length) :: particle_type_name
integer           :: i, n, hdferr

type(c_ptr) :: p_ptr
integer*8, dimension(1:1) :: tmp, maxdims

! Preparation
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)
allocate(particles_per_proc(0:n_cpu-1))

! Create file property list for parallel access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdferr)
call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

! Open the file
call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr, access_prp=plist)
call h5pclose_f(plist, hdferr)

! TODO read time and type of geometry, give an error if geometry differs

call h5gopen_f(file, '/groups/', group_id, hdferr)
call h5gget_info_f(group_id, storage_type, n, max_corder, hdferr)
call h5gclose_f(group_id, hdferr)

allocate(sim%groups(n))
do i=1,n
  ! Open the dataset
  write(dataset_name,'(A,i0.3)') 'groups/', i
  call h5dopen_f(file, dataset_name, dset, hdferr)

  ! Open the file dataspace
  call h5dget_space_f(dset, file_space, hdferr)

  ! Get the number of particles (fails if dataset is not 1-dimensional!)
  call h5sget_simple_extent_ndims_f(file_space, rank, hdferr)
  if (rank .gt. 1) then
    write(*,*) "Reading >1-dimensional data not supported yet!"
    call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
  endif
  call h5sget_simple_extent_dims_f(file_space, tmp, maxdims, hdferr)
  n_particles = int(tmp(1),4)

  ! Divide particles over processors
  particles_per_proc(0)         = n_particles/n_cpu + modulo(n_particles, n_cpu)
  particles_per_proc(1:n_cpu-1) = n_particles/n_cpu

  ! Get the particle type from the attribute
  call h5aopen_f(dset, particle_type_name_field_name, attr_id, hdferr)
  call h5aget_type_f(attr_id, atype_id, hdferr)
  call h5aread_f(attr_id, atype_id, particle_type_name, [1_HSIZE_T], hdferr)
  call h5tclose_f(atype_id, hdferr)
  call h5aclose_f(attr_id, hdferr)


  ierr = 0
  select case (particle_type_name)
  case ('particle_boris')
    allocate(particle_boris::sim%groups(i)%particles(particles_per_proc(my_id)), stat=ierr)
  case default
    write(*,*) "error: missing type name declaration for read"
    call exit(1)
  end select
  if (ierr .gt. 0) write(*,"(i3,a,i12,a)") my_id, &
      "unable to allocate particles(", particles_per_proc(my_id), ")"

  call h5screate_simple_f(1, (/int(particles_per_proc(my_id),kind=HSIZE_T)/), mem_space, hdferr)

  ! Select hyperslab in the file (offset only, no stride)
  call h5sselect_hyperslab_f(file_space, H5S_SELECT_SET_F, &
      start=(/int(sum(particles_per_proc(0:my_id-1)),HSIZE_T)/), &
      count=(/int(particles_per_proc(my_id),kind=HSIZE_T)/), &
      hdferr=hdferr, stride=(/1_HSIZE_T/), block=(/1_HSIZE_T/))

  ! Read the dataset independently
  ! ugly workaround, remove as soon as gfortran supports the 2012 interop TS
  ! for now, copy this code for each particle type expected
  select type (p => sim%groups(i)%particles)
  type is (particle_boris)
    p_ptr = C_LOC(p(1))
  class default
    write(*,*) "ERROR: missing type declaration for read"
  end select

  ! Get the data type for this kind of particle
  data_type = get_hdf5_particle_data_type(sim%groups(i)%particles)
  call h5dread_f(dset, data_type, p_ptr, &
    hdferr, mem_space_id=mem_space, file_space_id=file_space)

  call h5dclose_f(dset, hdferr)
  call h5sclose_f(file_space, hdferr)
  call h5sclose_f(mem_space, hdferr)
end do

! Close everything else
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)

end subroutine read_simulation_hdf5
end module mod_particle_io
