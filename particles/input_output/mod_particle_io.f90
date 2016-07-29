!> Particle input-output module, containing hdf5 data_type and writing routines
module mod_particle_io
  use mod_particle_base
  use hdf5, only: HID_T

  type :: particle_hdf5_io
    class(particle_base), public, allocatable, dimension(:) :: particle_type !< Example of this particle
    integer(HID_T), private                                 :: data_type !< HDF5 data_type
  contains
    procedure, pass, public :: set_data_type
    procedure, pass, non_overridable, public :: write => export_particles
    procedure, pass, non_overridable, public :: read  => import_particles
  end type
contains

!> Create hdf5 data type
subroutine set_data_type(this)
use iso_c_binding
use hdf5
implicit none
class(particle_hdf5_io), intent(inout) :: this !< Object-bound argument
integer                     :: hdferr
integer(HID_T)              :: st_array, x_array
integer, parameter :: n_dim = 2 !< JOREK integration

! Create the compound this%data_type
call h5tcreate_f(H5T_COMPOUND_F, H5OFFSETOF(C_LOC(this%particle_type(1)), &
                                            C_LOC(this%particle_type(2))), &
                                 this%data_type, hdferr)

call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, (/int(n_dim,HSIZE_T)/), st_array, hdferr)
call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, (/3_HSIZE_T/), x_array, hdferr)


! Fill type
call h5tinsert_f(this%data_type, "x [m] at time t", &
     H5OFFSETOF(C_LOC(this%particle_type(1)),C_LOC(this%particle_type(1)%x)), x_array, hdferr)

call h5tinsert_f(this%data_type, "mass [atomic mass units]", &
     H5OFFSETOF(C_LOC(this%particle_type(1)),C_LOC(this%particle_type(1)%mass)), &
     H5T_NATIVE_REAL, hdferr)

call h5tinsert_f(this%data_type, "weight (number of particles)", &
     H5OFFSETOF(C_LOC(this%particle_type(1)),C_LOC(this%particle_type(1)%weight)), &
     H5T_NATIVE_REAL, hdferr)


call h5tinsert_f(this%data_type, "st", &
     H5OFFSETOF(C_LOC(this%particle_type(1)),C_LOC(this%particle_type(1)%st)), st_array, hdferr)

call h5tinsert_f(this%data_type, "i_elm", &
     H5OFFSETOF(C_LOC(this%particle_type(1)),C_LOC(this%particle_type(1)%i_elm)), &
     H5T_NATIVE_INTEGER, hdferr)


call h5tinsert_f(this%data_type, "q [electron charges]", &
     H5OFFSETOF(C_LOC(this%particle_type(1)),C_LOC(this%particle_type(1)%q)), &
     h5kind_to_type(kind(this%particle_type(1)%q),H5_INTEGER_KIND), hdferr)

call h5tinsert_f(this%data_type, "label", &
     H5OFFSETOF(C_LOC(this%particle_type(1)),C_LOC(this%particle_type(1)%label)), &
     h5kind_to_type(kind(this%particle_type(1)%label),H5_INTEGER_KIND), hdferr)

call h5tinsert_f(this%data_type, "lost", &
     H5OFFSETOF(C_LOC(this%particle_type(1)),C_LOC(this%particle_type(1)%lost)), &
     h5kind_to_type(kind(this%particle_type(1)%lost),H5_INTEGER_KIND), hdferr)
end subroutine set_data_type




!> Export all particles using HDF5 Parallel File io
subroutine export_particles(this, filename, particles)
use mpi
use hdf5
implicit none

class(particle_hdf5_io), intent(in)                          :: this
character*(*)          , intent(in)                          :: filename !< File to dump particles in
class(particle_base)   , intent(inout), target, dimension(:) :: particles

integer :: my_id, n_cpu, info, ierr
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 writing
integer(HID_T)                :: file, file_space, mem_space, dset, plist ! handles
character*(*), parameter      :: dataset_name = 'particles' ! maybe make this an argument?
integer                       :: hdferr

! Preparatthisn
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)

allocate(particles_per_proc(0:n_cpu-1))

! Find the number of particles on each node
call MPI_AllGather(size(particles,1),1,MPI_INTEGER,&
    particles_per_proc,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

! Create file property list for parallel access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdferr)
call h5pset_fapl_mpthis_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

! Create file
call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr, access_prp=plist)
if (hdferr .gt. 0) then
  write(*,*) "file open failed:", hdferr
  call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
endif
call h5pclose_f(plist, hdferr)

! Create dataspace for file and memory separately
call h5screate_simple_f(1, (/int(sum(particles_per_proc),HSIZE_T)/), file_space, hdferr)
call h5screate_simple_f(1, (/size(particles,dim=1,kind=HSIZE_T)/), mem_space, hdferr)

! Create dataset for file
call h5dcreate_f(file, dataset_name, this%data_type, file_space, dset, hdferr)
call h5sclose_f(file_space, hdferr)

! Select hyperslab in the file (offset only, no stride)
call h5dget_space_f(dset, file_space, hdferr)
call h5sselect_hyperslab_f(file_space, H5S_SELECT_SET_F, &
    start=(/int(sum(particles_per_proc(0:my_id-1)),HSIZE_T)/), &
    count=(/size(particles,dim=1,kind=HSIZE_T)/), &
    hdferr=hdferr, stride=(/1_HSIZE_T/), block=(/1_HSIZE_T/))

! Write the dataset independently
call h5dwrite_f(dset, this%data_type, C_LOC(particles(1)), &
     hdferr, file_space_id = file_space, mem_space_id = mem_space)

! Close everything
call h5sclose_f(file_space, hdferr)
call h5sclose_f(mem_space, hdferr)
call h5dclose_f(dset, hdferr)
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)
end subroutine export_particles


!> Import all particles using MPI File this
!> Reads the number of particles from a file, determines a
!> particle distributthisn over all processors and read this many
!> particles per processor.
subroutine import_particles(this, filename, particles)
use mpi
use hdf5
implicit none

class(particle_hdf5_io), intent(in)                                     :: this
character*(*)          , intent(in)                                     :: filename
class(particle_base)   , intent(out), dimension(:), allocatable, target :: particles

integer                            :: my_id, n_cpu, ierr, rank, n_particles, info
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 reading
integer(HID_T)                :: file, file_space, mem_space, dset, plist ! handles
character*(*), parameter      :: dataset_name = 'particles' ! maybe make this an argument?
integer                       :: hdferr

type(c_ptr) :: f_ptr
integer*8, dimension(1:1) :: tmp, maxdims

! Preparatthisn
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)

allocate(particles_per_proc(0:n_cpu-1))

! Create access property list for parallel access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdferr)
call h5pset_fapl_mpthis_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

! Open the file
call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr, access_prp=plist)
call h5pclose_f(plist, hdferr)

! Open the dataset
call h5dopen_f(file, dataset_name, dset, hdferr)

! Open the file dataspace
call h5dget_space_f(dset, file_space, hdferr)

! Get the number of particles (fails if dataset is not 1-dimensthisnal!)
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

! And allocate the required space
allocate(particles(particles_per_proc(my_id)), source=this%particle_type, stat=ierr)
if (ierr .gt. 0) write(*,"(i3,a,i12,a)") my_id, &
    "unable to allocate particles(", particles_per_proc(my_id), ")"
call h5screate_simple_f(1, (/size(particles,dim=1,kind=HSIZE_T)/), mem_space, hdferr)

! Select hyperslab in the file (offset only, no stride)
call h5sselect_hyperslab_f(file_space, H5S_SELECT_SET_F, &
    start=(/int(sum(particles_per_proc(0:my_id-1)),HSIZE_T)/), &
    count=(/size(particles,dim=1,kind=HSIZE_T)/), &
    hdferr=hdferr, stride=(/1_HSIZE_T/), block=(/1_HSIZE_T/))

! Read the dataset independently
f_ptr = C_LOC(particles(1))
call h5dread_f(dset, this%data_type, f_ptr, &
    hdferr, mem_space_id=mem_space, file_space_id=file_space)

! Close everything
call h5sclose_f(file_space, hdferr)
call h5sclose_f(mem_space, hdferr)
call h5dclose_f(dset, hdferr)
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)

end subroutine import_particles
end module mod_particle_io
