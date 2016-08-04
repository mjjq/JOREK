!> Particle input-output module, containing hdf5 data_type and writing routines
module mod_particle_io
  use mod_particle_base
  use hdf5, only: HID_T

  type :: particle_hdf5_io
    class(particle_base), public, allocatable, dimension(:) :: particle_type !< Example of this particle
    integer(HID_T), public                                  :: data_type !< HDF5 data_type
  contains
    procedure, pass, public :: set_data_type
    procedure, pass, non_overridable, public :: write => export_particles
    procedure, pass, non_overridable, public :: read  => import_particles
  end type
contains

!> Create hdf5 data type.
!> The code below ([[import_particles]] and [[export_particles]]) is slightly 
!> illegal according to the standard.
!> - interoperability between C and fortran is not supported for polymorphism
!> - we don't know if all of the particles will follow eachother in memory
!> but it seems to work in ifort and in gfortran with a workaround for C_LOC.
!> This can be removed as soon as gfortran relaxes the restrictions on C_LOC.
subroutine set_data_type(this)
use iso_c_binding
use hdf5
use mod_particle_boris !< Gfortran workaround for C_LOC not allowing polymorphism
implicit none
class(particle_hdf5_io), intent(inout) :: this !< Object-bound argument
integer                     :: hdferr
integer(HID_T)              :: st_array, x_array
integer(HSIZE_T), parameter :: st_dim(1) = (/2/) !< JOREK integration
integer(HSIZE_T), parameter :: x_dim(1) = (/3/) !< JOREK integration
!class(particle), pointer, dimension(:) :: p
integer(HSIZE_T), dimension(0:8) :: offsets

! ugly workaround, remove as soon as gfortran supports the 2012 interop TS
! for now, copy this code for each particle type expected
select type(p => this%particle_type)
type is (particle_boris)
  offsets(0) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(2))) ! full size
  offsets(1) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%x))
  offsets(2) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%mass))
  offsets(3) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%weight))
  offsets(4) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%st))
  offsets(5) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%i_elm))
  offsets(6) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%q))
  offsets(7) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%label))
  offsets(8) = H5OFFSETOF(C_LOC(p(1)), C_LOC(p(1)%lost))
end select

! Reinitialize the library
call h5open_f(hdferr)

! Create the compound this%data_type
call h5tcreate_f(H5T_COMPOUND_F, offsets(0), &
                                 this%data_type, hdferr)

call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, st_dim, st_array, hdferr)
call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, x_dim, x_array, hdferr)

! Fill type
call h5tinsert_f(this%data_type, "x [m] at time t", &
     offsets(1), x_array, hdferr)
call h5tinsert_f(this%data_type, "mass [atomic mass units]", &
     offsets(2), H5T_NATIVE_REAL, hdferr)
call h5tinsert_f(this%data_type, "weight (number of particles)", &
     offsets(3), H5T_NATIVE_REAL, hdferr)


call h5tinsert_f(this%data_type, "st", offsets(4), st_array, hdferr)
call h5tinsert_f(this%data_type, "i_elm", offsets(5), H5T_NATIVE_INTEGER, hdferr)


call h5tinsert_f(this%data_type, "q [electron charges]", offsets(6), &
     h5kind_to_type(kind(this%particle_type(1)%q),H5_INTEGER_KIND), hdferr)

call h5tinsert_f(this%data_type, "label", offsets(7), &
     h5kind_to_type(kind(this%particle_type(1)%label),H5_INTEGER_KIND), hdferr)

call h5tinsert_f(this%data_type, "lost", offsets(8), &
     h5kind_to_type(kind(this%particle_type(1)%lost),H5_INTEGER_KIND), hdferr)
end subroutine set_data_type




!> Export all particles using HDF5 Parallel File io
subroutine export_particles(this, filename, particles)
use mpi
use hdf5
use mod_particle_boris !< Gfortran workaround for C_LOC not allowing polymorphism
implicit none

class(particle_hdf5_io), intent(in)                       :: this
character*(*)          , intent(in)                       :: filename !< File to dump particles in
class(particle_base)   , intent(in), target, dimension(:) :: particles

integer :: my_id, n_cpu, ierr
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 writing
integer(HID_T)                :: file, file_space, mem_space, dset, plist ! handles
character*(*), parameter      :: dataset_name = 'particles' ! maybe make this an argument?
integer                       :: hdferr
type(c_ptr) :: p_ptr

! Preparation
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)

allocate(particles_per_proc(0:n_cpu-1))

! Find the number of particles on each node
call MPI_AllGather(size(particles,1),1,MPI_INTEGER,&
    particles_per_proc,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

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

! ugly workaround, remove as soon as gfortran supports the 2012 interop TS
! for now, copy this code for each particle type expected
select type (p => particles)
type is (particle_boris)
  p_ptr = C_LOC(p(1))
end select
! Write the dataset independently
call h5dwrite_f(dset, this%data_type, p_ptr, &
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
!> particle distribution over all processors and read this many
!> particles per processor.
subroutine import_particles(this, filename, particles)
use mpi
use hdf5
use mod_particle_boris !< Gfortran workaround for C_LOC not allowing polymorphism
implicit none

class(particle_hdf5_io), intent(in)                                     :: this
character*(*)          , intent(in)                                     :: filename
class(particle_base)   , intent(out), dimension(:), allocatable, target :: particles

integer                            :: my_id, n_cpu, ierr, rank, n_particles
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 reading
integer(HID_T)                :: file, file_space, mem_space, dset, plist ! handles
character*(*), parameter      :: dataset_name = 'particles' ! maybe make this an argument?
integer                       :: hdferr

type(c_ptr) :: p_ptr
integer*8, dimension(1:1) :: tmp, maxdims

! Preparation
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)

allocate(particles_per_proc(0:n_cpu-1))

! Create access property list for parallel access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdferr)
call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

! Open the file
call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr, access_prp=plist)
call h5pclose_f(plist, hdferr)

! Open the dataset
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

! And allocate the required space
! allocate(particles(particles_per_proc(my_id)), source=this%particle_type, stat=ierr)
! the above fails in gfortran (allocates to size of this%particle_type), workaround below
select type(p => this%particle_type)
type is (particle_boris)
  allocate(particle_boris::particles(particles_per_proc(my_id)), stat=ierr)
end select
if (ierr .gt. 0) write(*,"(i3,a,i12,a)") my_id, &
    "unable to allocate particles(", particles_per_proc(my_id), ")"
call h5screate_simple_f(1, (/size(particles,dim=1,kind=HSIZE_T)/), mem_space, hdferr)

! Select hyperslab in the file (offset only, no stride)
call h5sselect_hyperslab_f(file_space, H5S_SELECT_SET_F, &
    start=(/int(sum(particles_per_proc(0:my_id-1)),HSIZE_T)/), &
    count=(/size(particles,dim=1,kind=HSIZE_T)/), &
    hdferr=hdferr, stride=(/1_HSIZE_T/), block=(/1_HSIZE_T/))

! Read the dataset independently
! ugly workaround, remove as soon as gfortran supports the 2012 interop TS
! for now, copy this code for each particle type expected
select type (p => particles)
type is (particle_boris)
  p_ptr = C_LOC(p(1))
end select
call h5dread_f(dset, this%data_type, p_ptr, &
    hdferr, mem_space_id=mem_space, file_space_id=file_space)

! Close everything
call h5sclose_f(file_space, hdferr)
call h5sclose_f(mem_space, hdferr)
call h5dclose_f(dset, hdferr)
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)

end subroutine import_particles
end module mod_particle_io
