module mod_import_export_particles
  use mod_particles
contains

function info_romio_locking_disabled() result(info)
  use mpi_mod
  implicit none
  integer :: info, ierr
  call mpi_info_create(info, ierr)
  call mpi_info_set(info, "romio_ds_write", "disable", ierr)
  call mpi_info_set(info, "romio_ds_read", "disable", ierr)
end function info_romio_locking_disabled


subroutine get_particle_hdf5_datatypes(memtype, filetype)
use iso_c_binding
use hdf5
integer(HID_T), intent(out) :: memtype, filetype

integer                     :: hdferr
integer(HID_T)              :: st_array_mem, x_array_mem, &
                               st_array_file, x_array_file

type(type_particle_list), target :: p


p%n_particles = 2
allocate(p%particle(p%n_particles))

!
! Create the compound datatype for file and memory
! We must calculate the offsets manually
!
call h5tcreate_f(H5T_COMPOUND_F, H5OFFSETOF(C_LOC(p%particle(1)), &
                                            C_LOC(p%particle(2))), &
                                 memtype, hdferr)
call h5tcreate_f(H5T_COMPOUND_F, INT(8*n_dim + 2*3*8 + 3*4 + 3*1 + 1, size_t), &
                                 filetype, hdferr)
! Extra + 1 is for aligning in file

call h5tarray_create_f(H5T_NATIVE_DOUBLE, &
    1, (/int(n_dim,HSIZE_T)/), st_array_mem, hdferr)
call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, (/int(n_dim,HSIZE_T)/), st_array_file, hdferr)

call h5tarray_create_f(H5T_NATIVE_DOUBLE, &
    1, (/3_HSIZE_T/), x_array_mem, hdferr)
call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, (/3_HSIZE_T/), x_array_file, hdferr)


! Fill memory and file types
call h5tinsert_f(memtype, "st", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%st)), st_array_mem, hdferr)
call h5tinsert_f(filetype, "st", &
     0_size_t, st_array_file, hdferr)

call h5tinsert_f(memtype, "x (R,Z,phi) [m] at time t", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%x)), x_array_mem, hdferr)
call h5tinsert_f(filetype, "x (R,Z,phi) [m] at time t", &
     int(n_dim*8, size_t), x_array_file, hdferr)

call h5tinsert_f(memtype, "v (R,Z,phi) [m/s * sqrt(mu0 rho0)] at time t-dt/2", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%v)), x_array_mem, hdferr)
call h5tinsert_f(filetype, "v (R,Z,phi) [m/s * sqrt(mu0 rho0)] at time t-dt/2", &
     int(n_dim*8 + 3*8, size_t), x_array_file, hdferr)

call h5tinsert_f(memtype, "mass [atomic mass units]", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%mass)), &
     H5T_NATIVE_REAL, hdferr)
call h5tinsert_f(filetype, "mass [atomic mass units]", &
     int(n_dim*8 + 2*3*8, size_t), H5T_NATIVE_REAL, hdferr)

call h5tinsert_f(memtype, "weight (number of particles)", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%weight)), &
     H5T_NATIVE_REAL, hdferr)
call h5tinsert_f(filetype, "weight (number of particles)", &
     int(n_dim*8 + 2*3*8 + 4, size_t), H5T_NATIVE_REAL, hdferr)

call h5tinsert_f(memtype, "i_elm", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%i_elm)), &
     H5T_NATIVE_INTEGER, hdferr)
call h5tinsert_f(filetype, "i_elm", &
     int(n_dim*8 + 2*3*8 + 4 + 4, size_t), H5T_NATIVE_INTEGER, hdferr)

call h5tinsert_f(memtype, "q [electron charges]", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%q)), &
     h5kind_to_type(kind(p%particle(1)%q),H5_INTEGER_KIND), hdferr)
call h5tinsert_f(filetype, "q [electron charges]", &
     int(n_dim*8 + 2*3*8 + 4 + 4 + 4, size_t), H5T_STD_I8LE, hdferr)

call h5tinsert_f(memtype, "label", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%label)), &
     h5kind_to_type(kind(p%particle(1)%label),H5_INTEGER_KIND), hdferr)
call h5tinsert_f(filetype, "label", &
     int(n_dim*8 + 2*3*8 + 4 + 4 + 4 + 1, size_t), H5T_STD_I8LE, hdferr)

call h5tinsert_f(memtype, "lost", &
     H5OFFSETOF(C_LOC(p%particle(1)),C_LOC(p%particle(1)%lost)), &
     h5kind_to_type(kind(p%particle(1)%lost),H5_INTEGER_KIND), hdferr)
call h5tinsert_f(filetype, "lost", &
     int(n_dim*8 + 2*3*8 + 4 + 4 + 4 + 1 + 1, size_t), H5T_STD_I8LE, hdferr)

! Pack the filetype to make it more efficient (does not help much)
!call h5tpack_f(filetype, hdferr)
end subroutine get_particle_hdf5_datatypes


!> Export all particles using HDF5 Parallel File IO
subroutine export_particles(particle_file, particle_list)
use mpi_mod
use hdf5
use iso_c_binding
implicit none

character*(*)            , intent(in)         :: particle_file
type (type_particle_list), intent(inout), target :: particle_list

integer              :: my_id, n_cpu, info, ierr
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 writing
integer(HID_T)                :: file, filetype, memtype, filespace, memspace, dset, plist ! handles
character*(*), parameter      :: dataset_name = 'particles'
integer                       :: hdferr

! Debugging output
real*8 :: t_start, t_end

! Preparation
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)
call get_particle_hdf5_datatypes(memtype, filetype)

allocate(particles_per_proc(0:n_cpu-1))

if (my_id .eq. 0) then
  write(*,*) '***********************************'
  write(*,*) '*       export particles          *'
  call cpu_time(t_start)
endif

! Find the number of particles on each node
call MPI_AllGather(particle_list%n_particles,1,MPI_INTEGER,particles_per_proc,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

! Create file property list for parallel access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdferr)
call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

! Create file
call h5fcreate_f(particle_file, H5F_ACC_TRUNC_F, file, hdferr, access_prp=plist)
if (hdferr .gt. 0) then
  write(*,*) "file open failed:", hdferr
  call MPI_Abort(MPI_COMM_WORLD)
endif
call h5pclose_f(plist, hdferr)

! Create dataspace for file and memory separately
call h5screate_simple_f(1, (/int(sum(particles_per_proc),HSIZE_T)/), filespace, hdferr)
call h5screate_simple_f(1, (/int(particle_list%n_particles,HSIZE_T)/), memspace, hdferr)

! Create dataset for file
call h5dcreate_f(file, dataset_name, filetype, filespace, dset, hdferr)
call h5sclose_f(filespace, hdferr)

! Select hyperslab in the file (offset only, no stride) (maybe use other dimension for particle type?)
call h5dget_space_f(dset, filespace, hdferr)
call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
    start=(/int(sum(particles_per_proc(0:my_id-1)),HSIZE_T)/), count=(/int(particle_list%n_particles,HSIZE_T)/), &
    hdferr=hdferr, stride=(/1_HSIZE_T/), block=(/1_HSIZE_T/))

! Write the dataset independently
call h5dwrite_f(dset, filetype, C_LOC(particle_list%particle(1)), &
     hdferr, file_space_id = filespace, mem_space_id = memspace)

! Close everything
call h5tclose_f(filetype, hdferr)
call h5tclose_f(memtype, hdferr)
call h5sclose_f(filespace, hdferr)
call h5sclose_f(memspace, hdferr)
call h5dclose_f(dset, hdferr)
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)

if (my_id .eq. 0) then
  call cpu_time(t_end)
  write(*,"(A,f8.4,A)") ' *       elapsed time: ',t_end-t_start, '    *'
  write(*,*) '***********************************'
endif
end subroutine export_particles


!> Import all particles using MPI File IO
!! Reads the number of particles from a file, determines a
!! particle distribution over all processors and read this many
!! particles per processor.
subroutine import_particles(particle_file, particle_list)
use mpi_mod
use hdf5
use iso_c_binding
implicit none

character*(*)            , intent(in)          :: particle_file
type (type_particle_list), intent(out), target :: particle_list

integer              :: my_id, n_cpu, ierr, rank, n_particles, info
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 reading
integer(HID_T)                :: file, filetype, memtype, filespace, memspace, dset, plist ! handles
character*(*), parameter      :: dataset_name = 'particles'
integer                       :: hdferr
real*8 :: t_start, t_end
type(c_ptr) :: f_ptr
integer*8, dimension(1:1) :: tmp, maxdims

! Preparation
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)
call get_particle_hdf5_datatypes(memtype, filetype)

allocate(particles_per_proc(0:n_cpu-1))

if (my_id .eq. 0) then
  write(*,*) '***********************************'
  write(*,*) '*       import particles          *'
  call cpu_time(t_start)
endif

! Create access property list for parallel access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdferr)
call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

! Open the file
call h5fopen_f(particle_file, H5F_ACC_RDONLY_F, file, hdferr, access_prp=plist)
call h5pclose_f(plist, hdferr)

! Open the dataset
call h5dopen_f(file, dataset_name, dset, hdferr)

! Open the file dataspace
call h5dget_space_f(dset, filespace, hdferr)

! Get the number of particles (fails if dataset is not 1-dimensional!)
call h5sget_simple_extent_ndims_f(filespace, rank, hdferr)
if (rank .gt. 1) then
  write(*,*) "Reading >1-dimensional data not supported yet!"
  call MPI_Abort(MPI_COMM_WORLD)
endif
call h5sget_simple_extent_dims_f(filespace, tmp, maxdims, hdferr)
n_particles = int(tmp(1),4)

! Divide particles over processors
particles_per_proc(0)         = n_particles/n_cpu + modulo(n_particles, n_cpu)
particles_per_proc(1:n_cpu-1) = n_particles/n_cpu

particle_list%n_particles = particles_per_proc(my_id)
! And allocate the required space
allocate(particle_list%particle(particle_list%n_particles), stat=ierr)
if (ierr .gt. 0) write(*,"(i3,a,i12,a)") my_id, &
    "unable to allocate particle_list%particle(", particle_list%n_particles, ")"
call h5screate_simple_f(1, (/int(particle_list%n_particles,HSIZE_T)/), memspace, hdferr)

! Select hyperslab in the file (offset only, no stride)
call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
    start=(/int(sum(particles_per_proc(0:my_id-1)),HSIZE_T)/), &
    count=(/int(particle_list%n_particles,HSIZE_T)/), &
    hdferr=hdferr, stride=(/1_HSIZE_T/), block=(/1_HSIZE_T/))

! Read the dataset independently
f_ptr = C_LOC(particle_list%particle(1))
call h5dread_f(dset, memtype, f_ptr, &
    hdferr, mem_space_id=memspace, file_space_id=filespace)

! Close everything
call h5tclose_f(filetype, hdferr)
call h5tclose_f(memtype, hdferr)
call h5sclose_f(filespace, hdferr)
call h5sclose_f(memspace, hdferr)
call h5dclose_f(dset, hdferr)
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)

if (my_id .eq. 0) then
  call cpu_time(t_end)
  write(*,"(A,f8.4,A)") ' *       elapsed time: ',t_end-t_start, '    *'
  write(*,*) '*       particles imported        *'
  write(*,*) '***********************************'
endif
end subroutine import_particles
end module mod_import_export_particles
