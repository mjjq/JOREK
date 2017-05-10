!> Particle input-output module, containing hdf5 data_type and writing routines
!> TODO: add metadata and/or use H5MD format (http://nongnu.org/h5md/h5md.html)
module mod_particle_io
use hdf5_io_module
use hdf5
use mpi
use mod_particle_types
use mod_particle_sim
implicit none
private
public write_simulation_hdf5, read_simulation_hdf5

integer(HSIZE_T), parameter :: particle_type_name_length = 40 !< length of the string used to identify a specific type of particle
character(len=*), parameter :: particle_type_name_field_name = 'particle_type' !< Name of the field containing the particle_type_name
contains

!> Export all particles using HDF5 Parallel File IO
subroutine write_simulation_hdf5(sim, filename)
type(particle_sim)   , intent(in) :: sim
character*(*)        , intent(in) :: filename

integer :: my_id, n_cpu, ierr
integer :: n_here, n_total
integer(HSIZE_T) :: i_here
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 writing
integer(HID_T)                :: file, create_file_space, write_file_space, dset, plist ! handles
integer(HID_T)                :: group_id
integer(HID_T)                :: data_type
integer(HID_T)                :: time_space_id, time_set_id
character(len=12)             :: group_name
character(len=particle_type_name_length) :: particle_type_name
integer                       :: i, j, hdferr
type(c_ptr) :: p_ptr
real*8, dimension(:,:), allocatable :: real8_2D
real*8, dimension(:), allocatable   :: real8_1D
real*4, dimension(:), allocatable   :: real4_1D
integer, dimension(:), allocatable  :: int4_1D

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


! Write the time
call HDF5_real_saving(file,sim%time,'/time')


if (allocated(sim%groups)) then
  do i=1,size(sim%groups,1)
    if (.not. allocated(sim%groups(i)%particles)) then
      write(*,*) "WARNING: group ", i, " not allocated, exiting"
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    ! Find the number of particles on each node
    n_here = size(sim%groups(i)%particles,1)
    call MPI_AllGather(n_here,1,MPI_INTEGER,&
        particles_per_proc,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    n_total = sum(particles_per_proc)
    i_here = sum(particles_per_proc(0:my_id-1))

    ! Create group to write in
    write(group_name,"(A,i0.3,A)") "/groups/", i, "/"
    call h5gcreate_f(file, group_name, group_id, hdferr)
    call h5gclose_f(group_id, hdferr)

    ! particle_base properties
    ! x
    allocate(real8_2D(3,n_here))
    do j=1,n_here
      real8_2D(:,j) = sim%groups(i)%particles(j)%x
    end do
    call HDF5_array2D_saving(file,real8_2D,3,n_total,group_name//"x",start=[0_HSIZE_T,i_here])
    deallocate(real8_2D)

    ! st
    allocate(real8_2D(2,n_here))
    do j=1,n_here
      real8_2D(:,j) = sim%groups(i)%particles(j)%st
    end do
    call HDF5_array2D_saving(file,real8_2D,2,n_total,group_name//"st",start=[0_HSIZE_T,i_here])
    deallocate(real8_2D)

    ! weight
    allocate(real4_1D(n_here))
    do j=1,n_here
      real4_1D(j) = sim%groups(i)%particles(j)%weight
    end do
    call HDF5_array1D_saving_r4(file,real4_1D,n_total,group_name//"weight",start=[i_here])
    deallocate(real4_1D)

    ! i_elm
    allocate(int4_1D(n_here))
    do j=1,n_here
      int4_1D(j) = sim%groups(i)%particles(j)%i_elm
    end do
    call HDF5_array1D_saving_int(file,int4_1D,n_total,group_name//"i_elm",start=[i_here])
    deallocate(int4_1D)

    ! Write out stuff depending on particle type
    select type (p => sim%groups(i)%particles)
    type is (particle_kinetic)
      particle_type_name = 'particle_kinetic'

      ! v
      allocate(real8_2D(3,n_here))
      do j=1,n_here
        real8_2D(:,j) = p(j)%v
      end do
      call HDF5_array2D_saving(file,real8_2D,3,n_total,group_name//"v",start=[0_HSIZE_T,i_here])
      deallocate(real8_2D)

      ! q
      allocate(int4_1D(n_here))
      do j=1,n_here
        int4_1D(j) = p(j)%q
      end do
      call HDF5_array1D_saving_int(file,int4_1D,n_total,group_name//"q",start=[i_here])
      deallocate(int4_1D)
    type is (particle_kinetic_leapfrog)
      particle_type_name = 'particle_kinetic_leapfrog'

      ! v
      allocate(real8_2D(3,n_here))
      do j=1,n_here
        real8_2D(:,j) = p(j)%v
      end do
      call HDF5_array2D_saving(file,real8_2D,3,n_total,group_name//"v",start=[0_HSIZE_T,i_here])
      deallocate(real8_2D)

      ! q
      allocate(int4_1D(n_here))
      do j=1,n_here
        int4_1D(j) = p(j)%q
      end do
      call HDF5_array1D_saving_int(file,int4_1D,n_total,group_name//"q",start=[i_here])
      deallocate(int4_1D)
    type is (particle_gc)
      particle_type_name = 'particle_gc'

      ! E
      allocate(real8_1D(n_here))
      do j=1,n_here
        real8_1D(j) = p(j)%E
      end do
      call HDF5_array1D_saving(file,real8_1D,n_total,group_name//"E",start=[i_here])
      deallocate(real8_1D)

      ! mu
      allocate(real8_1D(n_here))
      do j=1,n_here
        real8_1D(j) = p(j)%mu
      end do
      call HDF5_array1D_saving(file,real8_1D,n_total,group_name//"mu",start=[i_here])
      deallocate(real8_1D)

      ! q
      allocate(int4_1D(n_here))
      do j=1,n_here
        int4_1D(j) = p(j)%q
      end do
      call HDF5_array1D_saving_int(file,int4_1D,n_total,group_name//"q",start=[i_here])
      deallocate(int4_1D)
    class default
      write(*,*) "error: missing type name declaration for write"
      call exit(1)
    end select

    call HDF5_char_saving(file,particle_type_name,group_name//"type")
  end do
end if

! Close everything
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)

if (my_id .eq. 0) write(*,*) "Writing particle output file to ", filename, " succeeded"
end subroutine write_simulation_hdf5




!> Import all particles.
!> Reads the number of particles from a file, determines a
!> particle distribution over all processors and read this many
!> particles per processor.
subroutine read_simulation_hdf5(sim, filename)
type(particle_sim) , intent(out) :: sim
character*(*)      , intent(in)  :: filename

integer                            :: my_id, n_cpu, ierr, rank, n_particles
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 reading
integer(HID_T)    :: file, file_space, mem_space, dset, plist ! handles
integer(HID_T)    :: data_type
integer(HID_T)    :: group_id
integer(HID_T)    :: time_set_id
integer(HSIZE_T)  :: i_here
integer           :: n_here
integer           :: storage_type, max_corder
character(len=12) :: group_name
character(len=particle_type_name_length) :: particle_type_name
integer           :: i, j, n, hdferr

type(c_ptr) :: p_ptr
integer*8, dimension(1:2) :: tmp, maxdims
real*8, dimension(:,:), allocatable :: real8_2D
integer*4, dimension(:), allocatable :: int4_1D
real*4, dimension(:), allocatable :: real4_1D
real*8, dimension(:), allocatable :: real8_1D

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

! read the time
call HDF5_real_reading(file,sim%time,'/time')

! Get the number of groups
call h5gopen_f(file, '/groups/', group_id, hdferr)
call h5gget_info_f(group_id, storage_type, n, max_corder, hdferr)
call h5gclose_f(group_id, hdferr)

allocate(sim%groups(n))
do i=1,n
  ! Open the dataset for x
  write(group_name,'(A,i0.3,A)') '/groups/', i, '/'
  call h5dopen_f(file, group_name//"x", dset, hdferr)

  ! Open the file dataspace
  call h5dget_space_f(dset, file_space, hdferr)

  ! Get the number of particles
  call h5sget_simple_extent_ndims_f(file_space, rank, hdferr)
  call h5sget_simple_extent_dims_f(file_space, tmp, maxdims, hdferr)
  n_particles = int(tmp(2),4)
  call h5sclose_f(file_space, hdferr)
  call h5dclose_f(dset, hdferr)

  ! Divide particles over processors
  particles_per_proc(0)         = n_particles/n_cpu + modulo(n_particles, n_cpu)
  particles_per_proc(1:n_cpu-1) = n_particles/n_cpu
  i_here = sum(particles_per_proc(0:my_id-1))
  n_here = particles_per_proc(my_id)

  ! Get the particle type from the attribute
  call HDF5_char_reading(file,particle_type_name,group_name//"type")

  ierr = 0
  select case (trim(particle_type_name))
  case ('particle_kinetic')
    allocate(particle_kinetic::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_kinetic_leapfrog')
    allocate(particle_kinetic_leapfrog::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_gc')
    allocate(particle_gc::sim%groups(i)%particles(n_here), stat=ierr)
  case default
    write(*,*) "error: missing type name declaration for read"
    call exit(1)
  end select
  if (ierr .gt. 0) write(*,"(i3,a,i12,a)") my_id, &
      "unable to allocate particles(", particles_per_proc(my_id), ")"

  ! Read base particle attributes
  ! x
  allocate(real8_2D(3,n_here))
  call HDF5_array2D_reading(file, real8_2D, group_name//"x",start=[0_HSIZE_T,i_here])
  do j=1,n_here
    sim%groups(i)%particles(j)%x = real8_2D(:,j)
  end do
  deallocate(real8_2D)

  ! st
  allocate(real8_2D(2,n_here))
  call HDF5_array2D_reading(file, real8_2D, group_name//"st",start=[0_HSIZE_T,i_here])
  do j=1,n_here
    sim%groups(i)%particles(j)%st = real8_2D(:,j)
  end do
  deallocate(real8_2D)

  ! weight
  allocate(real4_1D(n_here))
  call HDF5_array1D_reading_r4(file, real4_1D, group_name//"weight",start=[i_here])
  do j=1,n_here
    sim%groups(i)%particles(j)%weight = real4_1D(j)
  end do
  deallocate(real4_1D)
  ! i_elm
  allocate(int4_1D(n_here))
  call HDF5_array1D_reading_int(file, int4_1D, group_name//"i_elm", start=[i_here])
  do j=1,n_here
    sim%groups(i)%particles(j)%i_elm = int4_1D(j)
  end do
  deallocate(int4_1D)

  select type (p => sim%groups(i)%particles)
  type is (particle_kinetic)
    ! v
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%v = real8_2D(:,j)
    end do
    deallocate(real8_2D)

    ! q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)
  type is (particle_kinetic_leapfrog)
    ! v
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%v = real8_2D(:,j)
    end do
    deallocate(real8_2D)

    ! q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)
  type is (particle_gc)
    ! E
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"E",start=[i_here])
    do j=1,n_here
      p(j)%E = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! mu
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"mu",start=[i_here])
    do j=1,n_here
      p(j)%mu = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)
  end select

end do

! Close everything else
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)

end subroutine read_simulation_hdf5
end module mod_particle_io
