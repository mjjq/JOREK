!> This module contains some routines for calculating diagnostics on particles
module mod_particle_diagnostics
use mod_io_actions
use data_structure
use mod_particle_sim
use hdf5
implicit none
private
public write_constants_of_motion

!> Action to calculate pphi_H_mu and write this to an HDF5 file
!> in an extensible (in the time-dimension) dataset
type, extends(io_action) :: write_constants_of_motion
  integer(HID_T) :: file_id !< file identifier
contains
  procedure :: do => do_write_constants_of_motion
end type write_constants_of_motion
interface write_constants_of_motion
  module procedure new_write_constants_of_motion
end interface write_constants_of_motion


contains
!> Constructor. Must use this or open the HDF5 file manually.
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_write_constants_of_motion(filename) result(new)
  use hdf5_io_module
  use mpi
  type(write_constants_of_motion) :: new
  character(len=*), intent(in)    :: filename
  integer :: my_id, n_cpu, ierr
  integer(HID_T) :: group_id
  logical :: link_exists
  new%filename = filename
  new%name = "WriteConstantsOfMotion"
  new%log = .true.
  ! Open existing HDF5 file
  call HDF5_open_or_create(filename, parallel=.true., file_id=new%file_id, ierr=ierr)
  if (ierr .ne. 0) then
    write(*,*) "ERROR: cannot open HDF5 file"
    call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
  end if

  ! Create group to write particle groups in if it does not exist yet
  call h5lexists_f(new%file_id, "/groups", link_exists, ierr)
  if (.not. link_exists) then
    call h5gcreate_f(new%file_id, "/groups", group_id, ierr)
    call h5gclose_f(group_id, ierr)
  end if
  ! assume that if it exists it's a group
end function new_write_constants_of_motion

!> Action to calculate all of these values and write them to an HDF5 file
subroutine do_write_constants_of_motion(this, sim)
  use mpi
  !$ use omp_lib
  class(write_constants_of_motion), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  integer :: i, my_id, n_cpu, ierr, rank
  integer(HSIZE_T), parameter :: n_time = 1_HSIZE_T, n_var = 3_HSIZE_T
  integer(HSIZE_T)  :: obj_count, data_dims(3), time_dims(1), data_maxdims(3), time_maxdims(1) ! equality with parameters above is coincidence
  integer(HID_T)    :: dspace, dset, plist, mem_space
  integer(HID_T)    :: group_id
  integer(HID_T)    :: tspace, tset
  real*8 :: last_time
  logical :: link_exists
  character(len=80) :: dataset_name, timeset_name
  integer, dimension(:), allocatable :: particles_per_proc
  real*8, dimension(:,:,:), allocatable, target :: stats ! Data storage order: particle index, variable number, time (in fortran)

  ! Safety checks
  if (.not. allocated(sim%groups)) return
  ! test if handle is open, give error about using constructor
  call h5fget_obj_count_f(this%file_id, H5F_OBJ_FILE_F, obj_count, ierr)
  if (ierr .ne. 0 .or. obj_count .eq. 0) then
    write(*,*) "file does not seem to be opened, did you call the constructor?"
    call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
  endif

  ! Preparation
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
  allocate(particles_per_proc(0:n_cpu-1))

  ! For each of the groups
  do i=lbound(sim%groups,1),ubound(sim%groups,1)
    ! Find the number of particles on each node
    call MPI_AllGather(size(sim%groups(i)%particles,1),1,MPI_INTEGER,&
        particles_per_proc,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    if (allocated(stats)) deallocate(stats)
    allocate(stats(size(sim%groups(i)%particles,1),3,1))

    ! Check the dataset existence and properties
    write(dataset_name,'(A,i0.3)') 'groups/', i
    call h5lexists_f(this%file_id, trim(dataset_name), link_exists, ierr)
    write(*,*) "Dataset exists? ", link_exists

    if (link_exists) then
      call h5dopen_f(this%file_id, trim(dataset_name), dset, ierr)
      write(*,*) "trying to open ierr=", ierr
      if (ierr .ne. 0) then
        write(*,*) "Error opening dataset", i
        call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
      elseif (.not. dataset_is_in_diag_format(dset)) then
        write(*,*) "ERROR: dataset " // trim(dataset_name) // " in wrong format"
        call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
      else
        call h5dget_space_f(dset, dspace, ierr)
      end if
    else
      call create_constants_dataset(this%file_id, dataset_name, &
          int(sum(particles_per_proc,1),HSIZE_T), dset, dspace)
    end if

    ! Check the timeset existence and properties
    write(timeset_name,'(A,i0.3,A)') 'groups/', i, '_t'
    call h5lexists_f(this%file_id, trim(timeset_name), link_exists, ierr)
    write(*,*) "Timeset exists? ", link_exists

    if (link_exists) then
      call h5dopen_f(this%file_id, trim(timeset_name), tset, ierr)
      write(*,*) "trying to open ierr=", ierr
      if (ierr .ne. 0) then
        write(*,*) "Error opening timeset", i
        call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
      elseif (.not. timeset_is_in_diag_format(dset)) then
        write(*,*) "ERROR: timeset " // trim(timeset_name) // " in wrong format"
        call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
      else
        call h5dget_space_f(tset, tspace, ierr)
      end if
    else
      call create_constants_time_dataset(this%file_id, trim(timeset_name), &
          tset, tspace)
    end if

    ! Get the current sizes
    call h5sget_simple_extent_dims_f(dspace, data_dims, data_maxdims, ierr)
    call h5sget_simple_extent_dims_f(tspace, time_dims, time_maxdims, ierr)
    if (time_dims(1) .ne. data_dims(3)) then
      write(*,*) "ERROR: data and time series are not of equal length"
      call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
    end if

    ! Check the last value in the time dataset
    ! We assume that time is monotonous (in this dataset)
    ! so that the last value is the largest.
    ! If this is not true raise an error and exit
    call h5soffset_simple_f(tspace, time_dims(1), ierr)
    last_time = -1d99
    call h5dread_f(tset, H5T_NATIVE_DOUBLE, last_time, [n_time], ierr, tspace)
    if (last_time .gt. sim%time) then
      write(*,*) "ERROR: overwriting existing timesteps not supported yet"
      call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
    end if

    ! Extend the dataset by 1 in the time-dimension
    data_dims(3) = data_dims(3) + 1_HSIZE_T
    time_dims(1) = time_dims(1) + 1_HSIZE_T
    data_dims(1) = max(data_dims(1), int(sum(particles_per_proc,1),HSIZE_T)) ! only grow this set
    call h5dset_extent_f(dset, data_dims, ierr)
    call h5dset_extent_f(tset, time_dims, ierr)
    ! Apparently after this tspace and dspace are invalid!

    ! Create dataspace for file and memory separately
    call h5screate_simple_f(3, [size(sim%groups(i)%particles,dim=1,kind=HSIZE_T),n_var,n_time], mem_space, ierr)

    call h5dget_space_f(dset, dspace, ierr)
    call h5sselect_none_f(dspace, ierr)
    ! Select hyperslab in the file (offset only, no stride)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, &
        start=[int(sum(particles_per_proc(0:my_id-1)),HSIZE_T), 0_HSIZE_T, data_dims(3)-1_HSIZE_T], &
        count=[size(sim%groups(i)%particles,dim=1,kind=HSIZE_T), n_var, n_time], &
        hdferr=ierr, stride=[1_HSIZE_T,1_HSIZE_T,1_HSIZE_T], block=[1_HSIZE_T])

    ! Calculate the statistics
    call calculate_pphi_H_mu(sim%fields, sim%time, sim%groups(i)%particles, sim%groups(i)%mass, stats(:,:,1))

    ! Write the dataset independently
    data_dims = size(stats,kind=HSIZE_T)
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, stats, data_dims, &
         ierr, file_space_id = dspace, mem_space_id = mem_space)
    call h5sclose_f(mem_space, ierr)
    call h5sclose_f(dspace, ierr)
    call h5dclose_f(dset, ierr)
    ! close everything

    ! Add the current time to the timeset
    ! TODO

  end do
end subroutine do_write_constants_of_motion

!> Create a new dataset for diagnostics with the right dimensions in file_id
subroutine create_constants_dataset(file_id, dataset_name, n_particles, dset, dspace)
  integer(HID_T), intent(in)   :: file_id
  character(len=*), intent(in) :: dataset_name
  integer(HID_T), intent(out)  :: dset, dspace
  integer(HSIZE_T), intent(in) :: n_particles
  integer :: ierr
  integer(HID_T) :: crp_list
  integer(HSIZE_T), parameter :: chunk_size(3) = [100,3,1]

  ! Create a dataspace with unlimited dimensions, 
  call h5screate_simple_f(3, [n_particles,3_HSIZE_T,0_HSIZE_T], dspace, ierr, &
      maxdims=[H5S_UNLIMITED_F,3_HSIZE_T,H5S_UNLIMITED_F])
  ! Create a dataset property list, enable chunking
  call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, ierr)
  call h5pset_chunk_f(crp_list, 3, chunk_size, ierr)
  ! Create a dataset with initial dimensions 0,n_particles,n_var (c-style)
  call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_DOUBLE, dspace, dset, ierr, crp_list)

  ! TODO: test options we have for compression
  ! ZFP filter?
  ! bitshuffle?
  ! shuffle filter + zlib (usually present)
  ! For now run without, test difference after creating the dataset
end subroutine create_constants_dataset

!> Create a new dataset for time data (1-d extensible, chunked (required for extensibility))
subroutine create_constants_time_dataset(file_id, dataset_name, dset, dspace)
  integer(HID_T), intent(in)   :: file_id
  character(len=*), intent(in) :: dataset_name
  integer(HID_T), intent(out)  :: dset, dspace
  integer :: ierr
  integer(HID_T) :: crp_list
  integer(HSIZE_T), parameter :: chunk_size(1) = [100]

  ! Create a dataspace with unlimited dimensions, 
  call h5screate_simple_f(1, [0_HSIZE_T], dspace, ierr, &
      maxdims=[H5S_UNLIMITED_F])
  ! Create a dataset property list, enable chunking
  call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, ierr)
  call h5pset_chunk_f(crp_list, 1, chunk_size, ierr)
  call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_DOUBLE, dspace, dset, ierr, crp_list)
end subroutine create_constants_time_dataset

!> Check whether a dataset has the right format
function dataset_is_in_diag_format(dset) result(correct)
  integer(HID_T), intent(in) :: dset
  integer(HID_T) :: dspace
  integer :: rank, ierr
  integer(HSIZE_T) :: dims(3), max_dims(3)
  logical :: correct
  ! Get the dataspace for this group
  ! Get the number of dimensions for this space
  call h5dget_space_f(dset, dspace, ierr)
  if (ierr .ne. 0) call h5sget_simple_extent_ndims_f(dspace, rank, ierr)
  if (rank .ne. 3 .or. ierr .ne. 0) then ! if dataspace could not be found or has wrong dimensionality
    write(*,*) "rank", rank, ierr
    correct = .false.
  else
    call h5sget_simple_extent_dims_f(dspace, dims, max_dims, ierr)
    ! Require at least 3 dimensions for variables (last dimension)
    ! particles, conserved quantity, time
    write(*,*) "dims", dims, "max_dims", max_dims
    if (max_dims(1) .ne. H5S_UNLIMITED_F .or. max_dims(3) .ne. H5S_UNLIMITED_F .or. max_dims(2) .lt. 3) then
      correct = .false.
    else
      correct = .true.
    end if
  endif
end function dataset_is_in_diag_format

!> Check whether a timeset has the right format
function timeset_is_in_diag_format(dset) result(correct)
  integer(HID_T), intent(in) :: dset
  integer(HID_T) :: dspace
  integer :: rank, ierr
  integer(HSIZE_T) :: dims(1), max_dims(1)
  logical :: correct
  ! Get the dataspace for this group
  ! Get the number of dimensions for this space
  call h5dget_space_f(dset, dspace, ierr)
  if (ierr .ne. 0) call h5sget_simple_extent_ndims_f(dspace, rank, ierr)
  if (rank .ne. 1 .or. ierr .ne. 0) then ! if dataspace could not be found or has wrong dimensionality
    correct = .false.
  else
    call h5sget_simple_extent_dims_f(dspace, dims, max_dims, ierr)
    correct = max_dims(1) .eq. H5S_UNLIMITED_F
  endif
end function timeset_is_in_diag_format

!> Calculate P_phi, H and mu for a list of particles.
!> mask is .f. if a particle is lost. These values in out are 0.d0
subroutine calculate_pphi_H_mu(fields, time, particles, mass, out, mask)
  use mod_particle_types
  use phys_module, only: F0
  use constants
  class(fields_base), intent(in)                      :: fields
  real*8, intent(in)                                  :: time
  class(particle_base), intent(in), dimension(:)      :: particles
  real*8, intent(in)                                  :: mass
  real*8, dimension(3,size(particles,1)), intent(out) :: out !< List of values
  logical, dimension(size(particles,1)), intent(out), optional :: mask !< Mask containing .f. if particle is lost
  real*8, dimension(1) :: P, P_s, P_t, P_phi, P_time
  real*8               :: inv_st_jac, psi_R, psi_Z, B(3), B_hat(3), B_norm
  real*8               :: R, R_s, R_t, Z, Z_s, Z_t

  integer :: i

  if (present(mask)) mask = .true.
  out  = 0.d0
  !$omp parallel do default(none) &
  !$omp shared(particles, fields, out, mask, time, f0, mass) &
  !$omp private(P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t, &
  !$omp inv_st_jac, psi_R, psi_Z, B, B_hat)
  do i=1,size(particles,1)
    if (particles(i)%i_elm .lt. 1) then
      if (present(mask)) mask(i) = .false.
    else
      select type (pa => particles(i))
      !> Should we calculate this for the guiding-centers?
      type is (particle_kinetic_leapfrog)
        ! Calculate psi and B
        call fields%interp_PRZ(time, pa%i_elm, &
                      [1], 1, pa%st(1),pa%st(2), &
                      pa%x(3), P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
        inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
        psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
        psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
        ! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
        B        = [+psi_Z, -psi_R, F0] / R
        B_hat = B/norm2(B)

        ! Calculate output variables
        out(1,i) = real(pa%q,8) * EL_CHG * P(1) + mass * ATOMIC_MASS_UNIT * R * pa%v(3)
        out(2,i) = mass * ATOMIC_MASS_UNIT * 0.5d0 * dot_product(pa%v,pa%v)
        out(3,i) = mass * ATOMIC_MASS_UNIT * 0.5d0 * dot_product(&
            pa%v - dot_product(pa%v,B_hat)*B_hat, &
            pa%v - dot_product(pa%v,B_hat)*B_hat &
            )/norm2(B)
      class default
        write(*,*) "ERROR: calculate_pphi_H_mu not implemented for this particle type"
      end select
    endif
  enddo
  !$omp end parallel do
end subroutine calculate_pphi_H_mu



!> Calculate particles present in specific regions on all particles
!> Performs MPI communications to sum values, returns the value
!> corresponding to all particles on node 0, and the value for each node on this node
!> Regions are: DOMAIN_PLASMA, DOMAIN_SOL, DOMAIN_OUTER_SOL,
!> DOMAIN_UPPER_PRIVATE, DOMAIN_LOWER_PRIVATE
function particles_in_regions(node_list, element_list, particles)
  use data_structure
  use phys_module, only: DOMAIN_PLASMA, DOMAIN_SOL, DOMAIN_OUTER_SOL, DOMAIN_UPPER_PRIVATE,        &
      DOMAIN_LOWER_PRIVATE, xpoint, xcase
  use mod_particle_types
  use domains
  use mpi
  implicit none

  type(type_node_list), intent(in)     :: node_list
  type(type_element_list), intent(in)  :: element_list
  class(particle_base), intent(in), dimension(:) :: particles

  integer, dimension(DOMAIN_PLASMA:DOMAIN_LOWER_PRIVATE) :: particles_in_regions, tmp
  integer :: i, ifail, my_id
  integer :: domain, i_elm_axis, i_elm_xpoint(2)
  real*8  :: psi, psi_s, psi_t, psi_st, psi_ss, psi_tt
  real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis, psi_limit
  real*8, dimension(2) :: psi_xpoint, R_xpoint, Z_xpoint, s_xpoint, t_xpoint


  !! Preparation (force my_id to 1 to suppress message)
  call find_axis(1,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

  if (xpoint) then
    call find_xpoint(1,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
    psi_limit  = psi_xpoint(1)
    if( (xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
      psi_limit = psi_xpoint(2)
    endif
  else
    psi_limit = 0.d0
  endif

  ! Call which_domain once to setup saved values
  domain = which_domain(node_list, element_list, &
      0.d0, 0.d0, &
      0.d0, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
      R_axis, Z_axis, psi_axis)

  tmp = 0
  !$omp parallel do default(none) &
  !$omp shared(node_list, element_list, particles, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
  !$omp     R_axis, Z_axis, psi_axis) &
  !$omp private(domain, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt) &
  !$omp reduction(+:tmp)
  do i=1,size(particles,1)
    associate (p => particles(i))
    if (p%i_elm .le. 0 .or. p%i_elm .gt. element_list%n_elements) cycle
    call interp(node_list, element_list, p%i_elm, 1, 1, & ! force i_harm to 1
        p%st(1), p%st(2), psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)

    domain = which_domain(node_list, element_list, &
        p%x(1), p%x(2), &
        psi, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
        R_axis, Z_axis, psi_axis)
    end associate

    tmp(domain) = tmp(domain) + 1
  end do
  !$omp end parallel do

  ! Save values on nodes
  particles_in_regions = tmp
  ! Mpi communication to get the total answer on node 0
  call MPI_Reduce(tmp, particles_in_regions, DOMAIN_LOWER_PRIVATE-DOMAIN_PLASMA, &
    MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ifail)
  call MPI_Comm_Rank(MPI_COMM_WORLD, my_id, ifail)
end function particles_in_regions
end module mod_particle_diagnostics
