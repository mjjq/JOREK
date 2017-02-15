!> This module contains some routines for calculating diagnostics on particles.
!> Outputs are:
!> 1. Energy
!> 2. Magnetic moment
!> 3. P_phi (generalized toroidal momentum)
!> 4. Psi_bar (P_phi/q)
!> 5. Psi
!> 6. q (charge)
module mod_particle_diagnostics
use mod_io_actions
use mod_particle_sim
use hdf5
implicit none
private
public write_particle_diagnostics, calculate_particle_diagnostics

integer(HSIZE_T), parameter :: n_var = 6

!> Action to calculate pphi_H_mu and write this to an HDF5 file
!> in an extensible (in the time-dimension) dataset
type, extends(io_action) :: write_particle_diagnostics
  integer(HID_T) :: file_id !< file identifier
contains
  procedure :: do => do_write_particle_diagnostics
end type write_particle_diagnostics
interface write_particle_diagnostics
  module procedure new_write_particle_diagnostics
end interface write_particle_diagnostics


contains
!> Constructor. Must use this or open the HDF5 file manually.
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_write_particle_diagnostics(filename) result(new)
  type(write_particle_diagnostics) :: new
  character(len=*), intent(in)    :: filename
  new%filename = filename
  new%name = "WriteConstantsOfMotion"
  new%log = .true.
end function new_write_particle_diagnostics

!> Action to calculate all of these values and write them to an HDF5 file
subroutine do_write_particle_diagnostics(this, sim, ev)
  use hdf5_io_module
  use mpi
  use mod_event
  class(write_particle_diagnostics), intent(inout) :: this
  type(particle_sim), intent(inout)               :: sim
  type(event), intent(inout), optional            :: ev

  integer(HSIZE_T), parameter :: n_time = 1_HSIZE_T

  integer(HSIZE_T)  :: data_dims(3), time_dims(1), data_maxdims(3), time_maxdims(1)
  integer(HSIZE_T)  :: npoints_mem, npoints_file
  integer           :: i, my_id, n_cpu, ierr, rank
  integer(HID_T)    :: dspace, dset, mem_space
  integer(HID_T)    :: tspace, t_mem_space, tset, group_id, plist
  logical           :: link_exists
  character(len=80) :: dataset_name, timeset_name
  integer, dimension(:), allocatable            :: particles_per_proc
  real*4, dimension(:,:,:), allocatable, target :: stats ! Data storage order: particle index, variable number, time (in fortran)

  call h5open_f(ierr)

  ! Create file property list for parallel access
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist, ierr)
  call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)

  call HDF5_open_or_create(this%filename, plist, file_id=this%file_id, ierr=ierr)
  if (ierr .ne. 0) then
    write(*,*) "ERROR: cannot open HDF5 file"
    call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
  end if
  call h5pclose_f(plist, ierr)

  ! Create group to write particle groups in if it does not exist yet
  call h5lexists_f(this%file_id, "/groups", link_exists, ierr)
  if (.not. link_exists) then
    call h5gcreate_f(this%file_id, "/groups", group_id, ierr)
    call h5gclose_f(group_id, ierr)
  end if
  ! assume that if it exists it's a group

  ! Safety checks
  if (.not. allocated(sim%groups)) return

  ! Preparation
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
  allocate(particles_per_proc(0:n_cpu-1))

  ! For each of the groups
  do i=lbound(sim%groups,1),ubound(sim%groups,1)
    ! Find the number of particles on each node
    call MPI_AllGather(size(sim%groups(i)%particles,1),1,MPI_INTEGER,&
        particles_per_proc,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    allocate(stats(size(sim%groups(i)%particles,1),n_var,1))

    ! Check the dataset existence and properties
    write(dataset_name,'(A,i0.3)') 'groups/', i
    call h5lexists_f(this%file_id, trim(dataset_name), link_exists, ierr)

    if (link_exists) then
      !write(*,*) "DEBUG: link to ", trim(dataset_name), " exists, trying to open"
      call h5dopen_f(this%file_id, trim(dataset_name), dset, ierr)
      if (ierr .ne. 0) then
        write(*,*) "Error opening dataset", i
        call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
      else
        call h5dget_space_f(dset, dspace, ierr)
        !write(*,*) "DEBUG: dataset opened, getting space (", ierr, ")"
      end if
    else
      call create_constants_dataset(this%file_id, dataset_name, &
          int(sum(particles_per_proc,1),HSIZE_T), dset, dspace)
      !write(*,*) "DEBUG: created new dataset ", trim(dataset_name)
    end if

    ! Check the timeset existence and properties
    write(timeset_name,'(A,i0.3,A)') 'groups/', i, '_t'
    call h5lexists_f(this%file_id, trim(timeset_name), link_exists, ierr)

    if (link_exists) then
      write(*,*) "DEBUG: link to ", trim(timeset_name), " exists, trying to open"
      call h5dopen_f(this%file_id, trim(timeset_name), tset, ierr)
      if (ierr .ne. 0) then
        write(*,*) "Error opening timeset", i
        call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
      else
        call h5dget_space_f(tset, tspace, ierr)
        !write(*,*) "DEBUG: dataset opened, getting space: ", ierr
      end if
    else
      call create_constants_time_dataset(this%file_id, trim(timeset_name), &
          tset, tspace)
      !write(*,*) "DEBUG: created new timeset ", trim(timeset_name)
    end if

    ! Get the current sizes
    call h5sget_simple_extent_dims_f(dspace, data_dims, data_maxdims, ierr)
    call h5sget_simple_extent_dims_f(tspace, time_dims, time_maxdims, ierr)
    !write(*,"(A,3i6,A,3i6)") " DEBUG: dspace sizes(3): ", data_dims, " maxdims(3): ", data_maxdims
    !write(*,"(A,i6,A,i6)") " DEBUG: tspace sizes(1): ", time_dims, " maxdims(1): ", time_maxdims
    if (time_dims(1) .ne. data_dims(3)) then
      write(*,*) "ERROR: data and time series are not of equal length"
      call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
    end if

    ! Extend the dataset by 1 in the time-dimension
    ! After extending, dspace and tspace are invalid. Close them already
    call h5sclose_f(dspace, ierr)
    call h5sclose_f(tspace, ierr)
    time_dims(1) = time_dims(1) + 1_HSIZE_T
    data_dims(1) = max(data_dims(1), int(sum(particles_per_proc,1),HSIZE_T)) ! only grow this set
    data_dims(2) = n_var
    data_dims(3) = data_dims(3) + 1_HSIZE_T
    !write(*,"(A,3i6)") " DEBUG: extending dset to ", data_dims
    call h5dset_extent_f(dset, data_dims, ierr)
    !write(*,"(A,1i6)") " DEBUG: extending tset to ", time_dims
    call h5dset_extent_f(tset, time_dims, ierr)

    ! Create dataspace for file and memory separately
    !write(*,*) "DEBUG: creating dataspace for memory and file"
    call h5screate_simple_f(3, [size(sim%groups(i)%particles,dim=1,kind=HSIZE_T), n_var, n_time], mem_space, ierr)
    ! We must create a hyperslab space here because simple does not play well
    ! with HDF5 (but does not tell you this!).
    ! It is also important to get the space belonging to this dataset, instead of
    ! creating a new one.
    call h5dget_space_f(dset, dspace, ierr)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, &
        start=[int(sum(particles_per_proc(0:my_id-1)),HSIZE_T), 0_HSIZE_T, data_dims(3)-1_HSIZE_T], &
        count=[1_HSIZE_T,1_HSIZE_T,1_HSIZE_T], &
        block=[size(sim%groups(i)%particles,dim=1,kind=HSIZE_T),n_var,n_time], &
        hdferr=ierr)

    ! Calculate the statistics
    call calculate_particle_diagnostics(sim%fields, sim%time, sim%groups(i)%particles, sim%groups(i)%mass, stats(:,:,1))
    !write(*,*) "DEBUG: calculated statistics"

    ! Write the dataset independently
    data_dims(1) = size(stats,dim=1,kind=HSIZE_T)
    data_dims(2) = n_var
    data_dims(3) = n_time
    !write(*,"(A,3i6)") " DEBUG: writing statistics to block. dims=", data_dims
    call h5dwrite_f(dset, H5T_NATIVE_REAL, stats, [1_HSIZE_T,1_HSIZE_T,1_HSIZE_T], &
         ierr, file_space_id = dspace, mem_space_id = mem_space)
    !write(*,*) "EXTRA_DEBUG: out(8,2)", stats(8,2,1)
    !write(*,*) "EXTRA_DEBUG: out(2,2)", stats(2,2,1)
    write(*,*) "Done writing, sum=", sum(stats) ! output here is to stop gfortran (tried with 6.2.1) optimizing away the result

    ! Add the current time to the timeset
    call h5dget_space_f(tset, tspace, ierr)
    call h5sselect_hyperslab_f(tspace, H5S_SELECT_SET_F, &
        start=[time_dims(1)-1_HSIZE_T], count=[n_time], hdferr=ierr)
    call h5screate_f(H5S_SCALAR_F, t_mem_space, ierr)
    !write(*,*) "DEBUG: writing single time value", sim%time
    call h5dwrite_f(tset, H5T_NATIVE_REAL, sim%time, [n_time], ierr, file_space_id=tspace, mem_space_id=t_mem_space)
    !write(*,*) "Done"

    call h5sclose_f(t_mem_space, ierr)
    call h5sclose_f(mem_space, ierr)
    call h5sclose_f(dspace, ierr)
    call h5sclose_f(tspace, ierr)
    call h5dclose_f(dset, ierr)
    call h5dclose_f(tset, ierr)
    deallocate(stats)
  end do
  call h5fclose_f(this%file_id, ierr)
  call h5close_f(ierr)
end subroutine do_write_particle_diagnostics

!> Create a new dataset for diagnostics with the right dimensions in file_id
subroutine create_constants_dataset(file_id, dataset_name, n_particles, dset, dspace)
  integer(HID_T), intent(in)   :: file_id
  character(len=*), intent(in) :: dataset_name
  integer(HID_T), intent(out)  :: dset, dspace
  integer(HSIZE_T), intent(in) :: n_particles
  integer :: ierr
  integer(HID_T) :: crp_list
  integer(HSIZE_T), parameter :: chunk_size(3) = [10000_HSIZE_T,1_HSIZE_T,1_HSIZE_T]

  ! Create a dataspace with unlimited dimensions, 
  call h5screate_simple_f(3, [n_particles,n_var,0_HSIZE_T], dspace, ierr, &
      maxdims=[H5S_UNLIMITED_F,n_var,H5S_UNLIMITED_F])
  ! Create a dataset property list, enable chunking
  call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, ierr)
  call h5pset_chunk_f(crp_list, 3, chunk_size, ierr)
  call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_REAL, dspace, dset, ierr, crp_list)
  call h5pclose_f(crp_list, ierr)

  ! shuffle filter + zlib (present on most system) (give 20% compression, don't use for)
end subroutine create_constants_dataset

!> Create a new dataset for time data (1-d extensible, chunked (required for extensibility))
subroutine create_constants_time_dataset(file_id, dataset_name, dset, dspace)
  integer(HID_T), intent(in)   :: file_id
  character(len=*), intent(in) :: dataset_name
  integer(HID_T), intent(out)  :: dset, dspace
  integer :: ierr
  integer(HID_T) :: crp_list
  integer(HSIZE_T), parameter :: chunk_size(1) = [1000]

  ! Create a dataspace with unlimited dimensions, 
  call h5screate_simple_f(1, [0_HSIZE_T], dspace, ierr, &
      maxdims=[H5S_UNLIMITED_F])
  ! Create a dataset property list, enable chunking
  call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, ierr)
  call h5pset_chunk_f(crp_list, 1, chunk_size, ierr)
  call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_REAL, dspace, dset, ierr, crp_list)
  call h5pclose_f(crp_list, ierr)
end subroutine create_constants_time_dataset


!> Calculate H, mu, P_phi, Psi_bar, Psi, q for a list of particles.
!> mask is .f. if a particle is lost. These values in out are 0.d0
subroutine calculate_particle_diagnostics(fields, time, particles, mass, out, mask)
  use mod_particle_types
  use phys_module, only: F0
  use constants
  use mod_boris
  use mod_fields_linear
  class(fields_base), intent(in)                               :: fields
  real*8, intent(in)                                           :: time
  class(particle_base), intent(in), dimension(:)               :: particles
  real*8, intent(in)                                           :: mass
  real*4, dimension(:,:), intent(out)                          :: out !< List of values (is actually size(particles,1),n_var big, but gfortran doesn't like that
  logical, dimension(:), intent(out), optional :: mask !< Mask containing .f. if particle is lost
  real*8, dimension(1) :: P, P_s, P_t, P_phi, P_time
  real*8               :: inv_st_jac, psi_R, psi_Z, B(3), B_norm, B_hat(3), v_par, v_perp(3)
  real*8               :: R, R_s, R_t, Z, Z_s, Z_t
  type(particle_gc)    :: particle

  integer :: i

  if (present(mask)) mask = .true.
  out  = 0.d0
  !$omp parallel do default(none) &
  !$omp shared(particles, fields, out, mask, time, f0, mass) &
  !$omp private(P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t, &
  !$omp inv_st_jac, psi_R, psi_Z, B, B_norm, B_hat, particle, v_par, v_perp)
  do i=1,size(particles,1)
    if (particles(i)%i_elm .lt. 1) then
      if (present(mask)) mask(i) = .false.
    else
      ! Calculate psi and B in the current particle location (either GC or kinetic)
      call fields%interp_PRZ(time, particles(i)%i_elm, &
                    [1], 1, particles(i)%st(1),particles(i)%st(2), &
                    particles(i)%x(3), P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
      inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
      psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
      psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
      ! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
      B        = [+psi_Z, -psi_R, F0] / R

      select type (particle_in => particles(i))
      type is (particle_kinetic_leapfrog)
        out(i,3) = particle_in%q * EL_CHG * P(1) + mass * ATOMIC_MASS_UNIT * R * particle_in%v(3)
        ! Let the conversion calculate the conserved quantities
        particle = kinetic_leapfrog_to_gc(fields%node_list, fields%element_list, particle_in, B, mass)

        if (particle%i_elm .eq. 0) cycle

        ! This recalculates P in the gc position also
        call fields%interp_PRZ(time, particle%i_elm, &
                      [1], 1, particle%st(1),particle%st(2), &
                      particle%x(3), P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)

      type is (particle_gc)
        v_par    = sign(sqrt(2*(particle%E-particle%mu*norm2(B))*EL_CHG/(mass*ATOMIC_MASS_UNIT)),particle%mu)
        particle = particle_in
        out(i,3) = real(particle_in%q,8) * EL_CHG * P(1) + mass * ATOMIC_MASS_UNIT * R * v_par * B(3)/norm2(B)
      class default
        write(*,*) "ERROR: calculate_pphi_H_mu not implemented for this particle type"
        cycle ! skip this iteration
      end select


      ! Calculate output variables
      ! 1. Energy
      out(i,1) = particle%E
      ! 2. Magnetic moment
      out(i,2) = particle%mu
      ! 3. P_phi (generalized toroidal momentum)
      out(i,3) = out(i,3) / EL_CHG ! normalize to ZPsi
      ! 4. Psi_bar (P_phi/q)
      out(i,4) = out(i,3) / real(particle%q)
      ! 5. Psi (at GC position)
      out(i,5) = P(1)
      ! 6. q (charge)
      out(i,6) = real(particle%q)
    endif
  enddo
  !$omp end parallel do
end subroutine calculate_particle_diagnostics


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
