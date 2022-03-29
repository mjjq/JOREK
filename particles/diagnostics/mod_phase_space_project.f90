!> Module containing functions for phase-space diagnostics
!> In usage, quite similar to particle projections, with the addition
!> of not only value projection functions but also grid projection functions
!> as to calculate the position of the particle on the grid & coordinates.
module mod_phase_space_project
  use mod_io_actions
  use data_structure
  use mod_particle_sim
  use mod_particle_types
  use mod_fields
  use mod_project_particles
  use hdf5
  implicit none
  private
  public phase_space_projection
  public new_phase_space_projection, output_phase_project, calc_index_shaped_part_x, project_single_particle_x

  ! A type to project particle quantities on arbitrary coordinates. Examples include projecting on
  ! invariant phase-space coordinates (E, \lambda, P_\phi) or real-space 1D histograms of initialization
  ! quantities.
  type, extends(io_action) :: phase_space_projection

    ! As the amount of dimensions is not specified and the resolution is not necessarily the same in all directions
    ! one large 1D array will be used both for the n-dimensional grid and the meshgrid. The indices then can be calculated seperately.
    real*8,       dimension(:),   allocatable         :: grids,dx
    real*8,       dimension(:),   allocatable         :: values

    ! Dimensions and size of values array.
    integer                                           :: ndim, val_size
    ! Resolution of each dimensions, index in grids array where the grid starts for each dimension
    ! and increments of the indices of each dimension in the values array.
    integer,      dimension(:),   allocatable         :: res(:), previndex(:), mult_index(:)

    ! Functions used to calculate the projected quantity and the grids/coordinates on which
    ! it is projected
    type(proj_f)                                      :: f_proj
    type(proj_f), dimension(:),   allocatable         :: f_grids

    ! Quantities for defined particle shapes.
    integer,     dimension(:),    allocatable         :: support,prevsupp, multsupp(:)
    real*8,      dimension(:),    allocatable         :: bandwidths

    integer,     dimension(:,:),  allocatable         :: totsupp_to_nD
    integer                                           :: totsupport,sumsupport, my_id
    logical                                           :: support_present

contains
procedure :: do => project_phase_space
  end type phase_space_projection

contains
! Constructor of the phase-space projection.
!  > ndim             : amount of dimensions
!  > res(ndim)        : resolution (amount of points) of the grids in each dimension
!  > start(ndim)      : starting points of each grid
!  > end(ndim)        : end points of each grid
!  > bandwidths(ndim) : ndim bandwidths (shape extent) in each dimension (optional, but highly recommended. Otherwise only nearest-neighbour)
!  > f_proj           : function that calculates the projected quantity (defined in mod_particle_projection)
!                      (only needed if you want the option to call with(sim,phase_space_proj), but this is most of the time unncessary and slow)
!  > f_grids(ndim)    : ndim functions that calculates the positions in each dimension
!                      (only needed if you want the option to call with(sim,phase_space_proj), but this is most of the time unncessary and slow)
! The reasons to use the shaped particles are threefold:
!   1.) This will smooth the result somewhat.
!   2.) The particles are then kernels, which allows for integration to the exact projected quantity (e.g. for power exchange, the integrated projection is the total exchanged power).
!   3.) Decoupling grid and particle projection.
! Nearest-neighbour projection is not very accurate and smooth, but if you have large amount of particles and not too many gridpoints it can be used.

function new_phase_space_projection(sim,ndim,res,start,end,bandwidths,f_proj,f_grids) result(new)
  type(phase_space_projection)                          :: new
  type(particle_sim),          intent(in)               :: sim
  integer,                     intent(in)               :: ndim
  integer,                     intent(in)               :: res(ndim)
  real*8,                      intent(in)               :: start(ndim), end(ndim)
  real*8,                      intent(in), optional     :: bandwidths(ndim)
  type(proj_f),                intent(in), optional     :: f_proj
  type(proj_f),                intent(in), optional     :: f_grids(ndim)
  integer                                               :: it, jt, gridpoints, valuepoints,j,mult_tmp, my_id
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,it)
  new%my_id = my_id
  ! Allocating the members of the phase-space projection. Have to allocate here due to varying ndim.
  allocate(new%res(ndim))
  allocate(new%previndex(ndim))
  allocate(new%mult_index(ndim))
  allocate(new%dx(ndim))
  allocate(new%f_grids(ndim))
  if(sim%my_id .eq. 0) then
    write(*,"(A,I1,A)") " PARTICLES: Phase space projection with ",ndim, " dimensions."
    if(ndim>7) then
      ! I cannot foresee a situation where this would come up, but for completeness sake.
      write(*,*) "PARTICLES: ", ndim, "dimensions is not supported by some compilers. "
    endif
  endif



  new%support_present=.false.
  new%ndim = ndim
  new%res = res
  if(present(f_proj) .and. present(f_grids)) then
    new%f_proj  = f_proj
    new%f_grids = f_grids
    if(sim%my_id.eq. 0) write(*,*)"PARTICLES: Initialised projection functions"
  else
    if(sim%my_id.eq. 0) write(*,*)"PARTICLES: No projection function initialised. "
  endif

  ! Calculate the total values needed for each dimension
  gridpoints=sum(res)
  ! Calculate the total grid points (on a meshgrid)
  valuepoints=product(res)

  ! Array size of 1D array containing all values.
  new%val_size=valuepoints

  ! This is the array which contains the starting point of each grid.
  new%previndex=0

  ! Allocated on all mpi processes as this is not very large most of the time, reduction can be done later
  allocate(new%grids(gridpoints))
  allocate(new%values(valuepoints))
  new%values = 0.d0
  new%grids  = 0.d0

  ! Pre-calculate dimensions multiplication factor. First dimension varies the slowest,
  ! last the quickest. Therefore, the index in the values array that is incremented by
  ! the first dimension increasing one is the multiplication of all other dimension
  ! and so on.
  mult_tmp=valuepoints
  do it=1,ndim
    mult_tmp=mult_tmp/res(it)
    new%mult_index(it)=mult_tmp
    if(it>1) then
      new%previndex(it)=sum(res(1:it-1)) ! Starting point of each of the grids in the grids array
    endif

    ! Fill the grids uniformly
    do j=1,res(it)
      new%grids(new%previndex(it)+j)=start(it)+(end(it)-start(it))/(res(it)-1)*(j-1)

    enddo ! Grids loop
    new%dx(it)=new%grids(new%previndex(it)+2)-new%grids(new%previndex(it)+1)
  enddo   ! Dimension loop


  if (present(bandwidths)) then
    write(*,*) "PARTICLES: Bandwidths given, particles will now have finite shapes "
    allocate(new%bandwidths(new%ndim))
    allocate(new%support(new%ndim))
    allocate(new%prevsupp(new%ndim))
    allocate(new%multsupp(new%ndim))
    new%support_present = .true.
    new%bandwidths=bandwidths

    new%prevsupp=0
    do it=1,ndim
      new%support(it)=floor(bandwidths(it)/new%dx(it))+1!




      if (it > 1 ) then
        new%prevsupp(it)=sum(new%support(1:(it-1)))
      endif

      if(new%support(it) < 3 .and. sim%my_id .eq. 0) then
        write(*,"(A,I1,A)") " PARTICLES: Support for dimension ", it, " is very low. Will not integrate or project correctly."
      endif
      if(sim%my_id .eq. 0) then
        write(*,"(A,I1,A,E12.4,A,I4)") " PARTICLES: In dimension ",it," the particle width is ",bandwidths(it), " giving a support of",new%support(it)
      endif
    enddo
    new%totsupport = product(new%support)
    new%sumsupport = sum(new%support)
    mult_tmp = new%totsupport
    do it=1,ndim
      mult_tmp=mult_tmp/new%support(it)
      new%multsupp(it)=mult_tmp
    enddo

    allocate(new%totsupp_to_nD(new%ndim,new%totsupport))

    do it=1,new%totsupport
      new%totsupp_to_nD(:, it) = support_to_nD(new,it)
    enddo
    
  endif
  if((.not. new%support_present ).and. (sim%my_id .eq. 0)) then
    write(*,*)"PARTICLES: Trying to use shaped particles without initializing WILL result in undefined behaviour."
    write(*,*) "PARTICLES: Not checked for every projection, use carefully."
  endif

  if(sim%my_id .eq. 0 ) write(*,*) "PARTICLES: Constructed phase space projection"

end function new_phase_space_projection

! > Subroutine for projecting the whole particle distribution
! > on the grids. Not suitable for averaging over single particle
! > motion. Usage can include projecting initialization quantities directly.
! > As the do member of the io-action phase_space_projection points to this
! > function, it is called if calling with(sim, phase_space_event).
subroutine project_phase_space(this, sim, ev)
  use mod_event
  class(phase_space_projection), intent(inout)     :: this
  type(particle_sim), intent(inout)                :: sim
  type(event), intent(inout), optional             :: ev
  integer                                          :: iR
  integer,dimension(:,:),allocatable               :: index_arr
  real*8 ,dimension(:,:),allocatable               :: weight_arr
  real*8 ,dimension(:),allocatable                 :: val_arr

end subroutine project_phase_space

! Output to h5 file. Structure: /values for the meshgrid-evaluated values. /grids/grid_i for 1D grids /grids/mgrid_i for ndim meshgrids
! (can be immediately plotted using these meshgrids)
subroutine output_phase_project(this,ino,output_grids_in)
  class(phase_space_projection), intent(inout)     :: this
  integer, intent(in)                              :: ino
  logical, intent(in), optional                    :: output_grids_in
  real*8, dimension(:), allocatable                :: val_output
  real*8, dimension(:,:), allocatable                :: grid_mesh
  character(len=1024)                              :: filename
  integer                                          :: my_id, ierr,i,index_arr_tmp(this%ndim)
  integer(HID_T)                                   :: file_id, group_id_grid,dspace,dset_id
  integer                                          :: ierrhdf5
  integer                                          :: it,j
  CHARACTER(LEN=8)                                 :: tmp_name
  integer                                          :: res_tmp(7),order_tmp(7),res_tmp2(this%ndim) ! Change this if more dimensions needed
  logical                                          :: output_grids
  !  but this should never come up as particles live in 7D at most.

  if(present(output_grids_in)) then
    output_grids = output_grids_in
  else
    output_grids=.false.
  endif

  ! Output check for too large dimensional arrays for Fortran to handle.
  if (this%ndim > 7) then
    write(*,*) "Can't output for ndim > 7, exiting output"
    return
  endif

  ! For the reshaping of the arrays into ndim arrays. E.g. if 2D, it would look like [res1,res2,0,0,0,0,0]
  ! with HDF5 ignoring the 0 dimensions.
  res_tmp=0
  res_tmp(1:this%ndim)=(/ (this%res(this%ndim-I),I=0,this%ndim-1) /)

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

  ! Allocate only properly on root process
  if (my_id .eq. 0) then
    allocate(val_output, source=this%values)
    allocate(grid_mesh(size(this%values),this%ndim))
  else
    allocate(val_output(0))
    allocate(grid_mesh(0,0))
  endif

  ! Reduce output to the root process for output
  val_output=0.d0
  call MPI_Reduce(this%values,val_output,size(this%values), MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  ! Output to HDF5
  if(my_id .eq. 0 ) then

    !Calculate meshgrids
    do it=1, size(this%values)
      index_arr_tmp = calc_reverse_index_phase_proj(this,it)
      do j=1, this%ndim
        grid_mesh(it,j) = this%grids(index_arr_tmp(j)+this%previndex(j))
      enddo
    enddo


    ! HDF5 file creation
    call h5open_f(ierrhdf5)
    write(filename,"(A5,i5.5,A)") "proj_" ,ino, ".h5"
    call H5Fcreate_f(filename,H5F_ACC_TRUNC_F, file_id, ierrhdf5)
    if(output_grids)then
      call h5gcreate_f(file_id, "grids", group_id_grid, ierrhdf5)
    endif


    res_tmp2=res_tmp(1:this%ndim)

    if(output_grids) then
      do it=1,this%ndim

        ! Output all the 1D grids in each dimensions under the /grids/ group
        write(tmp_name,fmt="(A,I1)") "grid_",it
        call h5screate_simple_f(1, [int(this%res(it),kind=HSIZE_T)], dspace, ierr)!, &
        call h5dcreate_f(group_id_grid, tmp_name, H5T_NATIVE_DOUBLE, dspace, &
                       dset_id, ierrhdf5)
        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,this%grids((this%previndex(it)+1):(this%previndex(it)+this%res(it))),[int(this%res(it),kind=HSIZE_T)],ierr)
        call h5dclose_f(dset_id, ierr)
        call h5sclose_f(dspace,ierr)

        !Same but meshgrid
        call h5screate_simple_f(this%ndim, int(res_tmp2,kind=HSIZE_T),dspace,ierr )
        write(tmp_name,fmt="(A,I1)") "mgrid_",it
        call h5dcreate_f(group_id_grid, tmp_name, H5T_NATIVE_DOUBLE, dspace, &
                       dset_id, ierrhdf5)

        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,RESHAPE(grid_mesh(:,it),res_tmp),int(res_tmp2,kind=HSIZE_T),ierr)
        call h5dclose_f(dset_id,ierr)
        call h5sclose_f(dspace,ierr)

      enddo ! dimension loop grids
    endif
    ! Output meshgrids of the value
    call h5screate_simple_f(this%ndim, int(res_tmp2,kind=HSIZE_T), dspace, ierr)!, &
    call h5dcreate_f(file_id,"values", H5T_NATIVE_DOUBLE, dspace, dset_id, ierrhdf5)

    ! The trick to gaining arbitrary dimensional arrays into HDF5 is by using reshape immediately in the function
    ! as this saves us a maximum-dimension array allocation in the code.
    ! Beware: stack limits for reshape! set "ulimit -s unlimited"
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,RESHAPE(val_output(1:size(val_output)),res_tmp),int(res_tmp2,kind=HSIZE_T),ierr)

    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(dspace,ierr)
    if(output_grids) then
      call h5gclose_f(group_id_grid, ierr)
    endif


    call h5fclose_f(file_id,  ierrhdf5)
    call h5close_f(ierrhdf5)
    write(*,* ) "PARTICLES: Written Phase Space Projection to ", trim(filename)
  endif  ! mpi
  deallocate(val_output)
  deallocate(grid_mesh)
end subroutine output_phase_project

! Subroutine for calculating an array of indices  index_val for a single particle
! corresponding to their support in all dimensions. Also, the shape is
! calculated here in val_val
! > this: phase space projection for the arrays needed to project
! > x_in: location of particle on the dimensions to project
! > value_arr: array to be projected on. Should be the same size as this%values
! > value: projection value (i.e. if only particle weight, then just particle%weight)
subroutine project_single_particle_x(this,x_in,value_arr,proj_value)
  class(phase_space_projection), intent(in)        :: this
  real*8,                        intent(in)        :: x_in(this%ndim)
  real*8,                        intent(inout)     :: value_arr(this%val_size) ! Large!
  real*8,                        intent(in)        :: proj_value
  integer                                          :: i, j, main_ind, min_ind,min_ind_cont(this%ndim), mesh_tmp,ind_mesh_tmp, i_phase
  integer                                          :: index_tmp, indices_phase_grids_tmp(this%ndim)
  real*8                                           :: dx, x, minx, maxx, val_phase_grids_tmp, supp_weight(this%sumsupport), distance

  index_tmp=-1

  supp_weight=0.d0
  
  ! whole loop ~ <5% execution time
  do i=1,this%ndim
    dx = this%dx(i)
    minx=this%grids(this%previndex(i)+1)
    maxx=this%grids(this%previndex(i)+this%res(i))
    x = x_in(i)
    ! Not test if main x is in grid yet in order for particles outside the considered grid to
    ! contribute with their support (i.e. considering R=9.5-10.0, particle at 9.45 with 0.1
    ! bandwidth will still contribute slightly)

    ! Calculate minimum index (will be < 0 if particle outside of grid, but will be tested later)
    ! 1-indexing, should be 1 higher to indicate a true minimum index in the array. But in the support loop this is taken care of.
    min_ind = ceiling( (x-this%bandwidths(i)/2.d0-minx)/dx)
    min_ind_cont(i)=min_ind

    !Now, we add all the grid points in the particles shape into the grids
    do j=1,this%support(i)
      ! Only calculate weight if it makes sense
      if(min_ind+j>0 .and. min_ind+j < this%res(i)+1) then

        ! Distance w.r.t particle centre.
        distance=abs((x-this%grids(this%previndex(i)+min_ind+j))/this%bandwidths(i)*2.d0)

        if (distance > 1.d0) then
          supp_weight(j+this%prevsupp(i))=0.d0
        else

          ! This is the kernel, i.e. https://en.wikipedia.org/wiki/Kernel_(statistics)
          ! but we can only use finite support ones naturally (i.e. no Gaussian)
          ! Will be scaled later to ensure integration to 1.
          !Epachnikov
          !supp_weight(j+this%prevsupp(i))=3.d0/4.d0*(1-distance**2.d0)
          !Linear
          !supp_weight(j+this%prevsupp(i))=1-abs(distance)
          !Uniform
          !supp_weight(j+this%prevsupp(i))=0.5d0
          !Quartic
          !supp_weight(j+this%prevsupp(i))=15.d0/16.d0*(1-distance**2)**2
          !Triweight
          supp_weight(j+this%prevsupp(i))=35.d0/32.d0*(1-distance**2)**3

        endif ! distance < 1.0        
      endif ! min_ind+j > 0 && min_ind+j < this%res(i)+1
    enddo !support points
  enddo ! n_dim

  ! At this point we have the same situation as outputting grids (1D supports and values). For every point in totsupport we have to find
  ! the indices in the main value array.
  do i=1,this%totsupport

    val_phase_grids_tmp=1.d0
    ! Total weight of the support point i is calculated by the product of the weights in all dimensions. 
    do j=1,this%ndim
      val_phase_grids_tmp=val_phase_grids_tmp*supp_weight(this%totsupp_to_nD(j,i)+this%prevsupp(j))*1/this%bandwidths(j)*2.d0 !Bandwidth = 2*bandwidth_wiki
    enddo
  
    ! We need to check that this point is in fact on the grids
    if (all(this%totsupp_to_nD(:,i)+min_ind_cont >0) .and. all(this%totsupp_to_nD(:,i)+min_ind_cont < this%res+1)) then ! 13% execution time
      index_tmp=calc_index_phase_proj(this, this%totsupp_to_nD(:,i)+min_ind_cont)
        ! 43% execution time
        value_arr(index_tmp)=value_arr(index_tmp)+1.d0*val_phase_grids_tmp*proj_value
      endif
  enddo ! support points
end subroutine project_single_particle_x


pure function support_to_nD(this, index_in) result(index_out)
  class(phase_space_projection),  intent(in)       :: this
  integer,                        intent(in)       :: index_in
  integer                                          :: index_out(this%ndim)
  integer                                          :: index_tmp,it
  index_out = 0
  index_tmp=index_in-1
  do it=1,this%ndim
    index_out(it) = index_tmp/this%multsupp(it)+1
    index_tmp = modulo(index_tmp, this%multsupp(it))
  enddo
end function
! Same as previous subroutine, but with the option of providing the x_in on the grids.
! This can be quite useful if the x calculation involves interpolating some quantity that is
! already known at the moment this function is called
subroutine calc_index_shaped_part_x(this, particle_in,index_val,val_val,sim,x_in)
  type(particle_kinetic_leapfrog),intent(in)       :: particle_in
  class(phase_space_projection),  intent(in)       :: this
  real*8,                         intent(in)       :: x_in(this%ndim)
  type(particle_sim)                               :: sim
  real*8                                           :: value
  real*8                                           :: x,dx, minx, maxx, xmin,xmain,distance
  real*8                                           :: bandwidths(this%ndim)
  integer                                          :: index_arr(this%ndim),totsupp,index_supp(this%sumsupport),i,j, main_ind, min_ind,ind_mesh_tmp(this%ndim)
  integer                                          :: mesh_tmp,mesh_tmp2, indices_phase_grids_tmp(this%ndim)
  integer,                        intent(out)      :: index_val(this%totsupport)
  real*8,                         intent(out)      :: val_val(this%totsupport)
  real*8                                           :: supp_weight(this%sumsupport),val_phase_grids_tmp

  ! Initialize with 0 and -1 values and indices respectively.
  index_val=-1
  val_val= 0.d0
  supp_weight=0.d0
  index_supp=-1

  do i=1, this%ndim
    dx=this%grids(this%previndex(i)+2)-this%grids(this%previndex(i)+1)

    minx=this%grids(this%previndex(i)+1)

    maxx=this%grids(this%previndex(i)+this%res(i))
    if (particle_in%i_elm > 0) then
      x = x_in(i)

      ! Not test if main x is in grid yet in order for particles outside the considered grid to
      ! contribute with their support (i.e. considering R=9.5-10.0, particle at 9.45 with 0.1
      ! bandwidth will still contribute slightly)

      ! Calculate minimum index (will be < 0 if particle outside of grid, but will be tested later)
      ! Correct for even/uneven support
      if (mod(this%support(i),2).eq. 0 ) then
        main_ind=nint((x-minx)/dx)+1
        min_ind = main_ind - (this%support(i)-1)/2-1
      else
        main_ind=ceiling((x-minx)/dx)+1
        min_ind = main_ind - this%support(i)/2-1
      endif

      !Now, we add all the grid points in the particles shape into the grids
      do j=1,this%support(i)
        ! Only calculate weight if it makes sense wrt the grid
        if(min_ind+j>0 .and. min_ind+j < this%res(i)+1) then

          ! Distance w.r.t particle centre.
          distance=abs((x-this%grids(this%previndex(i)+min_ind+j))/this%bandwidths(i)*2.d0)

          if (distance > 1.d0) then
            supp_weight(j+this%prevsupp(i))=0.d0
          else

            ! This is the kernel, i.e. https://en.wikipedia.org/wiki/Kernel_(statistics)
            ! but we can only use finite support ones naturally (i.e. no Gaussian)

            !Epachnikov
            !supp_weight(j+this%prevsupp(i))=3.d0/4.d0*(1-distance**2)

            !Triweight (preffered by me due to C1 continuity)
            supp_weight(j+this%prevsupp(i))=35.d0/32.d0*(1-distance**2)**3

          endif ! distance < 1.0

          index_supp(j+this%prevsupp(i))=min_ind+j

        endif ! min_ind+j > 0 && min_ind+j < this%res(i)+1
      enddo !support points


    endif ! particle elm > 0

  enddo   ! ndim
  ! At this point we have the same situation as outputting grids (1D supports and values). For every point in totsupport we have to find
  ! the indices in the main value array.
  do i=1,this%totsupport

    !Calculate meshgrid indices
    mesh_tmp=i-1

    do j=1,this%ndim

      ind_mesh_tmp(j)=floor(real(mesh_tmp)/real(this%multsupp(j)))+1
      mesh_tmp= modulo(mesh_tmp, this%multsupp(j)) !Remainder after subtracting first stride.

    enddo
    val_phase_grids_tmp=1.d0

    do j=1,this%ndim

      indices_phase_grids_tmp(j)=index_supp(ind_mesh_tmp(j)+this%prevsupp(j))
      ! Scale the values by 1/bandwidth to be able to compare units w/ different bandwidths in each dimension
      val_phase_grids_tmp=val_phase_grids_tmp*supp_weight(ind_mesh_tmp(j)+this%prevsupp(j))*1/this%bandwidths(j)*2.d0 !Bandwidth = 2*bandwidth_wiki

    enddo

    ! Again, if any index of the grids of this specific support point is less than 0, the particle is not considered.
    if (minval(indices_phase_grids_tmp)> 0) then
      index_val(i)=calc_index_phase_proj(this,indices_phase_grids_tmp)

      val_val(i)=val_phase_grids_tmp

    endif
  enddo ! support points


end subroutine calc_index_shaped_part_x
! To go from 1D large meshgrid to ndim meshgrids function
function calc_index_phase_proj(this, index_arr) result(index_values)
  class(phase_space_projection), intent(in) :: this
  integer, intent(in)                       :: index_arr(this%ndim)
  integer                                   :: it
  integer                                   :: index_values,index_tmp

  index_tmp=1
  do it=1,this%ndim
    index_tmp=index_tmp+this%mult_index(it)*(index_arr(it)-1)
  enddo
  index_values=index_tmp
end function
! To go from index in 1D array to ndim indices.
function calc_reverse_index_phase_proj(this, index_values) result(index_arr)
  class(phase_space_projection), intent(in) :: this
  integer                                   :: index_arr(this%ndim)
  integer                                   :: it, tmp,tmp2
  integer                                   :: index_values,index_tmp

  index_arr=0
  index_tmp=index_values-1
  do it=1,this%ndim
    tmp = floor(real(index_tmp)/real(this%mult_index(it))) ! First stride. E.g. 10 -> 2,2,2,2 -> 0,
    tmp2 = modulo(index_tmp, this%mult_index(it)) !Remainder after subtracting first stride.
    index_arr(it)=tmp+1
    index_tmp = tmp2

  enddo

end function

end module mod_phase_space_project
