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
  public new_phase_space_projection, output_phase_project,calc_index_val_phaseproj, calc_index_shaped_part, calc_index_shaped_part_x

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
    integer                                           :: totsupport,sumsupport
    logical                                           :: support_present

    contains
      procedure :: do => project_phase_space
  end type phase_space_projection

contains
! Constructor of the phase-space projection.
!  > ndim            : amount of dimensions
!  > res(ndim)       : resolution (amount of points) of the grids in each dimension
!  > start(ndim)     : starting points of each grid
!  > end(ndim)       : end points of each grid
!  > f_proj          : function that calculates the projected quantity (defined in mod_particle_projection)
!  > f_grids(ndim)   : ndim functions that calculates the positions in each dimension
!  > bandwidths(ndim): ndim bandwidths (shape extent) in each dimension
function new_phase_space_projection(sim,ndim,res,start,end,f_proj,f_grids,bandwidths) result(new)
  type(phase_space_projection) :: new
  type(particle_sim), intent(inout)                :: sim
  integer, intent(in)          :: ndim
  integer, intent(in)          :: res(ndim)
  real*8, intent(in), OPTIONAL :: bandwidths(ndim)
  integer                      :: it, gridpoints, valuepoints,j,mult_tmp
  type(proj_f)                 :: f_proj
  type(proj_f)                 :: f_grids(ndim)
  real*8,  intent(in)          ::  start(ndim), end(ndim)

  ! Allocating the members of the phase-space projection.
  allocate(new%res(ndim))
  allocate(new%previndex(ndim))
  allocate(new%mult_index(ndim))
  allocate(new%dx(ndim))
  allocate(new%f_grids(ndim))
  if(sim%my_id .eq. 0) then
    write(*,"(A,I1,A)") "PARTICLES: init phase space projection with ",ndim, " dimensions."
    if(ndim>7) then
      ! I cannot foresee a situation where this would come up, but for completeness sake.
      write(*,*) "PARTICLES: ", ndim, "dimensions is not supported by some compilers. "
      write(*,*) "This is not a problem for calculating values, but output into HDF5 is not"
      write(*,*) "possible (can be manually changed if your compiler supports it)"
      write(*,*) "However, it should not be a useful statistic anyway."
    endif
  endif



  new%support_present=.false.
  new%ndim = ndim
  new%res = res
  new%f_proj  = f_proj
  new%f_grids = f_grids

  ! Calculate the total values needed for each dimension
  gridpoints=sum(res)

  ! Calculate the total grid points (on a meshgrid)
  valuepoints=product(res)

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

    !Fill the grids uniformly
    do j=1,res(it)
      new%grids(new%previndex(it)+j)=start(it)+(end(it)-start(it))/(res(it)-1)*(j-1)

    enddo ! Grids loop
    new%dx(it)=new%grids(new%previndex(it)+2)-new%grids(new%previndex(it)+1)
  enddo   ! Dimension loop

  
  if (present(bandwidths)) then
    write(*,*) "Phase space projection: bandwidths given, particles will now have defined "
    write(*,*) "shapes in the projected dimensions"
    allocate(new%bandwidths(new%ndim))
    allocate(new%support(new%ndim))
    allocate(new%prevsupp(new%ndim))
    allocate(new%multsupp(new%ndim))
    new%support_present = .true.
    new%bandwidths=bandwidths

    new%prevsupp=0
    do it=1,ndim
      new%support(it)=ceiling(bandwidths(it)/new%dx(it))!


      write(*,*) "BANDWIDTH: ",bandwidths(it)

      if (it > 1 ) then
        new%prevsupp(it)=sum(new%support(1:(it-1)))
      endif

      write(*,*) "SUPPORT: ", it,new%support(it)
    enddo
    new%totsupport = product(new%support)
    new%sumsupport = sum(new%support)
    mult_tmp = new%totsupport
    do it=1,ndim
      mult_tmp=mult_tmp/new%support(it)
      new%multsupp(it)=mult_tmp
    enddo
  endif
  if((.not. new%support_present ).and. (sim%my_id .eq. 0)) then 
    write(*,*)"WARNING: Trying to use shaped particles without initializing WILL result in undefined behaviour"
    write(*,*) "Not checked for every projection, use carefully"
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

  allocate(index_arr(size(sim%groups(1)%particles),this%ndim))
  allocate(weight_arr(size(sim%groups(1)%particles),this%ndim))
  allocate(val_arr(size(sim%groups(1)%particles)))

  val_arr=0.d0

  call find_indices(this,sim,index_arr,weight_arr,val_arr)



  call project_linear(sim,this,index_arr,weight_arr,val_arr)

  deallocate(index_arr)
  deallocate(weight_arr)
  deallocate(val_arr)

end subroutine project_phase_space

! Output to h5 file. Structure: /values for the meshgrid-evaluated values. /grids/grid_i for 1D grids /grids/mgrid_i for ndim meshgrids
! (can be immediately plotted using these meshgrids)
subroutine output_phase_project(this)
  class(phase_space_projection), intent(inout)     :: this
  real*8, dimension(:), allocatable                :: val_output
  real*8, dimension(:,:), allocatable                :: grid_mesh
  integer                                          :: my_id, ierr,i,index_arr_tmp(this%ndim)
  integer(HID_T)                                   :: file_id, group_id_grid,dspace,dset_id
  integer                                          :: ierrhdf5
  integer                                          :: it,j
  CHARACTER(LEN=8)                                 :: tmp_name
  integer                                          :: res_tmp(7),order_tmp(7),res_tmp2(this%ndim) ! Change this if more dimensions needed
  !  but this should never come up as particles live in 7D at most.

  ! Output check for too large dimensional arrays for Fortran to handle.
  if (this%ndim > 7) then
    write(*,*) "Can't output for ndim > 7, aborting."
    call  MPI_ABORT(MPI_COMM_WORLD,-1, ierr)
  endif
  open(24,file="ree.txt")
  do it=1, this%val_size
    write(24,"(E11.4)") this%values(it)
  enddo
  close(24)

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
    call H5Fcreate_f("test.h5",H5F_ACC_TRUNC_F, file_id, ierrhdf5)
    call h5gcreate_f(file_id, "grids", group_id_grid, ierrhdf5)


    res_tmp2=res_tmp(1:this%ndim)
    ! Grid output
    do it=1,7
      write(*,*) "RESOLUTIONS: ",res_tmp(it)
    enddo

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

    ! Output meshgrids of the value
    call h5screate_simple_f(this%ndim, int(res_tmp2,kind=HSIZE_T), dspace, ierr)!, &
    call h5dcreate_f(file_id,"values", H5T_NATIVE_DOUBLE, dspace, dset_id, ierrhdf5)

    ! The trick to gaining arbitrary dimensional arrays into HDF5 is by using reshape immediately in the function
    ! as this saves us a maximum-dimension array allocation in the code.
    ! Compiler does not like using RESHAPE(val_output, res_tmp), make bounds explicit for it to work (although it seems rather pointless)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,RESHAPE(val_output(1:size(val_output)),res_tmp),int(res_tmp2,kind=HSIZE_T),ierr)

    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(dspace,ierr)

    call h5gclose_f(group_id_grid, ierr)


    call h5fclose_f(file_id,  ierrhdf5)
    call h5close_f(ierrhdf5)
  endif  ! mpi
  deallocate(val_output)
  deallocate(grid_mesh)
end subroutine output_phase_project


!> Subroutine for projecting a single a value for a single particle
!> in Nearest-Neighbour interpolation
subroutine calc_index_val_phaseproj(this, particle_in,index_val,sim)
  type(particle_kinetic_leapfrog), intent(in) :: particle_in
  class(phase_space_projection), intent(in)     :: this
  type(particle_sim)                               :: sim
  real*8                                           :: value
  real*8                                           :: x,dx, minx, maxx
  integer                                          :: index_arr(this%ndim),i
  integer,                        intent(out)      :: index_val

  !Calculate indices of the particle in the grids.
  index_arr=-1
  index_val=-1
  


  do i=1, this%ndim
    dx=this%grids(this%previndex(i)+2)-this%grids(this%previndex(i)+1)

    minx=this%grids(this%previndex(i)+1)

    maxx=this%grids(this%previndex(i)+this%res(i))
    if (particle_in%i_elm > 0) then
      x = this%f_grids(i)%f(sim,1,particle_in)
      if (x>maxx .or. x< minx) then
        ! by design indices are -1 if not in grids.
        index_arr(i)=-1
      else
        ! Nearest-Neighbour interpolation
        index_arr(i)=nint((x-minx)/dx)+1     
      endif
    else
      index_arr(i)=-1
    endif
  enddo

  ! If any of the indices are -1, particle is not considered.
  if ( minval(index_arr)> 0) then
    index_val=calc_index_phase_proj(this,index_arr)  ! ndim
  else
    index_val = -1
  endif

end subroutine calc_index_val_phaseproj

! Subroutine for calculating an array of indices  index_val for a single particle
! corresponding to their support in all dimensions. Also, the shape is 
! calculated here in val_val
subroutine calc_index_shaped_part(this, particle_in,index_val,val_val,sim)
  type(particle_kinetic_leapfrog), intent(in)      :: particle_in
  class(phase_space_projection), intent(in)        :: this
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
      x = this%f_grids(i)%f(sim,1,particle_in)
      
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
       val_phase_grids_tmp=val_phase_grids_tmp*supp_weight(ind_mesh_tmp(j)+this%prevsupp(j))*1/this%bandwidths(j) 
       
    enddo
    
    ! Again, if any index of the grids of this specific support point is less than 0, the particle is not considered.
    if (minval(indices_phase_grids_tmp)> 0) then
      index_val(i)=calc_index_phase_proj(this,indices_phase_grids_tmp)
      
      val_val(i)=val_phase_grids_tmp
      
    endif
  enddo ! support points


end subroutine calc_index_shaped_part

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
       val_phase_grids_tmp=val_phase_grids_tmp*supp_weight(ind_mesh_tmp(j)+this%prevsupp(j))*1/this%bandwidths(j) 
       
    enddo
    
    ! Again, if any index of the grids of this specific support point is less than 0, the particle is not considered.
    if (minval(indices_phase_grids_tmp)> 0) then
      index_val(i)=calc_index_phase_proj(this,indices_phase_grids_tmp)
      
      val_val(i)=val_phase_grids_tmp
      
    endif
  enddo ! support points


end subroutine calc_index_shaped_part_x


! Nearest neighbour projection of whole particle sim
subroutine project_linear(sim,this,index_arr,weight_arr,val_arr)
  class(phase_space_projection), intent(inout)     :: this

  type(particle_sim)                               :: sim
  integer,dimension(:,:),allocatable, intent(in)   :: index_arr
  real*8 ,dimension(:,:),allocatable, intent(in)   :: weight_arr
  real*8 ,dimension(:),allocatable, intent(in)     :: val_arr
  real*8 ,dimension(:), allocatable                :: val_tmp
  integer                                          :: i,j,index_val

  allocate(val_tmp,source=this%values)
  val_tmp=0.d0

  !$omp parallel do default(none)&
  !$omp shared(this,sim,index_arr,weight_arr,val_arr)&
  !$omp private(i,index_val)&
  !$omp reduction(+:val_tmp) schedule(dynamic,10)
  do i=1,size(index_arr,1)
    if (MINVAL(index_arr(i,1:this%ndim)) < 1) cycle !Using previous way of not including particle if not in the grid
    index_val = calc_index_phase_proj(this,index_arr(i,1:this%ndim))
    val_tmp(index_val) = val_tmp(index_val)+val_arr(i)
  enddo
  !$omp end parallel do
  this%values=this%values+val_tmp

  deallocate(val_tmp)

end subroutine project_linear

! Routine to find the indices on all the 1D grids for every particle for projection
subroutine find_indices(this, sim, index_arr,weight_arr,val_arr)!this,sim,index_arr,weight_arr)
  class(phase_space_projection),intent(in)         :: this
  type(particle_sim),intent(in)                    :: sim
  integer,  intent(inout)                          :: index_arr(:,:)
  real*8 ,  intent(inout)                          :: weight_arr(:,:)
  real*8,   intent(inout)                          :: val_arr(:)
  type(particle_kinetic_leapfrog)                  :: particle_tmp
  real*8                                           :: minx,maxx,dx, x
  integer                                          :: j,i

  do i=1, this%ndim
    dx=this%grids(this%previndex(i)+2)-this%grids(this%previndex(i)+1)

    minx=this%grids(this%previndex(i)+1)

    maxx=this%grids(this%previndex(i)+this%res(i))


    !$omp parallel do default(none)&
    !$omp shared(dx,minx,maxx,this,sim,weight_arr,index_arr,i,val_arr)&
    !$omp private(j,x)
    do j=1,size(sim%groups(1)%particles)
      !Grid spacing and minimal and maximal values to calculate element



      !Calculate the value of this grid dimension of the particle

      if (sim%groups(1)%particles(j)%i_elm < 1) then
        index_arr(j,i)=-1
        weight_arr(j,i)=-1.d0
        val_arr(j)=0.d0
        cycle
      endif

      x = this%f_grids(i)%f(sim,1,sim%groups(1)%particles(j))

      if (x>maxx .or. x< minx) then
        index_arr(j,i)=-1
        weight_arr(j,i)=-1.d0
        val_arr(j)=0.d0
      else
        ! Calculate projected quantity
        val_arr(j) = this%f_proj%f(sim,1,sim%groups(1)%particles(j))*sim%groups(1)%particles(j)%weight

        ! Grid element of this dimension and the corresponding weight to the left (lower) grid point
        ! Nearest neighbour interpolation
        index_arr(j,i)=nint((x-minx)/dx)+1

        weight_arr(j,i)=1-modulo(x-minx,dx)/dx
      endif ! Include only those particles fitting on the grid




    enddo  ! particles
    !$omp end parallel do
    write(*,*) "PARTICLES: calculated projections succesfully"
  enddo   ! ndim

end subroutine find_indices

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
