module mod_catalyst_adaptor
#ifdef USE_CATALYST

  use, intrinsic :: iso_c_binding
  use nodes_elements
  use mod_interp

  implicit none

  public :: catalyst_adaptor_initialise
  public :: catalyst_adaptor_execute
  public :: catalyst_adaptor_finalise
  public :: catalyst_adaptor

  private

  character(kind=c_char, len=4096), public :: catalyst_script !< Path to Catalyst pipeline script
  integer(c_int), public                   :: catalyst_nsub = 5 !< Number of subdivisions of each JOREK element
  integer                                  :: nnos !< Number of nodes in the Catalyst grid
  integer                                  :: nel !< Number of cells in the Catalyst grid
  integer                                  :: nnoel = 4 !< Number of points to define a cell
  real(c_float), pointer                   :: coords_RZst(:) !< pointer to coords array in catalyst_adaptor.cpp
                                                             !< only valid in catalyst_adaptor_execute

  interface

    subroutine catalyst_adaptor_initialise(a_catalyst_script) bind(C)
      use, intrinsic :: iso_c_binding
      character(kind=c_char), intent(in) :: a_catalyst_script 
    end subroutine catalyst_adaptor_initialise

    subroutine catalyst_adaptor_execute(a_step_index, a_time) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in) :: a_step_index
      real(c_double), intent(in) :: a_time
    end subroutine catalyst_adaptor_execute

    subroutine catalyst_adaptor_finalise() bind(C)
    end subroutine catalyst_adaptor_finalise

    ! empty function to get the dependency generator to work with catalyst_adaptor.cpp
    subroutine catalyst_adaptor() bind(C)
    end subroutine catalyst_adaptor

  end interface

  contains

    subroutine catalyst_get_grid_params(a_nsub, a_n_elements) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(out) :: a_nsub
      integer(c_int), intent(out) :: a_n_elements

      a_nsub = catalyst_nsub
      a_n_elements = element_list%n_elements
    end subroutine catalyst_get_grid_params

    subroutine catalyst_interp_grid(a_nnos, a_nel, a_coords_RZst, a_cell_points) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in) :: a_nnos
      integer(c_int), intent(in) :: a_nel
      real(c_float), intent(inout), dimension(4 * a_nnos), target :: a_coords_RZst
      integer(c_int), intent(inout), dimension(nnoel * a_nel) :: a_cell_points

      integer :: i, j, ielm, inode, k
      real*8 :: s, t, R, R_s, R_t, Z, Z_s, Z_t

      nnos = a_nnos
      nel = a_nel
      inode = 0
      ielm = 0

      ! Store a pointer to the coords in this module
      coords_RZst => a_coords_RZst

      ! Create points for each element
      do i=1,element_list%n_elements
        do j=1,catalyst_nsub
          s = float(j-1)/float(catalyst_nsub-1)
          ! Create nsub^2 points per element at regularly spaced intervals
          do k=1,catalyst_nsub
            t = float(k-1)/float(catalyst_nsub-1)
            call interp_RZ(node_list,element_list,i,s,t,R,R_s,R_t,Z,Z_s,Z_t)
            inode = inode+1
            a_coords_RZst(4*inode-3:4*inode) = real([R, Z, s, t], c_float)
          enddo
        enddo

        ! Calculate connectivity of each subelement
        do j=1,catalyst_nsub-1
          do k=1,catalyst_nsub-1
            ielm = ielm+1
            ! Hopefully Catalyst expects the same convention as VTK
            ! Conduit doesn't specify a convention
            a_cell_points(4*ielm-3) = inode - catalyst_nsub*catalyst_nsub + catalyst_nsub*(j-1) + k-1 ! 0 based indices as for VTK
            a_cell_points(4*ielm-2) = inode - catalyst_nsub*catalyst_nsub + catalyst_nsub*(j  ) + k-1
            a_cell_points(4*ielm-1) = inode - catalyst_nsub*catalyst_nsub + catalyst_nsub*(j  ) + k
            a_cell_points(4*ielm  ) = inode - catalyst_nsub*catalyst_nsub + catalyst_nsub*(j-1) + k
          enddo
        enddo
      enddo
    end subroutine catalyst_interp_grid


#endif 
end module mod_catalyst_adaptor