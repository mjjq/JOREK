module mod_catalyst_adaptor
#ifdef USE_CATALYST

  use, intrinsic :: iso_c_binding
  use basis_at_gaussian
  use nodes_elements
  use mod_interp
  use mod_parameters, only: n_var, variable_names

  implicit none

  public :: catalyst_adaptor_initialise
  public :: catalyst_adaptor_execute
  public :: catalyst_adaptor_finalise
  public :: catalyst_adaptor

  private

  character(kind=c_char, len=4096), public :: catalyst_script !< Path to Catalyst pipeline script
  integer(c_int), public                   :: catalyst_nsub = 5 !< Number of subdivisions of each JOREK element
  integer                                  :: n_scalars !< The number of scalar variables passed to Catalyst
  integer                                  :: nnos !< Number of nodes in the Catalyst grid
  integer                                  :: nel !< Number of cells in the Catalyst grid
  integer                                  :: nnoel = 4 !< Number of points to define a cell
  real(c_float), pointer                   :: coords_RZst(:) !< pointer to coords array in catalyst_adaptor.cpp
                                                             !< only valid in catalyst_adaptor_execute
  integer                                  :: i_plane = 1 !< which plane to interpolate at

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

    subroutine catalyst_get_params(a_nsub, a_n_elements, a_n_scalars) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(out) :: a_nsub
      integer(c_int), intent(out) :: a_n_elements
      integer(c_int), intent(out) :: a_n_scalars

      a_nsub = catalyst_nsub
      a_n_elements = element_list%n_elements
      ! Just pass all variables as scalars for now (can change this later)
      n_scalars = n_var
      a_n_scalars = n_scalars
    end subroutine catalyst_get_params

    subroutine catalyst_get_scalar_name(a_scalar_name, a_iscalar) bind(C)
      use, intrinsic :: iso_c_binding
      character(kind=c_char), intent(out), dimension(36) :: a_scalar_name
      integer(c_int), intent(in) :: a_iscalar

      integer :: ichar
      character(kind=c_char,len=12) :: trimmed_name
      trimmed_name = trim(variable_names(a_iscalar)) // c_null_char
      do ichar=1,len(trimmed_name)
        a_scalar_name(ichar) = trimmed_name(ichar:ichar)
      enddo
    end subroutine catalyst_get_scalar_name

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
            inode = inode + 1
            a_coords_RZst(4*inode-3:4*inode) = real([R, Z, s, t], c_float)
          enddo
        enddo

        ! Calculate connectivity of each subelement
        do j=1,catalyst_nsub-1
          do k=1,catalyst_nsub-1
            ielm = ielm + 1
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

    subroutine catalyst_interp_scalars(a_scalars) bind(C)
      use, intrinsic :: iso_c_binding
      real(c_float), intent(inout), dimension(nnos*n_scalars) :: a_scalars

      integer :: i !< JOREK element index
      integer :: iflat !< index in the flattened array a_scalars
      integer :: iscalar, inode, i_tor
      real*8 :: s, t, P, P_s, P_t, P_ss, P_st, P_tt

      i = 0
      a_scalars = 0.0
      do inode=1,nnos
        ! Increment the JOREK element index by 1 every nsub^2 nodes
        if (modulo(inode - 1, catalyst_nsub * catalyst_nsub) .eq. 0) then
          i = i + 1
        endif
        s = coords_RZst(4*inode-1)
        t = coords_RZst(4*inode  )
        do iscalar=1,n_scalars
          iflat = n_scalars * (inode - 1) + iscalar
          do i_tor=1,n_tor
            call interp(node_list,element_list,i,iscalar,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            a_scalars(iflat) = a_scalars(iflat) + real(P * HZ(i_tor,i_plane), c_float)
          enddo ! i_tor
        enddo ! iscalar
      enddo ! inode

    end subroutine catalyst_interp_scalars

#endif 
end module mod_catalyst_adaptor