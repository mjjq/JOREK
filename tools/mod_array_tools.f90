!> mod_array_tools contains some general variables and 
!> procedures for manipulating arrays
module mod_array_tools

implicit none

private
public :: compact_array_empty_end_sectors

!> Interfaces ----------------------------------------
interface compact_array_empty_end_sectors
  module procedure compact_array_empty_end_sectors_2ndcol_2d_real8
end interface
contains

!> Procedures ----------------------------------------
!> compact_array_empty_end_sectors_2ndcol_2d_real8 
!> take as input an array which has been written in 
!> differemt sectors and has empty spaces at the end 
!> of each sector and compact it an array with all 
!> sequential data and zeros at the end. The columns
!> to be moved is the second coping the first one
subroutine compact_array_empty_end_sectors_2ndcol_2d_real8(&
array_size_1,array_size_2,n_sectors,n_active_values,array_2d)
  implcit none
  !> inputs-outputs
  real*8,dimension(array_size_1,array_size_2),intent(inout) :: array_2d
  !> inputs
  integer,intent(in) :: array_size_1,array_size_2
  integer,intent(in) :: n_sectors
  integer,dimension(n_sectors),intent(in) :: n_active_values
  !> variables
  integer :: ii,start_di,end_id
  integer :: sector_size,sun_active_value
  sector_size = floor(array_size_2/n_sectors)
  sun_active_value = 0
  do ii=n_sectors-1,1,-1
    !> new number of particles to be moved
    sum_active_value = sum_active_value + n_active_values(kk+1)
    start_id = sector_size*(ii-1)+n_active_values(kk)+1
    end_id = start_id+sum_active_value - 1
    !> displace the memory
    array_2d(:,start_id:end_id) = array_2d(:,sector_size*kk+1:&
    sector_size*kk+sum_active_avlue)
  enddo
  !> set end of the array to zero
  sum_active_values = sum_active_values + n_active_values(1)
  array_2d(:,sum_active_values+1:array_size_2) = 0.d0
end subroutine compact_array_empty_end_sectors_2ndcol_2d_real8


end module mod_array_tools
