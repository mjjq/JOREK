subroutine ideal_wall_starwall(my_id,node_list,bnd_elm_list,bnd_node_list)
!-------------------------------------------------------------------
! Reads the ideal wall vacuum response matrices written out by STARWALL
!-------------------------------------------------------------------

  use data_structure
  use vacuum_response_module
  use phys_module
  
  implicit none


  integer,                     intent(in) :: my_id            ! MPI thread number of current thread
  type(type_node_list),        intent(in) :: node_list        ! List of boundary nodes
  type(type_bnd_element_list), intent(in) :: bnd_elm_list     ! List of boundary elements
  type(type_bnd_node_list),    intent(in) :: bnd_node_list    ! List of boundary nodes


  integer :: dim
  integer :: itor, jtor
  integer :: inode, jnode
  integer :: ibas, jbas
  integer :: iindex, jindex
  integer :: iindex2, jindex2
  real*8  :: response
  real*8  :: TWOPI
  integer :: ierr, itor_start, itor_end
  

  TWOPI=8.d0*atan(1.d0)
  
  write(*,*) '************************************'
  write(*,*) '*        ideal_wall_starwall       *'
  write(*,*) '************************************'
  
  
  ! --- Read response matrix.
  write(*,*) 'Reading response matrix from file "vacuum_response_starwall".'
  open(11, FILE='vacuum_response_starwall', status='old', action='read', iostat=ierr)
  if ( ierr /= 0 ) then
    write(*,*) 'FATAL ERROR: Could not open file.'
    stop
  end if

  read(11,*) dim
  write(*,*) 'dim=',dim
  
  if (n_tor .eq. 1) then
    itor_start = 1
    itor_end   = 1
  elseif (n_tor .eq. 3) then
    itor_start = 2
    itor_end   = 3
  else
    stop 'multiple harmonics in starwall not yet implemented!'
  endif
  
  ! Outer loops (response index)
  do inode = 1, bnd_elm_list%n_bnd_elements              ! loop over nodes
     do itor = itor_start,itor_end                                         ! loop over toroidal modes
      do ibas = 1, 2                                      ! loop over basis functions
        
        iindex = 2*(inode-1) +   (ibas-1) + 1             ! first index in response matrix
        
        ! Inner loops (perturbation index)
        do jnode = 1, bnd_elm_list%n_bnd_elements        ! loop over nodes
          do jtor = itor_start,itor_end                                  ! loop over toroidal harmonics
            do jbas = 1, 2                                ! loop over basis functions
              
              jindex = 2*(jnode-1) +   (jbas-1) + 1       ! second index in response matrix
                
              read(11,*) iindex2, jindex2, response
              
              if ( itor == jtor ) then
                vacuum_response(iindex,jindex,itor) = - TWOPI * response
              endif
                  
            end do
          end do
        end do
 
      end do
    end do
  end do

!  vacuum_response(:,:,3) = vacuum_response(:,:,2)

  close(11)
  
  write(*,*) 'END ideal_wall_starwall'
  
end subroutine ideal_wall_starwall
