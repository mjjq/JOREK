subroutine ideal_wall_starwall(my_id,node_list,boundary_list,bnd_node_list)
!-------------------------------------------------------------------
! Reads the ideal wall vacuum response matrices written out by STARWALL
!-------------------------------------------------------------------

  use data_structure
  use vacuum_response_module
  use phys_module
  
  implicit none


  integer,                     intent(in) :: my_id            ! MPI thread number of current thread
  type(type_node_list),        intent(in) :: node_list        ! List of boundary nodes
  type(type_bnd_element_list), intent(in) :: boundary_list    ! List of boundary elements
  type(type_bnd_node_list),    intent(in) :: bnd_node_list    ! List of boundary nodes


  integer :: dim
  integer :: itor, jtor
  integer :: inode, jnode
  integer :: ibas, jbas
  integer :: iindex, jindex
  integer :: iindex2, jindex2
  real*8  :: response
  real*8  :: TWOPI
  integer :: ierr
  

  TWOPI=8.*atan(1.)
  
  write(*,*) '************************************'
  write(*,*) '*        ideal_wall_starwall       *'
  write(*,*) '************************************'
  
  
  ! --- Read response matrix.
  write(*,*) 'Reading response matrix from file "vacuum_response_starwall".'
  open(42, FILE='vacuum_response_starwall', status='old', action='read', iostat=ierr)
  if ( ierr /= 0 ) then
    write(*,*) 'FATAL ERROR: Could not open file.'
    stop
  end if

  read(42,*) dim
  write(*,*) 'dim=',dim
  

  ! Outer loops (response index)
  do inode = 1, boundary_list%n_bnd_elements              ! loop over nodes
    do itor = 2,3                                         ! loop over toroidal modes
      do ibas = 1, 2                                      ! loop over basis functions
        
        iindex = 2*(inode-1) +   (ibas-1) + 1             ! first index in response matrix
        
        ! Inner loops (perturbation index)
        do jnode = 1, boundary_list%n_bnd_elements        ! loop over nodes
          do jtor = 2, 3                                  ! loop over toroidal harmonics
            do jbas = 1, 2                                ! loop over basis functions
              
              jindex = 2*(jnode-1) +   (jbas-1) + 1       ! second index in response matrix
                
              read(42,*) iindex2, jindex2, response
              
              if ( itor == jtor ) then
                vacuum_response(iindex,jindex,itor) = - TWOPI * response
              endif
                  
            end do
          end do
        end do
 
      end do
    end do
  end do

  vacuum_response(:,:,3) = vacuum_response(:,:,2)

  close(42)
  
  

  write(*,*) 'END ideal_wall_starwall'
  
  

end subroutine ideal_wall_starwall
