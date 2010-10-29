subroutine resistive_wall_starwall(my_id, node_list, boundary_list,bnd_node_list)
!-------------------------------------------------------------------
! Reads the resistive wall vacuum response matrices written out by STARWALL
!-------------------------------------------------------------------

!### TODO: Let only one MPI thread read the files and broadcast them to the others.

  use data_structure
  use vacuum_response_module
  use phys_module
  
  
  implicit none


  integer,                     intent(in) :: my_id            ! MPI thread number of current thread
  type(type_node_list),        intent(in) :: node_list        ! List of boundary nodes
  type(type_bnd_element_list), intent(in) :: boundary_list    ! List of boundary elements
  type(type_bnd_node_list),    intent(in) :: bnd_node_list    ! List of boundary nodes


  integer :: dim, dim1, dim2
  integer :: itor, jtor
  integer :: inode, jnode
  integer :: ibas, jbas
  integer :: iindex, jindex
  integer :: iindex2, jindex2
  integer :: iwall, jwall
  real*8  :: response
  real*8  :: TWOPI
  integer :: ierr
  
  TWOPI=8.d0*atan(1.d0)
  
  write(*,*) '************************************'
  write(*,*) '*     resistive_wall_starwall      *'
  write(*,*) '************************************'
  
  
  
  ! --- Read 'EE' matrix.
  write(*,*) 'Reading "matrix_ee".'
  open(42, FILE='matrix_ee', status='old', action='read', iostat=ierr)
  if ( ierr /= 0 ) then
    write(*,*) 'FATAL ERROR: Could not open file.'
    stop
  end if
  
  read(42,*) dim
  write(*,*) 'dim=',dim
  
  if ( allocated(matrix_ee) ) deallocate(matrix_ee)
  allocate( matrix_ee(n_dof_bnd,n_dof_bnd,n_tor) )
  matrix_ee = 0.d0
  
  ! Outer loops (response index)
  do inode = 1, boundary_list%n_bnd_elements          ! loop over nodes
    do itor = 2, 3                                    ! loop over toroidal modes
      do ibas = 1, 2                                  ! loop over basis functions
        
        iindex = 2*(inode-1) +   (ibas-1) + 1         ! first index in response matrix
        
        ! Inner loops (perturbation index)
        do jnode = 1, boundary_list%n_bnd_elements    ! loop over nodes
          do jtor = 2, 3                              ! loop over toroidal harmonics
            do jbas = 1, 2                            ! loop over basis functions
              
              jindex = 2*(jnode-1) +   (jbas-1) + 1   ! second index in response matrix
                
              read(42,*) iindex2, jindex2, response   ! read one response matrix entry

              if ( itor == jtor ) then
                matrix_ee(iindex,jindex,itor) = - TWOPI * response
              endif
                  
            end do
          end do
        end do
 
      end do
    end do
  end do

  matrix_ee(:,:,3) = matrix_ee(:,:,2)
  
  close(42)
  
  
  
  ! --- Read 'EY' matrix.
  write(*,*) 'Reading "matrix_ey".'
  open(42, FILE='matrix_ee', status='old', action='read', iostat=ierr)
  if ( ierr /= 0 ) then
    write(*,*) 'FATAL ERROR: Could not open file.'
    stop
  end if

  read(42,*) dim1, dim2
  write(*,*) 'dim1=',dim1,', dim2=',dim2
  
  n_wall_curr = dim2 + 1
  
  if ( allocated(matrix_ey) ) deallocate(matrix_ey)
  allocate( matrix_ey(n_dof_bnd,n_wall_curr,n_tor) )
  
  matrix_ey = 0.d0
  
  ! Outer loops (response index)
  do inode = 1, boundary_list%n_bnd_elements          ! loop over nodes
    do itor = 2, 3                                    ! loop over toroidal modes
      do ibas = 1, 2                                  ! loop over basis functions
        
        iindex = 2*(inode-1) +   (ibas-1) + 1         ! first index in response matrix
        
        ! Inner loop (wall current index)
        do jwall = 1, n_wall_curr - 1 ! (no response for last wall current potential)
          
          read(42,*) iindex2, jindex2, response       ! read one response matrix entry

          matrix_ey(iindex,jwall,itor) = - TWOPI * response    !### really -TWOPI???
          
        end do
 
      end do
    end do
  end do

  matrix_ey(:,:,3) = matrix_ey(:,:,2)
  
  close(42)
  
  

  ! --- Read 'YE' matrix.
  write(*,*) 'Reading "matrix_ye".'
  open(42, FILE='matrix_ee', status='old', action='read', iostat=ierr)
  if ( ierr /= 0 ) then
    write(*,*) 'FATAL ERROR: Could not open file.'
    stop
  end if

  read(42,*) dim1, dim2
  write(*,*) 'dim1=',dim1,', dim2=',dim2
  
  if ( allocated(matrix_ye) ) deallocate(matrix_ye)
  allocate( matrix_ye(n_wall_curr,n_dof_bnd,n_tor) )
  
  matrix_ye = 0.d0
  
  ! Outer loop (wall current index)
  do iwall = 1, n_wall_curr - 1 ! (no response for last wall current potential)
    
    ! Inner loops (perturbation index)
    do jnode = 1, boundary_list%n_bnd_elements    ! loop over nodes
      do jtor = 2, 3                              ! loop over toroidal harmonics
        do jbas = 1, 2                            ! loop over basis functions
          
          jindex = 2*(jnode-1) +   (jbas-1) + 1   ! second index in response matrix
            
          read(42,*) iindex2, jindex2, response   ! read one response matrix entry
          
          matrix_ye(iwall,jindex,jtor) = - TWOPI * response    !### really -TWOPI???
                  
        end do
      end do
    end do
 
  end do

  matrix_ye(:,:,3) = matrix_ye(:,:,2)
  
  close(42)
  
  

  ! --- Read 'YY' matrix.
  write(*,*) 'Reading "matrix_yy".'
  open(42, FILE='matrix_ee', status='old', action='read', iostat=ierr)
  if ( ierr /= 0 ) then
    write(*,*) 'FATAL ERROR: Could not open file.'
    stop
  end if

  read(42,*) dim
  write(*,*) 'dim=',dim
  
  if ( allocated(diagonal_yy) ) deallocate(diagonal_yy)
  allocate( diagonal_yy(n_wall_curr) )
  
  diagonal_yy = 0.d0
  
  ! (wall current index)
  do iwall = 1, n_wall_curr - 1 ! (no response for last wall current potential)
    
    read(42,*) iindex2, response   ! read one response matrix entry
    
    diagonal_yy(iwall) = - TWOPI * response    !### really -TWOPI???
    
  end do
  
  close(42)
  
  
  
  ! --- Prepare the wall current array.
  if ( allocated(wall_curr) ) deallocate(wall_curr)
  allocate( wall_curr(n_wall_curr) )
  wall_curr = 0.d0
  
  
  
  
  write(*,*) 'END resistive_wall_starwall'
  
  

end subroutine resistive_wall_starwall

