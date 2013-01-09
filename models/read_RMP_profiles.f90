subroutine read_RMP_profiles(bnd_node_list)

  !> Read RMP input profiles:
  !> - psi_RMP_cos from 'RMP_psi_cos_file'
  !> - psi_RMP_sin from 'RMP_psi_sin_file'
 
  use phys_module
  use data_structure, only: type_bnd_node_list
  use tr_module
  
  implicit none

  TYPE (type_bnd_node_list):: bnd_node_list
  integer :: j, err, k, ierr
  logical :: exist


  open (unit = 87, file=TRIM(RMP_psi_cos_file), &
       FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
  if ( err /= 0 ) then
     write(*,*) 'ERROR in boundary conditions: Cannot open file ''RMP_psi_cos_file'' '
     return
  end if
  
  open (unit = 88, file =TRIM(RMP_psi_sin_file), &
       FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
  if ( err /= 0 ) then
     write(*,*) 'ERROR in boundary conditions: Cannot open file ''RMP_psi_sin_file'' '
     return
  end if

  write (*,*) 'RMP_psi_cos_file = ', trim(RMP_psi_cos_file)
  write (*,*) 'RMP_psi_sin_file = ', trim(RMP_psi_sin_file)
  write (*,*) 'bnd_node_list%n_bnd_nodes = ', bnd_node_list%n_bnd_nodes

  k=0
  do
     k=k+1
     read(87,'()',iostat=ierr)
     if(ierr .ne. 0) exit
  end do
  k=k-1
  !write(*,*) 'k = ', k
  rewind(87)

  if ( k .ne. bnd_node_list%n_bnd_nodes ) then
     write(*,*) 'ERROR in read_RMP_profiles:  ''RMP_psi_cos_file'' has a wrong dimension'
     write(*,*) 'k_err = ', k
     return
  end if
  k=0
  do
     k=k+1
     read(88,'()',iostat=ierr)
     if(ierr .ne. 0) exit
  end do
  k=k-1
  !write(*,*) 'k1 = ', k
  rewind(88)
  if ( k .ne. bnd_node_list%n_bnd_nodes ) then
     write(*,*) 'ERROR in read_RMP_profiles:  ''RMP_psi_sin_file'' has a wrong dimension'
     write(*,*) 'k1_err = ', k
     return
  end if

  if (allocated(psi_RMP_cos))         call tr_deallocate(psi_RMP_cos,"psi_RMP_cos",CAT_UNKNOWN)
  if (allocated(dpsi_RMP_cos_dR))     call tr_deallocate(dpsi_RMP_cos_dR,"dpsi_RMP_cos_dR",CAT_UNKNOWN)
  if (allocated(dpsi_RMP_cos_dZ))     call tr_deallocate(dpsi_RMP_cos_dZ,"dpsi_RMP_cos_dZ",CAT_UNKNOWN)
  if (allocated(psi_RMP_sin))         call tr_deallocate(psi_RMP_sin,"psi_RMP_sin",CAT_UNKNOWN)
  if (allocated(dpsi_RMP_sin_dR))     call tr_deallocate(dpsi_RMP_sin_dR,"dpsi_RMP_sin_dR",CAT_UNKNOWN)
  if (allocated(dpsi_RMP_sin_dZ))     call tr_deallocate(dpsi_RMP_sin_dZ,"dpsi_RMP_sin_dZ",CAT_UNKNOWN)
  
  call tr_allocate(psi_RMP_cos,1, bnd_node_list%n_bnd_nodes,"psi_RMP_cos",CAT_UNKNOWN)
  call tr_allocate(dpsi_RMP_cos_dR,1,bnd_node_list%n_bnd_nodes,"dpsi_RMP_cos_dR",CAT_UNKNOWN)
  call tr_allocate(dpsi_RMP_cos_dZ,1,bnd_node_list%n_bnd_nodes,"dpsi_RMP_cos_dZ",CAT_UNKNOWN)
  call tr_allocate(psi_RMP_sin,1,bnd_node_list%n_bnd_nodes,"psi_RMP_sin",CAT_UNKNOWN)
  call tr_allocate(dpsi_RMP_sin_dR,1,bnd_node_list%n_bnd_nodes,"dpsi_RMP_sin_dR",CAT_UNKNOWN)
  call tr_allocate(dpsi_RMP_sin_dZ,1,bnd_node_list%n_bnd_nodes,"dpsi_RMP_sin_dZ",CAT_UNKNOWN)
                             
  do j=1, bnd_node_list%n_bnd_nodes  
     read (87, *) psi_RMP_cos(j),  dpsi_RMP_cos_dR(j), dpsi_RMP_cos_dZ(j)
     read (88, *) psi_RMP_sin(j),  dpsi_RMP_sin_dR(j), dpsi_RMP_sin_dZ(j)
  end do
  write (*,*) '**** before multiplying by time variation: ****'
  write (*,*) 'psi_RMP_cos(1) = ', psi_RMP_cos(1)
  write (*,*) 'dpsi_RMP_cos_dR(1) = ', dpsi_RMP_cos_dR(1)
  write (*,*) 'dpsi_RMP_cos_dZ(1) = ', dpsi_RMP_cos_dZ(1)
  write (*,*) 'psi_RMP_sin(1) = ', psi_RMP_sin(1)
  write (*,*) 'dpsi_RMP_sin_dR(1) = ', dpsi_RMP_sin_dR(1)
  write (*,*) 'dpsi_RMP_sin_dZ(1) = ', dpsi_RMP_sin_dZ(1)




  close (unit=87)
  close (unit=88)

  INQUIRE ( FILE ='RMP_start_time.dat', EXIST = exist )
  write (*,*) 'Does RMP_start_time file exist?', exist

  if (exist) then
     open (unit=97, file='RMP_start_time.dat', IOSTAT=err)
     if ( err /= 0 ) then
        write(*,*) 'ERROR in boundary_conditions.f90: Cannot open file RMP_start_time.dat'
        stop 'error file'
     end if

     read (97,*) RMP_start_time
     close(97)     
     !        RMP_time = RMP_time+tstep

        
     
  else 
     open (unit=97, file='RMP_start_time.dat', STATUS='NEW', IOSTAT=err)
     RMP_start_time = t_start
     write (97,*) RMP_start_time
     close(97)
     

  end if

     write(*,*) 'RMP_start_time', RMP_start_time
     write (*,*) 't_start', t_start
     write (*,*) 't_start - RMP_start_time - tset =', t_start - RMP_start_time - tset

     return
   end subroutine read_RMP_profiles
