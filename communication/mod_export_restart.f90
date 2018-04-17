module mod_export_restart
implicit none
contains
!> Export the current simulation state as a restart file that can be read back into JOREK or into
!! a diagnostic program by the routine import_restart.
subroutine export_restart(node_list,element_list,filename)

  use mod_parameters
  use data_structure
  use phys_module
  use pellet_module

  implicit none

  ! --- Routine parameters
  type(type_node_list),    intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  character(len=*)       , intent(in) :: filename

  character*17 :: fileout

  if ( rst_hdf5 == 0 ) then
    ! --- Write restart binary file
    fileout = trim(filename)//".rst"
    write (6,*) " =============>, jorek2, filename = ", fileout
    call export_binary_restart(node_list, element_list, fileout)
  elseif ( rst_hdf5 == 1 ) then
    ! --- Write restart HDF5 file
    fileout = trim(filename)//".h5"
    write (6,*) " =============>, jorek2, filename = ", fileout
    call export_hdf5_restart(node_list, element_list, fileout)
  end if

end subroutine export_restart

!
! Export in a binary restart file
subroutine export_binary_restart(node_list,element_list,filename)

  use mod_parameters
  use data_structure
  use phys_module
  use pellet_module
  use vacuum, only: export_restart_vacuum

  implicit none

#include "version.h"

  ! --- Routine parameters
  type(type_node_list),    intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  character(len=*),        intent(in) :: filename

  ! --- Local variables
  integer :: i
  character*50 :: version_control

  ! -> Write binary restart file
  open(21, file=filename, form='unformatted', status='replace', action='write')

  write(21) n_tor
  write(21) node_list%n_nodes,element_list%n_elements
  write(21) node_list%n_dof

  do i=1,node_list%n_nodes
     write(21) node_list%node(i)%x
     write(21) node_list%node(i)%values
     write(21) node_list%node(i)%deltas
#ifdef fullmhd
     write(21) node_list%node(i)%psi_eq               !< equilibrium flux at the nodes
     write(21) node_list%node(i)%Fprof_eq             !< equilibrium profile R*B_phi at the nodes
#elif altcs
     write(21) node_list%node(i)%psi_eq               !< equilibrium flux at the nodes
#endif
     write(21) node_list%node(i)%index
     write(21) node_list%node(i)%boundary
     write(21) node_list%node(i)%parents
     write(21) node_list%node(i)%parent_elem
     write(21) node_list%node(i)%ref_lambda
     write(21) node_list%node(i)%ref_mu
     write(21) node_list%node(i)%constrained
  enddo

  write(21) element_list%element(1:element_list%n_elements)
  write(21) tstep,eta,visco,visco_par
  write(21) index_now
  write(21) t_now

  if (index_now .gt. 0) then
     write(21) xtime(1:index_now)
     write(21) energies(:,:,1:index_now)
#ifdef JECCD
     write(21) energies2(:,:,1:index_now)
     write(21) energies3(:,:,1:index_now)
#ifdef JEC2DIAG
     write(21) energies4(:,:,1:index_now)
#endif
#endif
  endif

  call export_restart_vacuum(21, freeboundary, resistive_wall)

  if (use_pellet) then
     if (index_now .gt. 0) then
        write(21) xtime_pellet_R(1:index_now)
        write(21) xtime_pellet_Z(1:index_now)
        write(21) xtime_pellet_psi(1:index_now)
        write(21) xtime_pellet_particles(1:index_now)
        write(21) xtime_phys_ablation(1:index_now)
     endif
     write(21) pellet_particles, pellet_R, pellet_Z
  endif

   
  ! save Revision control
  write(version_control,'(A)') trim(adjustl(RCS_VERSION))
  write(21) version_control

  ! save parameters
  write(21) jorek_model
  
  write(21) n_var
  write(21) n_order
  write(21) n_tor
  write(21) n_period
  write(21) n_plane
  write(21) n_vertex_max
  write(21) n_nodes_max
  write(21) n_elements_max
  write(21) n_boundary_max
  write(21) n_pieces_max
  write(21) n_degrees
  write(21) nref_max
  write(21) n_ref_list

  close(21)

  return
end subroutine export_binary_restart

 ! 
 ! Export in a HDF5 binary restart file
subroutine export_hdf5_restart(node_list,element_list,filename)
 
  use data_structure
  use phys_module
  use pellet_module
  use vacuum, only : export_HDF5_restart_vacuum
  
#ifdef USE_HDF5
  use hdf5
  use hdf5_io_module
  use tr_module
  use mod_parameters
#endif
 
  implicit none
 
#include "version.h"
 
  ! --- Routine parameters
  type(type_node_list),    intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  character*(*),           intent(in) :: filename

  ! --- Local variables
  integer :: i
  character(len=50)        :: version_control

#ifdef USE_HDF5
  integer(HID_T)     :: file_id
  integer            :: ind, ierr

  ! type_node, node_list%n_nodes
  real(RKIND), allocatable :: t_x(:,:,:)                   ! n_order+1, n_dim
  real(RKIND), allocatable :: t_values(:,:,:,:)            ! n_tor, n_order+1, n_var
  real(RKIND), allocatable :: t_deltas(:,:,:,:)            ! n_tor, n_order+1, n_var

  real(RKIND), allocatable :: t_psi_eq(:,:)                ! n_order+1
  real(RKIND), allocatable :: t_Fprof_eq(:,:)              ! n_order+1

  integer,     allocatable :: t_index(:,:)                 ! n_order+1
  integer,     allocatable :: t_boundary(:)                ! 
  integer,     allocatable :: t_parents(:,:)               ! 2
  integer,     allocatable :: t_parent_elem(:)             ! 
  real(RKIND), allocatable :: t_ref_lambda(:)
  real(RKIND), allocatable :: t_ref_mu(:)
  character,   allocatable :: t_constrained(:)     

  ! element, element_list%n_elements
  integer,     allocatable :: t_vertex(:,:)                ! n_vertex_max
  integer,     allocatable :: t_neighbours(:,:)            ! n_vertex_max
  real(RKIND), allocatable :: t_size(:,:,:)                ! n_vertex_max,n_order+1
  integer,     allocatable :: t_father(:)
  integer,     allocatable :: t_n_sons(:)
  integer,     allocatable :: t_n_gen(:)
  integer,     allocatable :: t_sons(:,:)                  ! 4
  integer,     allocatable :: t_contain_node(:,:)          ! 5
  integer,     allocatable :: t_nref(:)

  ! index_now+nstep
  real(RKIND), allocatable :: t_xtime(:)                   ! nstep
  real(RKIND), allocatable :: t_energies(:,:,:)            ! n_tor,2,index_start+nstep
#ifdef JECCD                                          
  real(RKIND), allocatable :: t_energies2(:,:,:)           ! n_tor,2,index_start+nstep
  real(RKIND), allocatable :: t_energies3(:,:,:)           ! n_tor,2,index_start+nstep
#ifdef JEC2DIAG                                       
  real(RKIND), allocatable :: t_energies4(:,:,:)           ! n_tor,2,index_start+nstep
#endif
#endif

  ! type_node, node_list%n_nodes
  call tr_allocate(t_x,1,node_list%n_nodes,1,n_order+1,1,n_dim, &
       "node_list%x",CAT_UNKNOWN)
  call tr_allocate(t_values,1,node_list%n_nodes,1,n_tor,1,n_order+1,1,n_var, &
       "node_list%values",CAT_UNKNOWN)
  call tr_allocate(t_deltas,1,node_list%n_nodes,1,n_tor,1,n_order+1,1,n_var, &
       "node_list%deltas",CAT_UNKNOWN)

#ifdef fullmhd
  call tr_allocate(t_psi_eq,1,node_list%n_nodes,1,n_order+1, &
       "node_list%psi_eq",CAT_UNKNOWN)
  call tr_allocate(t_Fprof_eq,1,node_list%n_nodes,1,n_order+1, &
       "node_list%Fprof_eq",CAT_UNKNOWN)
#elif altcs
  call tr_allocate(t_psi_eq,1,node_list%n_nodes,1,n_order+1, &
       "node_list%psi_eq",CAT_UNKNOWN)
#endif

  call tr_allocate(t_index,1,node_list%n_nodes,1,n_order+1,"index",CAT_UNKNOWN)
  call tr_allocate(t_boundary,1,node_list%n_nodes,"boundary",CAT_UNKNOWN)
  call tr_allocate(t_parents,1,node_list%n_nodes,1,2,"parent",CAT_UNKNOWN)
  call tr_allocate(t_parent_elem,1,node_list%n_nodes,"parent_elem",CAT_UNKNOWN)
  call tr_allocate(t_ref_lambda,1,node_list%n_nodes,"ref_lambade",CAT_UNKNOWN)
  call tr_allocate(t_ref_mu,1,node_list%n_nodes,"ref_mu",CAT_UNKNOWN)
  call tr_allocate(t_constrained,1,node_list%n_nodes,"constrained",CAT_UNKNOWN)

  ! element_list%n_elements
  call tr_allocate(t_vertex,1,element_list%n_elements,1,n_vertex_max,"vertex",CAT_UNKNOWN)
  call tr_allocate(t_neighbours,1,element_list%n_elements,1,n_vertex_max,"neighbours",CAT_UNKNOWN)
  call tr_allocate(t_size,1,element_list%n_elements,1,n_vertex_max,1,n_order+1,"size",CAT_UNKNOWN)
  call tr_allocate(t_father,1,element_list%n_elements,"father",CAT_UNKNOWN)
  call tr_allocate(t_n_sons,1,element_list%n_elements,"n_sons",CAT_UNKNOWN)
  call tr_allocate(t_n_gen,1,element_list%n_elements,"n_gen",CAT_UNKNOWN)
  call tr_allocate(t_sons,1,element_list%n_elements,1,4,"sons",CAT_UNKNOWN)
  call tr_allocate(t_contain_node,1,element_list%n_elements,1,5,"contain_node",CAT_UNKNOWN)
  call tr_allocate(t_nref,1,element_list%n_elements,"nref",CAT_UNKNOWN)

  ! index_now+nstep
  if (index_now .gt. 0) then
     if (allocated(t_xtime)) call tr_deallocate(t_xtime,"xtime",CAT_UNKNOWN)
     call tr_allocate(t_xtime,1,index_now,"xtime",CAT_UNKNOWN)
     t_xtime(:) = xtime(1:index_now)

     if (allocated(t_energies)) call tr_deallocate(t_energies,"energies",CAT_UNKNOWN)
     call tr_allocate(t_energies,1,n_tor,1,2,1,index_now,"energies",CAT_UNKNOWN)
     t_energies(:,:,:) = energies(:,:,1:index_now)
#ifdef JECCD
     if (allocated(t_energies2)) call tr_deallocate(t_energies2,"energies2",CAT_UNKNOWN)
     call tr_allocate(t_energies2,1,n_tor,1,2,1,index_now,"energies2",CAT_UNKNOWN)
     t_energies2(:,:,:) = t_energies2(:,:,1:index_now)

     if (allocated(t_energies3)) call tr_deallocate(t_energies3,"energies3",CAT_UNKNOWN)
     call tr_allocate(t_energies3,1,n_tor,1,2,1,index_now, "energies3",CAT_UNKNOWN)
     t_energies3(:,:,:) = t_energies3(:,:,1:index_now)

#ifdef JEC2DIAG
     if (allocated(t_energies4)) call tr_deallocate(t_energies4,"energies4",CAT_UNKNOWN)
     call tr_allocate(t_energies4,1,n_tor,1,2,1,index_now, "energies4",CAT_UNKNOWN)
     t_energies4(:,:,:) = t_energies4(:,:,1:index_now)
#endif
#endif

  end if

  !
  do i=1,node_list%n_nodes
     t_x(i,:,:)        = node_list%node(i)%x
     t_values(i,:,:,:) = node_list%node(i)%values
     t_deltas(i,:,:,:) = node_list%node(i)%deltas

#ifdef fullmhd
     t_psi_eq(i,:)     = node_list%node(i)%psi_eq
     t_Fprof_eq(i,:)   = node_list%node(i)%Fprof_eq
#elif altcs
     t_psi_eq(i,:)     = node_list%node(i)%psi_eq
#endif

     t_index(i,:)      = node_list%node(i)%index
     t_boundary(i)     = node_list%node(i)%boundary
     t_parents(i,:)    = node_list%node(i)%parents(1:2)
     t_parent_elem(i)  = node_list%node(i)%parent_elem
     t_ref_lambda(i)   = node_list%node(i)%ref_lambda
     t_ref_mu(i)       = node_list%node(i)%ref_mu
     if (node_list%node(i)%constrained) then
        t_constrained(i)  = 'T'
     else
        t_constrained(i)  = 'F'
     end if
  end do

  do i=1,element_list%n_elements
     t_vertex(i,:)       = element_list%element(i)%vertex
     t_neighbours(i,:)   = element_list%element(i)%neighbours
     t_size(i,:,:)       = element_list%element(i)%size
     t_father(i)         = element_list%element(i)%father
     t_n_sons(i)         = element_list%element(i)%n_sons
     t_n_gen(i)          = element_list%element(i)%n_gen
     t_sons(i,:)         = element_list%element(i)%sons
     t_contain_node(i,:) = element_list%element(i)%contain_node
     t_nref(i)           = element_list%element(i)%nref
  end do

  ! -> Create and open HDF5 file
  write (6,*) " HDF5 file ", filename
  call HDF5_create(trim(filename),file_id,ierr)
  if (ierr.ne.0) then
     print*,'pglobal_id = ',pglobal_id, &
          ' ==> error for opening of HDF5 file',filename
  end if
  
  ! Store the hdf5 restart file version
  write(*,*) 'Exporting HDF5 restart file with rst_hdf5_version=', rst_hdf5_version
  if ( rst_hdf5_version > rst_hdf5_version_supported ) then
    write(*,*) 'ERROR: Cannot write HDF5 restart file with rst_hdf5_version=', rst_hdf5_version
    write(*,*) '  This code version supports rst_hdf5_version_supported=', rst_hdf5_version_supported
    write(*,*) '  Please select an appropriate value in the namelist or use the default value.'
    stop
  end if
  call HDF5_integer_saving(file_id,rst_hdf5_version,'rst_hdf5_version'//char(0)) 

  ! -> Save version of revision control system
  write(version_control,'(A)') trim(adjustl(RCS_VERSION))
  version_control = trim(adjustl(version_control))
  call HDF5_char_saving(file_id,version_control,"RCS_version"//char(0))

  ! -> Save parameters
  call HDF5_integer_saving(file_id,jorek_model,'jorek_model'//char(0))
  call HDF5_integer_saving(file_id,n_var,'n_var'//char(0))
  call HDF5_integer_saving(file_id,n_dim,'n_dim'//char(0))
  call HDF5_integer_saving(file_id,n_order,'n_order'//char(0))
  call HDF5_integer_saving(file_id,n_tor,'n_tor'//char(0))
  call HDF5_integer_saving(file_id,n_period,'n_period'//char(0))
  call HDF5_integer_saving(file_id,n_plane,'n_plane'//char(0))
  call HDF5_integer_saving(file_id,n_vertex_max,'n_vertex_max'//char(0))
  call HDF5_integer_saving(file_id,n_nodes_max,'n_nodes_max'//char(0))
  call HDF5_integer_saving(file_id,n_elements_max,'n_elements_max'//char(0))
  call HDF5_integer_saving(file_id,n_boundary_max,'n_boundary_max'//char(0))
  call HDF5_integer_saving(file_id,n_pieces_max,'n_pieces_max'//char(0))
  call HDF5_integer_saving(file_id,n_degrees,'n_degrees'//char(0))
  call HDF5_integer_saving(file_id,nref_max,'nref_max'//char(0))
  call HDF5_integer_saving(file_id,n_ref_list,'n_ref_list'//char(0))

  ! -> 
  call HDF5_integer_saving(file_id,node_list%n_nodes,'n_nodes'//char(0))
  call HDF5_integer_saving(file_id,element_list%n_elements,'n_elements'//char(0))
  call HDF5_integer_saving(file_id,node_list%n_dof,'n_dof'//char(0))

  call HDF5_array3D_saving(file_id,t_x, &
       node_list%n_nodes,n_order+1,n_dim,'x'//char(0))
  call HDF5_array4D_saving(file_id,t_values, &
       node_list%n_nodes,n_tor,n_order+1,n_var,'values'//char(0))
  call HDF5_array4D_saving(file_id,t_deltas, &
       node_list%n_nodes,n_tor,n_order+1,n_var,'deltas'//char(0))

#ifdef fullmhd
  call HDF5_array2D_saving(file_id,t_psi_eq, &
       node_list%n_nodes,n_order+1,'psi_eq'//char(0))
  call HDF5_array2D_saving(file_id,t_Fprof_eq, &
       node_list%n_nodes,n_order+1,'Fprof_eq'//char(0))
#elif altcs
  call HDF5_array2D_saving(file_id,t_psi_eq, &
       node_list%n_nodes,n_order+1,'psi_eq'//char(0))
#endif

  call HDF5_array2D_saving_int(file_id,t_index, &
       node_list%n_nodes,n_order+1,'index'//char(0))
  call HDF5_array1D_saving_int(file_id,t_boundary, &
       node_list%n_nodes,'boundary'//char(0))
  call HDF5_array2D_saving_int(file_id,t_parents, &
       node_list%n_nodes,2,'parents'//char(0))
  call HDF5_array1D_saving_int(file_id,t_parent_elem, &
       node_list%n_nodes,'parent_elem'//char(0))
  call HDF5_array1D_saving(file_id,t_ref_lambda, &
       node_list%n_nodes,'ref_lambda'//char(0))
  call HDF5_array1D_saving(file_id,t_ref_mu, &
       node_list%n_nodes,'ref_mu'//char(0))
  call HDF5_array1D_saving_char(file_id,t_constrained, &
       node_list%n_nodes,'constrained'//char(0))

  call HDF5_array2D_saving_int(file_id,t_vertex, &
       element_list%n_elements,n_vertex_max,'vertex'//char(0))
  call HDF5_array2D_saving_int(file_id,t_neighbours, &
       element_list%n_elements,n_vertex_max,'neighbours'//char(0))
  call HDF5_array3D_saving(file_id,t_size, &
       element_list%n_elements,n_vertex_max,n_order+1,'size'//char(0))
  call HDF5_array1D_saving_int(file_id,t_father, &
       element_list%n_elements,'father'//char(0))
  call HDF5_array1D_saving_int(file_id,t_n_sons, &
       element_list%n_elements,'n_sons'//char(0))
  call HDF5_array1D_saving_int(file_id,t_n_gen, &
       element_list%n_elements,'n_gen'//char(0))
  call HDF5_array2D_saving_int(file_id,t_sons, &
       element_list%n_elements,4,'sons'//char(0))
  call HDF5_array2D_saving_int(file_id,t_contain_node, &
       element_list%n_elements,5,'contain_node'//char(0))
  call HDF5_array1D_saving_int(file_id,t_nref, &
       element_list%n_elements,'nref'//char(0))
  call HDF5_real_saving(file_id,tstep,'tstep'//char(0))
  call HDF5_real_saving(file_id,eta,'eta'//char(0))
  call HDF5_real_saving(file_id,visco,'visco'//char(0))
  call HDF5_real_saving(file_id,visco_par,'visco_par'//char(0))
  call HDF5_integer_saving(file_id,index_now,'index_now'//char(0)) 
  call HDF5_real_saving(file_id,t_now,'t_now'//char(0))
  call HDF5_real_saving(file_id,central_density,'central_density'//char(0))
  call HDF5_real_saving(file_id,central_mass,'central_mass'//char(0))
  call HDF5_real_saving(file_id,F0,'F0'//char(0))
  call HDF5_real_saving(file_id,sqrt_mu0_rho0,'sqrt_mu0_rho0'//char(0))
  call HDF5_real_saving(file_id,sqrt_mu0_rho0,'t_norm'//char(0))
  call HDF5_real_saving(file_id,sqrt_mu0_over_rho0,'sqrt_mu0_over_rho0'//char(0))

  if (index_now .gt. 0) then
     call HDF5_array1D_saving(file_id,t_xtime,index_now,'xtime'//char(0))
     call HDF5_array3D_saving(file_id,t_energies, &
          n_tor,2,index_now,'energies'//char(0))
     !           n_tor,2,index_now,'energies'//char(0))

     call HDF5_array1D_saving(file_id,R_axis_t(1:index_now),index_now,'R_axis_t'//char(0))
     call HDF5_array1D_saving(file_id,Z_axis_t(1:index_now),index_now,'Z_axis_t'//char(0))
     call HDF5_array1D_saving(file_id,psi_axis_t(1:index_now),index_now,'psi_axis_t'//char(0))
     call HDF5_array1D_saving(file_id,current_t(1:index_now),index_now,'current_t'//char(0))
     call HDF5_array1D_saving(file_id,beta_p_t(1:index_now),index_now,'beta_p_t'//char(0))
     call HDF5_array1D_saving(file_id,beta_t_t(1:index_now),index_now,'beta_t_t'//char(0))
     call HDF5_array1D_saving(file_id,beta_n_t(1:index_now),index_now,'beta_n_t'//char(0))
     call HDF5_array1D_saving(file_id,density_in_t(1:index_now),index_now,'density_in_t'//char(0))
     call HDF5_array1D_saving(file_id,density_out_t(1:index_now),index_now,'density_out_t'//char(0))
     call HDF5_array1D_saving(file_id,pressure_in_t(1:index_now),index_now,'pressure_in_t'//char(0))
     call HDF5_array1D_saving(file_id,pressure_out_t(1:index_now),index_now,'pressure_out_t'//char(0))
     call HDF5_array1D_saving(file_id,heat_src_in_t(1:index_now),index_now,'heat_src_in_t'//char(0))
     call HDF5_array1D_saving(file_id,heat_src_out_t(1:index_now),index_now,'heat_src_out_t'//char(0))
     call HDF5_array1D_saving(file_id,part_src_in_t(1:index_now),index_now,'part_src_in_t'//char(0))
     call HDF5_array1D_saving(file_id,part_src_out_t(1:index_now),index_now,'part_src_out_t'//char(0))

#ifdef JECCD                   
     call HDF5_array3D_saving(file_id,t_energies2, &
          n_tor,2,index_now,'energies2'//char(0))
     !           n_tor,2,index_now,'energies2'//char(0))
     call HDF5_array3D_saving(file_id,t_energies3, &
          n_tor,2,index_now,'energies3'//char(0))
     !           n_tor,2,index_now,'energies3'//char(0))
#ifdef JEC2DIAG
     call HDF5_array3D_saving(file_id,t_energies4, &
          n_tor,2,index_now,'energies4'//char(0))
     !           n_tor,2,index_now,'energies4'//char(0))
#endif
#endif
  end if

  if (use_pellet) then
     if (index_now .gt. 0) then
        call HDF5_array1D_saving(file_id,xtime_pellet_R, &
             index_now,'xtime_pellet_R'//char(0))
        call HDF5_array1D_saving(file_id,xtime_pellet_Z, &
             index_now,'xtime_pellet_Z'//char(0))
        call HDF5_array1D_saving(file_id,xtime_pellet_psi, &
             index_now,'xtime_pellet_psi'//char(0))
        call HDF5_array1D_saving(file_id,xtime_pellet_particles, &
             index_now,'xtime_pellet_particles'//char(0))
        call HDF5_array1D_saving(file_id,xtime_phys_ablation, &
             index_now,'xtime_phys_ablation'//char(0))
     end if
     call HDF5_real_saving(file_id,pellet_particles,"pellet_particles"//char(0))
     call HDF5_real_saving(file_id,pellet_R,"pellet_R"//char(0))
     call HDF5_real_saving(file_id,pellet_Z,"pellet_Z"//char(0))
  end if


  ! Export restart vacuum 
  call export_HDF5_restart_vacuum(file_id, freeboundary, resistive_wall)

  ! -> clode file
  call HDF5_close(file_id)

  ! -> Deallocate arrays
  call tr_deallocate(t_x,"x",CAT_UNKNOWN)
  call tr_deallocate(t_values,"values",CAT_UNKNOWN)
  call tr_deallocate(t_deltas,"deltas",CAT_UNKNOWN)
#ifdef fullmhd
  call tr_deallocate(t_psi_eq,"psi_eq",CAT_UNKNOWN)
  call tr_deallocate(t_Fprof_eq,"Fprof_eq",CAT_UNKNOWN)
#elif altcs
  call tr_deallocate(t_psi_eq,"psi_eq",CAT_UNKNOWN)
#endif
  call tr_deallocate(t_index,"index",CAT_UNKNOWN)
  call tr_deallocate(t_boundary,"boundary",CAT_UNKNOWN)
  call tr_deallocate(t_parents,"parents",CAT_UNKNOWN)
  call tr_deallocate(t_parent_elem,"parent_elem",CAT_UNKNOWN)
  call tr_deallocate(t_ref_lambda,"ref_lambda",CAT_UNKNOWN)
  call tr_deallocate(t_ref_mu,"ref_mu",CAT_UNKNOWN)
  call tr_deallocate(t_constrained,"constrained",CAT_UNKNOWN)

  call tr_deallocate(t_vertex,"vertex",CAT_UNKNOWN)
  call tr_deallocate(t_neighbours,"neighbours",CAT_UNKNOWN)
  call tr_deallocate(t_size,"size",CAT_UNKNOWN)
  call tr_deallocate(t_father,"father",CAT_UNKNOWN)
  call tr_deallocate(t_n_sons,"n_sons",CAT_UNKNOWN)
  call tr_deallocate(t_n_gen,"n_gen",CAT_UNKNOWN)
  call tr_deallocate(t_sons,"sons",CAT_UNKNOWN)
  call tr_deallocate(t_contain_node,"contain_node",CAT_UNKNOWN)
  call tr_deallocate(t_nref,"nref",CAT_UNKNOWN)

  if (index_now .gt. 0) then
     call tr_deallocate(t_xtime,"xtime",CAT_UNKNOWN)
     call tr_deallocate(t_energies,"energies",CAT_UNKNOWN)
#ifdef JECCD
     call tr_deallocate(t_energies2,"energies2",CAT_UNKNOWN)
     call tr_deallocate(t_energies3,"energies3",CAT_UNKNOWN)
#ifdef JEC2DIAG
     call tr_deallocate(t_energies4,"energies4",CAT_UNKNOWN)
#endif
#endif
  end if

#else
  write (6,*) " ERROR: trying to export with hdf5 but USE_HDF5 was not set at compile-time"
#endif

  return
end subroutine export_hdf5_restart
end module mod_export_restart
