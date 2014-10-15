!> Export the current simulation state as a restart file that can be read back into JOREK or into
!! a diagnostic program by the routine import_restart.
subroutine export_restart(node_list,element_list,filename)

   use data_structure
   use phys_module
   use vacuum, only: export_restart_vacuum
#ifdef USE_HDF5
   use hdf5
   use HDF5_io_module
   use parameters, ONLY : n_tor, n_var, n_order
#endif

   implicit none

   ! --- Routine parameters
   type(type_node_list),    intent(in)    :: node_list
   type(type_element_list), intent(in)    :: element_list
   character*(*),           intent(in)    :: filename
   
   ! --- Local variables
   integer :: i
   
#ifdef USE_HDF5
   real*8 :: t_values(node_list%n_nodes,n_tor,n_order+1,n_var)
   real*8 :: t_psi_eq(node_list%n_nodes,n_tor,n_order+1,n_var)

   character(LEN=100) :: fileh5
   integer(HID_T)     :: file_id
   integer            :: ind, ierr
   integer            :: Nbr, Nbc
#endif
   
   open(21, file=filename, form='unformatted', status='replace', action='write')
   
   write(21) n_tor
   write(21) node_list%n_nodes,element_list%n_elements
   write(21) node_list%n_dof
   
   do i=1,node_list%n_nodes
      write(21) node_list%node(i)%x
      write(21) node_list%node(i)%values
#ifdef USE_HDF5
      t_values(i,:,:,:) = node_list%node(i)%values
#endif
      write(21) node_list%node(i)%deltas
#ifdef fullmhd
      write(21) node_list%node(i)%psi_eq               !< equilibrium flux at the nodes
#ifdef USE_HDF5
      t_psi_eq(i,:,:,:) = node_list%node(i)%psi_eq
#endif
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
   
#ifdef USE_HDF5
   write(21) h5_nbsave_all
   
   fileh5  = ' '
   ind     = index(filename,".rst")
   write(fileh5(1:),'(A,A)') trim(filename(1:ind)),'h5'
   write (6,*) ' file HDF5 = ',trim(fileh5)

   ! -> Create and open HDF5 file
   call HDF5_create(trim(fileh5),file_id,ierr)
   if (ierr.ne.0) then
      print*,'pglobal_id = ',pglobal_id, &
           ' ==> error for opening of HDF5 file',fileh5
   end if

   ! -> save : 'n_nodes'
   call HDF5_integer_saving(file_id,node_list%n_nodes,'n_nodes'//char(0))
   ! -> save : 'n_tor, n_order, n_var' 
   call HDF5_integer_saving(file_id,n_tor,'n_tor'//char(0))
   call HDF5_integer_saving(file_id,n_order,'n_order'//char(0))
   call HDF5_integer_saving(file_id,n_var,'n_var'//char(0))
   ! -> save : 'h5_nbsave_all'
   call HDF5_integer_saving(file_id,h5_nbsave_all,'h5_nbsave_all'//char(0))
   ! -> save : 't_now'
   call HDF5_real_saving(file_id,t_now,'t_now'//char(0))

   ! -> save : 'node_list%node(i)%values' 
   call HDF5_array4D_saving(file_id,t_values, &
        node_list%n_nodes,n_tor,n_order+1,n_var,'values')

   ! -> save : 'node_list%node(i)%psi_eq' 
   call HDF5_array4D_saving(file_id,t_psi_eq, &
        node_list%n_nodes,n_tor,n_order+1,n_var,'psi_eq')
   ! -> clode file
   call HDF5_close(file_id)
#endif
   
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
   
   if (use_pellet) then
      write(21) pellet_particles, pellet_R, pellet_Z
   endif
   
   call export_restart_vacuum(21, freeboundary, resistive_wall)
   
   close(21)
   
   return
end subroutine export_restart

