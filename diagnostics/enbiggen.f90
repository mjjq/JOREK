!> Program to convert a JOREK2 restart file into binary VTK format
program enbiggen

use parameters, only: n_var, variable_names
use data_structure
use phys_module
!use vacuum, only: export_restart_vacuum

implicit none

type (type_node_list)   , pointer :: node_list
type (type_element_list), pointer :: element_list

character *(256)      :: fname,fname2,fname3,fname4,fname5,fname6
integer               :: jj,ii,fct,snum, my_id
integer               :: ierr
character(Len=100)    :: buf
integer               :: jlen
!**********************************
!*** Get Command Line Arguments ***
!**********************************

        jj = getoargc(buf)
        If (jj == 0) Then
                Write(0,*) "Need number of restart file!"
                Write(0,*) ""
                Stop
        Else
                Read(buf,*) snum
        End If

write(6,*) 'enbiggen:  expand a jorek restart file from model303 to model305'
call flush_it(6)
allocate(node_list)
allocate(element_list)

! new do-loop for more vtk files
Write(fname,"('jorek',i5.5)") snum
Write(fname2,"('jorek',i5.5,'.rst_mod305')") snum
call import_restart(node_list,element_list, fname, rst_format, ierr)

write(6,*) 'enbiggen:  read frame successfully'
write(6,*) 'enbiggen:  original n_var, n_tor = ', n_var, n_tor, index_start

write(6,*) 'enbiggen:  initial frame exported successfully'

open(25, file=fname2, form='unformatted', status='replace', action='write')

write(25) n_tor
write(25) node_list%n_nodes,element_list%n_elements
write(25) node_list%n_dof

do ii=1,node_list%n_nodes
  write(25) node_list%node(ii)%x
  write(25) node_list%node(ii)%values,node_list%node(ii)%values(:,:,1)*0.
! ad values for new equation
!!
  write(25) node_list%node(ii)%deltas, node_list%node(ii)%deltas(:,:,1)*0.
! ad values for new equation
!!
  write(25) node_list%node(ii)%index
  write(25) node_list%node(ii)%boundary
  write(25) node_list%node(ii)%parents
  write(25) node_list%node(ii)%parent_elem
  write(25) node_list%node(ii)%ref_lambda
  write(25) node_list%node(ii)%ref_mu
  write(25) node_list%node(ii)%constrained
enddo
!! tests 

  !write(6,*) "values of values:"
  !write(6,*) node_list%node(1)%values(:,:,1)

  !write(6,*) "values of values:"
  !write(6,*) node_list%node(1)%values(:,:,2)

  !write(6,*) "values of values:"
  !write(6,*) node_list%node(1)%values
  !write(6,*) node_list%node(1)%values(:,:,2)*0.

  !write(6,*) "times!!! ",index_start,t_start
  !write(6,*) "xtimes!!! ",xtime(index_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(25) element_list%element(1:element_list%n_elements)
write(25) tstep,eta,visco,visco_par
write(25) index_start
write(25) t_start

write(25) 1

if (index_start .gt. 0) then
  write(25) xtime(1:index_start)
  write(25) energies(:,:,1:index_start)
  write(25) energies(:,:,1:index_start)*0.
  write(25) energies(:,:,1:index_start)*0.
endif

!call export_restart_vacuum(25, freeboundary, resistive_wall, index_start)

    if(freeboundary .and. (index_start > 0)) then
      write(25) resistive_wall
    endif

close(25)

write(*,*) 'enbiggen:  finished.'


Contains
Integer Function getoargc (buf)
        Implicit None
        Character(Len=*)::      buf
        Integer::               iargc
        Integer, Save:: ncall = 0
        ncall = ncall + 1
        If (ncall > iargc()) Then
                getoargc = 0
                Return
        Endif
        Call getarg (ncall, buf)
        getoargc = 1
End Function getoargc
 
end program enbiggen
