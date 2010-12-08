subroutine import_external_fields
!--------------------------------------------------------
! reads the external fields and the poloidal field coils
! from the STARWALL output
!--------------------------------------------------------

use vacuum_response_module
implicit none

integer  :: n_dof_bnd_tmp, i, j, n_coils_tmp, itmp, jtmp

write(*,*) '**************************'
write(*,*) '* import external fields *'
write(*,*) '**************************'
write(*,*) ' reading from file : response__bext_par'

open(11,file='response_bext_par')                         ! external fields form starwall

read(11,*)
read(11,*) n_dof_bnd_tmp, n_coils
              
n_dof_bnd_tmp = 2*n_dof_bnd_tmp
       
if (n_dof_bnd .ne. n_dof_bnd_tmp) write(*,'(A,2i4)') ' ERROR : wrong number of boundary points ;',n_dof_bnd,n_dof_bnd_tmp
       
allocate(external_field(n_dof_bnd_tmp,n_coils))           ! allocate external_field matrix
allocate(I_coils(n_coils))                                 ! allocate the coil currents
       
do i=1,n_dof_bnd_tmp
  do j=1,n_coils
    read(11,*) itmp,jtmp,external_field(i,j)
  enddo
enddo

close(11)

open(11,file='coils.txt')                                ! read poloidal field coil geometry (only for plotting)

read(11,*);read(11,*);read(11,*);read(11,*)

read(11,*) n_coils_tmp

if (n_coils_tmp .ne. n_coils) write(*,*) 'INCONSISTENT NUMBER OF COILS (coil.txt)'

read(11,*);read(11,*)

allocate(R_coils(n_coils_tmp),Z_coils(n_coils_tmp),dR_coils(n_coils_tmp),dZ_coils(n_coils_tmp))

do i=1,n_coils_tmp
  read(11,*) R_coils(i), Z_coils(i), dR_coils(i), dZ_coils(i)
enddo

close(11)

return
end
