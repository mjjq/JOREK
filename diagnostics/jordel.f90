program jordel
!> Small program to extract a radial profile of values from jorek restart files
!> This program takes as input in the command line the number of a beginning
!> restart file and a final restart file.  An appropriate call is thus:
!> ./jordel 10 300 <inputfile
!> Data is produced in a file "deltas" which the user must create.
!> The radial cut can be manipulated by defining a single point
!> with the variables rpos and zpos. The linear cut will pass through that
!> point and also the origin.

use parameters, only: n_var, variable_names
use data_structure
use phys_module
use basis_at_gaussian
use import_restart

implicit none

type (type_node_list)   , pointer :: node_list
type (type_element_list), pointer :: element_list

integer               :: inode, nnos, nsub
real*8,allocatable    :: outs(:,:)
integer               :: i, j, k, m, etype, irst, int, i_var, i_tor, i_tor_old, i_plane, index, index_node, my_id
character             :: buffer*80, lf*1, str1*12, str2*12
real*8                :: s, t,rpos,zpos,tanth
real*8                :: P,P_s,P_t,P_st,P_ss,P_tt,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
real*8                :: Psi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt, ZJ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt
real*8                :: ps_x,ps_y,ps_xx,ps_yy
real*8                :: psi_bnd,psi_axis,R_axis,Z_axis,s_axis,t_axis
real*8                :: psi_xpoint(2),R_xpoint(2),Z_xpoint(2),s_xpoint(2),t_xpoint(2)
real*8                :: xjac_x, xjac_y,xjac
integer               :: i_elm_axis, i_elm_xpoint(2), k_tor, ifail, ierr
!====================== --- add the diagnostics Er, Vtheta and [not yet Vneo]
real*8                :: Er, psi_abs, Vtheta, Btheta
! alterations
character *(256)      :: fname,jname1,jname2,jname3,jname4,jname5,jname6,jname7
integer               :: jj,ii,kk,ij,jk,fct,snum
character(Len=100)    :: buf

real*8                :: ztop,zbot

!**********************************
!*** Get Command Line Arguments ***
!**********************************

        jj = getoargc(buf)
        If (jj == 0) Then
                Write(0,*) "Need number of first file!"
                Write(0,*) ""
                Stop
        Else
                Read(buf,*) snum
        End If

        jj = getoargc(buf)
        If (jj == 0) Then
                Write(0,*) "Need number of final file!"
                Write(0,*) ""
                Stop
        Else
                Read(buf,*) fct
        End If

write(*,*) 'jor2del'

allocate(node_list)
allocate(element_list)

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters
nsub      = 5             ! Number of subdivisions of the cubic finite elements into linear pieces
i_tor     = -1            ! If i_tor > 0, only this mode will be included in the vtk file...
i_plane   = 1             ! ... otherwise, all modes will be summed up at the toroidal plane i_plane

! here you choose a r and z position.  the variable is printed out along this
! line.
zpos=0.45
rpos=-0.24
tanth=(rpos/zpos)  ! tangent of an angle

! do-loop for multiple restart files
do ii=snum,fct
Write(fname,"('jorek',i5.5,'.rst')") ii
if(ii.eq.fct) Write(fname,"('jorek_restart.rst')")
call import_binary_restart(node_list,element_list, fname, rst_format, ierr)

do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo

call initialise_basis                              ! define the basis functions at the Gaussian points

nnos = nsub*nsub*node_list%n_nodes
write(6,*) "nnos = ", nnos
allocate(outs(1:nnos,1:3))

  call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
  psi_bnd = 0.d0
inode=0
do i=1,element_list%n_elements
   do j=1,nsub
      s = float(j-1)/float(nsub-1)
      do k=1,nsub
         t = float(k-1)/float(nsub-1)
         call interp_RZ(node_list,element_list,i,s,t,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
         xjac  = R_s * Z_t - R_t * Z_s
         if ( xjac == 0.d0 ) xjac = 1.d-8
         xjac_x  = (R_ss*Z_t**2 - Z_ss*R_t*Z_t - 2.d0*R_st*Z_s*Z_t   & 
              + Z_st*(R_s*Z_t + R_t*Z_s) + R_tt*Z_s**2 - Z_tt*R_s*Z_s) / xjac
         xjac_y  = (Z_tt*R_s**2 - R_tt*Z_s*R_s - 2.d0*Z_st*R_t*R_s   &
              + R_st*(Z_t*Z_s + Z_s*R_t) + Z_ss*R_t**2 - R_ss*Z_t*R_t) / xjac
         inode = inode+1

         ! here change the variable number (8) to the one you want to plot in the
         ! poloidal plane.
         call interp(node_list,element_list,i,8,i_tor,s,t,Psi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt)

        outs(inode,1) = R-10.  ! subtract R_0, in this case 10
        outs(inode,2) = Z
        outs(inode,3) = Psi

      enddo  !k loop
   enddo  !j loop
enddo  ! i loop

!output to data files
Write(jname1,"('deltas/dr',i5.5,'.jdat')") ii
Write(jname2,"('deltas/dz',i5.5,'.jdat')") ii
Write(jname3,"('deltas/psi',i5.5,'.jdat')") ii
open(unit=121,file=jname1,action='write',position='append')
open(unit=122,file=jname2,action='write',position='append')
open(unit=123,file=jname3,action='write',position='append')

do jk=1,nnos
  zbot=outs(jk,2)*tanth-.013
  ztop=outs(jk,2)*tanth+.013
  if(outs(jk,1)>0.d0) then
  if(outs(jk,2)<0.d0) then
  if(outs(jk,1)<=ztop.AND.outs(jk,1)>=zbot) then
     write(121,*) outs(jk,1)  ! R
     write(122,*) outs(jk,2)  ! Z
     write(123,*) outs(jk,3)  ! psi (or other variable chosen above)
  endif
  endif
  endif
enddo

do jk=1,nnos
  zbot=outs(jk,2)*tanth-.013
  ztop=outs(jk,2)*tanth+.013
  if(outs(jk,1)<0.d0) then
  if(outs(jk,2)>0.d0) then
  if(outs(jk,1)<=ztop.AND.outs(jk,1)>=zbot) then
     write(121,*) outs(jk,1)  ! R
     write(122,*) outs(jk,2)  ! Z
     write(123,*) outs(jk,3)  ! psi
  endif
  endif
  endif
enddo

close(121)
close(122)
close(123)

! write out stuff
deallocate(outs)

write(*,*) 'finished.'

enddo  ! loop over all files

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
 
end program jordel
