program jordel
!> Program to dump radial psi values extracted from restart file

use parameters, only: n_var, variable_names
use data_structure
use phys_module
use basis_at_gaussian

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
character *(256)      :: fname,jname1,jname2,jname3,jname4,jname5,jname6
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

write(*,*) 'jorek2del'

allocate(node_list)
allocate(element_list)

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters
nsub      = 5             ! Number of subdivisions of the cubic finite elements into linear pieces
i_tor     = 2  !-1        ! If i_tor > 0, only this mode will be included in the vtk file...
i_plane   = 1             ! ... otherwise, all modes will be summed up at the toroidal plane i_plane

rpos=.17
zpos=.44
tanth=tan(zpos/rpos)


! do-loop for multiple restart files
do ii=snum,fct
Write(fname,"('jorek',i5.5,'.rst')") ii
if(ii.eq.fct) Write(fname,"('jorek_restart.rst')")
call import_restart(node_list,element_list, fname, rst_format, ierr)

do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo

call initialise_basis                              ! define the basis functions at the Gaussian points

nnos = nsub*nsub*node_list%n_nodes
write(6,*) "nnos = ", nnos
allocate(outs(1:nnos,1:5))

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

         call interp(node_list,element_list,i,1,i_tor,s,t,Psi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt)
         !ps_x  = (   Z_t * PS_s - Z_s * PS_t ) / xjac
         !ps_y  = ( - R_t * PS_s + R_s * PS_t ) / xjac

         !ps_xx = (ps_ss * Z_t**2 - 2.d0*ps_st * Z_s*Z_t + ps_tt * Z_s**2  &
         !        + ps_s * (Z_st*Z_t - Z_tt*Z_s )                                   & 
         !        + ps_t * (Z_st*Z_s - Z_ss*Z_t ) )     / xjac**2                   & 
         !        - xjac_x * (ps_s* Z_t - ps_t * Z_s)  / xjac**2

         !ps_yy = (ps_ss * R_t**2 - 2.d0*ps_st * R_s*R_t + ps_tt * R_s**2  &
         !       + ps_s * (R_st*R_t - R_tt*R_s )                                    &
         !       + ps_t * (R_st*R_s - R_ss*R_t ) )         / xjac**2                &
         !       - xjac_y * (- ps_s * R_t + ps_t * R_s )  / xjac**2

        !outs(inode,1) = ps_x
        !outs(inode,2) = ps_y
        !outs(inode,3) = ps_xx
        !outs(inode,4) = ps_yy

        outs(inode,1) = R-10.  ! subtract R_0
        outs(inode,2) = Z
        outs(inode,3) = Psi
        outs(inode,4) = Ps_s
        outs(inode,5) = Ps_ss

      enddo  !k loop
   enddo  !j loop
enddo  ! i loop

!output to data files
!Write(jname1,"('deltas/psidr1',i5.5,'.jdat')") ii
!Write(jname2,"('deltas/psidr2',i5.5,'.jdat')") ii
Write(jname1,"('deltas/dr',i5.5,'.jdat')") ii
Write(jname2,"('deltas/dz',i5.5,'.jdat')") ii
Write(jname3,"('deltas/psi',i5.5,'.jdat')") ii
Write(jname4,"('deltas/psi1',i5.5,'.jdat')") ii
Write(jname5,"('deltas/psi2',i5.5,'.jdat')") ii
open(unit=121,file=jname1,action='write',position='append')
open(unit=122,file=jname2,action='write',position='append')
open(unit=123,file=jname3,action='write',position='append')
open(unit=124,file=jname4,action='write',position='append')
open(unit=125,file=jname5,action='write',position='append')

do jk=1,nnos
  zbot=outs(jk,1)*tanth-.013
  ztop=outs(jk,1)*tanth+.013
  if(outs(jk,1)>0.AND.outs(jk,2)<=ztop.AND.outs(jk,2)>=zbot) then
     write(121,*) outs(jk,1)  ! R
     write(122,*) outs(jk,2)  ! Z
     write(123,*) outs(jk,3)  ! psi
     write(124,*) outs(jk,4)  ! psi/dr
     write(125,*) outs(jk,5)  ! psi/dr^2
  endif
enddo

close(121)
close(122)
close(123)
close(124)
close(125)

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
