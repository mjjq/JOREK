subroutine respline_flux_surfaces(node_list,element_list,flux_list)

  use data_structure
  use grid_xpoint_data
  use py_plots_grids
  
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),         intent(in)     :: node_list
  type (type_element_list),      intent(in)     :: element_list
  type (type_surface_list),      intent(inout)  :: flux_list
  
  ! --- Internal parameters
  integer :: i_surf, i_piece
  real*8  :: s_int(2)
  real*8  :: s_cub1d(4), t_cub1d(4)
  integer :: failed
  real*8, parameter :: tolerance = 1.d-14
  
  logical, parameter :: debug = .false.
  real*8             :: progress
  
  
  ! --- Some printout
  if (debug) write(*,*)'Resplining surfaces (new)'
    
  ! --- Loop over each surface
  do i_surf = 1,flux_list%n_psi
    
    progress = 1.d2 * float(i_surf) / float(flux_list%n_psi)
    progress = max(0.d0,progress)
    progress = min(1.d2,progress)
    if (debug) write(*,"(' Processing  ... ',I0,'%',A,$)") int(progress),CHAR(13)
    
    ! --- Loop over each piece
    do i_piece = 1,flux_list%flux_surfaces(i_surf)%n_pieces
      
      s_cub1d(1:4) = flux_list%flux_surfaces(i_surf)%s(1:4,i_piece)
      t_cub1d(1:4) = flux_list%flux_surfaces(i_surf)%t(1:4,i_piece)
      s_int(1) = 1.d0/3.d0 ; s_int(2) = 2.d0/3.d0
      call respline_cub1d_piece(node_list,element_list,s_cub1d,t_cub1d,s_int,tolerance,failed)
      if (failed .eq. 2) then
        s_cub1d(1:4) = flux_list%flux_surfaces(i_surf)%s(1:4,i_piece)
        t_cub1d(1:4) = flux_list%flux_surfaces(i_surf)%t(1:4,i_piece)
        s_int(1) = 8.d0/10.d0 ; s_int(2) = 9.d0/10.d0
        call respline_cub1d_piece(node_list,element_list,s_cub1d,t_cub1d,s_int,tolerance,failed)
        !if (failed) write(*,*)'Warning: failed respline twice'
      endif
      
      if (failed .eq. 0) then
        flux_list%flux_surfaces(i_surf)%s(1:4,i_piece) = s_cub1d(1:4)
        flux_list%flux_surfaces(i_surf)%t(1:4,i_piece) = t_cub1d(1:4)
      else
        ! --- If the piece is completely symmetric, introduce artificial error (for the denominator of find_R_surface and find_Z_surface)
        if (failed .eq. 2) then
          if (abs(s_cub1d(2)-s_cub1d(4)) .le. tolerance) then
            s_cub1d(4) = s_cub1d(4) + 10.d0*tolerance
          else
            t_cub1d(4) = t_cub1d(4) + 10.d0*tolerance
          endif
          flux_list%flux_surfaces(i_surf)%s(1:4,i_piece) = s_cub1d(1:4)
          flux_list%flux_surfaces(i_surf)%t(1:4,i_piece) = t_cub1d(1:4)
        endif
      endif
!flux_list%flux_surfaces(i_surf)%s(2,i_piece) = 0.d0
!flux_list%flux_surfaces(i_surf)%s(4,i_piece) = 0.d0
!flux_list%flux_surfaces(i_surf)%t(2,i_piece) = 0.d0
!flux_list%flux_surfaces(i_surf)%t(4,i_piece) = 0.d0
      
    enddo
  
  enddo
  
  
  if (debug) write(*,*) 'Processing  ... 100'
  if (debug) write(*,*) 'finished resplining'
  
  
  return

end subroutine respline_flux_surfaces








subroutine respline_cub1d_piece(node_list,element_list,s_cub1d,t_cub1d,s_int,tolerance,failed)

  use data_structure
  use grid_xpoint_data
  use py_plots_grids
  
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),         intent(in)     :: node_list
  type (type_element_list),      intent(in)     :: element_list
  real*8,                        intent(inout)  :: s_cub1d(4),t_cub1d(4)
  real*8,                        intent(in)     :: s_int(2)
  real*8,                        intent(in)     :: tolerance
  integer,                       intent(out)    :: failed
  
  ! --- Internal parameters
  integer :: i_surf, i_piece, i
  real*8  :: s_beg, s_end
  real*8  :: t_beg, t_end
  real*8  :: s_find(2),  t_find(2)
  logical :: found
  real*8  :: st, st_find, new_st(2)
  real*8  :: ds_tmp
  real*8  :: dt_tmp
  real*8  ::  H1_s1 ,H1_s2,  H2_s1 ,H2_s2
  real*8  :: dH1_s1,dH1_s2, dH2_s1,dH2_s2
  real*8  :: cub_H1 , cub_H2 
  real*8  :: cub_dH1, cub_dH2
  real*8  :: Hjac
  real*8  ::  sX1,  sX2 
  real*8  :: dsX1, dsX2
  real*8  ::  tX1,  tX2 
  real*8  :: dtX1, dtX2
  
  failed = 0
  
  s_beg = s_cub1d(1)
  s_end = s_cub1d(3)
  t_beg = t_cub1d(1)
  t_end = t_cub1d(3)
  if (abs(s_end-s_beg) .gt. abs(t_end-t_beg)) then
    do i=1,2
      s_find(i) = s_beg + s_int(i) * (s_end-s_beg)
      call find_polar_intersection(s_cub1d, s_find(i), st_find, found)
      if (found) then
        call CUB1D(t_cub1d(1),t_cub1d(2),t_cub1d(3),t_cub1d(4), st_find, t_find(i), dt_tmp)
      else
        call find_polar_intersection_slow(s_cub1d, s_find(i), st_find, found)
        if (found) then
          call CUB1D(t_cub1d(1),t_cub1d(2),t_cub1d(3),t_cub1d(4), st_find, t_find(i), dt_tmp)
        else
          write(*,'(A)')'Warning: didnt find s'
          return
        endif
      endif
    enddo
  else
    do i=1,2
      t_find(i) = t_beg + s_int(i) * (t_end-t_beg)
      call find_polar_intersection(t_cub1d, t_find(i), st_find, found)
      if (found) then
        call CUB1D(s_cub1d(1),s_cub1d(2),s_cub1d(3),s_cub1d(4), st_find, s_find(i), ds_tmp)
      else
        call find_polar_intersection_slow(t_cub1d, t_find(i), st_find, found)
        if (found) then
          call CUB1D(s_cub1d(1),s_cub1d(2),s_cub1d(3),s_cub1d(4), st_find, s_find(i), ds_tmp)
        else
          write(*,*)'Warning: didnt find s',s_beg,s_end,s_find(i)
          return
        endif
      endif
    enddo
  endif
  
  do i=1,2
    new_st(i) = -1.d0 + 2.d0 * s_int(i)
  enddo
  H1_s1  = cub_H1 (new_st(1)) ; H1_s2  = cub_H1 (new_st(2))
  H2_s1  = cub_H2 (new_st(1)) ; H2_s2  = cub_H2 (new_st(2))
  dH1_s1 = cub_dH1(new_st(1)) ; dH1_s2 = cub_dH1(new_st(2))
  dH2_s1 = cub_dH2(new_st(1)) ; dH2_s2 = cub_dH2(new_st(2))
  
  Hjac = dH1_s1*dH2_s2 - dH1_s2*dH2_s1
  if (abs(Hjac) .lt. 1.d-12) then
    !write(*,*)'Problem: Hjac too small',Hjac
    failed = 1
    return
  endif
  
  ! --- This is the inverse solution of the system with 4 unknown of cubic splines
  ! --- Two of the solutions are at the edges, for X1 and X2, and the other two are in the middle (given by above) for dX1 and dX2
  sX1  = s_beg
  sX2  = s_end
  dsX1 = + dH2_s2/Hjac * (s_find(1) - sX1*H1_s1 - sX2*H2_s1) &
         - dH2_s1/Hjac * (s_find(2) - sX1*H1_s2 - sX2*H2_s2)
  dsX2 = - dH1_s2/Hjac * (s_find(1) - sX1*H1_s1 - sX2*H2_s1) &
         + dH1_s1/Hjac * (s_find(2) - sX1*H1_s2 - sX2*H2_s2)
  
  tX1  = t_beg
  tX2  = t_end
  dtX1 = + dH2_s2/Hjac * (t_find(1) - tX1*H1_s1 - tX2*H2_s1) &
         - dH2_s1/Hjac * (t_find(2) - tX1*H1_s2 - tX2*H2_s2)
  dtX2 = - dH1_s2/Hjac * (t_find(1) - tX1*H1_s1 - tX2*H2_s1) &
         + dH1_s1/Hjac * (t_find(2) - tX1*H1_s2 - tX2*H2_s2)
  
  if ( (abs(dsX1-dsX2) .le. tolerance) .or. (abs(dtX1-dtX2) .le. tolerance) ) then
    !write(*,*)'Problem: deriv too close',abs(dsX1-dsX2),abs(dtX1-dtX2)
    failed = 2
  endif
      
  ! --- Copy back into spline
  s_cub1d(2) = dsX1
  s_cub1d(4) = dsX2
  t_cub1d(2) = dtX1
  t_cub1d(4) = dtX2
  
  
  return

end subroutine respline_cub1d_piece







subroutine from_polar_to_cubic(s_polar,s_cubic)

  implicit none

  ! --- Routine variables
  real*8, intent(in)  :: s_polar(4)
  real*8, intent(out) :: s_cubic(4)
  
  s_cubic = (/ s_polar(1), 3.d0/2.d0 *(s_polar(2)-s_polar(1)), &
               s_polar(4), 3.d0/2.d0 *(s_polar(4)-s_polar(3))  /)
  return
end subroutine from_polar_to_cubic








subroutine find_polar_intersection(s_polar, s_find, st_find, found)

  implicit none
  
  ! --- Routine parameters
  real*8,  intent(in)  :: s_polar(4)
  real*8,  intent(in)  :: s_find
  real*8,  intent(out) :: st_find
  logical, intent(out) :: found
  
  ! --- Internal parameters
  real*8,  parameter :: accuracy = 1.d-6
  integer, parameter :: n_step = 50
  real*8  :: st_left, st_mid, st_right
  real*8  :: s_left,  s_mid,  s_right, ds_tmp
  integer :: step
  
  st_left  = -1.d0 - accuracy
  st_mid   = 0.d0
  st_right = +1.d0 + accuracy
  
  st_find  = 0.d0
  found    = .false.
  
  do step = 1,n_step
  
    call CUB1D(s_polar(1), s_polar(2), s_polar(3), s_polar(4), st_left,  s_left,  ds_tmp)
    call CUB1D(s_polar(1), s_polar(2), s_polar(3), s_polar(4), st_mid,   s_mid,   ds_tmp)
    call CUB1D(s_polar(1), s_polar(2), s_polar(3), s_polar(4), st_right, s_right, ds_tmp)
    
    if (abs(s_mid-s_find) .lt. accuracy) then
      st_find = st_mid
      found   = .true.
      exit
    endif
    
    if (    ( (s_mid .le. s_find) .and. (s_find .le. s_right) ) &
        .or.( (s_mid .ge. s_find) .and. (s_find .ge. s_right) ) ) then
      st_left = st_mid
    elseif (    ( (s_left .le. s_find) .and. (s_find .le. s_mid) ) &
            .or.( (s_left .ge. s_find) .and. (s_find .ge. s_mid) ) ) then
      st_right = st_mid
    endif
    st_mid = 0.5d0 * (st_left+st_right)
  
  enddo
  
  ! --- Make sure we are inside
  if (found) then
    if (st_find .lt. -1.d0) st_find = -1.d0
    if (st_find .gt. +1.d0) st_find = +1.d0
  endif
  
  return

end subroutine find_polar_intersection






subroutine find_polar_intersection_slow(s_polar, s_find, st_find, found)

  implicit none
  
  ! --- Routine parameters
  real*8,  intent(in)  :: s_polar(4)
  real*8,  intent(in)  :: s_find
  real*8,  intent(out) :: st_find
  logical, intent(out) :: found
  
  ! --- Internal parameters
  real*8,  parameter :: accuracy = 1.d-6
  integer, parameter :: n_step = 50
  real*8  :: st_left, st_mid, st_right, st, st_min
  real*8  :: s_left,  s_mid,  s_right, s_tmp, ds_tmp
  real*8  :: diff_min, diff
  integer :: step
  
  st_left  = -1.d0 - accuracy
  st_mid   = 0.d0
  st_right = +1.d0 + accuracy
  
  st_find  = 0.d0
  found    = .false.
  
  diff_min = 1.d10
  do step = 1,n_step
  
    st = -1.d0 + 2.d0 * real(step-1)/real(n_step-1)
    call CUB1D(s_polar(1), s_polar(2), s_polar(3), s_polar(4), st,  s_tmp,  ds_tmp)
    diff = abs(s_tmp-s_find)
    if (diff .lt. diff_min) then
      diff_min = diff
      st_min   = st
    endif
  enddo
  
  st_left  = st_min - 2.d0/real(n_step-1)
  st_mid   = st_min
  st_right = st_min + 2.d0/real(n_step-1)
  
  do step = 1,n_step
  
    call CUB1D(s_polar(1), s_polar(2), s_polar(3), s_polar(4), st_left,  s_left,  ds_tmp)
    call CUB1D(s_polar(1), s_polar(2), s_polar(3), s_polar(4), st_mid,   s_mid,   ds_tmp)
    call CUB1D(s_polar(1), s_polar(2), s_polar(3), s_polar(4), st_right, s_right, ds_tmp)
    
    if (abs(s_mid-s_find) .lt. accuracy) then
      st_find = st_mid
      found   = .true.
      exit
    endif
    
    if (    ( (s_mid .le. s_find) .and. (s_find .le. s_right) ) &
        .or.( (s_mid .ge. s_find) .and. (s_find .ge. s_right) ) ) then
      st_left = st_mid
    elseif (    ( (s_left .le. s_find) .and. (s_find .le. s_mid) ) &
            .or.( (s_left .ge. s_find) .and. (s_find .ge. s_mid) ) ) then
      st_right = st_mid
    endif
    st_mid = 0.5d0 * (st_left+st_right)
  
  enddo
  
  ! --- Make sure we are inside
  if (found) then
    if (st_find .lt. -1.d0) st_find = -1.d0
    if (st_find .gt. +1.d0) st_find = +1.d0
  endif
  
  return

end subroutine find_polar_intersection_slow

















real*8 function cub_H1(ss)
  implicit none
  real*8, intent(in) :: ss
  cub_H1 = 0.25d0 * (ss-1.d0)**2 * (ss+2.d0)
end function cub_H1

real*8 function cub_H2(ss)
  implicit none
  real*8, intent(in) :: ss
  cub_H2 = - 0.25d0 * (ss+1.d0)**2 * (ss-2.d0)
end function cub_H2

real*8 function cub_dH1(ss)
  implicit none
  real*8, intent(in) :: ss
  cub_dH1 = 0.25d0 * (ss-1.d0)**2 * (ss+1.d0)
end function cub_dH1

real*8 function cub_dH2(ss)
  implicit none
  real*8, intent(in) :: ss
  cub_dH2 = 0.25d0 * (ss+1.d0)**2 * (ss-1.d0)
end function cub_dH2

