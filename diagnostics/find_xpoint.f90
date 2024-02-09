!> Routine determines the position(s) of the xpoint(s).
subroutine find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,proper_xpoint,ifail)

use constants
use data_structure
use gauss
use basis_at_gaussian
use phys_module, only: tokamak_device, Z_xpoint_limit
use equil_info, only : ES
use mod_interp

implicit none

! --- Routine parameters
integer,                  intent(in)    :: my_id
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
real*8,                   intent(out)   :: psi_xpoint(2)
real*8,                   intent(out)   :: R_xpoint(2)
real*8,                   intent(out)   :: Z_xpoint(2)
integer,                  intent(out)   :: i_elm_xpoint(2)
real*8,                   intent(out)   :: s_xpoint(2)
real*8,                   intent(out)   :: t_xpoint(2)
integer,                  intent(in)    :: xcase
logical,                  intent(out)   :: proper_xpoint(2)        
integer,                  intent(out)   :: ifail

! --- Local variables
real*8  :: ps_s, ps_t, ps_x, ps_y, xjac
real*8  :: R, R_s, R_t, Z, Z_s, Z_t, P, P_s, P_t, P_st, P_ss, P_tt
real*8  :: x(2), s, t, xerr, ferr, s_xp_init(2), t_xp_init(2)
real*8  :: R_axis0, Z_axis0, R_xpoint0, Z_xpoint0, r_margin, s_axis, t_axis, psi_axis, fac_axis_xpoint       
integer :: ij_xpoint(2,2), i, iv, ms, mt, kf, kv, i_tries, n_tries, i_init       
integer :: i_elm_xp_init(2), min_indices_lw(3), min_indices_up(3)
integer :: i_elm_axis, ifail_axis   
logical :: found_upper, found_lower
real*8,  allocatable :: grad_psi(:,:,:)
logical, allocatable :: include_pt_lw(:,:,:), include_pt_up(:,:,:)

if (my_id .eq. 0) then
  write(*,*) '*********************************'
  write(*,*) '*     find_xpoint               *'
  write(*,*) '*********************************'
endif

ifail   = 1
n_tries = 500             
r_margin = 0.05          ! X-point found in sqrt((R-R_axis)^2 + (Z-Z_axis)^2) < r_margin will be dismissed and excluded from next loop. ! Grids in this circle must < n_tries                
fac_axis_xpoint = 9      ! If the min(|grad_psi|) point fulfilling the previous comment is closer to the axis than (fac_axis_xpoint * r_margin), assume that x-point is not found properly.
                         ! X-point where |grad_psi|=0 has no root, but where |grad_psi| ~< r_margin * div_psi(axis) can be accepted. 
proper_xpoint(:) = .false.   !Is the X-point(s) farther than (fac_axis_xpoint * r_margin) from the axis?

psi_xpoint = 0.
R_xpoint   = 0.;    Z_xpoint = 0.
s_xpoint   = 0.;    t_xpoint = 0.
i_elm_xpoint = 0

allocate(grad_psi      (element_list%n_elements,n_gauss,n_gauss))            ! --- vector storing |grad_psi| at gaussian poitns
allocate(include_pt_lw (element_list%n_elements,n_gauss,n_gauss))
allocate(include_pt_up (element_list%n_elements,n_gauss,n_gauss))
grad_psi    = 0.d0
include_pt_lw = .false.
include_pt_up = .false.

found_upper = .false. 
found_lower = .false.




if (.not. ES%initialized) then    
  call find_axis(99, node_list, element_list, psi_axis, R_axis0, Z_axis0, i_elm_axis, s_axis, &
  t_axis, ifail_axis)
else
  R_axis0 = ES%R_axis
  Z_axis0 = ES%Z_axis
endif


do i=1,element_list%n_elements    ! --- loop over elements
  
  do ms = 1, n_gauss           ! Gaussian points
    do mt = 1, n_gauss         ! Gaussian points

      ps_s = 0.d0
      ps_t = 0.d0
      R_s  = 0.d0 
      Z_s  = 0.d0
      R_t  = 0.d0 
      Z_t  = 0.d0
      R    = 0.d0
      Z    = 0.d0

      do kf = 1, n_degrees ! basis functions
        do kv = 1, 4       ! 4 vertices

          iv = element_list%element(i)%vertex(kv)

          ps_s = ps_s + node_list%node(iv)%values(1,kf,1) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          ps_t = ps_t + node_list%node(iv)%values(1,kf,1) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)

          R   = R   + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)
          Z   = Z   + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)

          R_s = R_s + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          Z_s = Z_s + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          R_t = R_t + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)
          Z_t = Z_t + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)

        enddo
      enddo

      xjac = R_s * Z_t - R_t * Z_s
      ps_x = (  ps_s * Z_t - ps_t * Z_s)/ xjac
      ps_y = (- ps_s * R_t + ps_t * R_s)/ xjac

      grad_psi(i,ms,mt) = sqrt(ps_x*ps_x + ps_y*ps_y)

      
      ! --- Look for the lower Xpoint
      if (xcase .ne. UPPER_XPOINT) then
        if (  (Z .lt. (Z_axis0 + 0.2)) .and. ((tokamak_device(1:4) .ne. 'MAST')  &
            .or. ((tokamak_device(1:4) .eq. 'MAST') .and. (R .gt. 0.45d0) .and. (R .lt. 1.d0))) ) then
          include_pt_lw(i,ms,mt) = .true.        
        endif
      endif
      
      ! --- And for the upper Xpoint
      if (xcase .ne. LOWER_XPOINT) then
        if (   (Z .gt. (Z_axis0 - 0.2)) .and. ((tokamak_device(1:4) .ne. 'MAST') &
            .or. ((tokamak_device(1:4) .eq. 'MAST') .and. (R .gt. 0.45d0) .and. (R .lt. 1.d0)) ) ) then
          include_pt_up(i,ms,mt) = .true.
        endif
      endif

    enddo
  enddo

enddo    ! --- end loop over elements


if(xcase .ne. UPPER_XPOINT) then

  i_init = 0

  do i_tries=1,  n_tries  ! --- start attempts to find the lower x-point
    
    ! --- min_indices = indices for gaussian point with min |grad_psi|,   (1) = element index, (2) = s-gaussian point index, (3) = t-gaussian point index
    min_indices_lw(:) = minloc(grad_psi, mask=include_pt_lw)
    if (.not. any(include_pt_lw)) min_indices_lw = 0

    if ((min_indices_lw(1) == 0) .and. (i_tries == 1)) then     ! --- if all elements are initially excluded, stop search and initialize values
      found_lower      = .false.
      s_xp_init(1)     = 0.d0
      t_xp_init(1)     = 0.d0
      i_elm_xp_init(1) = 1
      exit
    elseif  (min_indices_lw(1) == 0) then   ! --- if all elements have been excluded, exit search
      found_lower = .false.
      exit
    endif
    
    i_elm_xpoint(1) = min_indices_lw(1)    ! --- element with minimum |grad_psi|
    s = Xgauss(min_indices_lw(2)) 
    t = Xgauss(min_indices_lw(3))
    
    call mnewtax(node_list,element_list,i_elm_xpoint(1),s,t,xerr,ferr,ifail)
    if (ifail .ne. 0 ) then      ! --- if Newton's method failed, exclude element in next search
      include_pt_lw(i_elm_xpoint(1),:,:) = .false.
    endif
    call interp_RZ(node_list,element_list,i_elm_xpoint(1),s,t,R_xpoint0,R_s,R_t,Z_xpoint0,Z_s,Z_t)
    if (sqrt((R_axis0-R_xpoint0)**2 + (Z_xpoint0-Z_axis0)**2) .lt. r_margin)  then
      include_pt_lw(i_elm_xpoint(1),:,:) = .false.                                  ! If the point is within the r=r_margin circle around axis, exclude it
    elseif (include_pt_lw(i_elm_xpoint(1),1,1)) then 
      found_lower   = .true.
      s_xpoint(1)   = s
      t_xpoint(1)   = t
      exit
    elseif (i_init == 0) then           ! --- save first attempt outside axis region in case all the attempts fail. Aka. the min|grad\psi| point not excluded
      s_xp_init(1)     = s              ! ---possibly a x-point where \psi map is nosiy and |grad\psi|=0 fails to be solved
      t_xp_init(1)     = t
      i_elm_xp_init(1) = i_elm_xpoint(1)       
      i_init = 1
    endif
     
  enddo
  
endif

if(xcase .ne. LOWER_XPOINT) then

  i_init = 0

  do i_tries=1,  n_tries  ! --- start attempts to find the upper x-point

    ! --- min_indices = indices for gaussian point with min |grad_psi|,   (1) = element index, (2) = s-gaussian point index, (3) = t-gaussian point index
    min_indices_up(:) = minloc(grad_psi, mask=include_pt_up)
    if (.not. any(include_pt_up)) min_indices_up = 0
    
    if ((min_indices_up(1) == 0) .and. (i_tries == 1)) then     ! --- if all elements are initially excluded, stop search and initialize values
      found_upper      = .false.
      s_xp_init(2)     = 0.d0                             
      t_xp_init(2)     = 0.d0
      i_elm_xp_init(2) = 1                
      exit
    elseif  (min_indices_up(1) == 0) then   ! --- if all elements have been excluded, exit search
      found_upper     = .false.
      exit
    endif

    i_elm_xpoint(2) = min_indices_up(1)    ! --- element with minimum |grad_psi|
    s = Xgauss(min_indices_up(2)) 
    t = Xgauss(min_indices_up(3))
    
    call mnewtax(node_list,element_list,i_elm_xpoint(2),s,t,xerr,ferr,ifail)
    if (ifail .ne. 0 ) then       ! --- if Newton's method failed, exclude element in next search
      include_pt_up(i_elm_xpoint(2),:,:) = .false.
    endif
    call interp_RZ(node_list,element_list,i_elm_xpoint(2),s,t,R_xpoint0,R_s,R_t,Z_xpoint0,Z_s,Z_t)
    if (sqrt((R_axis0-R_xpoint0)**2 + (Z_xpoint0-Z_axis0)**2) .lt. r_margin) then
      include_pt_up(i_elm_xpoint(2),:,:) = .false.                                                   ! If the point is within the r=r_margin circle around axis, exclude it
    elseif ( include_pt_up(i_elm_xpoint(2),1,1) ) then
      found_upper   = .true.
      s_xpoint(2)   = s
      t_xpoint(2)   = t
      exit
    elseif (i_init == 0) then        ! --- save first attempt outside axis region in case all the attempts fail. Aka. the min|grad\psi| point not excluded   
      s_xp_init(2)     = s
      t_xp_init(2)     = t
      i_elm_xp_init(2) = i_elm_xpoint(2)   ! ---possibly a x-point where \psi map is nosiy and |grad\psi|=0 fails to be solved 
      i_init = 1 
    endif    

  enddo ! --- end attempts     

endif  
  


if(xcase .ne. UPPER_XPOINT) then

  proper_xpoint(1) = .true.
  if (.not. found_lower) then    ! --- if all the attempts failed, take the initial solution
    s_xpoint(1)     = s_xp_init(1)     
    t_xpoint(1)     = t_xp_init(1)     
    i_elm_xpoint(1) = i_elm_xp_init(1) 
  endif

  call interp(node_list,element_list,i_elm_xpoint(1),1,1,s_xpoint(1),t_xpoint(1),psi_xpoint(1),P_s,P_t,P_st,P_ss,P_tt)
  call interp_RZ(node_list,element_list,i_elm_xpoint(1),s_xpoint(1),t_xpoint(1),R_xpoint(1),R_s,R_t,Z_xpoint(1),Z_s,Z_t)

  xjac = R_s * Z_t - R_t * Z_s
  ps_x = (  P_s * Z_t - P_t * Z_s)/ xjac
  ps_y = (- P_s * R_t + P_t * R_s)/ xjac
  
  if ((.not. found_lower) .and. (sqrt((R_axis0-R_xpoint(1))**2 + (Z_xpoint(1)-Z_axis0)**2) .lt. fac_axis_xpoint*r_margin))  proper_xpoint(1) = .false.              
 ! If d_{xpoint to axis}<fac_axis_xpoint*r_margin, lower xpoint is not at a proper position

  if (my_id .eq. 0) then
    write(*,'(A,i6,4f14.8)') ' Lower X-point : ',i_elm_xpoint(1),R_xpoint(1),Z_xpoint(1),psi_xpoint(1),sqrt(ps_x**2+ps_y**2)
  endif
  
  if ((.not. found_lower )) write(*,*) 'WARNING: lower X-point not properly found after ', n_tries, ' attempts'
  if (.not.  proper_xpoint(1)) write(*,*) 'WARNING: lower X-point may not exist!'
  
endif


if(xcase .ne. LOWER_XPOINT) then 

  proper_xpoint(2) = .true.

  if (.not. found_upper) then    ! --- if all the attempts failed, take the initial solution
    s_xpoint(2)     = s_xp_init(2)     
    t_xpoint(2)     = t_xp_init(2)     
    i_elm_xpoint(2) = i_elm_xp_init(2) 
  endif
  
  call interp(node_list,element_list,i_elm_xpoint(2),1,1,s_xpoint(2),t_xpoint(2),psi_xpoint(2),P_s,P_t,P_st,P_ss,P_tt)
  call interp_RZ(node_list,element_list,i_elm_xpoint(2),s_xpoint(2),t_xpoint(2),R_xpoint(2),R_s,R_t,Z_xpoint(2),Z_s,Z_t)

  xjac = R_s * Z_t - R_t * Z_s
  ps_x = (  P_s * Z_t - P_t * Z_s)/ xjac
  ps_y = (- P_s * R_t + P_t * R_s)/ xjac
  
  if ((.not. found_upper) .and. (sqrt((R_axis0-R_xpoint(2))**2 + (Z_xpoint(2)-Z_axis0)**2) .lt. fac_axis_xpoint*r_margin))  proper_xpoint(2) = .false.             
   ! If d_{xpoint to axis}<fac_axis_xpoint*r_margin, upper xpoint is not at a proper position


  if (my_id .eq. 0) then
    write(*,'(A,i6,4f14.8)') ' Upper X-point : ',i_elm_xpoint(2),R_xpoint(2),Z_xpoint(2),psi_xpoint(2),sqrt(ps_x**2+ps_y**2)
  endif
    
  if ((.not. found_upper )) write(*,*) 'WARNING: upper X-point not properly found after ', n_tries, ' attempts'
  if (.not.  proper_xpoint(2)) write(*,*) 'WARNING: upper X-point may not exist!'

endif


deallocate(include_pt_lw,include_pt_up, grad_psi)

return
end subroutine find_xpoint
