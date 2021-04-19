!< Calculates fields created by wall and PF coil currents at arbitrary points
!< Additionally wall forces are computed through the integral of the stress tensor
module mod_vacuum_fields

  use vacuum

  implicit none

  ! --- Scaled STARWALL wall triangle coordinates (n_tri, 3)
  real*8, allocatable :: xw_scaled(:,:), yw_scaled(:,:), zw_scaled(:,:)

   
  contains


  !< This routine re-calculates the STARWALL wall from wall Fourier harmonics and
  !! resizes it acording to a scale factor. This is useful to create a surface
  !! just outside the wall to perform wall forces integration. Note that using the
  !! same wall triangles would lead to singularities for the wall field calculation
  subroutine resize_starwall_wall(scale_fact)

    implicit none

    ! --- External parameters
    real*8, intent(in) :: scale_fact

    ! --- Local parameters
    real*8              :: pi2,fnu,alu,alv
    real*8              :: cm,cn,cov,siv,cou,siu,cop,sip,co,si
    real*8, allocatable :: r_w(:), x_w(:), y_w(:), z_w(:)
    integer             :: nwuv, i, j, kv, ku

    if (sr%file_version < 5) then
      write(*,*) 'STARWALL version not suitable for wall resizing'
      STOP
    endif
    if (sr%iwall /= 1) then
      write(*,*) 'Wall must be given in Fourier harmonics for resizing (iwall=1)'
      STOP
    endif

    nwuv  = sr%nwu*sr%nwv ! total number of wall nodes

    allocate (x_w(nwuv),y_w(nwuv),z_w(nwuv),r_w(nwuv))
    if (allocated(xw_scaled)) deallocate(xw_scaled)
    if (allocated(yw_scaled)) deallocate(yw_scaled)
    if (allocated(zw_scaled)) deallocate(zw_scaled)
    allocate (xw_scaled(sr%ntri_w,3),yw_scaled(sr%ntri_w,3),zw_scaled(sr%ntri_w,3))

    ! --- Resize the wall, wall construction copied from STARWALL
    pi2  =4.*asin(1.)
    fnu  = 1./float(sr%nwu) 
    alu  = pi2*fnu 
    alv  = pi2/float(sr%nwv)
    z_w = 0.
    r_w = 0.
    sr%rc_w(2:sr%mn_w) =sr%rc_w(2:sr%mn_w)*scale_fact
    sr%rs_w(2:sr%mn_w) =sr%rs_w(2:sr%mn_w)*scale_fact
    sr%zc_w(2:sr%mn_w) =sr%zc_w(2:sr%mn_w)*scale_fact
    sr%zs_w(2:sr%mn_w) =sr%zs_w(2:sr%mn_w)*scale_fact

    do  j =  1, sr%mn_w
      cm = sr%m_w(j)*pi2
      cn = sr%n_w_fourier(j)*pi2
      do  kv=1, sr%nwv
        cov = cos(alv*sr%n_w_fourier(j)*(kv-1))
        siv = sin(alv*sr%n_w_fourier(j)*(kv-1))
        do  ku=1,sr%nwu
          i = ku+sr%nwu*(kv-1)
          cou = cos(alu*sr%m_w(j)*(ku-1))
          siu = sin(alu*sr%m_w(j)*(ku-1))
    
          cop = cou*cov-siu*siv
          sip = siu*cov+cou*siv
    
          r_w(i) = r_w(i) + sr%rs_w(j)*sip + sr%rc_w(j)*cop
          z_w(i) = z_w(i) + sr%zs_w(j)*sip + sr%zc_w(j)*cop 
        end do
      end do
    end do
    
    do  kv = 1, sr%nwv
      co   = cos(alv*(kv-1)) 
      si   = sin(alv*(kv-1)) 
      do  ku = 1,sr%nwu
        i      = sr%nwu*(kv-1)+ku
        x_w(i)   =   co * r_w(i)
        y_w(i)   =   si * r_w(i)
      end do
    end do

    ! --- Reconstruct scaled wall triangle coordinates
    do i=1,2*nwuv
      xw_scaled(i,1) = x_w(sr%jpot_w(i,1))
      xw_scaled(i,2) = x_w(sr%jpot_w(i,2))
      xw_scaled(i,3) = x_w(sr%jpot_w(i,3))

      yw_scaled(i,1) = y_w(sr%jpot_w(i,1))
      yw_scaled(i,2) = y_w(sr%jpot_w(i,2))
      yw_scaled(i,3) = y_w(sr%jpot_w(i,3))

      zw_scaled(i,1) = z_w(sr%jpot_w(i,1))
      zw_scaled(i,2) = z_w(sr%jpot_w(i,2))
      zw_scaled(i,3) = z_w(sr%jpot_w(i,3))
    end do

    deallocate(x_w, y_w, z_w, r_w)

  end subroutine resize_starwall_wall


end module mod_vacuum_fields

