module divertor_desc
  implicit none
  ! identifier of one line on the divertor shape
  type divseg_num_t
     integer :: num_part, num_seg
  end type divseg_num_t
  

  type divertor_pos_t
     ! position on a divertor is a line identifier and a position on that line
     type(divseg_num_t) :: divseg_pos
     real*8 :: relpos
  end type divertor_pos_t

  type divpoint_t
     ! a point of the divertor shape
     real*8 absdist
     real*8 :: x(2)
  end type divpoint_t

  type divertor_part
     ! a divertor part is a contiguous sequence of segments defined by
     ! their endpoints
     integer npts
     real*8 :: zmin, zmax
     type(divpoint_t), pointer :: divpts(:)
  end type divertor_part
  
  type divertor_t
     ! the divertor shape is a sequence of discontiguous divertor parts
     type(divertor_part), pointer :: divparts(:)
  end type divertor_t

  type(divertor_t) :: divertor

  contains
    function divpos_r(divpos) result (p)
      real*8 p
      type(divertor_pos_t), intent(in) :: divpos

      ! between 0 and 1
      real*8 normpos
      integer i, j, index

      index = 1
      goto 1
      entry divpos_z(divpos) result (p)
      index = 2

1     continue
      i = divpos%divseg_pos%num_part
      j = divpos%divseg_pos%num_seg

      normpos = divpos%relpos/divseglen(divpos%divseg_pos)
      p = (1-normpos)* divertor%divparts(i)%divpts(j)%x(index) + &
           &normpos* divertor%divparts(i)%divpts(j+1)%x(index)
    end function divpos_r

    function divpos2divdist(divpos)
      real*8 divpos2divdist
      type(divertor_pos_t), intent(in) :: divpos

      integer i, j

      i = divpos%divseg_pos%num_part
      j = divpos%divseg_pos%num_seg
      divpos2divdist = divertor%divparts(i)%divpts(j)%absdist + divpos%relpos
    end function divpos2divdist

    function divdist2divpos(divdist)
      real*8 divdist(:)
      type(divertor_pos_t) ::  divdist2divpos(ubound(divdist, 1))
      
      integer r, i, j
      
      r = 1
      parts: do i = 1, ubound(divertor%divparts, 1)
         do j = 1, divertor%divparts(i)%npts-1
            startpts: do
               if (divdist(r) < divertor%divparts(i)%divpts(j+1)%absdist) then
                  divdist2divpos(r)%divseg_pos%num_part = i
                  divdist2divpos(r)%divseg_pos%num_seg = j
                  divdist2divpos(r)%relpos = &
                       divdist(r) - divertor%divparts(i)%divpts(j)%absdist
                  r = r + 1
                  if (r > ubound(divdist, 1)) then
                     exit parts
                  endif
               else
                  exit startpts
               endif
            enddo startpts
         enddo
      enddo parts
      if (r <= ubound(divdist, 1)) then
         write (0, *) 'r= ', r
         STOP 'requested point not in divertor'
      endif
    end function divdist2divpos

    function divseg_det(vector, divseg)
      ! positive if the vector endpoint is inside the divertor if the
      ! start point is on the divertor. Requires a defined divertor
      ! orientation which is that the inside of the plasma is at the
      ! left of the vectors of the divertor shape.
      real*8 divseg_det
      real*8 vector(:)
      type(divseg_num_t) :: divseg
      
      real*8 vd(2)
      vd = divvec(divseg)
      
      divseg_det = vd(1)*vector(2) - vd(2)*vector(1)
    end function divseg_det

    SUBROUTINE striked_wall(in_wall, R1, Z1, R2, Z2, &
         phi1, phi2, R_strike, Z_strike, phi_strike, pos_strike)
      real*8 R1, Z1, R2, Z2, &
           phi1, phi2, R_strike, Z_strike, phi_strike
      LOGICAL in_wall
      type(divertor_pos_t) :: pos_strike
      
      real*8 a, b, c, d, e, f, det, a2, b2, scalar_test1, scalar_ref1,&
           &scalar_test2, scalar_ref2 
      type(divseg_num_t) :: i_seg
      integer i, j

      in_wall = .false.

      parts: do i = 1, ubound(divertor%divparts, 1)
         i_seg%num_part=i
         if (divertor%divparts(i)%zmax < 0.) then
            ! lower divertor
            if (divertor%divparts(i)%zmax < Z1 .and. &
                 &divertor%divparts(i)%zmax < Z2) then
               ! we are certainly outside!
               return
            endif
         elseif (divertor%divparts(i)%zmin > 0.) then
            ! upper divertor
            if (divertor%divparts(i)%zmin > Z1 .and. &
                 &divertor%divparts(i)%zmin > Z2) then
               ! we are certainly outside!
               return
            endif
         else
            stop 'no support for divertor crossing the midplane!!!'
         endif
         do j = 1, divertor%divparts(i_seg%num_part)%npts-1
            i_seg%num_seg = j
            a = -(Z_div(i_seg,1)-Z_div(i_seg,2))
            b = (R_div(i_seg,1)-R_div(i_seg,2))
            c = -(Z1-Z2)
            d = (R1-R2)
            e = b*Z_div(i_seg,2) + a*R_div(i_seg,2)
            f = d*Z2 + c*R2
            det = a*d - b*c

            if (det<0.) then
               R_strike = (1./det)*(e*d-b*f)
               Z_strike = (1./det)*(a*f-c*e)

               a2 = (R_strike-R_div(i_seg,1)) * (R_div(i_seg,2)-R_div(i_seg,1))
               b2 = (Z_strike-Z_div(i_seg,1)) * (Z_div(i_seg,2)-Z_div(i_seg,1))

               scalar_test1 = a2 + b2
               scalar_ref1 = (R_div(i_seg,2)-R_div(i_seg,1))**2. &
                    + (Z_div(i_seg,2)-Z_div(i_seg,1))**2.

               scalar_test2 = (R_strike-R1)*(R2-R1) + (Z_strike-Z1)*(Z2-Z1)
               scalar_ref2 = (R2-R1)**2. + (Z2-Z1)**2.

               if ((0..le.scalar_test1).and.(scalar_test1.le.scalar_ref1).and.&
                  (scalar_test2.ge.0.).and.(scalar_test2.le.scalar_ref2)) then
                  in_wall = .true.
                  pos_strike%divseg_pos = i_seg
                  pos_strike%relpos= divseglen(i_seg)*scalar_test1/scalar_ref1
                  phi_strike = phi1 + (phi2-phi1)*scalar_test2/scalar_ref2
                  exit parts
               endif
            endif
         enddo
      enddo parts
    END SUBROUTINE striked_wall

    subroutine alloc_divertor(div, nparts)
      type(divertor_t) :: div
      integer nparts

      allocate(div%divparts(nparts))
    end subroutine alloc_divertor

    subroutine init_divertor(fdiv)
      integer fdiv

      integer nparts, i, j
      real*8 startpos
      logical explicit_coord
      read (fdiv, *) nparts, explicit_coord 

      call alloc_divertor(divertor, nparts)
      startpos = 0.
      do i = 1, nparts
         if (explicit_coord) then
            print *, "part ", i
            call init_divpart(fdiv, divertor%divparts(i))
         else
            print *, "part ", i, " from ", startpos
            call init_divpart(fdiv, divertor%divparts(i), startpos)
            do j = 1,divertor%divparts(i)%npts
               print *, divertor%divparts(i)%divpts(j)%absdist
            enddo
            print *, "to ", startpos
         endif
         print *, "points", divertor%divparts(i)%npts
      end do   
    end subroutine init_divertor

    subroutine save_divertor(fdiv, explicit_coord)
      integer, intent(in) :: fdiv
      logical, intent(in) :: explicit_coord

      integer nparts, i

      nparts = size(divertor%divparts)
      write (fdiv, *) nparts, explicit_coord 

      do i = 1, nparts
         call save_divpart(fdiv, divertor%divparts(i))
      end do   
    end subroutine save_divertor

    subroutine save_divpart(fdiv, divpart)
      integer fdiv
      type(divertor_part), intent(in) :: divpart
      
      integer i

      write (fdiv, *) divpart%npts

      do i = 1, divpart%npts
         write(fdiv, *) divpart%divpts(i)%absdist,divpart%divpts(i)%x
      end do
    end subroutine save_divpart

    subroutine init_divpart_data(divpart, rz, startpos)
      type(divertor_part) :: divpart
      real*8 :: rz(:, :)
      real*8, optional :: startpos

      integer i
      real*8 len

      divpart%npts = size(rz, 2)

      allocate(divpart%divpts(divpart%npts))
      do i = 1, divpart%npts
         divpart%divpts(i)%x = rz(:, i)
      end do
      if (present(startpos)) then
         divpart%divpts(1)%absdist = startpos
         do i = 2, divpart%npts
            len = sqrt(sum((divpart%divpts(i)%x - divpart%divpts(i-1)%x)**2))
            divpart%divpts(i)%absdist = divpart%divpts(i-1)%absdist + len
         end do
         startpos = divpart%divpts(divpart%npts)%absdist
      endif
      divpart%zmin = minval(divpart%divpts(:)%x(2))
      divpart%zmax = maxval(divpart%divpts(:)%x(2))
    end subroutine init_divpart_data

    subroutine init_divpart(fdiv, divpart, startpos)
      integer fdiv
      type(divertor_part) :: divpart
      real*8, optional :: startpos

      integer i
      real*8 len

      read (fdiv, *) divpart%npts
      allocate(divpart%divpts(divpart%npts))
      do i = 1, divpart%npts
         read(fdiv, *) divpart%divpts(i)%absdist,divpart%divpts(i)%x
      end do
      if (present(startpos)) then
         divpart%divpts(1)%absdist = startpos
         do i = 2, divpart%npts
            len = sqrt(sum((divpart%divpts(i)%x - divpart%divpts(i-1)%x)**2))
            divpart%divpts(i)%absdist = divpart%divpts(i-1)%absdist + len
         end do
         startpos = divpart%divpts(divpart%npts)%absdist
      endif
      divpart%zmin = minval(divpart%divpts(:)%x(2))
      divpart%zmax = maxval(divpart%divpts(:)%x(2))
    end subroutine init_divpart

    function R_div(divseg, pt) result (p)
      real*8 p
      type(divseg_num_t) :: divseg
      integer pt
      integer index, i, j

      index = 1
      goto 2
      entry Z_div(divseg, pt) result (p)
      index = 2

2     continue

      i = divseg%num_part
      j = divseg%num_seg

      p = divertor%divparts(i)%divpts(j+pt-1)%x(index)
    end function R_div
    

    function divseglen(divseg)
      real*8 divseglen
      type(divseg_num_t) :: divseg
      integer i, j

      i = divseg%num_part
      j = divseg%num_seg
      
      divseglen = divertor%divparts(i)%divpts(j+1)%absdist - &
           &divertor%divparts(i)%divpts(j)%absdist
    end function divseglen

    function divvec(divseg)
      real*8 divvec(2)
      type(divseg_num_t) :: divseg
      integer i, j

      i = divseg%num_part
      j = divseg%num_seg

      divvec = divertor%divparts(i)%divpts(j+1)%x - &
           &divertor%divparts(i)%divpts(j)%x
    end function divvec
end module divertor_desc
