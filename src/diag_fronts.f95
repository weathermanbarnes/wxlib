! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- Helper functions for front detection
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module adopted from code by Gareth Berry
! Module maintained by Clemens Spensberger
module diag_fronts
  use kind
  use config
  !
  implicit none
  ! The module is used (only) in the front_location subroutine in diag.f95
contains
  !
  ! Find locations where dat is zero by interpolating the 2-dim gridded data
  subroutine find_zeroloc(dat, nx,ny,nn, NaN, zeroloc, zerocnt)
    real(kind=nr), intent(in)  ::  dat(ny,nx)
    real(kind=nr), intent(out) :: zeroloc(nn,2_ni)
    real(kind=nr), intent(in) :: NaN
    integer(kind=ni), intent(in) :: nx,ny, nn
    integer(kind=ni), intent(out) :: zerocnt 
    !
    real   (kind=nr) :: latB, dist
    integer(kind=ni) :: i,j, ip1,gotflag
    ! -----------------------------------------------------------------
    !
    ! (1) scan along x
    zerocnt = 0_ni
    do i = 1_ni,nx-1_ni
       do j = 1_ni,ny-1_ni
          gotflag = 0_ni
          ! Zero line hits a grid point
          if (dat(j,i) == 0.0_nr) then
             zerocnt = zerocnt + 1_ni
             zeroloc(zerocnt,1_ni) = i
             zeroloc(zerocnt,2_ni) = j
             gotflag = 1_ni
          ! interpolate to find line, i direction first
          else   
             ip1 = i+1_ni
             if (dat(j,ip1) /= NaN .and. gotflag == 0_ni .and. & 
               &   dat(j,i) /= NaN) then
                if ((dat(j,i) > 0.0_nr .and. dat(j,ip1) > 0.0_nr) .or. &
                  & (dat(j,i) < 0.0_nr .and. dat(j,ip1) > 0.0_nr)) then
                   latB = dat(j,ip1)
                   if (dat(j,i) /= latB) then 
                      zerocnt = zerocnt + 1_ni
                      zeroloc(zerocnt,2_ni) = j
                      zeroloc(zerocnt,1_ni) = i + dat(j,i)/(dat(j,i) - dat(j,ip1))
                   ! points same magn, diff sign
                   else
                      zerocnt = zerocnt + 1_ni
                      zeroloc(zerocnt,2_ni) = j
                      zeroloc(zerocnt,1_ni) = i + 0.5_nr
                   end if    ! latB
                   gotflag = 1
                end if    ! diff signs
             end if    ! Missin ip1
          end if     ! zero exactly at grid point
       end do   
    end do
    ! 
    ! (2) scan along y
    do i = 1_ni,nx-1_ni
       do j = 1_ni,ny-1_ni
          gotflag = 0_ni
          ! Zero line hits a grid point
          if (dat(j,i) == 0) then
             zerocnt = zerocnt + 1_ni
             zeroloc(zerocnt,1_ni) = i
             zeroloc(zerocnt,2_ni) = j
             gotflag = 1_ni
          ! interpolate to find line, j direction first
          else   
             ip1 = j + 1_ni
             if (dat(ip1,i) /= NaN .and. gotflag == 0_ni .and. &
               & dat(j,i) /= NaN) then
                if ((dat(j,i) > 0.0_nr .and. dat(ip1,i) < 0.0_nr) .or. &
                  & (dat(j,i) < 0.0_nr .and. dat(ip1,i) > 0.0_nr)) then
                   latB = dat(ip1,i)
                   if (dat(j,i) /= latB) then
                      zerocnt = zerocnt + 1_ni
                      zeroloc(zerocnt,2_ni) = j + dat(j,i)/(dat(j,i) - dat(ip1,i))
                      zeroloc(zerocnt,1_ni) = i
                   ! points same magn, diff sign
                   else
                      zerocnt = zerocnt + 1_ni
                      zeroloc(zerocnt,2_ni) = j + 0.5_nr          
                      zeroloc(zerocnt,1_ni) = i
                   end if    ! latB
                   gotflag = 1
                end if    ! diff signs
             end if    ! Missin ip1
          end if    ! zero exactly at grid point
       end do   
    end do
    !
    return
  end  subroutine find_zeroloc
  !
  !
  ! Join a cloud of frontal points into frontal lines
  subroutine linejoin(cnt,lat,lon,searchrad,reclat,reclon)
    real(kind=nr), intent(in)  :: lat(cnt), lon(cnt), searchrad 
    real(kind=nr), intent(out) :: reclat(cnt,cnt), reclon(cnt,cnt)
    integer(kind=ni) :: cnt
    !
    real   (kind=nr) :: minlength, mindist, dist
    integer(kind=ni) :: n,m, na,nn, cur,next
    logical :: look(cnt)
    ! -----------------------------------------------------------------
    !
    look(:) = .false.
    na = 0_ni
    do n = 1_ni,cnt
       nn = 0_ni
       ! new line segment at point n
       if (.not. look(n)) then
          na = na + 1_ni
          !
          next = n
          mindist = 0.0_nr
          ! do as long as there is a continuation of the front
          do while (mindist <= searchrad)
             nn = nn + 1_ni
             look(next) = .true.
             reclat(na,nn) = lat(next)
             reclon(na,nn) = lon(next)
             cur = next
             !
             mindist = 181.0_nr
             ! find the closest point of all points not considered so far
             do m = 1_ni,cnt
                if (.not. look(m)) then
                   dist = sqrt((lat(cur)-lat(m))**2_ni + (lon(cur)-lon(m))**2_ni)
                   !
                   if (dist /= 0.0_nr .and. mindist > dist) then
                      mindist = dist 
                      next = m
                   end if
                end if
             end do
             !
          end do
       end if
    end do
    !print*,na,' Line segments'
    !
    return
  end subroutine linejoin
end module 
