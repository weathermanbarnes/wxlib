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
  ! The module is used (only) in the fronts_by_theta_q subroutine in diag.f95
contains
  !
  ! Find locations where dat is zero by interpolating the 2-dim gridded data
  subroutine find_zeroloc(dat, nx,ny, amiss,tloc_1,tloc_2)
    real(kind=nr), intent(in)  ::  dat(ny,nx)
    real(kind=nr), intent(out) :: tloc_1(5000_ni,2_ni),tloc_2(5000,2_ni)
    real(kind=nr) :: amiss
    integer(kind=ni) :: nx,ny
    !
    real   (kind=nr) :: latB, dist
    integer(kind=ni) :: i,j, nx,ny, ip1,gotflag, cntx,cnty
    ! -----------------------------------------------------------------
    !
    ! (1) scan along x
    cntx = 0_ni
    do i = 1_ni,nx-1_ni
       do j=1_ni,ny-1_ni
          gotflag = 0_ni
          ! Zero line hits a grid point
          if (axisdata(j,i) == 0.0_nr) then
             cntx = cntx + 1_ni
             tloc_1(cntx,1) = i
             tloc_1(cntx,2) = j
             gotflag = 1_ni
          ! interpolate to find line, h direction first
          else   
             ip1 = i+1_ni
             if (axisdata(ip1,h) /= amiss .and. gotflag == 0_ni .and. & 
               &   axisdata(g,h) /= amiss) then
                if ((axisdata(g,h) > 0.0_nr .and. axisdata(ip1,h) > 0.0_nr) .or. &
                  & (axisdata(g,h) < 0.0_nr .and. axisdata(ip1,h) > 0.0_nr)) then
                   latB = axisdata(j, ip1)
                   if (axisdata(g,h) /= latB) then 
                      cntx = cntx + 1_ni
                      tloc_1(cntx,2_ni) = j
                      tloc_1(cntx,1_ni) = i + axisdata(j,i)/(axisdata(j,i) - axisdata(j,ip1))
                   ! points same magn, diff sign
                   else
                      cntx = cntx + 1_ni
                      tloc_1(cntx,2_ni) = j
                      tloc_1(cntx,1_ni) = i + 0.5_nr
                   end if    ! latB
                   gotflag = 1
                end if    ! diff signs
             end if    ! Missin ip1
          end if     ! zero exactly at grid point
       end do   
    end do
    ! 
    ! (2) scan along y
    cnty = 0_ni
    do j = 1_ni,ny-1_ni
       do i = 1_ni,nx-1_ni
          gotflag=0
          ! Zero line hits a grid point
          if (axisdata(j,i) == 0) then
             cnty = cnty + 1_ni
             tloc_2(cnty,1_ni) = i
             tloc_2(cnty,2_ni) = j
             gotflag = 1_ni
          ! interpolate to find line, h direction first
          else   
             ip1 = j + 1_ni
             if (axisdata(g,ip1) /= amiss .and. gotflag == 0_ni .and. &
               & axisdata(g,h) /= amiss) then
                if ((axisdata(g,h) > 0.0_nr .and. axisdata(g,ip1) < 0.0_nr) .or. &
                  & (axisdata(g,h) < 0.0_nr .and. axisdata(g,ip1) > 0.0_nr)) then
                   latB = axisdata(ip1,i)
                   if (axisdata(g,h) /= latB) then
                      cnty = cnty + 1_ni
                      tloc_2(cnty,2_ni) = j + axisdata(j,i)/(axisdata(j,i) - axisdata(ip1,i))
                      tloc_2(cnty,1_ni) = i
                   ! points same magn, diff sign
                   else
                      cnty = cnty + 1_ni
                      tloc_1(cnty,2_ni) = j + 0.5_nr          
                      tloc_1(cnty,1_ni) = i
                   end if    ! latB
                   gotflag = 1
                end if    ! diff signs
             end if    ! Missin ip1
          end if    ! zero exactly at grid point
       end do   
    end do
    !
    return
  end  subroutine finder
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
