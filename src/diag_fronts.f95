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
  ! The module is used in the front/line/axis_location subroutines in diag.f95
contains
  !
  !
  subroutine line_locate(lines,lnoff,nx,ny,nz,no,nf,lnloc,lnint,searchrad,minlen,NaN)
    use consts
    !
    real(kind=nr), intent(in) :: lnloc(nz,ny,nx), lnint(nz,ny,nx), &
                 &               searchrad, NaN
    real(kind=nr), intent(out) :: lines(nz,no,3_ni), lnoff(nz,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf, minlen
    !f2py depend(nx,ny,nz) lnint
    !f2py depend(nz) lines, lnoff
    !
    integer(kind=ni), parameter :: nn = 30000_ni
    !
    real   (kind=nr), allocatable :: reci(:,:), recj(:,:)
    integer(kind=ni), allocatable :: linelen(:)
    !
    real(kind=nr) :: zeroloc(nn,2_ni)
    integer(kind=ni) :: k, m, n, zerocnt, ptcnt, linecnt, off
    ! -----------------------------------------------------------------
    !
    do k = 1_ni,nz
       write(*,'(I5,A4,I5,A)', advance='no') k, 'of', nz, cr
       !
       ! find fronts
       call find_zeroloc(lnloc(k,:,:), nx,ny,nn, NaN, zeroloc,zerocnt)
       ! 
       ! searchrad is in grid point indexes, as it is easier to transform one
       ! scalar instead of two arrays. For standard grids the two options are
       ! equivalent. In the long-term one should try to adapt the search radius and
       ! the minimum length to SI length scales (km)
       !
       ! todo why sort points?
       !
       allocate(recj(zerocnt,zerocnt), reci(zerocnt,zerocnt), linelen(zerocnt) )
       reci(:,:) = NaN
       recj(:,:) = NaN
       call linejoin(zerocnt, zeroloc(:,2_ni), zeroloc(:,1_ni), searchrad, recj, reci) 
       !
       ! filter by length
       linecnt    = 0_ni ! number of lines
       ptcnt      = 0_ni ! total numer of points
       linelen(:) = 0_ni ! number of points per line
       !
       off = 0_ni
       do n = 1_ni,zerocnt
          if (recj(n,1_ni) == NaN) then
             exit
          end if
          do m = 1_ni,zerocnt
             if (recj(n,m) == NaN) then
                exit
             end if
             linelen(n) = linelen(n) + 1_ni
          end do
          !
          ! filter fronts by length
          if (linelen(n) >= minlen) then
             linecnt = linecnt + 1_ni
             ptcnt = ptcnt + linelen(n)
             !
             ! check if results larger than output array
             if (ptcnt > no) then
                write(*,*) 'Found more points than output array allows: ', no
                stop 1
             end if
             if (linecnt > nf) then
                write(*,*) 'Found more fronts than output array allows: ', nf
                stop 1
             end if
             !
             ! write into output arrays fr and froff
             do m = 1_ni,linelen(n)
                lines(k,off+m,1_ni) = reci(n,m)
                lines(k,off+m,2_ni) = recj(n,m)
                lines(k,off+m,3_ni) = lnint(k,int(recj(n,m),ni),int(reci(n,m),ni))
             end do
             lnoff(k,linecnt) = off
             off = off + linelen(n)
          end if
       end do
       ! Save the ending of the last front by saving the beginning of the
       ! first non-existant
       lnoff(k,linecnt+1_ni) = off
       !
       deallocate(reci, recj, linelen)
       !
    end do ! loop over k
    !
    return
  end subroutine
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
                if ((dat(j,i) > 0.0_nr .and. dat(j,ip1) < 0.0_nr) .or. &
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
                   gotflag = 1_ni
                end if    ! diff signs
             end if    ! Missin ip1
          end if    ! zero exactly at grid point
       end do   
    end do
    !
    ! take into account periodicity in x
    if ( grid_cyclic_ew ) then
       do j = 1_ni,ny-1_ni
          gotflag = 0_ni
          ! Zero line hits a grid point
          if (dat(j,nx) == 0.0_nr) then
             zerocnt = zerocnt + 1_ni
             zeroloc(zerocnt,1_ni) = nx
             zeroloc(zerocnt,2_ni) = j
             gotflag = 1_ni
          ! interpolate to find line, i direction first
          else   
             ip1 = 1_ni
             if (dat(j,ip1) /= NaN .and. gotflag == 0_ni .and. & 
               &   dat(j,nx) /= NaN) then
                if ((dat(j,nx) > 0.0_nr .and. dat(j,ip1) < 0.0_nr) .or. &
                  & (dat(j,nx) < 0.0_nr .and. dat(j,ip1) > 0.0_nr)) then
                   latB = dat(j,ip1)
                   if (dat(j,nx) /= latB) then 
                      zerocnt = zerocnt + 1_ni
                      zeroloc(zerocnt,2_ni) = j
                      zeroloc(zerocnt,1_ni) = nx + dat(j,nx)/(dat(j,nx) - dat(j,ip1))
                   ! points same magn, diff sign
                   else
                      zerocnt = zerocnt + 1_ni
                      zeroloc(zerocnt,2_ni) = j
                      zeroloc(zerocnt,1_ni) = nx + 0.5_nr
                   end if    ! latB
                   gotflag = 1_ni
                end if    ! diff signs
             end if    ! Missin ip1
          end if     ! zero exactly at grid point
       end do   
    end if
    ! 
    ! (2) scan along y
    do i = 1_ni,nx-1_ni
       do j = 1_ni,ny-1_ni
          gotflag = 0_ni
          ! Zero line hits a grid point
          if (dat(j,i) == 0.0_nr) then
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
                   gotflag = 1_ni
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
             mindist = 362.0_nr
             ! find the closest point of all points not considered so far
             ! todo: generalize constants for different grids.
             do m = 1_ni,cnt
                if (.not. look(m)) then
                   dist = sqrt((lat(cur)-lat(m))**2_ni + (lon(cur)-lon(m))**2_ni)
                   ! take periodicity into account
                   if ( grid_cyclic_ew .and. dist > 705.0_nr ) then
                      dist = sqrt((lat(cur)-lat(m))**2_ni + (abs(lon(cur)-lon(m))-720.0_nr)**2_ni)
                   end if
                   !
                   if ( dist /= 0.0_nr .and. mindist > dist ) then
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
