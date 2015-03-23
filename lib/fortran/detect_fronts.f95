! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- Helper functions for front detection
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module adopted from code by Gareth Berry
! Module maintained by Clemens Spensberger
module detect_fronts
  use kind
  use config
  use derivatives
  !
  implicit none
  ! The module is used in the front/line/axis_location subroutines in diag.f95
contains
  !
  !
  subroutine line_locate(lines,lnoff,nx,ny,nz,no,nf,lnloc,lnint,searchrad,minlen,NaN,dx,dy)
    use consts
    !
    real(kind=nr), intent(in) :: lnloc(nz,ny,nx), lnint(nz,ny,nx), &
                 &               dx(ny,nx), dy(ny,nx), searchrad, NaN, minlen
    real(kind=nr), intent(out) :: lines(nz,no,3_ni), lnoff(nz,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf
    !f2py depend(nx,ny,nz) lnint
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) lines, lnoff
    !
    integer(kind=ni), parameter :: nn = 30000_ni
    !
    real   (kind=nr), allocatable :: reci(:,:), recj(:,:), linelen(:)
    integer(kind=ni), allocatable :: lineptcnt(:)
    !
    real(kind=nr) :: zeroloc(nn,2_ni)
    integer(kind=ni) :: i,j,k, di, m, n, zerocnt, ptcnt, linecnt, off
    ! -----------------------------------------------------------------
    !
    do k = 1_ni,nz
       write(*,'(I5,A4,I5,A)', advance='no') k, 'of', nz, cr
       !
       ! find fronts
       call find_zeroloc(lnloc(k,:,:), nx,ny,nn, NaN, zeroloc,zerocnt)
       ! 
       ! searchrad is in grid point indexes, as zero locations are found at grid
       ! resolution, hence the number of neighbours for a given searchrad does not
       ! depend on location within the grid
       !
       allocate(recj(zerocnt,zerocnt), reci(zerocnt,zerocnt), lineptcnt(zerocnt), linelen(zerocnt) )
       reci(:,:) = NaN
       recj(:,:) = NaN
       linecnt    = 0_ni ! number of lines
       ptcnt      = 0_ni ! total numer of points
       lineptcnt(:) = 0_ni ! number of points per line
       call linejoin(zerocnt, nf*3_ni, nx,ny, zeroloc(:,2_ni), zeroloc(:,1_ni), searchrad, &
               & recj, reci, linelen, lineptcnt, dx,dy) 
       !
       off = 0_ni
       do n = 1_ni,zerocnt
          if (recj(n,1_ni) == NaN) then
             exit
          end if
          !
          ! filter fronts by length
          if (linelen(n) >= minlen) then
          !if (.true.) then
             linecnt = linecnt + 1_ni
             ptcnt = ptcnt + lineptcnt(n)
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
             do m = 1_ni,lineptcnt(n)
                lines(k,off+m,1_ni) = reci(n,m)
                lines(k,off+m,2_ni) = recj(n,m)
                lines(k,off+m,3_ni) = lnint(k,int(recj(n,m),ni),int(reci(n,m),ni))
             end do
             lnoff(k,linecnt) = off
             off = off + lineptcnt(n)
          end if
       end do
       ! Save the ending of the last front by saving the beginning of the
       ! first non-existant
       lnoff(k,linecnt+1_ni) = off
       !
       deallocate(reci, recj, lineptcnt, linelen)
       !
    end do ! loop over k
    !
    return
  end subroutine
  !
  ! TODO: Reduce duplication with line_locate
  subroutine maxline_locate(lines,lnoff,nx,ny,nz,no,nf,dat,thres,searchrad,minlen,NaN,dx,dy)
    use consts
    !
    real(kind=nr), intent(in) :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx), thres, searchrad, NaN, minlen
    real(kind=nr), intent(out) :: lines(nz,no,3_ni), lnoff(nz,nf)
    integer(kind=ni) :: nx,ny,nz, no, nf
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) lines, lnoff
    !
    integer(kind=ni), parameter :: nn = 30000_ni
    !
    real   (kind=nr), allocatable :: reci(:,:), recj(:,:), linelen(:)
    integer(kind=ni), allocatable :: lineptcnt(:)
    !
    real(kind=nr) :: maxloc_(nz,nn,2_ni)
    integer(kind=ni) :: i,j,k, di, m, n, maxcnt(nz), ptcnt, linecnt, off
    ! -----------------------------------------------------------------
    !
    ! find lines by zero-criterion
    call find_maxloc(dat(:,:,:),thres, nx,ny,nz, nn, NaN, maxloc_(:,:,:),maxcnt(:), dx,dy)
    !
    do k = 1_ni,nz
       write(*,'(I5,A4,I5,A)', advance='no') k, 'of', nz, cr
       !
       ! searchrad is in grid point indexes, as zero locations are found at grid
       ! resolution, hence the number of neighbours for a given searchrad does not
       ! depend on location within the grid
       !
       allocate(recj(maxcnt(k),maxcnt(k)), reci(maxcnt(k),maxcnt(k)), lineptcnt(maxcnt(k)), linelen(maxcnt(k)) )
       reci(:,:) = NaN
       recj(:,:) = NaN
       linecnt    = 0_ni ! number of lines
       ptcnt      = 0_ni ! total numer of points
       lineptcnt(:) = 0_ni ! number of points per line
       call linejoin(maxcnt(k), nf*3_ni, nx,ny, maxloc_(k,:,2_ni), maxloc_(k,:,1_ni), searchrad, &
               & recj, reci, linelen, lineptcnt, dx,dy) 
       !
       off = 0_ni
       do n = 1_ni,maxcnt(k)
          if (recj(n,1_ni) == NaN) then
             exit
          end if
          !
          ! filter fronts by length
          if (linelen(n) >= minlen) then
          !if (.true.) then
             linecnt = linecnt + 1_ni
             ptcnt = ptcnt + lineptcnt(n)
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
             do m = 1_ni,lineptcnt(n)
                lines(k,off+m,1_ni) = reci(n,m)
                lines(k,off+m,2_ni) = recj(n,m)
                lines(k,off+m,3_ni) = dat(k,int(recj(n,m),ni),int(reci(n,m),ni))
             end do
             lnoff(k,linecnt) = off
             off = off + lineptcnt(n)
          end if
       end do
       ! Save the ending of the last front by saving the beginning of the
       ! first non-existant
       lnoff(k,linecnt+1_ni) = off
       !
       deallocate(reci, recj, lineptcnt, linelen)
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
    integer(kind=ni) :: i,j, ip1
    ! -----------------------------------------------------------------
    !
    ! (1) scan along x
    zerocnt = 0_ni
    do i = 1_ni,nx-1_ni
       ip1 = i+1_ni
       do j = 1_ni,ny
          ! Zero line hits a grid point
          if (dat(j,i) == 0.0_nr) then
             if (dat(j,ip1) == 0.0_nr) then
                zerocnt = zerocnt + 1_ni
                zeroloc(zerocnt,2_ni) = j
                zeroloc(zerocnt,1_ni) = i + 0.5_nr
             else
                zerocnt = zerocnt + 1_ni
                zeroloc(zerocnt,1_ni) = i
                zeroloc(zerocnt,2_ni) = j
             end if
          ! interpolate to find line, i direction first
          else   
             if (dat(j,ip1) /= NaN .and. dat(j,i) /= NaN) then
                if ((dat(j,i) > 0.0_nr .and. dat(j,ip1) < 0.0_nr) .or. &
                  & (dat(j,i) < 0.0_nr .and. dat(j,ip1) > 0.0_nr)) then
                   zerocnt = zerocnt + 1_ni
                   zeroloc(zerocnt,2_ni) = j
                   zeroloc(zerocnt,1_ni) = i + dat(j,i)/(dat(j,i) - dat(j,ip1))
                end if    ! diff signs
             end if    ! Missin ip1
          end if    ! zero exactly at grid point
       end do   
    end do
    !
    ! take into account periodicity in x
    if ( grid_cyclic_ew ) then
       i = nx
       ip1 = 1_ni
       do j = 1_ni,ny
          ! Zero line hits a grid point
          if (dat(j,i) == 0.0_nr) then
             if (dat(j,ip1) == 0.0_nr) then
                zerocnt = zerocnt + 1_ni
                zeroloc(zerocnt,2_ni) = j
                zeroloc(zerocnt,1_ni) = i + 0.5_nr
             else
                zerocnt = zerocnt + 1_ni
                zeroloc(zerocnt,1_ni) = i
                zeroloc(zerocnt,2_ni) = j
             end if
          ! interpolate to find line, i direction first
          else   
             if (dat(j,ip1) /= NaN .and. dat(j,i) /= NaN) then
                if ((dat(j,i) > 0.0_nr .and. dat(j,ip1) < 0.0_nr) .or. &
                  & (dat(j,i) < 0.0_nr .and. dat(j,ip1) > 0.0_nr)) then
                   zerocnt = zerocnt + 1_ni
                   zeroloc(zerocnt,2_ni) = j
                   zeroloc(zerocnt,1_ni) = i + dat(j,i)/(dat(j,i) - dat(j,ip1))
                end if    ! diff signs
             end if    ! Missin ip1
          end if    ! zero exactly at grid point
       end do
    end if
    ! 
    ! (2) scan along y
    do i = 1_ni,nx
       do j = 1_ni,ny-1_ni
          ip1 = j + 1_ni
          ! Zero line hits a grid point
          if (dat(j,i) == 0.0_nr) then
             if (dat(ip1,i) == 0.0_nr) then
                zerocnt = zerocnt + 1_ni
                zeroloc(zerocnt,2_ni) = j + 0.5_nr          
                zeroloc(zerocnt,1_ni) = i
             else
                zerocnt = zerocnt + 1_ni
                zeroloc(zerocnt,1_ni) = i
                zeroloc(zerocnt,2_ni) = j
             end if
          ! interpolate to find line, j direction first
          else   
             if (dat(ip1,i) /= NaN .and. dat(j,i) /= NaN) then
                if ((dat(j,i) > 0.0_nr .and. dat(ip1,i) < 0.0_nr) .or. &
                  & (dat(j,i) < 0.0_nr .and. dat(ip1,i) > 0.0_nr)) then
                   zerocnt = zerocnt + 1_ni
                   zeroloc(zerocnt,2_ni) = j + dat(j,i)/(dat(j,i) - dat(ip1,i))
                   zeroloc(zerocnt,1_ni) = i
                end if
             end if    ! Missin ip1
          end if    ! zero exactly at grid point
       end do   
    end do
    !
    return
  end  subroutine find_zeroloc
  !
  ! Find maximum axis locations by interpolating the 2-dim gridded data
  subroutine find_maxloc(dat,thres, nx,ny,nz, nn, NaN, maxloc_, maxcnt, dx,dy)
    real(kind=nr), intent(in)  ::  dat(nz,ny,nx), dx(ny,nx), dy(ny,nx), thres
    real(kind=nr), intent(out) :: maxloc_(nz,nn,2_ni)
    real(kind=nr), intent(in) :: NaN
    integer(kind=ni), intent(in) :: nx,ny,nz, nn
    integer(kind=ni), intent(out) :: maxcnt(nz)
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nz) maxloc_, maxcnt
    !
    real(kind=nr) :: datxx(nz,ny,nx), datyy(nz,ny,nx), datxy(nz,ny,nx), &
                &    datx (nz,ny,nx), daty (nz,ny,nx), dats, ds, di, dj
    real(kind=nr) :: hess(2_ni,2_ni), evalr(2_ni), evali(2_ni), evec(2_ni,2_ni), &
                &    dummy(2_ni,2_ni), work(68_ni)
    integer(kind=ni) :: i,j,k, minidx, info, cnt
    ! -----------------------------------------------------------------
    !
    call ddx(datx, nx,ny,nz, dat, dx,dy)
    call ddy(daty, nx,ny,nz, dat, dx,dy)
    call ddx2(datxx, nx,ny,nz, dat, dx,dy)
    call ddy2(datyy, nx,ny,nz, dat, dx,dy)
    call ddxy(datxy, nx,ny,nz, dat, dx,dy)
    !
    do k = 1_ni,nz
       cnt = 0_ni
       do i = 1_ni,nx
          do j = 1_ni,ny
             if ( dat(k,j,i) < thres ) cycle
             !
             hess(:,:) = reshape((/ datxx(k,j,i), datxy(k,j,i), datxy(k,j,i), datyy(k,j,i) /), (/ 2_ni, 2_ni /) )
             call dgeev('N','V',2_ni,hess,2_ni,evalr,evali,dummy,2_ni,evec,2_ni,work,68_ni,info)
             !
             ! Minimum eigenvalue: Largest negative curvature
             minidx = minloc(evalr, 1_ni)
             if ( evalr(minidx) > 0.0_nr ) cycle
             ! Enforce an elongated aspect ratio
             if ( evalr(minidx)/abs(evalr(3_ni-minidx)) > -1.0_nr ) cycle

             ! Slope in direction of the largest negative curvature
             dats = datx(k,j,i)*evec(1_ni,minidx) + daty(k,j,i)*evec(2_ni,minidx)
             !
             ! Extrapolate to find zero-slope
             ds = -dats/evalr(minidx)
             di = ds*evec(1_ni,minidx)/dx(j,i)
             dj = ds*evec(2_ni,minidx)/dy(j,i)
             if ( abs(di) > 0.5_nr .or. abs(dj) > 0.5_nr ) cycle
             !
             cnt = cnt + 1_ni
             maxloc_(k,cnt,2_ni) = j + 0.5_nr + dj
             maxloc_(k,cnt,1_ni) = i + 0.5_nr + di
             !
             !if ( i == 371_ni .and. j == 60_ni .and. k == 1_ni ) then
             !   hess(:,:) = reshape((/ datxx(k,j,i), datxy(k,j,i), datxy(k,j,i), datyy(k,j,i) /), (/ 2_ni, 2_ni /) )
             !   write(*,*) 'matrix', hess
             !   write(*,*) 'eval 1', evalr(1_ni), evali(1_ni)
             !   write(*,*) 'evec 1', evec(:,1_ni)
             !   write(*,*) 'eval 2', evalr(2_ni), evali(2_ni)
             !   write(*,*) 'evec 2', evec(:,2_ni)
             !   write(*,*) 'work len', work(1_ni)
             !   write(*,*) 'info', info
             !   write(*,*) 'min val', evalr(minidx)
             !   write(*,*) 'slope xys', datx(k,j,i), daty(k,j,i), dats
             !   write(*,*) 'zero dist', ds
             !   write(*,*) 'dpos', di, dj
             !   stop 1
             !end if
          end do
       end do
       maxcnt(k) = cnt
    end do
    !
  end subroutine
  !
  ! Join a cloud of frontal points into frontal lines
  subroutine linejoin(cnt,nstruct_max,nx,ny,jidx,iidx,searchrad,recj,reci,linelen,lineptcnt,dx,dy)
    real(kind=nr), intent(in)  :: jidx(cnt), iidx(cnt), searchrad, dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: recj(cnt,cnt), reci(cnt,cnt), linelen(cnt)
    integer(kind=ni), intent(in) :: cnt, nx,ny, nstruct_max
    integer(kind=ni), intent(out) :: lineptcnt(cnt)
    !
    real   (kind=nr) :: dists(cnt,cnt), longestdist
    integer(kind=ni) :: n,m,l, nn, inter(cnt,cnt), oidx(cnt), structcnt(nstruct_max), &
            &           donecnt,prevcnt, nstruct, startidx,endidx, &
            &           longestpath(2_ni), linecnt, pos,startpos,endpos
    logical :: used(cnt)
    ! -----------------------------------------------------------------
    !
    dists(:,:) = -1.0_nr ! distance metric between all pairs of points (sorted coherent structures);
    !                    ! -1.0_nr represents infinite distance ( = no connection)
    inter(:,:) = 0_ni    ! intermediate point defining the shortest path between each pair points
    used(:) = .false.    ! Used in the DFS? Is eqivalent to: Already assigned to a choherent structure?
    oidx(:) = 0_ni       ! Mapping of the new index (used in the dists matrix) to the old index used in jidx/iidx
    !
    !
    ! 1. Step: Find coherent structures via a Depth-First-Search
    nstruct = 0_ni
    donecnt = 0_ni
    do n = 1_ni,cnt
       if ( used(n) ) cycle
       !
       nstruct = nstruct + 1_ni
       if ( nstruct > nstruct_max ) then
          write(*,*) 'Found more than NSTRUCT_MAX=', nstruct_max, 'structures'
          stop 1
       end if
       prevcnt = donecnt
       call depth_first_search(cnt,donecnt, nx,ny, n,used,dists, oidx, jidx,iidx,searchrad, dx,dy)
       structcnt(nstruct) = donecnt - prevcnt
    end do
    !
    ! For each of the structures
    linecnt = 0_ni
    startidx = 1_ni
    do nn = 1_ni,nstruct
       endidx = startidx - 1_ni + structcnt(nn)
       !
       ! 2. Step: Find shortest paths between all pairs of points in the structures 
       !          using the Floyd-Warshall algorithm 
       do n = startidx,endidx
          do m = startidx,endidx
             do l = startidx,endidx
                if ( dists(n,m) > 0.0_nr .and. dists(n,l) > 0.0_nr ) then
                   if ( dists(n,m) + dists(n,l) < dists(m,l) .or. dists(m,l) < 0.0_nr ) then
                      dists(m,l) = dists(n,m) + dists(n,l)
                      inter(m,l) = n
                   end if
                end if
             end do
          end do
       end do
       !
       ! 3. Step: Reconstruct path for the longest shortest paths within each structure 
       !          using the "inter(mediate)" information from step 2
       longestpath = maxloc(dists(startidx:endidx,startidx:endidx))
       longestdist = maxval(dists(startidx:endidx,startidx:endidx))
       startpos = longestpath(1_ni) + startidx - 1_ni
       endpos   = longestpath(2_ni) + startidx - 1_ni
       !write(*,*) 'Debug struct ', nn, startidx, endidx, longestdist
       !write(*,*) 'Debug struct ', nn, jidx(startpos), iidx(startpos), jidx(endpos), iidx(endpos)
       !
       if ( longestdist >= minlen ) then
          ! Initialise and record start position
          linecnt = linecnt + 1_ni
          lineptcnt(linecnt) = 1_ni
          linelen(linecnt) = longestdist
          pos = startpos
          recj(linecnt,1_ni) = jidx(oidx(pos))
          reci(linecnt,1_ni) = iidx(oidx(pos))
          !
          ! Recursively Construct and record the path to the end position
          call construct_path(cnt, startpos,endpos, inter, oidx, jidx,iidx, &
                  & lineptcnt(linecnt), recj(linecnt,:), reci(linecnt,:))
          !
          !write(*,*) linecnt, linelen(linecnt), lineptcnt(linecnt)
       end if
       !
       startidx = endidx + 1_ni
    end do
    !
    return
  end subroutine linejoin
  !
  recursive subroutine depth_first_search(cnt,donecnt, nx,ny, n,used,dists, oidx, jidx,iidx,searchrad, dx,dy)
    real(kind=nr), intent(in)  :: jidx(cnt), iidx(cnt), searchrad, dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(inout) :: dists(cnt,cnt)
    integer(kind=ni), intent(in) :: cnt, nx,ny, n
    integer(kind=ni), intent(inout) :: oidx(cnt), donecnt
    logical, intent(inout) :: used(cnt)
    !
    real   (kind=nr) :: dist, distji
    integer(kind=ni) :: i,j, m, nn,nm, di
    ! -----------------------------------------------------------------
    !
    donecnt = donecnt + 1_ni
    nn = donecnt
    used(n) = .true.
    oidx(nn) = n
    dists(nn,nn) = 0.0_nr
    !
    i = iidx(n)
    j = jidx(n)
    !
    do m = 1_ni,cnt
       if ( used(m) ) cycle
       !
       di = iidx(m) - iidx(n)
       if ( di > size(dx,2_ni)/2_ni ) di = di - size(dx,2_ni)
       if ( di < -size(dx,2_ni)/2_ni ) di = di + size(dx,2_ni)
       dist = sqrt( (dx(j,i)*(di)/2.0_nr)**2.0_nr + &
                  & (dy(j,i)*(jidx(m)-jidx(n))/2.0_nr)**2.0_nr   )
       distji = sqrt(di**2_ni + (jidx(m)-jidx(n))**2_ni)
       !
       ! if we found a pair of close points
       if ( distji <= searchrad ) then
          nm = donecnt + 1_ni
          !write(*,*) n, m, nn, nm
          dists(nn,nm) = dist
          dists(nm,nn) = dist
          call depth_first_search(cnt,donecnt, nx,ny, m,used,dists, oidx, jidx,iidx,searchrad, dx,dy)
       end if
    end do
  end subroutine
  ! 
  recursive subroutine construct_path(cnt, startpos,endpos, inter, oidx,jidx,iidx, lineptcnt, recj, reci)
    real(kind=nr), intent(in)  :: jidx(cnt), iidx(cnt)
    real(kind=nr), intent(out) :: recj(cnt), reci(cnt)
    integer(kind=ni), intent(in) :: cnt, inter(cnt,cnt), oidx(cnt), startpos, endpos
    integer(kind=ni), intent(inout) :: lineptcnt
    !
    integer(kind=ni) :: interpos
    ! -----------------------------------------------------------------
    !
    interpos = inter(startpos,endpos)
    ! Nothing in between: directly connected
    if ( interpos == 0_ni ) then
       lineptcnt = lineptcnt + 1_ni
       recj(lineptcnt) = jidx(oidx(endpos))
       reci(lineptcnt) = iidx(oidx(endpos))
    ! No connection at all
    else if ( inter(startpos,endpos) < 0_ni ) then
       write(*,*) 'ERROR: Trying to join unconnected points'
       stop 1
    ! At least one point in between start and end
    else
       call construct_path(cnt, startpos,interpos, inter, oidx, jidx,iidx, lineptcnt, recj, reci)
       call construct_path(cnt, interpos,endpos, inter, oidx, jidx,iidx, lineptcnt, recj, reci)
    end if
  end subroutine
  !
end module 
