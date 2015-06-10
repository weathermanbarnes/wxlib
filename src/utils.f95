! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- utility functions
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module utils
  use kind
  use config
  !
  implicit none
contains
  !
  ! Prepare global data for FFT: 
  !
  ! Returns the data extended along complementary meridians (for fft)
  ! For each lon, the reflected (lon+180) is attached below
  ! so that data is periodic in x and y.
  ! NOTE: Input data must be lats -90 to 90!!! and nx must be even
  subroutine mirror_y_domain(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,2_ni*ny-2_ni,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz,iextra
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    do k=1_ni,nz
       do i=1_ni,nx
          iextra=i+nx/2_ni
          if (iextra.GT.nx) then
             iextra=iextra-nx          
          end if
          do j=2_ni,ny-1_ni
             res(k,j,i) = dat(k,j,i)
             res(k,2*ny-j,i)=dat(k,j,iextra)
          end do
          res(k,1_ni,i)=dat(k,1_ni,i)
          res(k,ny,i)=dat(k,ny,i) 
       end do
    end do
  end subroutine
  !
  ! Smooth in 2D. Applies the same filter for all (x,y)-slices in a 3D (Z/t,y,x) field
  !
  ! The routine is taken from NCL 6.1.2, where it is called [d]filter2d and used in wrf_smooth_2d
  subroutine smooth_xy(res,nx,ny,nz,dat,niter)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni), intent(in) :: niter
    integer(kind=ni) :: i,j,k, n, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !
    real(kind=nr) :: tmp(nz,ny,nx) ! temporary helper array to avoid changing dat
    ! -----------------------------------------------------------------
    !
    res(:,:,:) = dat(:,:,:)
    !
    do n = 1_ni,niter
       tmp(:,:,:) = res(:,:,:)
       do i = 1_ni,nx
          do j = 2_ni,ny-1_ni
             do k = 1_ni,nz
                res(k,j,i) = res(k,j,i) + smooth_coef * (tmp(k,j-1_ni,i)-2_ni*tmp(k,j,i)+tmp(k,j+1_ni,i))
             end do
          end do
       end do
       do i = 2_ni,nx-1_ni
          do j = 1_ni,ny
             do k = 1_ni,nz
                res(k,j,i) = res(k,j,i) + smooth_coef * (tmp(k,j,i-1_ni)-2_ni*tmp(k,j,i)+tmp(k,j,i+1_ni))
             end do
          end do
       end do
       if ( grid_cyclic_ew ) then
          do j = 1_ni,ny
             do k = 1_ni,nz
                res(k,j,1_ni) = res(k,j,1_ni) + smooth_coef * (tmp(k,j,nx)-2_ni*tmp(k,j,1_ni)+tmp(k,j,2_ni))
                res(k,j,nx) = res(k,j,nx) + smooth_coef * (tmp(k,j,nx-1_ni)-2_ni*tmp(k,j,nx)+tmp(k,j,1_ni))
             end do
          end do
       end if
    end do
    !
    return
  end subroutine
  !
  ! Filter in 2D. Applies the same filter for all (x,y)-slices in a 3D (Z/t,y,x) field
  ! `filtr` must be a odd-length 1D array with entries normalised to 1
  subroutine filter_xy(res,nx,ny,nz,nf,dat,filtr)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), filtr(nf)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz, nf
    integer(kind=ni) :: i,j,k, n,m
    !f2py depend(nx,ny,nz) res
    !
    real(kind=nr) :: tmp(nz,ny,nx) ! temporary helper array to avoid changing dat
    ! -----------------------------------------------------------------
    !
    n = (nf - 1_ni)/2_ni
    !
    !res(:,:,:) = dat(:,:,:)
    !
    do i = 1_ni,nx
       do j = 1_ni+n,ny-n
          do k = 1_ni,nz
             res(k,j,i) = sum(filtr(:) * dat(k,j-n:j+n,i))
          end do
       end do
    end do
    !
    tmp(:,:,:) = res(:,:,:)
    do i = 1_ni+n,nx-n
       do j = 1_ni,ny
          do k = 1_ni,nz
             res(k,j,i) = sum(filtr(:) * tmp(k,j,i-n:i+n))
          end do
       end do
    end do
    if ( grid_cyclic_ew ) then
       do j = 1_ni,ny
          do k = 1_ni,nz
             do i = 1_ni,n
                m = n-i+1_ni
                res(k,j,i) = sum(filtr(1_ni:m)*tmp(k,j,nx-m+1_ni:nx)) &
                      &    + sum(filtr(m+1_ni:)*dat(k,j,:nf-m))
             end do
             do i = nx-n,nx
                m = nf - (i-nx+n + 1_ni)
                res(k,j,i) = sum(filtr(1_ni:m)*tmp(k,j,nx-m+1_ni:nx)) &
                      &    + sum(filtr(m+1_ni:)*tmp(k,j,:nf-m))
             end do
          end do
       end do
    end if
    !
    return
  end subroutine
  !
  ! Scaling and offsetting data
  ! [frequently used on netCDF input, see attributes "add_offset" and "scale"]
  subroutine scaleoff(res,nx,ny,nz,dat,scale,offset)
    integer(kind=2), intent(in) :: dat(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr), intent( in) :: scale, offset
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    ! -----------------------------------------------------------------
    !
    do i=1_ni,nx
       do j=1_ni,ny
          do k=1_ni,nz
             res(k,j,i) = dat(k,j,i)*scale + offset
          end do
       end do
    end do
  end subroutine
  !
  ! Take line data and fuzzily mark line locations on the map
  subroutine mask_lines(res,nx,ny,nz,np,nl,line,lineoff)
    real(kind=nr), intent(in) :: line(nz,np,3_ni)
    integer(kind=ni), intent(in) :: lineoff(nz,nl)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,np,nl
    !f2py depend(nz) res, lineoff
    !
    integer(kind=ni) :: i,j,k, n
    ! -----------------------------------------------------------------
    !
    res(:,:,:) = 0.0_nr
    !
    do k = 1_ni,nz
       do n = 1_ni,maxval(lineoff(k,:))
          i = nint(line(k,n,1_ni), ni)
          if (i > nx) i = i - nx
          j = nint(line(k,n,2_ni), ni)
          res(k,j,i) = 1.0_nr
       end do
    end do
  end subroutine
  !
  ! Regression 
  subroutine regress(res,ts,nx,nxp,nyp,ny,nz,nn,dat,datproj,pat,mask)
    real(kind=nr),intent(in) :: dat(nz,ny,nx), datproj(nz,nyp,nxp), pat(nn,nyp,nxp)
    real(kind=nr),intent(out) :: res(nn,ny,nx), ts(nn,nz)
    logical, intent(in) :: mask(nz)
    integer(kind=ni), intent(in) :: nx,nxp,ny,nyp,nz,nn
    !f2py depend(nyp,nxp) pat
    !f2py depend(nn,nz) ts
    !f2py depend(nn,ny,nx) res
    !f2py depend(nz) mask
    !
    integer(kind=ni) :: k,n
    ! -----------------------------------------------------------------
    !
    do k = 1_ni,nz
       if ( mask(k) ) then
          do n = 1_ni,nn
             ts(n,k) = sum(datproj(k,:,:)*pat(n,:,:))
             res(n,:,:) = res(n,:,:) + ts(n,k)*dat(k,:,:)
          end do
       end if
    end do
    !
  end subroutine
  !
  ! Make a 2d-histogram based on given linear ranges
  subroutine hist2d(res, incnt,outcnt, nx,ny,nz, dat1,scale1,ns1, dat2,scale2,ns2)
    real(kind=nr), intent(in) :: dat1(nz,ny,nx), dat2(nz,ny,nx), &
            &                    scale1(ns1), scale2(ns2)
    integer(kind=ni), intent(in) :: nx,ny,nz, ns1, ns2
    integer(kind=ni), intent(out) :: res(ns2,ns1), incnt, outcnt
    !f2py depend(nx,ny,nz) dat2
    !f2py depend(ns1,ns2) res
    !
    real(kind=nr) :: lower1, lower2, step1, step2
    integer(kind=ni) :: i,j,k, idx1, idx2
    ! -----------------------------------------------------------------
    !
    res(:,:) = 0_ni
    !
    lower1 = scale1(1_ni)
    lower2 = scale2(1_ni)
    step1  = scale1(2_ni) - scale1(1_ni)
    step2  = scale2(2_ni) - scale2(1_ni)
    !
    do k = 1_ni,nz
       do j = 1_ni,ny
          do i = 1_ni,nx
             idx1 = int((dat1(k,j,i)-lower1)/step1,ni) + 1_ni
             idx2 = int((dat2(k,j,i)-lower2)/step2,ni) + 1_ni
             !
             if ( idx1 < 1_ni .or. idx1 > ns1 .or. idx2 < 1_ni .or. idx2 > ns2 ) then
                outcnt = outcnt + 1_ni
             else 
                incnt = incnt + 1_ni
                res(idx2,idx1) = res(idx2,idx1) + 1_ni
             end if
          end do
       end do
    end do
    !
    return
  end subroutine
  !
  ! Mask areas above given threshold and conntected to list of points
  subroutine mask_minimum_connect(res, nx,ny, ns, dat, seeds,thres)
    real(kind=nr), intent(in) :: dat(ny,nx), thres
    integer(kind=ni), intent(in) :: seeds(ns, 2_ni), nx,ny,ns
    logical, intent(out) :: res(ny,nx)
    !f2py depend(nx,ny) res
    !
    integer(kind=ni) :: n, i,j
    ! -----------------------------------------------------------------
    !
    res(:,:) = .false.
    do n = 1_ni,ns
       j = seeds(n,1_ni)
       i = seeds(n,2_ni)
       if ( dat(j,i) < thres ) then
          write(*,'(A18,I4,A24)') 'FATAL Error: seed ', n, ' itself below threshold.'
          stop 1
       end if
       if ( .not. res(j,i) ) then
          res(j,i) = .true.
          call minimum_connect(res, nx,ny, dat, i,j, thres)
       end if
    end do
  end subroutine
  !
  ! Mask areas above given threshold and conntected to list of points
  recursive subroutine minimum_connect(res, nx,ny, dat, i,j, thres)
    real(kind=nr), intent(in) :: dat(ny,nx), thres
    integer(kind=ni), intent(in) :: nx,ny, i,j
    logical, intent(inout) :: res(ny,nx)
    !f2py depend(nx,ny) res
    ! -----------------------------------------------------------------
    !
    ! North
    if ( j > 1_ni ) then
       if ( dat(j-1_ni,i) >= thres .and. .not. res(j-1_ni,i) ) then
          res(j-1_ni,i) = .true.
          call minimum_connect(res, nx,ny, dat, i,j-1_ni, thres)
       end if
    end if
    ! South
    if ( j < ny ) then
       if ( dat(j+1_ni,i) >= thres .and. .not. res(j+1_ni,i) ) then
          res(j+1_ni,i) = .true.
          call minimum_connect(res, nx,ny, dat, i,j+1_ni, thres)
       end if
    end if
    ! West
    if ( i > 1_ni ) then
       if ( dat(j,i-1_ni) >= thres .and. .not. res(j,i-1_ni) ) then
          res(j,i-1_ni) = .true.
          call minimum_connect(res, nx,ny, dat, i-1_ni,j, thres)
       end if
    end if
    ! East
    if ( i < nx ) then
       if ( dat(j,i+1_ni) >= thres .and. .not. res(j,i+1_ni) ) then
          res(j,i+1_ni) = .true.
          call minimum_connect(res, nx,ny, dat, i+1_ni,j, thres)
       end if
    end if
    ! Periodic boundary West
    if ( grid_cyclic_ew .and. i == 1_ni ) then
       if ( dat(j,nx) >= thres .and. .not. res(j,nx) ) then
          res(j,nx) = .true.
          call minimum_connect(res, nx,ny, dat, nx,j, thres)
       end if
    end if
    ! Periodic boundary East
    if ( grid_cyclic_ew .and. i == nx ) then
       if ( dat(j,1_ni) >= thres .and. .not. res(j,1_ni) ) then
          res(j,1_ni) = .true.
          call minimum_connect(res, nx,ny, dat, 1_ni,j, thres)
       end if
    end if
  end subroutine
  !
  ! Mask areas above given threshold and conntected to list of points
  subroutine mask_grow(res, nx,ny,nz, mask, nn)
    logical, intent(in) :: mask(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz, nn
    logical, intent(out) :: res(nz,ny,nx)
    !f2py depend(nx,ny,nz) res
    !
    logical :: tmp(nz,ny,nx)
    integer(kind=ni) :: n, i,j,k
    ! -----------------------------------------------------------------
    !
    ! Temporary copy to avoid overwriting mask
    tmp(:,:,:) = mask(:,:,:)
    res(:,:,:) = mask(:,:,:)
    !
    ! Grow the mask, one layer at a time
    do n = 1_ni,nn
       do k = 1_ni,nz
          do j = 1_ni,ny
             do i = 1_ni,nx
                if ( .not. res(k,j,i) ) then
                   res(k,j,i) = tmp(k,min(j+1_ni,ny),i) .or. tmp(k,max(j-1_ni,1_ni),i) .or. &
                              & tmp(k,j,min(i+1_ni,nx)) .or. tmp(k,j,max(i-1_ni,1_ni))
                end if
             end do
             ! Respect periodic boundaries in West and East
             if ( grid_cyclic_ew ) then
                res(k,j,1_ni) = res(k,j,1_ni) .or. tmp(k,j,nx)
                res(k,j,nx) = res(k,j,nx) .or. tmp(k,j,1_ni)
             end if
          end do
       end do
       tmp(:,:,:) = res(:,:,:)
    end do
  end subroutine
end module
