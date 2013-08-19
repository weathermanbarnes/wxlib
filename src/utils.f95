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
    ! TODO: Make accessible via config
    real(kind=nr) :: tmp(nz,ny,nx) ! temporary helper array to avoid changing dat
    real(kind=nr), parameter :: coef = 0.25_nr 
    ! -----------------------------------------------------------------
    !
    res(:,:,:) = dat(:,:,:)
    !
    do n = 1_ni,niter
       tmp(:,:,:) = res(:,:,:)
       do i = 1_ni,nx
          do j = 2_ni,ny-1_ni
             do k = 1_ni,nz
                res(k,j,i) = res(k,j,i) + coef * (tmp(k,j-1_ni,i)-2_ni*tmp(k,j,i)+tmp(k,j+1_ni,i))
             end do
          end do
       end do
       do i = 2_ni,nx-1_ni
          do j = 1_ni,ny
             do k = 1_ni,nz
                res(k,j,i) = res(k,j,i) + coef * (tmp(k,j,i-1_ni)-2_ni*tmp(k,j,i)+tmp(k,j,i+1_ni))
             end do
          end do
       end do
       if ( grid_cyclic_ew ) then
          do j = 1_ni,ny
             do k = 1_ni,nz
                res(k,j,1_ni) = res(k,j,1_ni) + coef * (tmp(k,j,nx)-2_ni*tmp(k,j,1_ni)+tmp(k,j,2_ni))
                res(k,j,nx) = res(k,j,nx) + coef * (tmp(k,j,nx-1_ni)-2_ni*tmp(k,j,nx)+tmp(k,j,1_ni))
             end do
          end do
       end if
    end do
    !
    return
  end subroutine
  !
end module
