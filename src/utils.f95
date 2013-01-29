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
    !real(kind=nr), intent(out) :: res(nz,2*ny-2_ni,nx)
    real(kind=nr), intent(out) :: res(nz,2*ny-2,nx)
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
end module
