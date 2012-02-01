! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- data type conversions
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module conv
  use kind
  use config
  !
  implicit none
contains
  ! Scaling and offsetting data
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
end module
