! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- statistical functions
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module stat
  use kind
  use config
  !
  implicit none
contains
  ! Scaling and offsetting data
  subroutine basic(mean,stddev, nx,ny,nt,dat)
    real(kind=nr), intent(in) :: dat(nt,ny,nx)
    real(kind=nr), intent(out) :: mean(ny,nx), stddev(ny,nx)
    integer(kind=ni) :: i,j,n, nx,ny,nt
    !f2py depend(nx,ny) mean, stddev
    ! -----------------------------------------------------------------
    !
    do i=1_ni,nx
       do j=1_ni,ny
          mean  (j,i) = 0._nr
          stddev(j,i) = 0._nr
          do n=1_ni,nt
             mean  (j,i) = mean(j,i)   + dat(n,j,i)
             stddev(j,i) = stddev(j,i) + dat(n,j,i)**2._nr
          end do
          mean  (j,i) = mean(j,i)/nt
          stddev(j,i) = sqrt((stddev(j,i)-(2._nr*nt-1._nr)*mean(j,i)**2._nr)/(nt-1_ni))
       end do
    end do
  end subroutine
end module
