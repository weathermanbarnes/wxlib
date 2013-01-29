! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- partial derivatives
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module derivatives
  use kind
  use config
  !
  implicit none
contains
  !
  ! Calculates partial x derivative: b = partial(a)/partial(x)
  !  Returns 0 on first and last lon for non-cyclic grid
  subroutine ddx(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(i = 2_ni:nx-1_ni, j = 1_ni:ny, k = 1_ni:nz)
       res(k,j,i) = (dat(k,j,i+1_ni)-dat(k,j,i-1_ni))/dx(j,i)
    end forall
    if (grid_cyclic_ew) then
       forall(j = 1_ni:ny, k = 1_ni:nz)
          res(k,j,1_ni) = (dat(k,j,  2_ni)-dat(k,j     ,nx))/dx(j,1_ni)
          res(k,j,nx  ) = (dat(k,j  ,1_ni)-dat(k,j,nx-1_ni))/dx(j,nx)
       end forall
    else 
       forall(j = 1_ni:ny, k = 1_ni:nz)
          res(k,j,1_ni) = 0._nr
          res(k,j,nx  ) = 0._nr
       end forall
    end if
  end subroutine
  !
  ! Calculates partial y derivative: b = partial(a)/partial(y)
  !  Returns 0 on first and last lat
  subroutine ddy(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(i = 1_ni:nx, j = 2_ni:ny-1_ni, k = 1_ni:nz)
       res(k,j,i) = (dat(k,j+1_ni,i)-dat(k,j-1_ni,i))/dy(j,i)
    end forall
    forall(i = 1_ni:nx, k = 1_ni:nz)
       res(k,1_ni,i)=0._nr
       res(k,ny,i)=0._nr
    end forall
  end subroutine
  !
  ! Calculates gradient [bx,by] = grad(a)
  !  Returns 0 on edges
  subroutine grad(resx,resy,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: resx(nz,ny,nx), resy(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) resx, resy
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(resx,nx,ny,nz,dat,dx,dy)
    call ddy(resy,nx,ny,nz,dat,dx,dy)
    ! set to 0 where grad result not valid: j=1 and ny
    resx(:,1_ni,:) = 0._nr
    resy(:,1_ni,:) = 0._nr
    resx(:,ny,:  ) = 0._nr
    resy(:,ny,:  ) = 0._nr
    if (.not.(grid_cyclic_ew)) then
       ! set to 0 where grad result not valid: i=1 and nx
       resx(:,:,1_ni) = 0._nr 
       resy(:,:,1_ni) = 0._nr
       resx(:,:,nx  ) = 0._nr
       resy(:,:,nx  ) = 0._nr
    end if
  end subroutine
  !
  ! Calculates  2-D laplacian lap2(a)
  !  returns 0 on edges
  subroutine lap2(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    forall(i = 2_ni:nx-1_ni, j = 2_ni:ny-1_ni, k = 1_ni:nz)
       res(k,j,i) = (dat(k,j,i+1_ni)+dat(k,j,i-1_ni)-2_nr*dat(k,j,i))/dx(j,i)**2_nr + &
                    (dat(k,j+1_ni,i)+dat(k,j-1_ni,i)-2_nr*dat(k,j,i))/dy(j,i)**2_nr
    end forall
    if (grid_cyclic_ew) then
       forall(j = 2_ni:ny-1_ni, k = 1_ni:nz)
          res(k,j,1_ni) = (dat(k,j,2_ni)+dat(k,j,nx)-2_nr*dat(k,j,1_ni))/dx(j,1_ni)**2_nr + &
              (dat(k,j+1_ni,1_ni)+dat(k,j-1_ni,1_ni)-2_nr*dat(k,j,1_ni))/dy(j,1_ni)**2_nr
          res(k,j,nx) = (dat(k,j,1_ni)+dat(k,j,nx-1_ni)-2_nr*dat(k,j,nx))/dx(j,nx)**2_nr + &
                     (dat(k,j+1_ni,nx)+dat(k,j-1_ni,nx)-2_nr*dat(k,j,nx))/dy(j,nx)**2_nr
       end forall
    else
       forall(j = 2_ni:ny-1_ni, k = 1_ni:nz)
          res(k,j,1_ni) = 0._nr
          res(k,j,nx) = 0._nr
       end forall
    end if
    forall(i = 1_ni:nx,k = 1_ni:nz)
       res(k,1_ni,i) = 0._nr
       res(k,ny  ,i) = 0._nr
    end forall
  end subroutine
  !
end module
