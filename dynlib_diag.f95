! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- diagnostical functions
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module diag
  use kind
  use config
  !
  implicit none
contains
  ! Calculates vorticity b = vor(a)
  subroutine vor(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    do k=1_ni,nz
       do j=2_ni,ny-1_ni
          do i=2_ni,nx-1_ni
             res(k,j,i) = (v(k,j,i+1_ni)-v(k,j,i-1_ni))/dx(j,i) &
                  &     - (u(k,j+1_ni,i)-u(k,j-1_ni,i))/dy(j,i)
          end do
          if (grid_cyclic_ew) then
             res(k,j,1_ni) = (v(k,j,  2_ni)-u(k,j     ,nx))/dx(j,i) &
                     &     - (u(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i)
             res(k,j,nx  ) = (v(k,j  ,1_ni)-u(k,j,nx-1_ni))/dx(j,i) &
                     &     - (u(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i)
          end if
       end do
    end do
  end subroutine
  !
  ! Calculates divergence b = div(a)
  subroutine div(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    do k=1_ni,nz
       do j=2_ni,ny-1_ni
          do i=2_ni,nx-1_ni
             res(k,j,i) = (u(k,j,i+1_ni)-u(k,j,i-1_ni))/dx(j,i) &
                  &     + (v(k,j+1_ni,i)-v(k,j-1_ni,i))/dy(j,i)
          end do
          if (grid_cyclic_ew) then
             res(k,j,1_ni) = (u(k,j,  2_ni)-u(k,j     ,nx))/dx(j,i) &
                     &     + (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i)
             res(k,j,nx  ) = (u(k,j  ,1_ni)-u(k,j,nx-1_ni))/dx(j,i) &
                     &     + (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i)
          end if
       end do
    end do
  end subroutine
  !
  ! Calculates shear deformation b = def_shear(a)
  subroutine def_shear(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    do k=1_ni,nz
       do j=2_ni,ny-1_ni
          do i=2_ni,nx-1_ni
             res(k,j,i) = (u(k,j+1_ni,i)-u(k,j-1_ni,i))/dy(j,i) &
                  &     + (v(k,j,i+1_ni)-v(k,j,i-1_ni))/dx(j,i)
          end do
          if (grid_cyclic_ew) then
             res(k,j,1_ni) = (u(k,j+1_ni,i)-u(k,j-1_ni, i))/dy(j,i) &
                     &     + (v(k,j,  2_ni)-v(k,j     ,nx))/dx(j,i)
             res(k,j,nx  ) = (u(k,j+1_ni,i)-u(k,j-1_ni, i))/dy(j,i) &
                     &     + (v(k,j  ,1_ni)-v(k,j,nx-1_ni))/dx(j,i)
          end if
       end do
    end do
  end subroutine
  !
  ! Calculates shear deformation b = def_stretch(a)
  subroutine def_stretch(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    do k=1_ni,nz
       do j=2_ni,ny-1_ni
          do i=2_ni,nx-1_ni
             res(k,j,i) = (u(k,j,i+1_ni)-u(k,j,i-1_ni))/dx(j,i) &
                  &     - (v(k,j+1_ni,i)-v(k,j-1_ni,i))/dy(j,i)
          end do
          if (grid_cyclic_ew) then
             res(k,j,1_ni) = (u(k,j,  2_ni)-u(k,j     ,nx))/dx(j,i) &
                     &     - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i)
             res(k,j,nx  ) = (u(k,j  ,1_ni)-u(k,j,nx-1_ni))/dx(j,i) &
                     &     - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i)
          end if
       end do
    end do
  end subroutine
  !
  ! Calculates total [and rotation invariant] deformation b = def(a)
  subroutine def_total(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    do k=1_ni,nz
       do j=2_ni,ny-1_ni
          do i=2_ni,nx-1_ni
             res(k,j,i) = ((u(k,j+1_ni,i)-u(k,j-1_ni,i))/dy(j,i)         &
                  &     +  (v(k,j,i+1_ni)-v(k,j,i-1_ni))/dx(j,i))**2._nr &
                  &     + ((u(k,j,i+1_ni)-u(k,j,i-1_ni))/dx(j,i)         &
                  &     -  (v(k,j+1_ni,i)-v(k,j-1_ni,i))/dy(j,i))**2._nr
             res(k,j,i) = sqrt(res(k,j,i))
          end do
          if (grid_cyclic_ew) then
             res(k,j,1_ni) = ((u(k,j+1_ni,i)-u(k,j-1_ni, i))/dy(j,i)         &
                     &     +  (v(k,j,  2_ni)-v(k,j     ,nx))/dx(j,i))**2._nr &
                     &     + ((u(k,j,  2_ni)-u(k,j     ,nx))/dx(j,i)         &
                     &     -  (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i))**2._nr
             res(k,j,nx  ) = ((u(k,j+1_ni,i)-u(k,j-1_ni, i))/dy(j,i)         &
                     &     +  (v(k,j  ,1_ni)-v(k,j,nx-1_ni))/dx(j,i))**2._nr &
                     &     + ((u(k,j  ,1_ni)-u(k,j,nx-1_ni))/dx(j,i)         &
                     &     -  (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i))**2._nr
             res(k,j,1_ni) = sqrt(res(k,j,1_ni))
             res(k,j,  nx) = sqrt(res(k,j,  nx))
          end if
       end do
    end do
  end subroutine
  !
  ! Calculates rotation angle  b = def_angle(a)  to achieve pure shear rotation
  subroutine def_angle(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: def_shear, def_stretch
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    do k=1_ni,nz
       do j=2_ni,ny-1_ni
          do i=2_ni,nx-1_ni
             def_shear   = (u(k,j+1_ni,i)-u(k,j-1_ni,i))/dy(j,i) &
                  &      + (v(k,j,i+1_ni)-v(k,j,i-1_ni))/dx(j,i)
             def_stretch = (u(k,j,i+1_ni)-u(k,j,i-1_ni))/dx(j,i) &
                  &      - (v(k,j+1_ni,i)-v(k,j-1_ni,i))/dy(j,i)
             res(k,j,i)  = 0.5_nr*atan2(def_shear, def_stretch)
          end do
          if (grid_cyclic_ew) then
             def_shear   = (u(k,j+1_ni,i)-u(k,j-1_ni, i))/dy(j,i) &
                     &   + (v(k,j,  2_ni)-v(k,j     ,nx))/dx(j,i)
             def_stretch = (u(k,j,  2_ni)-u(k,j     ,nx))/dx(j,i) &
                     &   - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i)
             res(k,j,1_ni) = 0.5_nr*atan2(def_shear, def_stretch)
             def_shear   = (u(k,j+1_ni,i)-u(k,j-1_ni, i))/dy(j,i) &
                     &   + (v(k,j  ,1_ni)-v(k,j,nx-1_ni))/dx(j,i)
             def_stretch = (u(k,j  ,1_ni)-u(k,j,nx-1_ni))/dx(j,i) &
                     &   - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i)
             res(k,j,nx) = 0.5_nr*atan2(def_shear, def_stretch)
          end if
       end do
    end do
  end subroutine
end module
