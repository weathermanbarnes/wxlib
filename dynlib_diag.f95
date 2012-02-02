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
    do i=2_ni,nx-1_ni
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,i) = (v(k,j,i+1_ni)-v(k,j,i-1_ni))/(2.0_nr*dx(j,i)) &
                  &     - (u(k,j+1_ni,i)-u(k,j-1_ni,i))/(2.0_nr*dy(j,i))
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,1_ni) = (v(k,j,  2_ni)-u(k,j     ,nx))/(2.0_nr*dx(j,1_ni)) &
                     &     - (u(k,j+1_ni,i)-v(k,j-1_ni, i))/(2.0_nr*dy(j,1_ni))
             res(k,j,nx  ) = (v(k,j  ,1_ni)-u(k,j,nx-1_ni))/(2.0_nr*dx(j,nx)) &
                     &     - (u(k,j+1_ni,i)-v(k,j-1_ni, i))/(2.0_nr*dy(j,nx))
          end do
       end do
    end if
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
    do i=2_ni,nx-1_ni
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,i) = (u(k,j,i+1_ni)-u(k,j,i-1_ni))/(2.0_nr*dx(j,i)) &
                  &     + (v(k,j+1_ni,i)-v(k,j-1_ni,i))/(2.0_nr*dy(j,i))
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,1_ni) = (u(k,j,  2_ni)-u(k,j     ,nx))/(2.0_nr*dx(j,1_ni)) &
                     &     + (v(k,j+1_ni,i)-v(k,j-1_ni, i))/(2.0_nr*dy(j,1_ni))
             res(k,j,nx  ) = (u(k,j  ,1_ni)-u(k,j,nx-1_ni))/(2.0_nr*dx(j,nx)) &
                     &     + (v(k,j+1_ni,i)-v(k,j-1_ni, i))/(2.0_nr*dy(j,nx))
          end do
       end do
    end if
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
    do i=2_ni,nx-1_ni
        do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,i) = (u(k,j+1_ni,i)-u(k,j-1_ni,i))/(2.0_nr*dy(j,i)) &
                  &     + (v(k,j,i+1_ni)-v(k,j,i-1_ni))/(2.0_nr*dx(j,i))
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,1_ni) = (u(k,j+1_ni,i)-u(k,j-1_ni, i))/(2.0_nr*dy(j,1_ni)) &
                     &     + (v(k,j,  2_ni)-v(k,j     ,nx))/(2.0_nr*dx(j,1_ni))
             res(k,j,nx  ) = (u(k,j+1_ni,i)-u(k,j-1_ni, i))/(2.0_nr*dy(j,nx)) &
                     &     + (v(k,j  ,1_ni)-v(k,j,nx-1_ni))/(2.0_nr*dx(j,nx))
          end do
       end do
    end if
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
    do i=2_ni,nx-1_ni
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,i) = (u(k,j,i+1_ni)-u(k,j,i-1_ni))/(2.0_nr*dx(j,i)) &
                  &     - (v(k,j+1_ni,i)-v(k,j-1_ni,i))/(2.0_nr*dy(j,i))
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,1_ni) = (u(k,j,  2_ni)-u(k,j     ,nx))/(2.0_nr*dx(j,1_ni)) &
                     &     - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/(2.0_nr*dy(j,1_ni))
             res(k,j,nx  ) = (u(k,j  ,1_ni)-u(k,j,nx-1_ni))/(2.0_nr*dx(j,nx)) &
                     &     - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/(2.0_nr*dy(j,nx))
          end do
       end do
    end if
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
    do i=2_ni,nx-1_ni
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,i) = ((u(k,j+1_ni,i)-u(k,j-1_ni,i))/(2.0_nr*dy(j,i))         &
                  &     +  (v(k,j,i+1_ni)-v(k,j,i-1_ni))/(2.0_nr*dx(j,i)))**2._nr &
                  &     + ((u(k,j,i+1_ni)-u(k,j,i-1_ni))/(2.0_nr*dx(j,i))         &
                  &     -  (v(k,j+1_ni,i)-v(k,j-1_ni,i))/(2.0_nr*dy(j,i)))**2._nr
             res(k,j,i) = sqrt(res(k,j,i))
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,1_ni) = ((u(k,j+1_ni,i)-u(k,j-1_ni, i))/(2.0_nr*dy(j,1_ni))         &
                     &     +  (v(k,j,  2_ni)-v(k,j     ,nx))/(2.0_nr*dx(j,1_ni)))**2._nr &
                     &     + ((u(k,j,  2_ni)-u(k,j     ,nx))/(2.0_nr*dx(j,1_ni))         &
                     &     -  (v(k,j+1_ni,i)-v(k,j-1_ni, i))/(2.0_nr*dy(j,1_ni)))**2._nr
             res(k,j,nx  ) = ((u(k,j+1_ni,i)-u(k,j-1_ni, i))/(2.0_nr*dy(j,nx))         &
                     &     +  (v(k,j  ,1_ni)-v(k,j,nx-1_ni))/(2.0_nr*dx(j,nx)))**2._nr &
                     &     + ((u(k,j  ,1_ni)-u(k,j,nx-1_ni))/(2.0_nr*dx(j,nx))         &
                     &     -  (v(k,j+1_ni,i)-v(k,j-1_ni, i))/(2.0_nr*dy(j,nx)))**2._nr
             res(k,j,1_ni) = sqrt(res(k,j,1_ni))
             res(k,j,  nx) = sqrt(res(k,j,  nx))
          end do
       end do
    end if
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
    do i=2_ni,nx-1_ni
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             def_shear   = (u(k,j+1_ni,i)-u(k,j-1_ni,i))/(2.0_nr*dy(j,i)) &
                  &      + (v(k,j,i+1_ni)-v(k,j,i-1_ni))/(2.0_nr*dx(j,i))
             def_stretch = (u(k,j,i+1_ni)-u(k,j,i-1_ni))/(2.0_nr*dx(j,i)) &
                  &      - (v(k,j+1_ni,i)-v(k,j-1_ni,i))/(2.0_nr*dy(j,i))
             res(k,j,i)  = 0.5_nr*atan2(def_shear, def_stretch)
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             def_shear   = (u(k,j+1_ni,i)-u(k,j-1_ni, i))/(2.0_nr*dy(j,1_ni)) &
                     &   + (v(k,j,  2_ni)-v(k,j     ,nx))/(2.0_nr*dx(j,1_ni))
             def_stretch = (u(k,j,  2_ni)-u(k,j     ,nx))/(2.0_nr*dx(j,1_ni)) &
                     &   - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/(2.0_nr*dy(j,1_ni))
             res(k,j,1_ni) = 0.5_nr*atan2(def_shear, def_stretch)
             def_shear   = (u(k,j+1_ni,i)-u(k,j-1_ni, i))/(2.0_nr*dy(j,nx)) &
                     &   + (v(k,j  ,1_ni)-v(k,j,nx-1_ni))/(2.0_nr*dx(j,nx))
             def_stretch = (u(k,j  ,1_ni)-u(k,j,nx-1_ni))/(2.0_nr*dx(j,nx)) &
                     &   - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/(2.0_nr*dy(j,nx))
             res(k,j,nx) = 0.5_nr*atan2(def_shear, def_stretch)
          end do
       end do
    end if
  end subroutine
end module
