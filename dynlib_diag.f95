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
             res(k,j,i) = 0.5_nr*((v(k,j,i+1_ni)-v(k,j,i-1_ni))/dx(j,i) &
                  &     -         (u(k,j+1_ni,i)-u(k,j-1_ni,i))/dy(j,i))
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,1_ni) = 0.5_nr*((v(k,j,     2_ni)-u(k,j        ,nx))/dx(j,i) &
                     &     -         (u(k,j+1_ni,1_ni)-v(k,j-1_ni, 1_ni))/dy(j,i))
             res(k,j,nx  ) = 0.5_nr*((v(k,j,     1_ni)-u(k,j,   nx-1_ni))/dx(j,i) &
                     &     -         (u(k,j+1_ni,  nx)-v(k,j-1_ni,   nx))/dy(j,i))
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
             res(k,j,i) = 0.5_nr*((u(k,j,i+1_ni)-u(k,j,i-1_ni))/dx(j,i) &
                  &        +      (v(k,j+1_ni,i)-v(k,j-1_ni,i))/dy(j,i))
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,1_ni) = 0.5_nr*((u(k,j,  2_ni)-u(k,j     ,nx))/dx(j,i) &
                     &             + (v(k,j+1_ni,1_ni)-v(k,j-1_ni, 1_ni))/dy(j,i))
             res(k,j,nx  ) = 0.5_nr*((u(k,j  ,1_ni)-u(k,j,nx-1_ni))/dx(j,i) &
                     &             + (v(k,j+1_ni,nx)-v(k,j-1_ni, nx))/dy(j,i))
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
             res(k,j,i) = 0.5_nr*((u(k,j+1_ni,i)-u(k,j-1_ni,i))/dy(j,i) &
                  &             + (v(k,j,i+1_ni)-v(k,j,i-1_ni))/dx(j,i))
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,1_ni) = 0.5_nr*((u(k,j+1_ni,1_ni)-u(k,j-1_ni,1_ni))/dy(j,i) &
                     &            + (v(k,j,  2_ni)-v(k,j     ,nx))/dx(j,i))
             res(k,j,nx  ) = 0.5_nr*((u(k,j+1_ni,nx)-u(k,j-1_ni, nx))/dy(j,i) &
                     &            + (v(k,j,  1_ni)-v(k,j,nx-1_ni))/dx(j,i))
          end do
       end do
    end if
  end subroutine
  !
  ! Calculates stretch deformation b = def_stretch(a)
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
             res(k,j,i) = 0.5_nr*((u(k,j,i+1_ni)-u(k,j,i-1_ni))/dx(j,i) &
                  &             - (v(k,j+1_ni,i)-v(k,j-1_ni,i))/dy(j,i))
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,1_ni) = 0.5_nr*((u(k,j,  2_ni)-u(k,j     ,nx))/dx(j,i) &
                     &             - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i))
             res(k,j,nx  ) = 0.5_nr*((u(k,j  ,1_ni)-u(k,j,nx-1_ni))/dx(j,i) &
                     &             - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i))
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
             res(k,j,i) = ((u(k,j+1_ni,i)-u(k,j-1_ni,i))/dy(j,i)         &
                  &     +  (v(k,j,i+1_ni)-v(k,j,i-1_ni))/dx(j,i))**2._nr &
                  &     + ((u(k,j,i+1_ni)-u(k,j,i-1_ni))/dx(j,i)         &
                  &     -  (v(k,j+1_ni,i)-v(k,j-1_ni,i))/dy(j,i))**2._nr
             res(k,j,i) = 0.5_nr*sqrt(res(k,j,i))
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             res(k,j,1_ni) = ((u(k,j+1_ni,i)-u(k,j-1_ni, i))/dy(j,i)         &
                     &     +  (v(k,j,  2_ni)-v(k,j     ,nx))/dx(j,i))**2._nr &
                     &     + ((u(k,j,  2_ni)-u(k,j     ,nx))/dx(j,i)         &
                     &     -  (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i))**2._nr
             res(k,j,nx  ) = ((u(k,j+1_ni,i)-u(k,j-1_ni, i))/dy(j,i)         &
                     &     +  (v(k,j  ,1_ni)-v(k,j,nx-1_ni))/dx(j,i))**2._nr &
                     &     + ((u(k,j  ,1_ni)-u(k,j,nx-1_ni))/dx(j,i)         &
                     &     -  (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i))**2._nr
             res(k,j,1_ni) = 0.5_nr*sqrt(res(k,j,1_ni))
             res(k,j,  nx) = 0.5_nr*sqrt(res(k,j,  nx))
          end do
       end do
    end if
  end subroutine
  !
  ! Calculates rotation angle  b = def_angle(a)  to achieve pure shear rotation
  ! This is the angle of axis of dilatation psi_d (called alpha in Markowski Richardson)
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
             def_shear   = 0.5_nr*((u(k,j+1_ni,i)-u(k,j-1_ni,i))/dy(j,i) &
                  &              + (v(k,j,i+1_ni)-v(k,j,i-1_ni))/dx(j,i))
             def_stretch = 0.5_nr*((u(k,j,i+1_ni)-u(k,j,i-1_ni))/dx(j,i) &
                  &              - (v(k,j+1_ni,i)-v(k,j-1_ni,i))/dy(j,i))
             res(k,j,i)  = 0.5_nr*atan2(def_shear, def_stretch)
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             def_shear   = 0.5_nr*((u(k,j+1_ni,1)-u(k,j-1_ni,1))/dy(j,i) &
                     &          +  (v(k,j,  2_ni)-v(k,j    ,nx))/dx(j,i))
             def_stretch = 0.5_nr*((u(k,j,  2_ni)-u(k,j     ,nx))/dx(j,i) &
                     &           - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i))
             res(k,j,1_ni) = 0.5_nr*atan2(def_shear, def_stretch)
             def_shear   = 0.5_nr*((u(k,j+1_ni,i)-u(k,j-1_ni, i))/dy(j,i) &
                     &           + (v(k,j  ,1_ni)-v(k,j,nx-1_ni))/dx(j,i))
             def_stretch = 0.5_nr*((u(k,j  ,1_ni)-u(k,j,nx-1_ni))/dx(j,i) &
                     &           - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,i))
             res(k,j,nx) = 0.5_nr*atan2(def_shear, def_stretch)
          end do
       end do
    end if
  end subroutine
  !
  ! Calculates cos(2*beta)  b = cos2beta(a), cos double the angle between dilatation axis (u,v) and iso-lines of pv
  subroutine cos2beta(res,nx,ny,nz,u,v,pv,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), pv(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: def_shear, def_stretch, alpha, beta, gradx, grady, pi=3.141592 !pi should be a constant...
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res, v, pv
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    do i=2_ni,nx-1_ni
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             def_shear   = 0.5_nr*((u(k,j+1_ni,i)-u(k,j-1_ni,i))/dy(j,i) &
                  &              + (v(k,j,i+1_ni)-v(k,j,i-1_ni))/dx(j,i))
             def_stretch = 0.5_nr*((u(k,j,i+1_ni)-u(k,j,i-1_ni))/dx(j,i) &
                  &              - (v(k,j+1_ni,i)-v(k,j-1_ni,i))/dy(j,i))
             gradx = 0.5_nr*(pv(k,j,i+1_ni)-pv(k,j,i-1_ni))/dx(j,i)
             grady = 0.5_nr*(pv(k,j+1_ni,i)-pv(k,j-1_ni,i))/dy(j,i)
             alpha = 0.5_nr*atan2(def_shear, def_stretch)
             beta  = alpha-atan2(-grady,gradx)
             res(k,j,i) = cos(2.*beta)
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             def_shear   = 0.5_nr*(u(k,j+1_ni,1_ni)-u(k,j-1_ni,1_ni))/dy(j,i) &
                  &      + 0.5_nr*(v(k,j,  2_ni)-v(k,j,   nx))/dx(j,i)
             def_stretch = 0.5_nr*(u(k,j,  2_ni)-u(k,j,   nx))/dx(j,i) &
                  &      - 0.5_nr*(v(k,j+1_ni,1_ni)-v(k,j-1_ni,1_ni))/dy(j,i)
             gradx = 0.5_nr*(pv(k,j,2_ni)-pv(k,j,nx))/dx(j,i)
             grady = 0.5_nr*(pv(k,j+1_ni,1_ni)-pv(k,j-1_ni,1_ni))/dy(j,i)
             alpha = 0.5_nr*atan2(def_shear, def_stretch)
             beta  = alpha-atan2(-grady,gradx)
             res(k,j,1_ni) = cos(2.*beta)

             def_shear   = 0.5_nr*(u(k,j+1_ni,nx)-u(k,j-1_ni,nx))/dy(j,i) &
                  &      + 0.5_nr*(v(k,j,  1_ni)-v(k,j,nx-1_ni))/dx(j,i)
             def_stretch = 0.5_nr*(u(k,j,  1_ni)-u(k,j,nx-1_ni))/dx(j,i) &
                  &      - 0.5_nr*(v(k,j+1_ni,nx)-v(k,j-1_ni,nx))/dy(j,i)
             gradx = 0.5_nr*(pv(k,j,1_ni)-pv(k,j,nx-1_ni))/dx(j,i)
             grady = 0.5_nr*(pv(k,j+1_ni,nx)-pv(k,j-1_ni,nx))/dy(j,i)
             alpha = 0.5_nr*atan2(def_shear, def_stretch)
             beta  = alpha-atan2(-grady,gradx)
             res(k,j,nx) = cos(2.*beta)
          end do
       end do
    end if
  end subroutine
  !
  ! Calculates gradient [bx,by] = grad(a)
  subroutine grad(resx,resy,nx,ny,nz,pv,dx,dy)
    real(kind=nr), intent(in)  :: pv(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: resx(nz,ny,nx), resy(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) resx, resy
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    do i=2_ni,nx-1_ni
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             resx(k,j,i) = 0.5_nr*((pv(k,j,i+1_ni)-pv(k,j,i-1_ni))/dx(j,i))
             resy(k,j,i) = 0.5_nr*((pv(k,j+1_ni,i)-pv(k,j-1_ni,i))/dy(j,i))
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             resx(k,j,1_ni) = 0.5_nr*((pv(k,j,  2_ni)-pv(k,j     ,nx))/dx(j,i)) 
             resy(k,j,1_ni) = 0.5_nr*((pv(k,j+1_ni,i)-pv(k,j-1_ni, i))/dy(j,i))
             resx(k,j,nx  ) = 0.5_nr*((pv(k,j  ,1_ni)-pv(k,j,nx-1_ni))/dx(j,i))
             resy(k,j,nx  ) = 0.5_nr*((pv(k,j+1_ni,i)-pv(k,j-1_ni, i))/dy(j,i))
          end do
       end do
    end if
  end subroutine
!
! Calculates frontogenesis F = d/dt|grad(pv)|=0.5*|grad(pv)|(D*cos(2*beta)-delta)
  subroutine fgen(res,nx,ny,nz,u,v,pv,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), pv(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: beta, delta
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) v, pv, res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !

  end subroutine
!
! Binary gradient: Flags 0 if meridional (PV) gradient .GE. thres, or 1 if gradient .LT. thres
! NOTE: type returned is real, not integer

  subroutine bgra(res,nx,ny,nz,pv,dx,dy)
    real(kind=nr), intent(in)  :: pv(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: gra, thres
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    !do i=2_ni,nx-1_ni
                !
    thres=-5E-12  ! binary threshold for PV gradient
                !
    do i=1_ni,nx
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             gra = 0.5_nr*((pv(k,j+1_ni,i)-pv(k,j-1_ni,i))/dy(j,i)) 
             if (gra < thres) then
                res(k,j,i)=1.
             else 
                res(k,j,i)=0.
             end if     
          end do
       end do
    end do
!   (does not need the code below if calculating meridional gradient only)
!    if (grid_cyclic_ew) then
!       do j=2_ni,ny-1_ni
!          do k=1_ni,nz
!            resx(k,j,1_ni) = (pv(k,j,  2_ni)-pv(k,j     ,nx))/dx(j,i) 
!             resy(k,j,1_ni) = (pv(k,j+1_ni,i)-pv(k,j-1_ni, i))/dy(j,i)
!             resx(k,j,nx  ) = (pv(k,j  ,1_ni)-pv(k,j,nx-1_ni))/dx(j,i)
!             resy(k,j,nx  ) = (pv(k,j+1_ni,i)-pv(k,j-1_ni, i))/dy(j,i)
!          end do
!       end do
!    end if
  end subroutine

end module
