! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- diagnostical functions
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!@ Diagnostic functions for dynamic variables. 
!@ 
!@ This module contains diagnostic functions useful for analysing the wind
!@ field. For diagnostics relating to thermodynamic variables, refer to the
!@ :mod:`dynlib.humidity` module.
!@
!@ This module contains those diagnostics that depend on tendencies
module diag_tend
  use kind
  use config
  use consts
  use derivatives
  use tend
  !
  implicit none
contains
  !
  !@ Calculate the ratio between effective rotation and strain rate
  !@ 
  !@ This ratio is the main result of Lapeyre, Klein and Hua (1999). They denote it "r".
  !@ Their effective rotation takes into account both vorticity and rotation of the axis 
  !@ of dilatation.
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ mont : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Montgomery potential on an isentropic surface or geopotential on pressure surface.
  !@ fcor : np.ndarray with shape (ny,nx) and dtype float64
  !@     Latitude of the grid point.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated effective-rotation-to-strain ratio.
  subroutine rotation_strain_ratio(res,nx,ny,nz,u,v,mont,fcor,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx),v(nz,ny,nx),mont(nz,ny,nx),fcor(ny,nx),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: ddtphi(nz,ny,nx), omega(nz,ny,nx), sig(nz,ny,nx), dump(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res,v,mont
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    ! -----------------------------------------------------------------
    !
    call def_prescor(dump,ddtphi, nx,ny,nz, u,v,mont,fcor, dx,dy)
    call vor(omega, nx,ny,nz, u,v, dx,dy)
    call def_total(sig, nx,ny,nz, u,v, dx,dy)
    !
    ! dot_def_angle is - dot(phi)
    res = (omega-2.*ddtphi)/sig
  end subroutine
  !
  !@ Calculate Lagrangian acceleration gradient tensor eigenvalues
  !@ 
  !@ For information on this diagnostic refer to Hua and Klein (1998).
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ mont : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Montgomery potential on an isentropic surface or geopotential in pressure surface.
  !@ lat : np.ndarray with shape (ny,nx) and dtype float64
  !@     Latitude of the grid point.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ 4-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The real and imaginary components of the two eigenvalues.
  subroutine accgrad_eigs(respr,respi,resmr,resmi,nx,ny,nz,u,v,mont,lat,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx),v(nz,ny,nx),mont(nz,ny,nx)
    real(kind=nr), intent(in)  :: lat(ny),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: respr(nz,ny,nx),respi(nz,ny,nx)
    real(kind=nr), intent(out) :: resmr(nz,ny,nx), resmi(nz,ny,nx)
    type::tensor
        real(kind=nr),dimension(2,2)::t
    end type
    type(tensor)::gamma(nz,ny,nx)
    real(kind=nr) :: accelx(nz,ny,nx),accely(nz,ny,nx)
    real(kind=nr) :: dxaccelx(nz,ny,nx),dyaccelx(nz,ny,nx)
    real(kind=nr) :: dxaccely(nz,ny,nx),dyaccely(nz,ny,nx)
    real(kind=nr) :: tem(2,2)
    real(kind=nr) :: eigensr(2),eigensi(2)
    real(kind=nr) :: dummy1(2,2),dummy2(2,2),dummy3(6)
    integer(kind=ni) :: i,j,k, nx,ny,nz, dummy4
    !f2py depend(nx,ny,nz) respr,resmr,respi,resmi,v,mont
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    ! -----------------------------------------------------------------
    !
    call uv_prescor(accelx,accely,nx,ny,nz,u,v,mont,lat,dx,dy)
    call grad(dxaccelx,dyaccelx,nx,ny,nz,accelx,dx,dy)
    call grad(dxaccely,dyaccely,nx,ny,nz,accely,dx,dy)
    forall(k = 1_ni:nz, j = 1_ni:ny, i = 1_ni:nx)
       gamma(k,j,i)%t(1,1) = dxaccelx(k,j,i)
       gamma(k,j,i)%t(1,2) = dyaccelx(k,j,i)
       gamma(k,j,i)%t(2,1) = dxaccely(k,j,i)
       gamma(k,j,i)%t(2,2) = dyaccely(k,j,i)
    end forall 
    do i=1_ni,nx
       do j=2_ni,ny-1_ni
          do k=1_ni,nz  
            tem=gamma(k,j,i)%t
            call dgeev('N','N',2,tem,2,eigensr,eigensi,dummy1,2,dummy2,2,dummy3,6,dummy4)
            !like eigens(1:2) = eig(gamma(k,j,i)%t)
            respr(k,j,i)=eigensr(1) 
            respi(k,j,i)=eigensi(1) 
            resmr(k,j,i)=eigensr(2) 
            resmi(k,j,i)=eigensi(2) 
          end do
       end do
    end do
    !
    respr(:,(/1,ny/),:)=nan
    resmr(:,(/1,ny/),:)=nan
    respi(:,(/1,ny/),:)=nan
    resmi(:,(/1,ny/),:)=nan
    if ( .not. grid_cyclic_ew ) then
       resmi(:,:,(/1,nx/))=nan
       respi(:,:,(/1,nx/))=nan
       resmr(:,:,(/1,nx/))=nan
       respr(:,:,(/1,nx/))=nan
    end if
  end subroutine
  !
end module
