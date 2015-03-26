! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- diagnostical functions
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!@ Diagnostic functions for dynamic variables. 
!@ 
!@ This module contains diagnostic functions useful for analysing the wind
!@ field. For diagnostics relating to thermodynamic variables, refer to the
!@ :mod:`dynlib.humidity` module.
module diag
  use kind
  use config
  use derivatives
  !
  implicit none
contains
  !
  !@ Calculate vorticity::
  !@
  !@     vor = dv/dx - du/dy 
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
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
  !@     Calculated vorticity.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`div`, :meth:`def_shear`, :meth:`def_stretch`
  subroutine vor(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dyu(nz,ny,nx),dxv(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(dxv,nx,ny,nz,v,dx,dy)
    call ddy(dyu,nx,ny,nz,u,dx,dy)
    res = dxv - dyu
  end subroutine
  !
  !@ Calculate divergence::
  !@
  !@     div = du/dx + dv/dy 
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
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
  !@     Calculated divergence.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`vor`, :meth:`def_shear`, :meth:`def_stretch`
  subroutine div(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dxu(nz,ny,nx),dyv(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(dxu,nx,ny,nz,u,dx,dy)
    call ddy(dyv,nx,ny,nz,v,dx,dy)
    res = dxu + dyv
  end subroutine
  !
  !@ Calculate shear deformation::
  !@
  !@     def_shear = du/dy + dv/dx
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
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
  !@     Calculated shear deformation.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_stretch`, :meth:`def_total`, :meth:`def_angle`, :meth:`vor`, :meth:`div`
  subroutine def_shear(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dyu(nz,ny,nx),dxv(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    call ddx(dxv,nx,ny,nz,v,dx,dy)
    call ddy(dyu,nx,ny,nz,u,dx,dy)
    res = dyu + dxv
  end subroutine
  !
  !@ Calculate stretch deformation::
  !@
  !@    def_stretch =  du/dx - dv/dy
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
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
  !@     Calculated stretching deformation.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_shear`, :meth:`def_total`, :meth:`def_angle`, :meth:`vor`, :meth:`div`
  subroutine def_stretch(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dxu(nz,ny,nx),dyv(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(dxu,nx,ny,nz,u,dx,dy)
    call ddy(dyv,nx,ny,nz,v,dx,dy)
    res = dxu - dyv
  end subroutine
  !
  !@ Calculate total deformation
  !@
  !@ Total deformation is a coordinate-system independent measure for the
  !@ strength of deformation::
  !@
  !@    def_total = sqrt(def_stretch^2 + def_shear^2)
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The u-wind velocity field.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The v-wind velocity field.
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
  !@     Calculated total deformation.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_angle`, :meth:`def_angle_nat`
  subroutine def_total(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)
    res = sqrt(sig_sh**2_ni + sig_st**2_ni)
  end subroutine
  !
  !@ Calculate total deformation tendency due the pressure term
  !@
  !@ TODO: Should the deformation tendencies be in the general library?
  !@ TODO: There seems to be a duplication with def_tend_prescor
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ z : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Geopotential on a pressure surface or equivalent potentials on other
  !@     surfaces.
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
  !@     Calculated total deformation tendency.
  subroutine def_total_tend_pres(res,nx,ny,nz,u,v,z,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), z(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx), zxxyy(nz,ny,nx), z2xy(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v, z
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)
    call antilap2(zxxyy,nx,ny,nz,z,dx,dy)
    call crosslap2(z2xy,nx,ny,nz,z,dx,dy)
    res = - (sig_st*zxxyy + sig_sh*z2xy)/sqrt(sig_sh**2_ni + sig_st**2_ni)
  end subroutine
  !
  !@ Calculate the angle between the x-axis and the axis of dilatation
  !@ 
  !@ This angle goes by many different symbols:
  !@
  !@  * Spensberger and Spengler (2014): ``gamma``
  !@  * Markowski and Richardson (2011): ``alpha``
  !@  * Keyser, Reeder and Reed (1988), Lapeyre, Klein and Hua (1999): ``gamma``
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
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
  !@     Calculated deformation angle.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_angle_nat`, :meth:`def_total`
  subroutine def_angle(res,nx,ny,nz,u,v,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)    
    res = 0.5_nr * atan2(sig_sh,sig_st)
  end subroutine
  !
  !@ Calculate the angle between the wind direction and the axis of dilatation
  !@ 
  !@ Yields the equivalent of def_angle in natural coordinates.
  !@
  !@ Parameters
  !@ ----------
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@
  !@ Other parameters
  !@ ----------------
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
  !@     Calculated deformation angle.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`def_angle`, :meth:`def_total`
  subroutine def_angle_nat(res,nx,ny,nz,u,v,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)    
    res = 0.5_nr * atan2(sig_sh,sig_st)
    !
    ! rotate along wind direction instead of x-axis
    res = res - atan2(v,u)
    where ( res >= pi/2.0_nr )
       res = res - pi
    end where
    where ( res < -pi/2.0_nr )
       res = res + pi
    end where
  end subroutine
  !
  !@ Calculate the local wind shear in natural coordinates
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
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
  !@     Calculated wind shear.
  subroutine shear_nat(res,nx,ny,nz,u,v,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: ff(nz,ny,nx), ffx(nz,ny,nx), ffy(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ff = sqrt(u(:,:,:)**2_ni + v(:,:,:)**2_ni)
    call grad(ffx, ffy, nx,ny,nz, ff, dx,dy)
    !
    res(:,:,:) = (u(:,:,:)*ffy(:,:,:) - v(:,:,:)*ffx(:,:,:))/ff(:,:,:)
  end subroutine
  !
  !@ Calculate the local wind shear gradient in natural coordinates
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
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
  !@     Calculated wind shear gradient.
  subroutine grad_shear_nat(res,nx,ny,nz,u,v,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: shear(nz,ny,nx), shearx(nz,ny,nx), sheary(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call shear_nat(shear, nx,ny,nz, u,v, dx,dy)
    call ddx(shearx, nx,ny,nz, shear, dx,dy)
    call ddy(sheary, nx,ny,nz, shear, dx,dy)
    !
    res(:,:,:) = (u(:,:,:)*sheary(:,:,:) - v(:,:,:)*shearx(:,:,:)) &
            &   /sqrt(u(:,:,:)**2_ni + v(:,:,:)**2_ni)
  end subroutine
  !
  !@ Calculate total deformation and the deformation angle with all shear removed
  !@
  !@ Wind shear is both reflected in vorticity and deformation. This function
  !@ calculates the part of deformation that cannot be interpreted as vorticity.
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
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
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated total deformation and deformation angle.
  subroutine def_nat_shearless(resabs,resang,nx,ny,nz,u,v,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: resabs(nz,ny,nx), resang(nz,ny,nx)
    real(kind=nr) :: ff2(nz,ny,nx), shear(nz,ny,nx), ux(nz,ny,nx), uy(nz,ny,nx), &
            &        vx(nz,ny,nx), vy(nz,ny,nx), def_sh(nz,ny,nx), def_st(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) resabs,resang, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call shear_nat(shear, nx,ny,nz, u,v, dx,dy)
    !
    ff2(:,:,:) = u(:,:,:)**2_ni + v(:,:,:)**2_ni
    !
    call grad(ux,uy, nx,ny,nz, u, dx, dy)
    call grad(vx,vy, nx,ny,nz, v, dx, dy)
    !
    ux(:,:,:) = ux(:,:,:) - shear(:,:,:)*(u(:,:,:)*v(:,:,:))/ff2(:,:,:)
    uy(:,:,:) = uy(:,:,:) - shear(:,:,:)*(u(:,:,:)**2_ni)/ff2(:,:,:)
    vx(:,:,:) = vx(:,:,:) - shear(:,:,:)*(v(:,:,:)**2_ni)/ff2(:,:,:)
    vy(:,:,:) = vy(:,:,:) - shear(:,:,:)*(u(:,:,:)*v(:,:,:))/ff2(:,:,:)
    !
    def_sh(:,:,:) = uy(:,:,:) + vx(:,:,:)
    def_st(:,:,:) = ux(:,:,:) - vy(:,:,:)
    !
    resabs = sqrt(def_sh**2_ni + def_st**2_ni)
    resang = 0.5_nr * atan2(def_sh, def_st)
    !
    ! rotate along wind direction instead of x-axis
    resang = resang - atan2(v,u)
    where ( resang >= pi/2.0_nr )
       resang = resang - pi
    end where
    where ( resang < -pi/2.0_nr )
       resang = resang + pi
    end where
  end subroutine
  !
  !@ Calculate the deformation tendencies due to horizontal advection
  !@ 
  !@ TODO: Should the deformation tendencies be in the general library?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
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
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated tendencies of total deformation and deformation angle.
  subroutine def_tend_adv(tendabs,tendang,nx,ny,nz,u,v,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: tendabs(nz,ny,nx), tendang(nz,ny,nx)
    real(kind=nr) :: defabs(nz,ny,nx), defabsx(nz,ny,nx), defabsy(nz,ny,nx), &
            &        defang(nz,ny,nx), defangx(nz,ny,nx), defangy(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) tendabs, tendang, v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_total(defabs, nx,ny,nz, u,v, dx,dy)
    call def_angle(defang, nx,ny,nz, u,v, dx,dy)
    !
    call grad(defabsx,defabsy, nx,ny,nz, defabs, dx,dy)
    call grad(defangx,defangy, nx,ny,nz, defang, dx,dy)
    !
    tendabs(:,:,:) = -(u(:,:,:)*defabsx(:,:,:) + v(:,:,:)*defabsy(:,:,:))/defabs(:,:,:)
    tendang(:,:,:) = -(u(:,:,:)*defangx(:,:,:) + v(:,:,:)*defangy(:,:,:))/defabs(:,:,:)
  end subroutine
  !
  !@ Calculate the deformation tendencies due to full 3d advection
  !@ 
  !@ TODO: Should the deformation tendencies be in the general library?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ uz : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Vertical shear of the u-wind velocity.
  !@ vz : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Vertical shear of the v-wind velocity.
  !@ w : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Vertical wind velocity.
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
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated tendencies of total deformation and deformation angle.
  subroutine def_tend_adv3d(tendabs,tendang,nx,ny,nz,u,v,uz,vz,w,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u (nz,ny,nx), v (nz,ny,nx), w(nz,ny,nx), &
            &                     uz(nz,ny,nx), vz(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: tendabs(nz,ny,nx), tendang(nz,ny,nx)
    real(kind=nr) :: defabs(nz,ny,nx), defabsx(nz,ny,nx), defabsy(nz,ny,nx), defabsz(nz,ny,nx), &
            &        defang(nz,ny,nx), defangx(nz,ny,nx), defangy(nz,ny,nx), defangz(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) tendabs, tendang, v, uz, vz, w
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_total(defabs, nx,ny,nz, u,v, dx,dy)
    call def_angle(defang, nx,ny,nz, u,v, dx,dy)
    !
    call def_total(defabsz, nx,ny,nz, uz,vz, dx,dy)
    call def_angle(defangz, nx,ny,nz, uz,vz, dx,dy)
    !
    call grad(defabsx,defabsy, nx,ny,nz, defabs, dx,dy)
    call grad(defangx,defangy, nx,ny,nz, defang, dx,dy)
    !
    tendabs(:,:,:) = -(u(:,:,:)*defabsx(:,:,:) + v(:,:,:)*defabsy(:,:,:) &
             &       + w(:,:,:)*defabsz(:,:,:))/defabs(:,:,:)
    tendang(:,:,:) = -(u(:,:,:)*defangx(:,:,:) + v(:,:,:)*defangy(:,:,:) & 
             &       + w(:,:,:)*defangz(:,:,:))/defabs(:,:,:)
  end subroutine
  !
  !@ Calculate the deformation tendencies due to the beta-effect
  !@ 
  !@ TODO: Should the deformation tendencies be in the general library?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ beta : np.ndarray with shape (ny,nx) and dtype float64
  !@     Locally varying beta parameter.
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
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated tendencies of total deformation and deformation angle.
  subroutine def_tend_beta(tendabs,tendang,nx,ny,nz,u,v,beta,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), beta(ny,nx), &
            &                     dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: tendabs(nz,ny,nx), tendang(nz,ny,nx)
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx), defabs2(ny,nx)
    integer(kind=ni) :: nx,ny,nz, k
    !f2py depend(nx,ny,nz) tendabs, tendang, v
    !f2py depend(nx,ny) beta, dx, dy
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)    
    !
    do k = 1_ni,nz
       defabs2(:,:) = sig_sh(k,:,:)**2_ni + sig_st(k,:,:)**2_ni
       tendabs(k,:,:) = beta(:,:)*(sig_st(k,:,:)*u(k,:,:) + sig_sh(k,:,:)*v(k,:,:))/defabs2
       tendang(k,:,:) = beta(:,:)*(sig_sh(k,:,:)*u(k,:,:) - sig_st(k,:,:)*v(k,:,:))/defabs2
    end do
  end subroutine
  !
  !@ Calculate the deformation tendencies due to the pressure terms and Coriolis
  !@ 
  !@ TODO: Should the deformation tendencies be in the general library?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ beta : np.ndarray with shape (ny,nx) and dtype float64
  !@     Locally varying beta parameter.
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
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated tendencies of total deformation and deformation angle.
  subroutine def_tend_prescor(tendabs,tendang,nx,ny,nz,u,v,geop,fcor,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), geop(nz,ny,nx), &
                &                 fcor(ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: tendabs(nz,ny,nx), tendang(nz,ny,nx)
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx), defabs2(ny,nx), &
                &    geop_xx(nz,ny,nx), geop_yy(nz,ny,nx), geop_xy(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz, k
    !f2py depend(nx,ny,nz) tendabs, tendang, v, geop
    !f2py depend(nx,ny) fcor, dx, dy
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh, nx,ny,nz, u,v, dx,dy)
    call def_stretch(sig_st, nx,ny,nz, u,v, dx,dy)
    call ddx2(geop_xx, nx,ny,nz, geop, dx,dy)
    call ddy2(geop_yy, nx,ny,nz, geop, dx,dy)
    call ddxy(geop_xy, nx,ny,nz, geop, dx,dy)
    !
    do k = 1_ni,nz
       defabs2(:,:) = sig_sh(k,:,:)**2_ni + sig_st(k,:,:)**2_ni
       tendabs(k,:,:) = -(sig_st(k,:,:)*(geop_xx(k,:,:) - geop_yy(k,:,:)) + &
               &     2_ni*sig_sh(k,:,:)* geop_xy(k,:,:))/defabs2
       tendang(k,:,:) = -(sig_sh(k,:,:)*(geop_xx(k,:,:) - geop_yy(k,:,:)) - &
               &     2_ni*sig_st(k,:,:)* geop_xy(k,:,:))/defabs2 - fcor
    end do
  end subroutine
  !
  !@ Calculate the deformation tendencies due to "tilting" 
  !@ 
  !@ TODO: Should the deformation tendencies be in the general library?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ uz : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Vertical shear of the u-wind velocity.
  !@ vz : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Vertical shear of the v-wind velocity.
  !@ w : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Vertical wind velocity.
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
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated tendencies of total deformation and deformation angle.
  subroutine def_tend_tilt(tendabs,tendang,nx,ny,nz,u,v,uz,vz,w,dx,dy)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), &
                &                 uz(nz,ny,nx), vz(nz,ny,nx), w(nz,ny,nx), &
                &                 dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: tendabs(nz,ny,nx), tendang(nz,ny,nx)
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx), defabs2(nz,ny,nx), &
                &    wx(nz,ny,nx), wy(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) tendabs, tendang, v, uz, vz, w
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh, nx,ny,nz, u,v, dx,dy)
    call def_stretch(sig_st, nx,ny,nz, u,v, dx,dy)
    call ddx(wx, nx,ny,nz, w, dx,dy)
    call ddy(wy, nx,ny,nz, w, dx,dy)
    !
    defabs2(:,:,:) = sig_sh(:,:,:)**2_ni + sig_st(:,:,:)**2_ni
    tendabs(:,:,:) = (sig_st(:,:,:)*(uz(:,:,:)*wx(:,:,:) - vz(:,:,:)*wy(:,:,:)) + &
            &         sig_sh(:,:,:)*(vz(:,:,:)*wx(:,:,:) + uz(:,:,:)*wy(:,:,:)) )/defabs2(:,:,:)
    tendang(:,:,:) = (sig_sh(:,:,:)*(uz(:,:,:)*wx(:,:,:) - vz(:,:,:)*wy(:,:,:)) - &
            &         sig_st(:,:,:)*(vz(:,:,:)*wx(:,:,:) + uz(:,:,:)*wy(:,:,:)) )/defabs2(:,:,:)
  end subroutine
  !
  !@ Calculate 3d deformation
  !@ 
  !@ **Note**: This subroutine is work-in-progress and likely to change in the
  !@ future. The generalisation of deformation to 3 dimensions is done using the
  !@ analogy to Lyapunov exponents.
  !@
  !@ TODO: Should work-in-progress routines live in a different name space?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ w : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Pressure vertical wind velocity.
  !@ rho : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Density of air, used to convert the pressure vertical velocity to a cartensian 
  !@     vertical velocity.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@ dz : np.ndarray with shape (nz-2,ny,nx) and dtype float64
  !@     The double grid spacing in z-direction to be directly for centered differences.
  !@     ``dy(k,j,i)`` is expected to contain the z-distance between ``(k+1,j,i)`` and ``(k-1,j,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Grid size in z-direction.
  !@ nt : int
  !@     Grid size in t-direction.
  !@
  !@ Returns
  !@ -------
  !@ 2-tuple of np.ndarray with shapes (nt,nz-2,ny,nx), (nt,nz-2,ny,nx,3) and dtype float64
  !@     Calculated 3D eigenvalues and eigenvectors showing the strength and
  !@     orientation of deformation in 3D.
  subroutine def_3d(res_eval,res_evec,nx,ny,nz,nt,u,v,w,rho,dx,dy,dz)
    use consts
    !
    real(kind=nr), intent(in)  :: u(nt,nz,ny,nx), v(nt,nz,ny,nx), w(nt,nz,ny,nx), & 
            &                     rho(nt,nz,ny,nx),                               &
            &                     dx(ny,nx), dy(ny,nx), dz(nt,2_ni:nz-1_ni,ny,nx)
    real(kind=nr), intent(out) :: res_eval(nt,2_ni:nz-1_ni,ny,nx,3_ni), &
                                  res_evec(nt,2_ni:nz-1_ni,ny,nx,3_ni,3_ni)
    integer(kind=ni) :: nx,ny,nz,nt
    !f2py depend(nx,ny,nz,nt) res, v, w
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: ux(nt,nz,ny,nx), uy(nt,nz,ny,nx), uz(nt,nz,ny,nx), &
            &        vx(nt,nz,ny,nx), vy(nt,nz,ny,nx), vz(nt,nz,ny,nx), &
            &        wx(nt,nz,ny,nx), wy(nt,nz,ny,nx), wz(nt,nz,ny,nx), &
            &        vor_x(nt,nz,ny,nx), vor_y(nt,nz,ny,nx), vor_z(nt,nz,ny,nx), &
            &        div(nt,nz,ny,nx), wm(nt,nz,ny,nx)
    real(kind=nr) :: jac(3_ni,3_ni), evalr(3_ni), dummy0(3_ni), dummy1(1_ni,3_ni), &
            &        evec(3_ni,3_ni), work(102_ni)
    integer(kind=ni) :: i,j,k,n,t, info
    ! -----------------------------------------------------------------
    !
    ! Crude transformation from Pa/s to m/s
    wm = -w(:,:,:,:) / (rho(:,:,:,:)*g)
    !
    call grad_3d(ux,uy,uz, nx,ny,nz,nt, u, dx,dy,dz)
    call grad_3d(vx,vy,vz, nx,ny,nz,nt, v, dx,dy,dz)
    call grad_3d(wx,wy,wz, nx,ny,nz,nt, wm, dx,dy,dz)
    !
    ! TEST: Disregard vertical wind shear
    uz = wx
    vz = wy
    !
    ! Subtracting divergence and vorticity
    div = ux + vy + wz
    ux = ux - 1.0_nr/3.0_nr * div
    vy = vy - 1.0_nr/3.0_nr * div
    wz = wz - 1.0_nr/3.0_nr * div
    !
    vor_x = wy - vz
    vor_y = uz - wx
    vor_z = vx - uy
    !
    wy = wy - 1.0_nr/2.0_nr * vor_x
    vz = vz + 1.0_nr/2.0_nr * vor_x
    uz = uz - 1.0_nr/2.0_nr * vor_y
    wx = wx + 1.0_nr/2.0_nr * vor_y
    vx = vx - 1.0_nr/2.0_nr * vor_z
    uy = uy + 1.0_nr/2.0_nr * vor_z
    !
    ! Find the eigenvectors of the remaining symmetric zero-trace matrixes
    do i = 1_ni,nx
       write(*,'(I5,A4,I5,A)', advance='no') i, 'of', nx, cr
       !
       do j = 1_ni,ny
          do k = 2_ni,nz-1_ni
             do t = 1_ni,nt
                ! jac is in major column order
                jac = reshape( (/ ux(t,k,j,i), vx(t,k,j,i), wx(t,k,j,i), &
                           &      uy(t,k,j,i), vy(t,k,j,i), wy(t,k,j,i), &
                           &      uz(t,k,j,i), vz(t,k,j,i), wz(t,k,j,i) /), (/ 3_ni, 3_ni /) )
                call dgeev('N', 'V', 3_ni, jac, 3_ni, evalr, dummy0, & 
                        & dummy1, 1_ni, evec, 3_ni, work, 102_ni, info)
                ! sort eigenvalues and eigenvectors by eigenvalue
                ! first: max
                n = maxloc(evalr, dim=1_ni)
                res_eval(t,k,j,i,1_ni)   = evalr(n)
                res_evec(t,k,j,i,1_ni,:) = evec(:,n)
                ! second: middle
                n = 6_ni - n - minloc(evalr, dim=1_ni)
                res_eval(t,k,j,i,2_ni)   = evalr(n)
                res_evec(t,k,j,i,2_ni,:) = evec(:,n)
                ! last and least
                n = minloc(evalr, dim=1_ni)
                res_eval(t,k,j,i,3_ni)   = evalr(n)
                res_evec(t,k,j,i,3_ni,:) = evec(:,n)
             end do
          end do
       end do
    end do
    !
    return
  end subroutine
  !
  !@ Calculate angle between the x-axis and isolines of a given field
  !@ 
  !@ The direction of isolines is calculated by::
  !@
  !@     k times grad(dat)
  !@
  !@ Keyser, Reeder and Reed (1988) call this angle "alpha".
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Input field for the isolines.
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
  !@     Calculated isoline angle.
  subroutine isoline_angle(res,nx,ny,nz,dat,dx,dy)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: gradx(nz,ny,nx), grady(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(gradx,nx,ny,nz,dat,dx,dy)
    call ddy(grady,nx,ny,nz,dat,dx,dy)  
    res = atan2(-gradx,grady)
  end subroutine
  !
  !@ Calculate angle between the axis of dilatation and isolines of a given field
  !@ 
  !@ This angle goes by different names in the literature:
  !@
  !@  * Keyser, Reeder and Reed (1988), Markowski and Richardson (2011): ``beta``
  !@  * Lapeyre, Klein and Hua (1999): ``theta + phi + pi/4``
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Input field for the isolines.
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
  !@     Calculated isoline angle.
  subroutine isoline_to_deformation_angle(res,nx,ny,nz,u,v,dat,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: alpha(nz,ny,nx), delta(nz,ny,nx) 
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v, dat
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ! beta = delta - alpha (Keyser Reeder Reed) 
    !  alpha = angle between x axis and iso-dat line (direction k X grad dat)
    !  delta = angle between x axis and dilatation axis (dynlib: def_angle)
    call isoline_angle(alpha,nx,ny,nz,dat,dx,dy)
    call def_angle(delta,nx,ny,nz,u,v,dx,dy)
    res = delta - alpha
  end subroutine
  !
  !@ Calculate the frontogenesis function
  !@
  !@ Calculates both the streching and stirring rates of gradients in a given field.
  !@ The stretching rate is defined as::
  !@
  !@      1/|grad(dat)| * d/dt(|grad(dat)|)
  !@
  !@ and called gamma in Lapeyre, Klein and Hua (1999). This quantity alone is often
  !@ referred to as the [Petterssen] frontogenesis function 
  !@ (e.g. in Markowski and Richardson 2011). 
  !@ 
  !@ Keyser, Reeder and Reed (1988) generalise the Petterssen frontogenesis function to
  !@ its vector form, introducing the stirring rate of the gradient as its second component.
  !@ The stirring rate is defined as::
  !@
  !@      grad(dat)/|grad(dat)|² * (k \times d/dt(grad(dat)))
  !@ 
  !@ In the nomenclature of Lapeyre, Klein and Hua, this quantity is simply::
  !@
  !@      d(theta)/dt
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Input field for the isolines.
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
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The streching and stirring component of the generalised frontogenesis function.
  subroutine frontogenesis(resstretch,resstir,nx,ny,nz,u,v,dat,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: resstretch(nz,ny,nx),resstir(nz,ny,nx)
    real(kind=nr) :: bet(nz,ny,nx),totdef(nz,ny,nx),divergence(nz,ny,nx),vorticity(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) resstretch,resstir,v,dat
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ! (Lapeyre Klein Hua):
    !  stretch = gamma = 1/rho * d(rho)/dt (rho=|gradPV|) = -sigma/2 * sin(2(theta+phi))
    ! (Keyser Reeder Reed):
    !  d(gradPV)/dt = (Fn,Fs), coordinates rotated to dat-contours (Keyser Reeder Reed) 
    !  stretch = -1/|gradPV| * Fn 
    !       Fn = 0.5*|gradPV|(D-E*cos(2*beta))
    !  stir = Fs/|gradPV|
    !    Fs = 0.5*|gradPV|(vort+E*sin(2*beta))
    call isoline_to_deformation_angle(bet,nx,ny,nz,u,v,dat,dx,dy)
    call def_total(totdef,nx,ny,nz,u,v,dx,dy)
    call div(divergence,nx,ny,nz,u,v,dx,dy)
    call vor(vorticity,nx,ny,nz,u,v,dx,dy)
    !
    resstretch = -0.5_nr * (divergence - totdef*cos(2_nr*bet))
    resstir = 0.5_nr * (vorticity + totdef*sin(2_nr*bet))
  end subroutine
  !
  !@ Calculate the Okubo-Weiss parameter
  !@ 
  !@ The parameter is defined as::
  !@
  !@     1/4 (total deformation^2 - vorticity^2)
  !@
  !@ the square of the eigenvalues in Okubo's paper (assuming small divergence).
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
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
  !@     Calculated Okubo-Weiss parameter.
  subroutine okuboweiss(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx),v(nz,ny,nx),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig(nz,ny,nx), omega(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res,v
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call def_total(sig,nx,ny,nz,u,v,dx,dy)
    call vor(omega,nx,ny,nz,u,v,dx,dy)
    res = 0.25 * (sig**2 - omega**2)
  end subroutine
  !
  !@ Calculate the wind acceleration due to Coriolis and pressure gradient forces
  !@ 
  !@ Assuming frictionless flow, and either
  !@
  !@  * z being the geopotential on a pressure surface, or
  !@  * z being the Montgomery potential on an isentropic surface, 
  !@
  !@ the result is the Lagrangian acceleration of fluid parcels.
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ z : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Geopotential on a pressure surface or equivalent potentials on other
  !@     surfaces.
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
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     The acceleration of the u and v velocity components.
  subroutine uv_tend_prescor(resx,resy,nx,ny,nz,u,v,z,lat,dx,dy)
    use consts!, only: pi, omega
    !
    real(kind=nr), intent(in)  :: u(nz,ny,nx),v(nz,ny,nx),z(nz,ny,nx),lat(ny),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: resx(nz,ny,nx),resy(nz,ny,nx)
    real(kind=nr) :: a_pressx(nz,ny,nx), a_pressy(nz,ny,nx), f(nz,ny,nx)
    integer(kind=ni) :: i,k,nx,ny,nz
    !f2py depend(nx,ny,nz) resx,resy,v,z
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    ! -----------------------------------------------------------------
    !
    forall(k = 1_ni:nz, i = 1_ni:nx)
       f(k,:,i) = 2.0_nr * omega * sin(lat*pi/180._nr)
    end forall
    !
    call grad(a_pressx,a_pressy,nx,ny,nz,z,dx,dy)
    resx = -a_pressx + f*v ! pressure force + coriolis force
    resy = -a_pressy - f*u ! pressure force + coriolis force
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
  !@     Montgomery potential on an isentropic surface.
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
    call uv_tend_prescor(accelx,accely,nx,ny,nz,u,v,mont,lat,dx,dy)
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
    respr(:,(/1,ny/),:)=0._nr
    !respr(:,:,(/1,nx/))=0._nr
    resmr(:,(/1,ny/),:)=0._nr
    !resmr(:,:,(/1,nx/))=0._nr
    respi(:,(/1,ny/),:)=0._nr
    !respi(:,:,(/1,nx/))=0._nr
    resmi(:,(/1,ny/),:)=0._nr
    !resmi(:,:,(/1,nx/))=0._nr
  end subroutine
  !
  !@ Calculate deformation angle tendency due the pressure and Coriolis terms
  !@
  !@ TODO: Should the deformation tendencies be in the general library?
  !@ TODO: There seems to be a duplication with def_tend_prescor
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ z : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Geopotential on a pressure surface or equivalent potentials on other
  !@     surfaces.
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
  !@     Calculated deformation angle tendency.
  subroutine def_angle_tend_prescor(res,nx,ny,nz,u,v,z,lat,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx),v(nz,ny,nx),z(nz,ny,nx),lat(ny),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig_st(nz,ny,nx),sig_sh(nz,ny,nx)
    real(kind=nr) :: accelx(nz,ny,nx),accely(nz,ny,nx)
    real(kind=nr) :: dxaccelx(nz,ny,nx),dyaccelx(nz,ny,nx)
    real(kind=nr) :: dxaccely(nz,ny,nx),dyaccely(nz,ny,nx)
    real(kind=nr) :: ddtsig_st(nz,ny,nx),ddtsig_sh(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res,v,z
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)
    call uv_tend_prescor(accelx,accely,nx,ny,nz,u,v,z,lat,dx,dy)
    ! calculate the lagrangian acceleration gradient tensor:
    call grad(dxaccelx,dyaccelx,nx,ny,nz,accelx,dx,dy)
    call grad(dxaccely,dyaccely,nx,ny,nz,accely,dx,dy)
    ! 1. Calculate d/dt(sig_st)
    ddtsig_st = -(dxaccelx - dyaccely)
    ! 2. calculate d/dt(sig_sh)
    ddtsig_sh = -(dyaccelx + dxaccely)
    ! 3. Calculate d/dt(phi)
    res = 0.5 * (sig_st*ddtsig_sh - sig_sh*ddtsig_st) / (sig_sh**2 + sig_st**2)
  end subroutine
  !
  !@ Calculate the ratio between effective rotation and strain rate
  !@ 
  !@ This ratio is the main result of Lapeyre, Klein and Hua (1999). They denote it "r".
  !@ Their effective rotation takes into account both vorticity and rotation of the axis 
  !@ of dilatation.
  !@ 
  !@ TODO: Do you need to assume z=mont, or does the same routine work also on pressure levels?
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ mont : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Montgomery potential on an isentropic surface.
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
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated effective-rotation-to-strain ratio.
  subroutine rotation_strain_ratio(res,nx,ny,nz,u,v,mont,lat,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx),v(nz,ny,nx),mont(nz,ny,nx),lat(ny),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: ddtphi(nz,ny,nx), omega(nz,ny,nx), sig(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res,v,mont
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    ! -----------------------------------------------------------------
    !
    call def_angle_tend_prescor(ddtphi,nx,ny,nz,u,v,mont,lat,dx,dy)
    call vor(omega,nx,ny,nz,u,v,dx,dy)
    call def_total(sig,nx,ny,nz,u,v,dx,dy)
    !
    ! dot_def_angle is - dot(phi)
    res = (omega-2.*ddtphi)/sig
  end subroutine
  !
  !@ Calculate geostrophic velocity
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ z : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Geopotential on a pressure surface or equivalent potentials on other
  !@     surfaces.
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
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Calculated geostrophic U- and V-velocity components.
  subroutine uv_geo_from_pot(resx,resy,nx,ny,nz,z,lat,dx,dy)
    use consts!, only: pi, omega
    !
    real(kind=nr), intent(in)  :: z(nz,ny,nx),lat(ny),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: resx(nz,ny,nx), resy(nz,ny,nx)
    real(kind=nr) :: zx(nz,ny,nx), zy(nz,ny,nx),f(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) resx, resy
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    ! -----------------------------------------------------------------
    !
    call grad(zx,zy,nx,ny,nz,z,dx,dy)
    !
    forall(k = 1_ni:nz, i = 1_ni:nx)
       f(k,:,i) = 2.0_nr * omega * sin(lat*pi/180._nr)
    end forall
    where (f == 0._nr) f = 9.E99_nr !avoid singularity calculating v_g at equator
    !
    resx = -zy/f
    resy =  zx/f
  end subroutine
  !
  !@ Calculate the geopotential from montgomery potential, potential temperature and pressure
  !@ 
  !@ Parameters
  !@ ----------
  !@
  !@ mont : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Montgomery potential on an isentropic surface.
  !@ theta : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Potential temperature.
  !@ p : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Pressure on the isentropic surface.
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
  !@     Calculated geopotential.
  subroutine geop_from_montgp(res,nx,ny,nz,mont,theta,p,dx,dy)
    use consts!, only: cp, Rl, p0
    !
    real(kind=nr), intent(in)  :: mont(nz,ny,nx), theta(nz,ny,nx), p(nz,ny,nx), &
                                & dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res, theta, p
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    do i=1_ni,nx
       do j=1_ni,ny
          do k=1_ni,nz
             res(k,j,i) = mont(k,j,i) - cp*theta(k,j,i)*(p(k,j,i)/p0)**(Rl/cp)
          end do
       end do
    end do
  end subroutine
  !
end module
