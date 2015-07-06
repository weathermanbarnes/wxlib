! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- tendencies of meteorological variables
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!@ Tendencies of meteorological variables. 
!@ 
!@ For each variable the adiabatic and diabatic tendencies are separate routines, 
!@ as far as available.
module tend
  use kind
  use config
  use derivatives
  use diag
  !
  implicit none
  !
  real(kind=nr), allocatable :: dzdth(:,:,:,:), slope(:,:,:,:), dzdx(:,:,:,:), dzdy(:,:,:,:), w(:,:,:,:)
  integer(kind=ni) :: slope_nt=-1_ni, slope_nx=-1_ni, slope_ny=-1_ni, slope_nz=-1_ni
contains
  !
  !@ Calculate the deformation tendencies due to horizontal advection
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
  subroutine def_adv(tendabs,tendang,nx,ny,nz,u,v,dx,dy)
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
  subroutine def_adv3d(tendabs,tendang,nx,ny,nz,u,v,uz,vz,w,dx,dy)
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
  subroutine def_beta(tendabs,tendang,nx,ny,nz,u,v,beta,dx,dy)
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
  !@ Parameters
  !@ ----------
  !@
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     U-wind velocity.
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     V-wind velocity.
  !@ geop : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Geopotential or a related potential.
  !@ fcor : np.ndarray with shape (ny,nx) and dtype float64
  !@     Locally varying Coriolis parameter.
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
  subroutine def_prescor(tendabs,tendang,nx,ny,nz,u,v,geop,fcor,dx,dy)
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
  subroutine def_tilt(tendabs,tendang,nx,ny,nz,u,v,uz,vz,w,dx,dy)
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
  subroutine uv_prescor(resx,resy,nx,ny,nz,u,v,z,lat,dx,dy)
    use consts!, only: pi, omega_rot
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
       f(k,:,i) = 2.0_nr * omega_rot * sin(lat*pi/180._nr)
    end forall
    !
    call grad(a_pressx,a_pressy,nx,ny,nz,z,dx,dy)
    resx = -a_pressx + f*v ! pressure force + coriolis force
    resy = -a_pressy - f*u ! pressure force + coriolis force
  end subroutine
  !
  !@ Tendency due to frontogenesis
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
  !@     Number of isentropic levels.
  !@ nt : int 
  !@     Grid size in t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Tendency of the slope due to frontogenesis.
  subroutine slope_fronto(tend, nx,ny,nz,nt, u,v, dx,dy)
    use consts !, only: p0, kappa
    !
    real(kind=nr), intent(in) :: u(nt,nz,ny,nx), v(nt,nz,ny,nx), &
            &                    dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: tend(nt,nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,nt
    !
    real(kind=nr) :: def_st(nz,ny,nx), def_sh(nz,ny,nx), diverg(nz,ny,nx)
    integer(kind=ni) :: n
    !f2py depend(nx,ny,nz,nt) v,z,tend
    !f2py depend(nx,ny,nz) dth
    !f2py depend(nx,ny) dx,dy
    ! -----------------------------------------------------------------
    !
    if ( slope_nt /= nt .or. slope_nz /= nz .or. slope_ny /= ny .or. slope_nx /= nx ) then
       write(*,*) 'Input does not match prepared shape, did you forget to call slope_prepare before?'
       stop 1
    end if
    !
    do n = 1_ni,nt
       call def_shear(def_sh, nx,ny,nz, u(n,:,:,:),v(n,:,:,:), dx,dy)
       call def_stretch(def_st, nx,ny,nz, u(n,:,:,:),v(n,:,:,:), dx,dy)
       call div(diverg, nx,ny,nz, u(n,:,:,:),v(n,:,:,:), dx,dy)
       !
       tend(n,:,:,:) = ( 0.5_nr * def_st * (dzdy(n,:,:,:)**2 - dzdx(n,:,:,:)**2)  & 
              &                 - def_sh *  dzdx(n,:,:,:)    * dzdy(n,:,:,:)      &
              &        - 0.5_nr * diverg * (dzdx(n,:,:,:)**2 + dzdy(n,:,:,:)**2) ) / slope(n,:,:,:)
    end do
    !
  end subroutine
  !
  !@ Tendency due to tilting
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
  !@     Number of isentropic levels.
  !@ nt : int 
  !@     Grid size in t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Tendency of the slope due to full vertical motion.
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Tendency of the slope due to cross-isentropic vertical motion.
  subroutine slope_tilt(tend,tend_ci, nx,ny,nz,nt, u,v, dx,dy)
    use consts !, only: p0, kappa
    !
    real(kind=nr), intent(in) :: u(nt,nz,ny,nx), v(nt,nz,ny,nx), &
            &                    dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: tend(nt,nz,ny,nx), tend_ci(nt,nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,nt
    !
    real(kind=nr) :: dwdx(nz,ny,nx), dwdy(nz,ny,nx), w_ci(nt,nz,ny,nx)
    integer(kind=ni) :: n
    !f2py depend(nx,ny,nz) u,v,w, z, tend, tend_ci
    !f2py depend(nx,ny,nz) dth
    !f2py depend(nx,ny) dx,dy
    ! -----------------------------------------------------------------
    !
    if ( slope_nt /= nt .or. slope_nz /= nz .or. slope_ny /= ny .or. slope_nx /= nx ) then
       write(*,*) 'Input does not match prepared shape, did you forget to call slope_prepare before?'
       stop 1
    end if
    !
    ! Determine cross isentropic vertical motion
    w_ci = w - u * dzdx - v * dzdy
    !
    do n = 1_ni,nt
       call grad(dwdx,dwdy, nx,ny,nz, w(n,:,:,:), dx,dy)
       tend(n,:,:,:) = (dwdx * dzdx(n,:,:,:) + dwdy * dzdy(n,:,:,:))/slope(n,:,:,:)
       !
       call grad(dwdx,dwdy, nx,ny,nz, w_ci(n,:,:,:), dx,dy)
       tend_ci(n,:,:,:) = (dwdx * dzdx(n,:,:,:) + dwdy * dzdy(n,:,:,:))/slope(n,:,:,:)
    end do
    !
  end subroutine
  !
  !@ Tendency due to horizontal advection
  !@ 
  !@ TODO: Clarify "residual" calculation in isal.f90
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
  !@     Number of isentropic levels.
  !@ nt : int 
  !@     Grid size in t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Tendency of the slope due to horizontal advection.
  subroutine slope_adv(tend, nx,ny,nz,nt, u,v, dx,dy)
    use consts !, only: p0, kappa
    !
    real(kind=nr), intent(in) :: u(nt,nz,ny,nx), v(nt,nz,ny,nx), &
            &                    dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: tend(nt,nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,nt
    !
    real(kind=nr) :: dslopedx(nz,ny,nx), dslopedy(nz,ny,nx), dslopedth(nz,ny,nx)
    integer(kind=ni) :: n
    !f2py depend(nx,ny,nz) u,v,w, z, tend, tend_ci
    !f2py depend(nx,ny,nz) dth
    !f2py depend(nx,ny) dx,dy
    ! -----------------------------------------------------------------
    !
    if ( slope_nt /= nt .or. slope_nz /= nz .or. slope_ny /= ny .or. slope_nx /= nx ) then
       write(*,*) 'Input does not match prepared shape, did you forget to call slope_prepare before?'
       stop 1
    end if
    !
    do n = 1_ni,nt
       call grad(dslopedx,dslopedy, nx,ny,nz, slope(n,:,:,:), dx,dy)
       tend(n,:,:,:) = u(n,:,:,:) * dslopedx + v(n,:,:,:) * dslopedy
    end do
    !
  end subroutine
  !
  !@ Tendency due to a given diabatic heating (physical temperature tendency)
  !@ 
  !@ Parameters
  !@ ----------
  !@ 
  !@ diab : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Diabatic heating at isentropic levels.
  !@ z : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Geopotential at isentropic levels.
  !@ p : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Pressure at isentropic levels.
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
  !@     Number of isentropic levels.
  !@ nt : int 
  !@     Grid size in t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Tendency of the slope due to the given diabatic heating field.
  subroutine slope_diab(tend, nx,ny,nz,nt, diab,p, dx,dy)
    use consts !, only: p0, kappa
    !
    real(kind=nr), intent(in) :: diab(nt,nz,ny,nx), p(nt,nz,ny,nx), &
            &                    dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: tend(nt,nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,nt
    !
    real(kind=nr) :: ddiabdx(nz,ny,nx), ddiabdy(nz,ny,nx)
    integer(kind=ni) :: n
    !f2py depend(nx,ny,nz) z,p, tend
    !f2py depend(nx,ny,nz) dth
    !f2py depend(nx,ny) dx,dy
    ! -----------------------------------------------------------------
    !
    if ( slope_nt /= nt .or. slope_nz /= nz .or. slope_ny /= ny .or. slope_nx /= nx ) then
       write(*,*) 'Input does not match prepared shape, did you forget to call slope_prepare before?'
       stop 1
    end if
    !
    tend = diab * (p0/p)**(Rl/cp)
    do n = 1_ni,nt
       call grad(ddiabdx,ddiabdy, nx,ny,nz, tend(n,:,:,:), dx,dy)
       tend(n,:,:,:) = - dzdth(n,:,:,:) * (ddiabdx * dzdx(n,:,:,:) + ddiabdy * dzdy(n,:,:,:))/slope(n,:,:,:)
    end do
    !
  end subroutine
  !
  !@ Prepare the calculation of the tendencies by setting some widely used variables
  !@ 
  !@ Parameters
  !@ ----------
  !@ 
  !@ diab : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Diabatic heating at isentropic levels.
  !@ z : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Geopotential at isentropic levels.
  !@ p : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Pressure at isentropic levels.
  !@ T : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Temperature at isentropic levels.
  !@ omega : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
  !@     Pressure vertical velocity at isentropic levels.
  !@ dx : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in x-direction to be directly for centered differences.
  !@     ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``.
  !@ dy : np.ndarray with shape (ny,nx) and dtype float64
  !@     The double grid spacing in y-direction to be directly for centered differences.
  !@     ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
  !@ dth : np.ndarray with shape (nz-2,ny,nx) and dtype float64
  !@     The double grid spacing in th-direction to be directly for centered differences.
  !@     ``dth(k,j,i)`` is expected to contain the th-distance between ``(k+1,j,i)`` and ``(k-1,j,i)``.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nz : int
  !@     Number of isentropic levels.
  !@ nt : int 
  !@     Grid size in t-direction.
  subroutine slope_prepare(nx,ny,nz,nt, z,p,T,omega, dx,dy,dth)
    real(kind=nr), intent(in) :: z(nt,nz,ny,nx), p(nt,nz,ny,nx), T(nt,nz,ny,nx), &
            &                    omega(nt,nz,ny,nx), &
            &                    dx(ny,nx), dy(ny,nx), dth(nz-2_ni,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz,nt
    !
    real(kind=nr) :: rho(nz,ny,nx)
    integer(kind=ni) :: n
    !f2py depend(nx,ny,nz,nt) p,T,omega
    !f2py depend(nx,ny,nz) dth
    !f2py depend(nx,ny) dx,dy
    ! -----------------------------------------------------------------
    !
    if ( allocated(slope) ) deallocate(slope)
    if ( allocated(dzdx) ) deallocate(dzdx)
    if ( allocated(dzdy) ) deallocate(dzdy)
    if ( allocated(dzdth) ) deallocate(dzdth)
    if ( allocated(w) ) deallocate(w)
    !
    allocate(slope(nt,nz,ny,nx), dzdx(nt,nz,ny,nx), dzdy(nt,nz,ny,nx), dzdth(nt,nz,ny,nx), w(nt,nz,ny,nx))
    !
    slope_nx = nx
    slope_ny = ny
    slope_nz = nz
    slope_nt = nt
    !
    call isen_geop_slope(slope,dzdx,dzdy,dzdth, nx,ny,nz,nt, z, dx,dy,dth)
    do n = 1_ni,nt
       call rho_from_T_p(rho, nx,ny,nz, T(n,:,:,:), p(n,:,:,:))
       call w_from_omega(w(n,:,:,:), nx,ny,nz, omega(n,:,:,:), rho)
    end do
    !
  end subroutine
end module
