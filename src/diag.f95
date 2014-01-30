! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- diagnostical functions
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module diag
  use kind
  use config
  use derivatives
  !
  implicit none
contains
  !
  ! Calculates vorticity vor=dxv-dyu
  subroutine vor(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dyu(nz,ny,nx),dxv(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    call ddx(dxv,nx,ny,nz,v,dx,dy)
    call ddy(dyu,nx,ny,nz,u,dx,dy)
    res=dxv-dyu
  end subroutine
  !
  ! Calculates divergence div=dxu+dyv
  subroutine div(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dxu(nz,ny,nx),dyv(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    call ddx(dxu,nx,ny,nz,u,dx,dy)
    call ddy(dyv,nx,ny,nz,v,dx,dy)
    res=dxu+dyv
  end subroutine
  !
  ! Calculates shear deformation def_shear = dyu+dxv
  subroutine def_shear(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dyu(nz,ny,nx),dxv(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    call ddx(dxv,nx,ny,nz,v,dx,dy)
    call ddy(dyu,nx,ny,nz,u,dx,dy)
    res=dyu+dxv
  end subroutine
  !
  ! Calculates stretch deformation def_stretch =  dxu-dyv
  subroutine def_stretch(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: dxu(nz,ny,nx),dyv(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    call ddx(dxu,nx,ny,nz,u,dx,dy)
    call ddy(dyv,nx,ny,nz,v,dx,dy)
    res=dxu-dyv
  end subroutine
  !
  ! Calculates total [coordinate independent] deformation 
  ! def_total = sqrt(def_stretch^2+def_shear^2)
  subroutine def_total(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v
    !f2py depend(nx,ny) dx, dy
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)
    res=sqrt(sig_sh**2+sig_st**2)
  end subroutine
  !
  ! Calculates total [coordinate independent] deformation tendency due to
  ! pressure
  subroutine tend_def_total_pres(res,nx,ny,nz,u,v,z,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), z(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig_sh(nz,ny,nx), sig_st(nz,ny,nx), zxxyy(nz,ny,nx), z2xy(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v, z
    !f2py depend(nx,ny) dx, dy
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)
    call antilap2(zxxyy,nx,ny,nz,z,dx,dy)
    call crosslap2(z2xy,nx,ny,nz,z,dx,dy)
    res = - (sig_st*zxxyy + sig_sh*z2xy)/sqrt(sig_sh**2+sig_st**2)
  end subroutine
  !
  ! Calculates angle from x-axis to axis of dilatation :
  !   = delta (Keyser Reeder Reed)
  !   =   (Lapeyre Klein Hua) 
  !   = alpha (Markowski Richardson)
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
  ! Calculates angle from wind direction to axis of dilatation 
  ! = the same as def_angle, but in natural coordinates
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
  ! Calculates the deformation tendencies due to the beta-effect
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
  ! Calculates the deformation tendencies due to the pressure terms and Coriolis
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
  ! Calculates the deformation tendencies due to "tilting" 
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
  ! Calculates 3d deformation 
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
  ! Calculates angle from x-axis to isolines of the input field dat (direction k X grad(dat))
  !   = alpha (Keyser Reeder Reed)
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
    res=atan2(-gradx,grady)
  end subroutine
  !
  ! Calculates beta = angle between dilatation axis (u,v) and iso-lines of the input field dat
  !   = beta (Keyser Reeder Reed)
  !   = beta (Markowski Richardson)
  !   = theta + phi + pi/4 (Lapeyre Klein Hua)
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
    res=delta-alpha
  end subroutine
  !
  ! FRONTOGENESIS: STRETCHING AND STIRRING RATES:
  ! Calculates (stretch, stir) where:
  ! stretch= fractional dat gradient stretching rate = 1/|grad(dat)| * d/dt(|grad(dat)|)
  !        = gamma, 'stretching rate' (Lapeyre Klein Hua)
  !        = -1/|grad(dat=PV)| * Fn (Keyser Reeder Reed)
  !        =  1/|grad(dat=PV)| * F  (Markowski Richardson)
  !   stir = angular rotation rate of grad(dat) (stirring rate)
  !        = d(theta)/dt      (Lapeyre Klein Hua)
  !        = 1/|grad(dat=PV)| * Fs  (Keyser Reeder Reed)
  subroutine stretch_stir(resstretch,resstir,nx,ny,nz,u,v,dat,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), dat(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: resstretch(nz,ny,nx),resstir(nz,ny,nx)
    real(kind=nr) :: bet(nz,ny,nx),totdef(nz,ny,nx),divergence(nz,ny,nx),vorticity(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) resstretch,resstir,v,dat
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
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
    resstretch=-0.5_nr*(divergence-totdef*cos(2_nr*bet))
    resstir=0.5_nr*(vorticity+totdef*sin(2_nr*bet))
  end subroutine
  !
  ! Calculates Okubo-Weiss criterion lambda_0=1/4 (sigma^2-omega^2)= 1/4 W
  ! This is the square of the eigenvalues in Okubo's paper (assumes small divergence)
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
    res=0.25*(sig**2-omega**2)
  end subroutine
  !
  ! Calculates Lagrangian acceleration on the isentropic surface
  subroutine dot_uv(resx,resy,nx,ny,nz,u,v,mont,lat,dx,dy)
    use consts!, only: pi, omega
    real(kind=nr), intent(in)  :: u(nz,ny,nx),v(nz,ny,nx),mont(nz,ny,nx),lat(ny),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: resx(nz,ny,nx),resy(nz,ny,nx)
    real(kind=nr) :: a_pressx(nz,ny,nx), a_pressy(nz,ny,nx), f(nz,ny,nx)
    integer(kind=ni) :: i,k,nx,ny,nz
    !f2py depend(nx,ny,nz) resx,resy,v,mont
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    ! -----------------------------------------------------------------
    !
    forall(k = 1_ni:nz, i = 1_ni:nx)
       f(k,:,i) = 2.*omega*sin(lat*pi/180._nr)
    end forall
    call grad(a_pressx,a_pressy,nx,ny,nz,mont,dx,dy)
    resx=-a_pressx+f*v ! pressure force + coriolis force
    resy=-a_pressy-f*u ! pressure force + coriolis force
  end subroutine
  !
  ! Calculates lagrangian acceleration gradient tensor eigenvalues
  ! ref: Hua and Klein (1998)
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
    call dot_uv(accelx,accely,nx,ny,nz,u,v,mont,lat,dx,dy)
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
  ! Calculates Lagrangian time derivative of compression axis angle
  !  d(gamma)/dt (ref Spensberger and Spengler 2013)
  ! -d(phi)/dt (ref Lapeyre et al 1999)
  subroutine dot_def_angle(res,nx,ny,nz,u,v,mont,lat,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx),v(nz,ny,nx),mont(nz,ny,nx),lat(ny),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig_st(nz,ny,nx),sig_sh(nz,ny,nx)
    real(kind=nr) :: accelx(nz,ny,nx),accely(nz,ny,nx)
    real(kind=nr) :: dxaccelx(nz,ny,nx),dyaccelx(nz,ny,nx)
    real(kind=nr) :: dxaccely(nz,ny,nx),dyaccely(nz,ny,nx)
    real(kind=nr) :: ddtsig_st(nz,ny,nx),ddtsig_sh(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res,v,mont
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    ! -----------------------------------------------------------------
    !
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)
    call dot_uv(accelx,accely,nx,ny,nz,u,v,mont,lat,dx,dy)
    ! calculate the lagrangian acceleration gradient tensor:
    call grad(dxaccelx,dyaccelx,nx,ny,nz,accelx,dx,dy)
    call grad(dxaccely,dyaccely,nx,ny,nz,accely,dx,dy)
    ! 1. Calculate d/dt(sig_st)
    ddtsig_st=-(dxaccelx-dyaccely)
    ! 2. calculate d/dt(sig_sh)
    ddtsig_sh=-(dyaccelx+dxaccely)
    ! 3. Calculate d/dt(phi)
    res=0.5*(sig_st*ddtsig_sh-sig_sh*ddtsig_st)/(sig_sh**2+sig_st**2)
  end subroutine
  !
  !Calculates the r diagnostic of Lapeyre et al (1999)
  ! r is the ratio of effective rotation to strain rate
  ! where effective rotation comprises both vorticity and strain-axes rotation
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
    call dot_def_angle(ddtphi,nx,ny,nz,u,v,mont,lat,dx,dy)
    call vor(omega,nx,ny,nz,u,v,dx,dy)
    call def_total(sig,nx,ny,nz,u,v,dx,dy)
    !
    ! dot_def_angle is - dot(phi)
    res=(omega-2.*ddtphi)/sig
  end subroutine
  !
  ! Calculates geostrophic velocity [ug,vg] = v_g(mont,lat)
  subroutine v_geo_from_montgp(resx,resy,nx,ny,nz,mont,lat,dx,dy)
    use consts!, only: pi, omega
    real(kind=nr), intent(in)  :: mont(nz,ny,nx),lat(ny),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: resx(nz,ny,nx), resy(nz,ny,nx)
    real(kind=nr) :: montx(nz,ny,nx), monty(nz,ny,nx),f(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) resx, resy
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    ! -----------------------------------------------------------------
    !
    call grad(montx,monty,nx,ny,nz,mont,dx,dy)
    forall(k = 1_ni:nz, i = 1_ni:nx)
       f(k,:,i) = 2*omega*sin(lat*pi/180._nr)
    end forall
    where (f==0._nr) f=9.E99_nr !avoid singularity calculating v_g at equator
    resx=-monty/f
    resy=montx/f
  end subroutine
  !
  ! Calculates the geopotential (res) from montgomery potential (m), potential
  ! temperature (theta) and pressure (p)
  subroutine geop_from_montgp(res,nx,ny,nz,m,theta,p,dx,dy)
    use consts!, only: cp, Rl, p0
    !
    real(kind=nr), intent(in)  :: m(nz,ny,nx), theta(nz,ny,nx), p(nz,ny,nx), &
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
             res(k,j,i) = m(k,j,i) - cp*theta(k,j,i)*(p(k,j,i)/p0)**(Rl/cp)
          end do
       end do
    end do
  end subroutine
  !
end module
