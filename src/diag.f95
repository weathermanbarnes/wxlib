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
  ! Calculates angle from x-axis to axis of dilatation :
  !   = delta (Keyser Reeder Reed)
  !   =   (Lapeyre Klein Hua) 
  !   = alpha (Markowski Richardson)
  subroutine def_angle(res,nx,ny,nz,u,v,dx,dy)
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
    res=0.5_nr*atan2(sig_sh,sig_st)
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
    forall(i = 1_ni:nx, k = 1_ni:nz)
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
    forall(i = 1_ni:nx, j = 1_ni:ny, k = 1_ni:nz)
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
  ! d(phi)/dt (ref Lapeyre et al 1999)
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
    res=0.5*(sig_sh*ddtsig_st-sig_st*ddtsig_sh)/(sig_sh**2+sig_st**2)
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
    res=(omega+2.*ddtphi)/sig
  end subroutine
  !
  ! Gradient reversal: At each (i,j,k) grid point, finds the reversals of pv y-gradient. 
  !                    and classes them as:
  !                      c (cyclonic)
  !                      a (anticyclonic) 
  !
  ! Arguments:         pv: Potential vorticity pv(k,j,i) on (time, lat, lon) grid.
  !            highenough: array of flags, highenough(k,j,i) = {0 or 1} (type int*1)
  !                        denoting whether to test the point for reversal. This is
  !                        typically the output of highenough() funtion, which returns
  !                        1 where the surface is sufficiently above ground level and
  !             .          0 elsewhere.
  !             latitudes: vector of latitudes of the pv array
  !              ddythres: Cutoff y-gradient for pv. The magnitude of (negative) 
  !                        d(pv)/dy must be above ddythres for reversal to be detected.
  !             
  !
  ! Returns: int*1:: revc,   reva   (reversal flag) (threshold test applied)
  !          real::  revci,  revai  (reversal absolute gradient) 
  !                                 (threshold test applied)		
  !          real::  revciy, revaiy (reversal absolute y-gradient) 
  !                                 (no threshold test applied)
  !          int*1:: tested (flag to 1 all tested points: where highenough==1 
  !                          and not on edge of grid)
  !
  subroutine grad_rev(resa,resc,resai,resci,resaiy,resciy,tested,nx,ny,nz, &
       pv,highenough,latitudes,ddythres,dx,dy)
    real(kind=nr), intent(in)  :: pv(nz,ny,nx), latitudes(ny), ddythres, dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: resai(nz,ny,nx), resci(nz,ny,nx), & 
       resaiy(nz,ny,nx), resciy(nz,ny,nx)
    integer(kind=1), intent(out) :: resa(nz,ny,nx), resc(nz,ny,nx), tested(nz,ny,nx)
    integer(kind=1), intent(in) :: highenough(nz,ny,nx)
    real(kind=nr) :: gradx, grady, gradmag, gradang, hemi, sig, neg, tempval,edge
    real(kind=nr) :: xgradient(nz,ny,nx), ygradient(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) resa,resc,resai,resci,resaiy,resciy,tested, highenough,xgradient,ygradient
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) latitudes
    ! -----------------------------------------------------------------
    !
    call grad(xgradient,ygradient,nx,ny,nz,pv,dx,dy) ! returns grad at all testable points, 0,0 elsewhere
    !
    do i=1_ni,nx
       do j=1_ni,ny
          if (latitudes(j) > 0._nr) then !NH
             hemi = 1._nr
          else if (latitudes(j) < 0._nr) then !SH
             hemi = -1._nr
          else !equator
             hemi = 0._nr
          end if
          if (abs(latitudes(j))==90._nr) then 
             edge = 1._nr !at edge of grid, must skip
          elseif ((.NOT.(grid_cyclic_ew)).AND.((i==1).OR.(i==nx))) then
             edge=1._nr !at edge of grid, must skip
          else
             edge=0._nr
          end if
          do k=1_ni,nz
             !initialise all to 0, change later if reversal
             tested(k,j,i) = 0 
             resa(k,j,i)   = 0
             resai(k,j,i)  = 0._nr
             resaiy(k,j,i) = 0._nr
             resc(k,j,i)   = 0
             resci(k,j,i)  = 0._nr
             resciy(k,j,i) = 0._nr
             if ((edge==0._nr).AND.(highenough(k,j,i)==1)) then ! not on grid edge, and high enough: test it
                tested(k,j,i) = 1
                gradx=xgradient(k,j,i)
                grady=ygradient(k,j,i)
                if (grady<0._nr) then 
                   neg=1._nr !y-gradient is negative
                else
                   neg=0._nr
                end if
                if (grady<(-ddythres)) then  
                   sig=1._nr ! ygradient is significantly negative 
                else
                   sig=0._nr
                end if
                gradmag=sqrt(gradx**2.+grady**2.)
                if (gradmag .NE. 0._nr) then
                   ! angle exists if gradmag not 0
                   gradang=atan2(grady,gradx) !NB slow! faster to do it by signs
                   if ((gradang+3.141592/2.)*hemi>0._nr) then ! CW NH or CCW SH -> ACYC
                      resa(k,j,i)=int(sig)                 ! flag significant ygradient reversal
                      resai(k,j,i)=gradmag*sig        ! |gradPV| of significant reversal
                      resaiy(k,j,i)=-grady*neg        ! sum |dPV/dy| of all reversals
                   else ! CCW NH or CCW SH -> CYC
                      resc(k,j,i)=int(sig)  
                      resci(k,j,i)=gradmag*sig
                      resciy(k,j,i)=-grady*neg
                   end if
                end if
             end if  ! end test
          end do !k
       end do !j
    end do !i
  end subroutine
  !
  ! RWB detection, adapted from algorithm by Riviere
  !
  ! Detects the occurrence of anticyclonic and cyclonic wave-breaking 
  ! events from a PV field on isentropic coordinates.
  !
  ! Reference: Rivière (2009, hereafter R09): Effect of latitudinal 
  ! variations in low-level baroclinicity on eddy life cycles and upper-
  ! tropospheric wave-breaking processes. J. Atmos. Sci., 66, 1569–1592.
  ! See the appendix C.
  !
  ! Inputs: 
  !      pv_in(nz,ny,nx) : isentropic pv. Should be on a regular lat-lon grid 
  !                         and 180W must be the first longitude. 
  !                        (If 180W is not the first longitude, the outputs
  !                        will have 180W as the first, so must be rearranged) 
  !        lonvalues(nx) : vector of longitudes
  !        latvalues(ny) : vector of latitudes
  !                 ncon : number of contours to test, normally 41 or 21
  !                  lev : potential temperature of the level
  ! Outputs:
  !  beta_a_out(nz,ny,nx) : flag array, =1 if anticyclonic wave breaking 
  !  beta_c_out(nz,ny,nx) : flag array, =1 if cyclonic wave breaking 
  subroutine contour_rwb(beta_a_out,beta_c_out,nx,ny,nz,pv_in,lonvalues,latvalues,ncon,lev,dx,dy)
    use diag_contour_rwb
    !
    implicit none
    !
    real(kind=nr), intent(in) :: pv_in(nz,ny,nx),lonvalues(nx),latvalues(ny),dx(ny,nx),dy(ny,nx),lev
    integer(kind=ni), intent(in) :: ncon
    integer(kind=ni), intent(out) :: beta_a_out(nz,ny,nx), beta_c_out(nz,ny,nx)
    integer(kind=ni) :: beta_a_outint(nz,ny,nx), beta_c_outint(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) beta_a_out, beta_c_out,beta_a_outint, beta_c_outint
    !f2py depend(nx,ny) dx, dy
    !f2py depend(nx) lonvalues
    !f2py depend(ny) latvalues
    ! -----------------------------------------------------------------
    !
    integer(kind=ni) :: nlon, nlat, ntime, nplev
    integer(kind=ni) :: i,it,sd,j,l,nn,idecal
    !
    real(kind=nr), dimension(:,:,:,:), allocatable :: pv,pvint
    real(kind=nr), dimension(:), allocatable :: rdate
    real(kind=nr) :: plev(1)
    !
    character :: grandeur*2
    character :: cx*18,cy*18
    integer(kind=ni), dimension(:,:,:,:), allocatable :: beta_a,beta_c,beta_cint
    real(kind=nr), dimension(:,:,:), allocatable :: xvalf,yvalf
    real(kind=nr) :: incr,incr2,incrlat
    integer(kind=ni) :: nlon2,nlat2,choicegrid
    real(kind=nr), dimension(:), allocatable :: lonval2,latval2,lonval,latval
    real(kind=nr), dimension(:,:,:,:), allocatable :: pv2
    real(kind=nr), dimension(:,:,:), allocatable :: gamma_a,gamma_c
    real(kind=nr), dimension(:,:), allocatable :: moy_ga,moy_gc
    real(kind=nr), dimension(:,:,:), allocatable :: gamma_a1,gamma_c1
    real(kind=nr), dimension(:,:), allocatable :: moy_ga1,moy_gc1
    real(kind=nr),dimension(:,:),allocatable :: zon_clim_u
    real(kind=nr),dimension(:,:,:),allocatable :: clim_u
    !
    real(kind=nr) :: lonleft,latleft,lonright,latright,xlength,ylength
    !
    nlon = nx
    nlat = ny
    ntime = nz
    nplev = 1
    plev(1) = lev
    incr=lonvalues(2)-lonvalues(1)
    incrlat=latvalues(2)-latvalues(1)
    !
    ! CHANGING THE ORDER OF THE VALUES TO GET AN INCREASE IN LATITUDES 
    ! AND TO GET -180._nr0 AS THE FIRST LONGITUDE    
    !
    ! In the case incrlat<0 that is if the latitudes (latvalues) decrease with the
    ! corresponding index, we have to reverse the order to orientate from the South
    ! to the North
    !
    allocate(lonval(nlon))
    allocate(latval(nlat))
    allocate(pv(nlon,nlat,nplev,ntime))
    allocate(pvint(nlon,nlat,nplev,ntime))
    !
    ! Arrange pv into format required for this routine
    do it=1,ntime
       nn=1 ! only 1 level
       do j=1,nlat
          do i=1,nlon
             pv(i,j,nn,it) = pv_in(it,j,i)
          end do
       end do
    end do
    !
    pvint=pv
    lonval=lonvalues
    latval=latvalues
    !
    if (incrlat.lt.0._nr) then
       do j=1,nlat
          pvint(:,j,:,:)=pv(:,nlat-j+1,:,:)
          latval(j)=latvalues(nlat-j+1)
       end do
    end if
    !
    ! In the case, lonvalues(1) different from -180.0, we translate the field
    ! to have the first longitude equal to -180.0
    pv=pvint
    idecal=int((180._nr+lonvalues(1))/incr)
    !
    if (lonvalues(1).ne.(-180._nr)) then
       do i=1,nlon
          if (i.le.idecal) then
             pv(i,:,:,:)=pvint(i+nlon-idecal,:,:,:)
             if (lonvalues(i+nlon-idecal).lt.180._nr) then
                lonval(i)=lonvalues(i+nlon-idecal)
             else
                lonval(i)=lonvalues(i+nlon-idecal)-360._nr
             end if
          else
             pv(i,:,:,:)=pvint(i-idecal,:,:,:)
             if (lonvalues(i-idecal).lt.180._nr) then
                lonval(i)=lonvalues(i-idecal)
             else
                lonval(i)=lonvalues(i-idecal)-360._nr
             end if
          end if
       end do
    end if
    !
    ! DEFINING THE SPATIAL GRID TO WHICH THE WAVE-BREAKING
    ! DETECTION ALGORITHM WILL BE APPLIED
    !
    ! If choicegrid=1 then INCR2=INCR (we keep all the points of the original grid)
    ! If choicegrid=2 then INCR2=2*INCR (we take one grid point over two)
    ! If choicegrid=3 then INCR2=3*INCR (we take one grid point over three) this is the 
    ! choice made in Michel and Riviere (2012). Starting from a 1.5x1.5 grid (ERA40),
    ! we applied our algorithm to a 4.5x4.5 grid
    choicegrid=1
    !
    nlon2=1  ! One more point to get the periodicity in longitudes 
    do i=1,nlon,choicegrid
       nlon2=nlon2+1
    end do
    nlat2=0
    do j=1,nlat,choicegrid
       nlat2=nlat2+1
    end do
    incr2=incr*choicegrid
    !
    allocate(lonval2(nlon2))
    lonval2=0._nr
    do i=1,nlon2
       lonval2(i)=lonval(1)+incr2*real(i-1)
    end do
    !
    allocate(latval2(nlat2))
    latval2=0._nr
    do j=1,nlat2
       latval2(j)=latval(1)+incr2*real(j-1)
    end do
    ! 
    allocate(pv2(nlon2,nlat2,nplev,ntime))
    do it=1,ntime
       do nn=1,nplev
          do i=1,nlon,choicegrid
             do j=1,nlat,choicegrid     
                if (choicegrid.eq.1) then
                   pv2(i,j,nn,it)=pv(i,j,nn,it)
                else
                   pv2(int(i/choicegrid)+1,int(j/choicegrid)+1,nn,it)=pv(i,j,nn,it)
                end if
             end do
          end do
       end do
    end do
    !
    ! periodicity in longitudes (the first and last longitudes are the same)
    do it=1,ntime
       do nn=1,nplev
          do j=1,nlat,choicegrid
             if (choicegrid.eq.1) then
                pv2(nlon2,j,nn,it)=pv(1,j,nn,it)
             else
                pv2(nlon2,int(j/choicegrid)+1,nn,it)=pv(1,j,nn,it)
             end if
          end do
       end do
    end do
    deallocate(pv)
    !
    ! Computation of the beta_a and beta_c as defined in Riviere (2009) and Riviere et al (2010)
    ! beta_a=1 if the grid point belongs to an anticyclonic breaking (AWB) and 0 otherwise
    ! beta_c=1 if the grid point belongs to a cyclonic breaking (CWB) and 0 otherwise
    !
    allocate(beta_a(nlon2,nlat2,nplev,ntime))
    allocate(beta_c(nlon2,nlat2,nplev,ntime))
    beta_a=0
    beta_c=0
    grandeur='PV'
    !
    allocate(xvalf(8000,nlat2,ncon))
    allocate(yvalf(8000,nlat2,ncon))
    xvalf=0._nr
    yvalf=0._nr
    !print *,'before detection'
    !print *,'variable=',grandeur
    !print *,'nlon2=',nlon2,'nlat2=',nlat2,'ncon=',ncon,' incr2=',incr2
    !print *,'lonval2(1)=',lonval2(1),'latval2(1)=',latval2(1)
    !
    do it=1,ntime
       !print *,'time=', it
       do nn=1,nplev
          xvalf=0._nr
          yvalf=0._nr
          call wb_detection(grandeur,pv2(:,:,nn,it),nlon2,nlat2,ncon,&
                       &    lonval2(1),latval2(1),incr2,lonval2,latval2,&
                       &    xvalf,yvalf,beta_a(:,:,nn,it),beta_c(:,:,nn,it))
       end do
    end do
    !
    ! In the Northern Hemisphere, the anticyclonic wave breaking is part
    ! of the contour which is locally overturned and oriented NE-SW while 
    ! in the Southern Hemisphere it is oriented SE-NW. 
    ! The output fields of the above routine WB_DETECTION beta_a and beta_c 
    ! correspond respectively to NE-SW and SE-NW oriented segments over the whole earth.
    ! Therefore, we have to reverse beta_a and beta_c in the Southern Hemisphere
    ! in order to have beta_a equivalent to an anticyclonic breaking and
    ! beta_c to a cyclonic wave-breaking over the whole earth.
    !
    allocate(beta_cint(nlon2,nlat2,nplev,ntime))
    beta_cint=beta_c
    do j=1,nlat2
       if (latval2(j).lt.0._nr) then
          beta_c(:,j,:,:)=beta_a(:,j,:,:)
          beta_a(:,j,:,:)=beta_cint(:,j,:,:)
       end if
    end do
    !
    ! allocate dynlib output arrays (NOTE!! ASSUMES CHOICEGRID=1!!)
    do it=1,ntime
       nn=1 ! only 1 level
       do j=1,nlat
          do i=1,nlon
             beta_a_out(it,j,i)=beta_a(i,j,nn,it)
             beta_c_out(it,j,i)=beta_c(i,j,nn,it)
          end do
       end do
    end do
    !
    ! If latitudes were reversed for calculation, reverse them back for outputs
    if (incrlat.lt.0._nr) then
       do j=1,nlat
          beta_a_outint(:,j,:)=beta_a_out(:,nlat-j+1,:)
          beta_c_outint(:,j,:)=beta_c_out(:,nlat-j+1,:)
       end do
    end if
    beta_a_out=beta_a_outint
    beta_c_out=beta_c_outint
    !
    ! Note: Need to rearrange lons if starting lon not -180 - not done yet!
    !
    deallocate(lonval) 
    deallocate(latval) 
    deallocate(pvint) 
    deallocate(lonval2) 
    deallocate(latval2) 
    deallocate(pv2) 
    deallocate(beta_a) 
    deallocate(beta_c) 
    deallocate(xvalf) 
    deallocate(yvalf) 
    deallocate(beta_cint) 
  end subroutine ! end of subroutine contour_rwb
  !
  ! Calculates geostrophic velocity [ug,vg] = v_g(mont,lat)
  subroutine v_g(resx,resy,nx,ny,nz,mont,lat,dx,dy)
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
    forall(i = 1_ni:nx, k = 1_ni:nz)
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
