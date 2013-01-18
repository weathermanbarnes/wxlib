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
             res(k,j,1_ni) = 0_nr
             res(k,j,nx  ) = 0_nr
       end forall
    end if
  end subroutine
  
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
          res(k,1_ni,i)=0_nr
          res(k,ny,i)=0_nr
    end forall
  end subroutine
 
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
    resx(:,1_ni,:) = 0_nr
    resy(:,1_ni,:) = 0_nr
    resx(:,ny,:  ) = 0_nr
    resy(:,ny,:  ) = 0_nr
    if (.not.(grid_cyclic_ew)) then
         ! set to 0 where grad result not valid: i=1 and nx
         resx(:,:,1_ni) = 0_nr 
         resy(:,:,1_ni) = 0_nr
         resx(:,:,nx  ) = 0_nr
         resy(:,:,nx  ) = 0_nr
    end if
  end subroutine

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
    end if
    forall(i = 1_ni:nx,k = 1_ni:nz)
        res(k,1_ni,i) = 0_nr
        res(k,ny  ,i) = 0_nr
    end forall
  end subroutine
  
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
  
  ! Calculates angle from x-axis to axis of dilatation :
  !   = delta (Keyser Reeder Reed)
  !   =   (Lapeyre Klein Hua) 
  !   = alpha (Markowski Richardson)
  subroutine dil_angle(res,nx,ny,nz,u,v,dx,dy)
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
             def_shear   = (u(k,j+1_ni,i)-u(k,j-1_ni,i))/dy(j,i) &
                  &      + (v(k,j,i+1_ni)-v(k,j,i-1_ni))/dx(j,i)
             def_stretch = (u(k,j,i+1_ni)-u(k,j,i-1_ni))/dx(j,i) &
                  &      - (v(k,j+1_ni,i)-v(k,j-1_ni,i))/dy(j,i)
             res(k,j,i)  = 0.5_nr*atan2(def_shear, def_stretch)
          end do
       end do
    end do
    if (grid_cyclic_ew) then
       do j=2_ni,ny-1_ni
          do k=1_ni,nz
             def_shear   = (u(k,j+1_ni,i)-u(k,j-1_ni, i))/dy(j,1_ni) &
                     &   + (v(k,j,  2_ni)-v(k,j     ,nx))/dx(j,1_ni)
             def_stretch = (u(k,j,  2_ni)-u(k,j     ,nx))/dx(j,1_ni) &
                     &   - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,1_ni)
             res(k,j,1_ni) = 0.5_nr*atan2(def_shear, def_stretch)
             def_shear   = (u(k,j+1_ni,i)-u(k,j-1_ni, i))/dy(j,nx) &
                     &   + (v(k,j  ,1_ni)-v(k,j,nx-1_ni))/dx(j,nx)
             def_stretch = (u(k,j  ,1_ni)-u(k,j,nx-1_ni))/dx(j,nx) &
                     &   - (v(k,j+1_ni,i)-v(k,j-1_ni, i))/dy(j,nx)
             res(k,j,nx) = 0.5_nr*atan2(def_shear, def_stretch)
          end do
       end do
    end if
  end subroutine
  
  ! Calculates angle from x-axis to iso-pv lines (direction k X gradPV)
  !   = alpha (Keyser Reeder Reed)
  subroutine isopv_angle(res,nx,ny,nz,pv,dx,dy)
    real(kind=nr), intent(in)  :: pv(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: gradx(nz,ny,nx), grady(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    call ddx(gradx,nx,ny,nz,pv,dx,dy)
    call ddy(grady,nx,ny,nz,pv,dx,dy)  
    res=atan2(-gradx,grady)
  end subroutine
  
  ! Calculates beta = angle between dilatation axis (u,v) and iso-lines of pv
  !   = beta (Keyser Reeder Reed)
  !   = beta (Markowski Richardson)
  !   = theta + phi + pi/4 (Lapeyre Klein Hua)
  subroutine beta(res,nx,ny,nz,u,v,pv,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), pv(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: alpha(nz,ny,nx), delta(nz,ny,nx) 
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, v, pv
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    ! beta = delta - alpha (Keyser Reeder Reed) 
    !  alpha = angle between x axis and iso-pv line (direction k X grad pv)
    !  delta = angle between x axis and dilatation axis (dynlib: dil_angle)
    call isopv_angle(alpha,nx,ny,nz,pv,dx,dy)
    call dil_angle(delta,nx,ny,nz,u,v,dx,dy)
    res=delta-alpha
  end subroutine
  
  ! PV FRONTOGENESIS: STRETCHING AND STIRRING RATES:
  ! Calculates (stretch, stir) where:
  ! stretch= fractional pv gradient stretching rate = 1/|gradPV| * d/dt(|gradPV|)
  !        = gamma, 'stretching rate' (Lapeyre Klein Hua)
  !        = -1/|gradPV| * Fn (Keyser Reeder Reed)
  !        =  1/|gradPV| * F  (Markowski Richardson)
  !   stir = angular rotation rate of gradPV (stirring rate)
  !        = d(theta)/dt      (Lapeyre Klein Hua)
  !        = 1/|gradPV| * Fs  (Keyser Reeder Reed)
  subroutine stretch_stir(resstretch,resstir,nx,ny,nz,u,v,pv,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx), v(nz,ny,nx), pv(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: resstretch(nz,ny,nx),resstir(nz,ny,nx)
    real(kind=nr) :: bet(nz,ny,nx),totdef(nz,ny,nx),divergence(nz,ny,nx),vorticity(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) resstretch,resstir,v,pv
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    ! (Lapeyre Klein Hua):
    !  stretch = gamma = 1/rho * d(rho)/dt (rho=|gradPV|) = -sigma/2 * sin(2(theta+phi))
    ! (Keyser Reeder Reed):
    !  d(gradPV)/dt = (Fn,Fs), coordinates rotated to PV-contours (Keyser Reeder Reed) 
    !  stretch = -1/|gradPV| * Fn 
    !       Fn = 0.5*|gradPV|(D-E*cos(2*beta))
    !  stir = Fs/|gradPV|
    !    Fs = 0.5*|gradPV|(vort+E*sin(2*beta))

    call beta(bet,nx,ny,nz,u,v,pv,dx,dy)
    call def_total(totdef,nx,ny,nz,u,v,dx,dy)
    call div(divergence,nx,ny,nz,u,v,dx,dy)
    call vor(vorticity,nx,ny,nz,u,v,dx,dy)
    resstretch=-0.5_nr*(divergence-totdef*cos(2_nr*bet))
    resstir=0.5_nr*(vorticity+totdef*sin(2_nr*bet))
  end subroutine
  
  ! Calculates the geopotential (res) from montgomery potential (m), potential
  ! temperature (theta) and pressure (p)
  subroutine geop_from_montgp(res,nx,ny,nz,m,theta,p,dx,dy)
    use consts, only: cp, Rl, p0
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
  subroutine rev(resa,resc,resai,resci,resaiy,resciy,tested,nx,ny,nz, &
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
    !call high_enough(highenough,nx,ny,nz,zdata,ztest,zthres,dx,dy)

    do i=1_ni,nx
       do j=1_ni,ny
          if (latitudes(j)>0.) then !NH
             hemi=1.
          else if (latitudes(j)<0.) then !SH
             hemi=-1.
          else !equator
             hemi=0.
          endif
          if (abs(latitudes(j))==90.) then 
             edge=1. !at edge of grid, must skip
          elseif ((.NOT.(grid_cyclic_ew)).AND.((i==1).OR.(i==nx))) then
             edge=1. !at edge of grid, must skip
          else
             edge=0.
          end if
          do k=1_ni,nz
            !initialise all to 0, change later if reversal
            tested(k,j,i) = 0 
            resa(k,j,i)   = 0
            resai(k,j,i)  = 0.
            resaiy(k,j,i) = 0.
            resc(k,j,i)   = 0
            resci(k,j,i)  = 0.
            resciy(k,j,i) = 0.
            if ((edge==0.).AND.(highenough(k,j,i)==1)) then ! not on grid edge, and high enough: test it
              tested(k,j,i) = 1
              gradx=xgradient(k,j,i)
              grady=ygradient(k,j,i)
              if (grady<0.) then 
                 neg=1. !y-gradient is negative
              else
                 neg=0.
              end if
              if (grady<(-ddythres)) then  
                  sig=1. ! ygradient is significantly negative 
              else
                  sig=0.
              end if	
              gradmag=sqrt(gradx**2.+grady**2.)
              if (gradmag .NE. 0.) then
                 ! angle exists if gradmag not 0
                 gradang=atan2(grady,gradx) !NB slow! faster to do it by signs
                 if ((gradang+3.141592/2.)*hemi>0.) then ! CW NH or CCW SH -> ACYC
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


  ! Prepare global data for FFT: 
  !
  ! Returns the data extended along complementary meridians (for fft)
  ! For each lon, the reflected (lon+180) is attached below
  ! so that data is periodic in x and y.
  ! NOTE: Input data must be lats -90 to 90!!! and nx must be even
  subroutine prepare_fft(res,nx,ny,nz,thedata,dx,dy)
    real(kind=nr), intent(in)  :: thedata(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    !real(kind=nr), intent(out) :: res(nz,2*ny-2_ni,nx)
    real(kind=nr), intent(out) :: res(nz,2*ny-2,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz,iextra
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    !
    do k=1_ni,nz
      do i=1_ni,nx
        iextra=i+nx/2_ni
        if (iextra.GT.nx) then
          iextra=iextra-nx          
        end if
        do j=2_ni,ny-1_ni
          res(k,j,i) = thedata(k,j,i)
            
          res(k,2*ny-j,i)=thedata(k,j,iextra)
        end do
        res(k,1_ni,i)=thedata(k,1_ni,i)
        res(k,ny,i)=thedata(k,ny,i) 
      end do
    end do
  end subroutine

  ! subroutine sum_kix: calculates sum along k dimension for k values 
  !                     which are flagged in kix vector
  ! returns:
  !    res(ny,nx) - thedata summed over k where kix==1
  !    nres       - sum(kix)
  subroutine sum_kix(res,nres,nx,ny,nz,thedata,kix,dx,dy)
    real(kind=nr), intent(in)  :: thedata(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    integer(kind=ni), intent(in) :: kix(nz)
    real(kind=nr), intent(out) :: res(ny,nx)
    integer(kind=ni), intent(out) :: nres
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny) res, dx, dy
    !f2py depend(nz) kix
    ! -----------------------------------------------------------------
    !
    do i=1_ni,nx
       do j=1_ni,ny
          res(j,i) = 0.
          do k=1_ni,nz
             if (kix(k)==1) then
                res(j,i) = res(j,i) + thedata(k,j,i)
             end if
          end do
       end do
    end do
    nres = 0
    do k=1_ni,nz
       if (kix(k)==1) then 
          nres = nres + 1
       end if
    end do
  end subroutine
  
  !
  ! Calculates sum of thedata(k,j,i) along k where:
  !  ( zdata(k,j,i) > ztest(1,j,i)+test_thres )
  !   zdata(k,j,i) is the geopotential of the isentropic surface.
  !   ztest(1,j,i) is the geopotential of earth's surface
  !   test_thres is minimum geopotential difference from earth's surface
  ! So it only sums values where the surface is sufficiently above ground
  ! OUTPUT ntimes IS NUMBER OF TIMES SUFFICIENTLY ABOVE GROUND AT EACH GRIDPOINT
  ! OUTPUT ntimesall IS TOTAL NUMBER OF TIMES BROADCAST OVER GRID
  subroutine summation(res,ntimes,ntimesall,nx,ny,nz,thedata,zdata,ztest,zthres,dx,dy)
    real(kind=nr), intent(in)  :: thedata(nz,ny,nx), zdata(nz,ny,nx), ztest(1,ny,nx), zthres, dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(ny,nx)
    integer(kind=ni), intent(out) :: ntimes(ny,nx),ntimesall(ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    integer(kind=1) :: high(nz,ny,nx)
    !f2py depend(nx,ny) res, ntimes, ntimesall, dx, dy
    !f2py depend(nx,ny,nz) high
    ! -----------------------------------------------------------------
    !
    call high_enough(high,nx,ny,nz,zdata,ztest,zthres,dx,dy)
    do i=1_ni,nx
       do j=1_ni,ny
          res(j,i) = 0.
          ntimes(j,i) = 0_ni
          ntimesall(j,i) = 0_ni
          do k=1_ni,nz
             if (high(k,j,i)==1) then
                res(j,i) = res(j,i) + thedata(k,j,i)
                ntimes(j,i) = ntimes(j,i)+1_ni
             end if
             ntimesall(j,i) = ntimesall(j,i)+1_ni
          end do
       end do
    end do
  end subroutine

  !
  ! Tests if level data (3-D input array zdata(nz,ny,nx)) geopotential is 
  !  sufficiently above test geopotential (eg, topography) ztest(1,ny,nx) 
  ! Returns 3-D flag array:
  !  = 1 if zdata(t,y,x) > (ztest(1,y,x) + zthres)
  !  = 0 otherwise
  !  test_thres is minimum geopotential height above ztest surface
 
  subroutine high_enough(res,nx,ny,nz,zdata,ztest,zthres,dx,dy)
    real(kind=nr), intent(in)  :: zdata(nz,ny,nx), ztest(1,ny,nx), zthres, dx(ny,nx), dy(ny,nx)
    integer(kind=1), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy, ztest
    ! -----------------------------------------------------------------
    !
    do i=1_ni,nx
       do j=1_ni,ny
          do k=1_ni,nz
             if (zdata(k,j,i).GT.(ztest(1,j,i)+zthres)) then
                res(k,j,i) = 1
             else
                res(k,j,i) = 0
             end if
          end do
       end do
    end do
  end subroutine

!=======================================================================
!=======================================================================
!=======================================================================
!=========RWB detection, adapted from algorithm by Riviere==============
!=======================================================================
!=======================================================================
!=======================================================================
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

 use netcdf
 implicit none

 real(kind=nr), intent(in) :: pv_in(nz,ny,nx),lonvalues(nx),latvalues(ny),dx(ny,nx),dy(ny,nx),lev
 integer(kind=ni), intent(in) :: ncon
 integer(kind=ni), intent(out) :: beta_a_out(nz,ny,nx), beta_c_out(nz,ny,nx)
 integer(kind=ni) :: beta_a_outint(nz,ny,nx), beta_c_outint(nz,ny,nx)
 integer(kind=ni) :: nx,ny,nz
 !f2py depend(nx,ny,nz) beta_a_out, beta_c_out,beta_a_outint, beta_c_outint
 !f2py depend(nx,ny) dx, dy
 !f2py depend(nx) lonvalues
 !f2py depend(ny) latvalues

 integer(kind=ni) :: nlon, nlat, ntime, nplev
 integer(kind=ni) :: i,it,sd,j,l,nn,idecal

 real(kind=nr), dimension(:,:,:,:), allocatable :: pv,pvint
 real(kind=nr), dimension(:), allocatable :: rdate
 real(kind=nr) :: plev(1)

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

 real(kind=nr) :: lonleft,latleft,lonright,latright,xlength,ylength

 nlon = nx
 nlat = ny
 ntime = nz
 nplev = 1
 plev(1) = lev
 incr=lonvalues(2)-lonvalues(1)
 incrlat=latvalues(2)-latvalues(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! CHANGING THE ORDER OF THE VALUES !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! TO GET AN INCREASE IN LATITUDES  !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! AND TO GET -180.0 AS THE FIRST   !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! LONGITUDE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! In the case incrlat<0 that is if the latitudes (latvalues) decrease with the
! corresponding index, we have to reverse the order to orientate from the South
! to the North

  allocate(lonval(nlon))
  allocate(latval(nlat))
  allocate(pv(nlon,nlat,nplev,ntime))
  allocate(pvint(nlon,nlat,nplev,ntime))
 
!Arrange pv into format required for this routine
  DO it=1,ntime
      nn=1 ! only 1 level
      DO j=1,nlat
        DO i=1,nlon
           pv(i,j,nn,it) = pv_in(it,j,i)
        ENDDO
      ENDDO
  ENDDO

  pvint=pv
 
  lonval=lonvalues
  latval=latvalues

  IF (incrlat.lt.0.) THEN
          DO j=1,nlat
          pvint(:,j,:,:)=pv(:,nlat-j+1,:,:)
          latval(j)=latvalues(nlat-j+1)
          ENDDO
  END IF

! In the case, lonvalues(1) different from -180.0, we translate the field
! to have the first longitude equal to -180.0
  pv=pvint
  idecal=int((180.0+lonvalues(1))/incr)

  IF (lonvalues(1).ne.(-180.0)) THEN
          DO i=1,nlon
             IF (i.le.idecal) THEN
                  pv(i,:,:,:)=pvint(i+nlon-idecal,:,:,:)
                  if (lonvalues(i+nlon-idecal).lt.180.0) then
                  lonval(i)=lonvalues(i+nlon-idecal)
                  else
                  lonval(i)=lonvalues(i+nlon-idecal)-360.0
                  end if
             ELSE
                  pv(i,:,:,:)=pvint(i-idecal,:,:,:)
                  if (lonvalues(i-idecal).lt.180.0) then
                  lonval(i)=lonvalues(i-idecal)
                  else
                  lonval(i)=lonvalues(i-idecal)-360.0
                  end if
             ENDIF
          ENDDO
  ENDIF
!  print *,'lonval=',lonval
          


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! DEFINING THE SPATIAL GRID TO WHICH THE WAVE-BREAKING !!!!!!
!!!!!!!!!!!! DETECTION ALGORITHM WILL BE APPLIED  !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      choicegrid=1 ! 
! If choicegrid=1 then INCR2=INCR (we keep all the points of the original grid)
! If choicegrid=2 then INCR2=2*INCR (we take one grid point over two)
! If choicegrid=3 then INCR2=3*INCR (we take one grid point over three) this is the 
! choice made in Michel and Riviere (2012). Starting from a 1.5x1.5 grid (ERA40),
! we applied our algorithm to a 4.5x4.5 grid
  nlon2=1  ! One more point to get the periodicity in longitudes 
  DO i=1,nlon,choicegrid
    nlon2=nlon2+1
  ENDDO
  nlat2=0
  DO j=1,nlat,choicegrid
    nlat2=nlat2+1
  ENDDO
  incr2=incr*choicegrid
  !!print *,'incr2, nlon2, nlat2=',incr2, nlon2,nlat2     

  allocate(lonval2(nlon2))
  lonval2=0.0
  DO i=1,nlon2
    lonval2(i)=lonval(1)+incr2*real(i-1)
  ENDDO
  !!print *,'lonval2=',lonval2
  allocate(latval2(nlat2))
  latval2=0.0
  DO j=1,nlat2
    latval2(j)=latval(1)+incr2*real(j-1)
  ENDDO
  !!print *,'latval2=',latval2
  
  allocate(pv2(nlon2,nlat2,nplev,ntime))
!! print *,nlon,nlat,nplev,ntime
!! print *,nlon2, nlat2, nplev,ntime
  DO it=1,ntime
    DO nn=1,nplev
      DO i=1,nlon,choicegrid
        DO j=1,nlat,choicegrid     
            IF (choicegrid.eq.1) THEN
               pv2(i,j,nn,it)=pv(i,j,nn,it)
            ELSE
               pv2(int(i/choicegrid)+1,int(j/choicegrid)+1,nn,it)=pv(i,j,nn,it)
            ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO

! periodicity in longitudes (the first and last longitudes are the same)
  DO it=1,ntime
    DO nn=1,nplev
      DO j=1,nlat,choicegrid
          IF (choicegrid.eq.1) THEN
             pv2(nlon2,j,nn,it)=pv(1,j,nn,it)
          ELSE
             pv2(nlon2,int(j/choicegrid)+1,nn,it)=pv(1,j,nn,it)
          ENDIF
      ENDDO
    ENDDO
  ENDDO
  deallocate(pv)
!!print *,pv2(1,1,1,1)!!'Min de pv2=',MINVAL(pv2)
!!  print *, 'Max de pv2=',MAXVAL(pv2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Computation of the beta_a and beta_c as defined in !!!!!!!!
!!!!!!!!!!!!! Riviere (2009) and Riviere et al (2010) !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! beta_a=1 if the grid point belongs to an !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! anticyclonic breaking (AWB) and 0 otherwise !!!!!!!!!!!!!!!
!!!!!!!!!!!!! beta_c=1 if the grid point belongs to a !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! cyclonic breaking (CWB) and 0 otherwise !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  allocate(beta_a(nlon2,nlat2,nplev,ntime))
  allocate(beta_c(nlon2,nlat2,nplev,ntime))
  beta_a=0
  beta_c=0
  grandeur='PV'

  allocate(xvalf(8000,nlat2,ncon))
  allocate(yvalf(8000,nlat2,ncon))
  xvalf=0.0
  yvalf=0.0
  print *,'before detection'
  print *,'variable=',grandeur
  print *,'nlon2=',nlon2,'nlat2=',nlat2,'ncon=',ncon,' incr2=',incr2
  print *,'lonval2(1)=',lonval2(1),'latval2(1)=',latval2(1)

  DO it=1,ntime
    print *,'time=', it
    DO nn=1,nplev
      xvalf=0.0
      yvalf=0.0
      call WB_DETECTION(grandeur,pv2(:,:,nn,it),nlon2,nlat2,ncon,&
 &    lonval2(1),latval2(1),incr2,lonval2,latval2,&
 &    xvalf,yvalf,beta_a(:,:,nn,it),beta_c(:,:,nn,it))
!      print *,'Max de beta_a=',MAXVAL(beta_a(:,:,nn,it))
!      print *,'Max de beta_c=',MAXVAL(beta_c(:,:,nn,it))
    ENDDO
!    print *,'fin wb', it
  ENDDO

!  print *,'Min de beta_a=',MINVAL(beta_a)
!  print *,'Max de beta_a=',MAXVAL(beta_a)
!  print *,'Min de beta_c=',MINVAL(beta_c)
!  print *,'Max de beta_c=',MAXVAL(beta_c)

! In the Northern Hemisphere, the anticyclonic wave breaking is part
! of the contour which is locally overturned and oriented NE-SW while 
! in the Southern Hemisphere it is oriented SE-NW. 
! The output fields of the above routine WB_DETECTION beta_a and beta_c 
! correspond respectively to NE-SW and SE-NW oriented segments over the whole earth.
! Therefore, we have to reverse beta_a and beta_c in the Southern Hemisphere
! in order to have beta_a equivalent to an anticyclonic breaking and
! beta_c to a cyclonic wave-breaking over the whole earth.

  allocate(beta_cint(nlon2,nlat2,nplev,ntime))
  beta_cint=beta_c
  Do j=1,nlat2
  if (latval2(j).lt.0.) then
  beta_c(:,j,:,:)=beta_a(:,j,:,:)
  beta_a(:,j,:,:)=beta_cint(:,j,:,:)
  End if
 Enddo

!allocate dynlib output arrays (NOTE!! ASSUMES CHOICEGRID=1!!)
 DO it=1,ntime
      nn=1 ! only 1 level
      DO j=1,nlat
        DO i=1,nlon
           beta_a_out(it,j,i)=beta_a(i,j,nn,it)
           beta_c_out(it,j,i)=beta_c(i,j,nn,it)
        ENDDO
      ENDDO
  ENDDO

! If latitudes were reversed for calculation, reverse them back for outputs
  IF (incrlat.lt.0.) THEN
          DO j=1,nlat
            !pvint(:,j,:,:)=pv(:,nlat-j+1,:,:) !leave these reversed- not needed
            !latval(j)=latvalues(nlat-j+1) !leave these reversed- not needed
            beta_a_outint(:,j,:)=beta_a_out(:,nlat-j+1,:)
            beta_c_outint(:,j,:)=beta_c_out(:,nlat-j+1,:)
          ENDDO
  END IF
beta_a_out=beta_a_outint
beta_c_out=beta_c_outint

!! Note !! Need to rearrange lons if starting lon not -180 - not done yet!

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

!======================================================================
!======================================================================
!======================================================================
!======================================================================
!======================================================================
! This section contains all subroutines needed to detect RWBs
!======================================================================
! This subroutine is that called in contour_rwb subroutine to detect RWBs
        SUBROUTINE WB_DETECTION(grd,PV2d,NX,NY,ncon,LONMIN,LATMIN,INCR,&
        lon,lat,XVALF,YVALF,beta_a2d,beta_c2d)
! The input variables of the subroutine (given at the subroutine)
! grd= PV (for potential vorticity) or VO (for absolute vorticity)
        CHARACTER(LEN=2),INTENT(IN) :: grd
! NX and NY are the numbers of longitude and latitude
! ncon is the number of contours we want to detect (generally =21)
        integer(kind=ni),INTENT(IN) :: NX,NY,ncon
! ndim is the highest number of points that can define a contour
        integer(kind=ni),PARAMETER :: ndim=8000
! lonmin is the first value of longitude (generally 180°W)
! latmin is the first value of latitude (either 90°S or 0°)
! incr is the grid step
        real(kind=nr),INTENT(IN) :: LONMIN,LATMIN,INCR
! lon is the table containing all longitudes (from 180°W to 180°E)
        real(kind=nr),DIMENSION(NX),INTENT(IN) :: lon
! lat is the table containing all latitudes
        real(kind=nr),DIMENSION(NY),INTENT(IN) :: lat
! PV2d is the initial field on which we want to detect RWBs (the first
! and last points must be equal because it is the same longitude 180W
! and 180E)
        real(kind=nr),DIMENSION(NX,NY),INTENT(IN) :: PV2d
! XVALF and YVALF are the longitudes and latitudes defining the contours
! It depends on the latitude NY and on the number of the contour
! detected ncon
        real(kind=nr),DIMENSION(ndim,NY,ncon),INTENT(INOUT) :: XVALF,YVALF
! beta_a2d is a table of integers which are equal to 1 when there is an
! AWB and 0 when there is not
! beta_c2d is the same table but for CWB
        integer(kind=ni),DIMENSION(NX,NY),INTENT(INOUT) :: beta_a2d,beta_c2d

        integer(kind=ni) :: i,j,pp,p,dd,iii,jjj
        integer(kind=ni) :: numc,icompt,ee,nbptx,nbc,icomptprecedent
        integer(kind=ni) :: il,jl,ilongl,ilongu,jlatl,jlatu
        integer(kind=ni) :: nbpta,nbptc,iv_err
        integer(kind=ni) :: nbpta_temp,nbptc_temp
        integer(kind=ni),PARAMETER :: npt=8
        real(kind=nr),DIMENSION(ncon) :: zz
        real(kind=nr),DIMENSION(ndim) :: XVAL,YVAL
        real(kind=nr),DIMENSION(npt) :: xseg,yseg
        integer(kind=ni),DIMENSION(NY) :: flag
        integer(kind=ni),DIMENSION(NX,NY) :: beta_a_temp,beta_c_temp
        LOGICAL :: LRead
        integer(kind=ni) :: sd1
        integer(kind=ni) :: nc1,sd

 !f2py depend(NX) lon
 !f2py depend(NY) lat,flag
 !f2py depend(NX,NY) PV2d,beta_a2d,beta_c2d,beta_a_temp,beta_c_temp
 !f2py depend(ndim,NY,ncon) XVALF, YVALF
 !f2py depend(ncon) zz
 !f2py depend(ndim) XVAL,YVAL
 !f2py depend(npt) xseg,yseg

! Initialisations
        zz=0.0
        xseg=0.0
        yseg=0.0
        XVAL=0.0
        YVAL=0.0
        flag=0
        beta_a_temp=0.0
        beta_c_temp=0.0
        XVALF=0.0
        YVALF=0.0
        sd1=0
        sd=0
        
! zz is the table containing contours values (with ncon=21 or 41)
! For PV: from 0 pvu to 10 pvu with an interval of 0.5 pvu
! For VO (absolute vorticity): form -20e-5 to 20e-5 with an interval of 2e-5 s-1
        nc1=int(ncon/2)+1
!        print '(A4,I2)','nc1=',nc1
!        print '(A3,I2)','ncon=',ncon
        DO pp=1,ncon
! if the initial field studied is PV on the northern hemisphere
          if (grd.eq.'PV') then
            zz(pp)=0.5e-6*(real(pp-nc1))
          endif
! if the initial field studied is the absolute vorticity on the sphere
          if (grd.eq.'VO') then
            zz(pp)=2.0e-5*(real(pp-nc1))
          endif
        ENDDO
!        print *,'MAXVAL(PV2d)=',MAXVAL(PV2d)

! Scheme of a box with the coordinates of the summits
!     j+1      j+1
!     i -------- i+1
!       |      |
!       |      |
!     i -------- i+1
!       j      j
!
! Initialisations
        nbc=0
        nbpta=0
        nbptc=0

!        print *,'NX=',NX
!        print *,'NY=',NY
!        print *,'lat(NY)=',lat(NY)
!        print *,'lat(1)=',lat(1)
        
! Loop on the contours' values
        DO p=1,ncon
!          print *,'zz(p)=',zz(p)
          i=1
          flag=1
          numc=1
! Loop on the latitudes
          DO j=1,NY-1
            icompt=0
            XVAL=0.0
            YVAL=0.0
            xseg=0.0
            yseg=0.0
            ee=0
            nbptx=0
            if (flag(j).ne.2) then
! conrec is the subroutine that detects segments of the contours in the
! box studied
              call conrec(PV2d(i:i+1,j:j+1),lon(i:i+1),&
              lat(j:j+1),zz(p),xseg,yseg,npt,nbc)
              if ((nbc.ne.0).and.(xseg(nbc).ne.(-180.0)).and.&
              (xseg(nbc-1).ne.(-180.0))) CYCLE
              if (nbc.ne.0) then
                flag(j)=2
! orgseg is the subroutine that orders the points found by conrec in the
! west-east direction
                call orgseg(xseg(1:nbc),yseg(1:nbc),XVAL,YVAL,nbc,ndim,icompt)
                if (YVAL(icompt).eq.lat(1)) then
                  ee=1
                  EXIT
                endif
                XVALF(1:icompt,numc,p)=XVAL(1:icompt)
                YVALF(1:icompt,numc,p)=YVAL(1:icompt)
                il=i
                jl=j
                LRead=.TRUE.
                if (sd1.eq.0) then
                DO WHILE (LRead)
                  icomptprecedent=icompt
! boitesuiv is the subroutine that defines the coordinate of the next
! box
                  call boitesuiv(il,jl,lon(il:il+1),lat(jl:jl+1),NX-1,NY,&
                  XVAL(icompt),YVAL(icompt),ilongl,ilongu,jlatl,jlatu)
                  xseg=0.0
                  yseg=0.0
! Again we use conrec to find segments in this new box
                  call conrec(PV2d(ilongl:ilongu,jlatl:jlatu),&
                  lon(ilongl:ilongl+1),lat(jlatl:jlatl+1),zz(p),xseg,yseg,npt,nbc)
                  if (nbc.eq.0) then
                    LRead=.FALSE.
                  endif
                  if (ilongl.eq.(il+1)) then
                    nbptx=nbptx+1
                  else
                    if (ilongl.eq.(il-1)) then
                      nbptx=nbptx-1
                    else
                      nbptx=nbptx
                    endif
                  endif
                  il=ilongl
                  jl=jlatl
                  if (nbc.ne.0) then
                    DO dd=1,nbc
                      if ((xseg(dd).eq.(-180.0)).or.(ilongl.eq.1)) then
                        flag(jl)=2
                      endif
                    ENDDO
                  endif
                  if (nbc.ne.0) then
! Again the orgseg subroutine order the pints detected by conrec
                    call orgseg(xseg(1:nbc),yseg(1:nbc),XVAL,YVAL,nbc,ndim,icompt)
                    if (icompt.eq.icomptprecedent) then
                      XVAL=0.0
                      LRead=.FALSE.
                      CYCLE
                    endif
                    XVALF(icomptprecedent+1:icompt,numc,p)=XVAL(icomptprecedent+1:icompt)
                    YVALF(icomptprecedent+1:icompt,numc,p)=YVAL(icomptprecedent+1:icompt)
                  else
                    LRead=.FALSE.
                  endif
                  if ((XVALF(1,numc,p).eq.(-XVALF(icompt,numc,p))).and.&
                  (YVALF(1,numc,p).eq.YVALF(icompt,numc,p))) then
                    LRead=.FALSE.
                  endif
                  DO dd=icomptprecedent+1,icompt
                    if (YVALF(dd,numc,p).ge.lat(NY)) then
                      LRead=.FALSE.
                      ee=1
                    endif
                    if (YVALF(dd,numc,p).le.lat(1)) then
                      LRead=.FALSE.
                      ee=1
                    endif
                  ENDDO
                  if (icompt.ge.6990) then
                    LRead=.FALSE.
                    ee=1
                  endif
                ENDDO
! End of the DO WHILE loop
                if (ee.eq.1) then
                  XVAL=0.0
                  YVAL=0.0
                  XVALF(:,numc,p)=0.0
                  YVALF(:,numc,p)=0.0
                  numc=numc-1
                else
                  if (nbptx.eq.0) then
                    XVAL=0.0
                    YVAL=0.0
                    XVALF(:,numc,p)=0.0
                    YVALF(:,numc,p)=0.0
                    numc=numc-1
                  endif
                endif
                numc=numc+1
                endif
              endif
            endif
            if (sd.eq.0) then
            if (XVAL(1).ne.0.0) then
! Once the whole contour is detected, the eventfunctions subroutine
! finds the WB areas and fills in the two beta tables
              call eventfunctions(icompt,NX,NY,LONMIN,LATMIN,INCR,&
              XVAL(1:icompt),YVAL(1:icompt),beta_a_temp(1:NX,1:NY),&
              beta_c_temp(1:NX,1:NY),nbpta_temp,nbptc_temp)
              nbpta=nbpta_temp
              nbptc=nbptc_temp
              DO iii=1,NX
                DO jjj=1,NY
                  if (beta_a_temp(iii,jjj).ne.0) then
                    beta_a2d(iii,jjj)=beta_a_temp(iii,jjj)
                  endif
                  if (beta_c_temp(iii,jjj).ne.0) then
                    beta_c2d(iii,jjj)=beta_c_temp(iii,jjj)
                  endif
                ENDDO
              ENDDO
            endif
            endif
! End of the latitudes loop
          ENDDO
          if (p.eq.ncon) then
!          print *,'fin contours'
          endif
! End of the loop on the contours we want to detect
        ENDDO
!        print *,'Fin WB_DETECTION'
!        If (MAXVAL(beta_a2d).ne.0) THEN
!        print *,'Max de beta_a=',MAXVAL(beta_a2d)
!        print *,'Max de beta_c=',MAXVAL(beta_c2d)
!        print *,'Max de XVALF=',MAXVAL(XVALF)
!        print *,'Max de YVALF=',MAXVAL(YVALF)
!        Endif

        END SUBROUTINE

!======================================================================
! This subroutine finds segments of the contour in a rectangle box 
! divided in four triangles
        SUBROUTINE conrec(d,x,y,z,XVAL,YVAL,ndim,nbc)

        real(kind=nr),DIMENSION(2,2),INTENT(IN) :: d
        real(kind=nr),DIMENSION(2),INTENT(IN) :: x,y
        integer(kind=ni),INTENT(IN) :: ndim
        integer(kind=ni) :: i,j,k,m,m1,m2,m3
        real(kind=nr),INTENT(IN) :: z
        real(kind=nr),DIMENSION(ndim),INTENT(OUT) :: XVAL,YVAL
        integer(kind=ni),INTENT(INOUT) :: nbc
        real(kind=nr) x1,y1,x2,y2
!
!     Local declarations
!
        real(kind=nr) h(0:4)
        integer(kind=ni) sh(0:4)
        real(kind=nr) xh(0:4),yh(0:4)
        integer(kind=ni) im(1:4),jm(1:4)
        integer(kind=ni) case
        integer(kind=ni) castab(-1:1,-1:1,-1:1)
        integer(kind=ni) p1,p2,ii
        real(kind=nr) l1,l3,pi,incr
        real(kind=nr) l11,l31,l41,l12,l32,l42
        real(kind=nr) dmin,dmax
        data im/0,1,1,0/
        data jm/0,0,1,1/
        data castab/0,0,9,  0,1,5,  7,4,8, &
       &            0,3,6,  2,3,2,  6,3,0, &
       &            8,4,7,  5,1,0,  9,0,0/


        pi=atan(1.)*4
        ii=0
        i=1
        j=1

        dmin = min(d(i,j),d(i,j+1),d(i+1,j),d(i+1,j+1))
        dmax = max(d(i,j),d(i,j+1),d(i+1,j),d(i+1,j+1))

        if ((dmax.ge.z).and.(dmin.le.z)) then
          do 22 m=4,0,-1
            if (m.gt.0) then
              h(m)=d(i+im(m),j+jm(m))-z
              xh(m)=x(i+im(m))
              yh(m)=y(j+jm(m))
            else
              l1=sqrt(((cos(y(j)*pi/180.)*x(i+1)-cos(y(j)*pi/180.)*x(i))/2)**2+((y(j+1)-y(j))/2)**2)
              l3=sqrt(((cos(y(j+1)*pi/180.)*x(i+1)-cos(y(j+1)*pi/180.)*x(i))/2)**2+((y(j+1)-y(j))/2)**2)
              h(0)=(h(1)/l1+h(2)/l1+h(3)/l3+h(4)/l3)/(2/l1+2/l3)
              xh(0)=0.5*(x(i)+x(i+1))
              yh(0)=0.5*(y(j)+y(j+1))
            endif
            if (h(m).gt.0.0) then
              sh(m)=+1
            else if (h(m).lt.0.0) then
              sh(m)=-1
            else
              sh(m)=0
            endif
22        continue
          do 60 m=1,4
            m1=m
            m2=0
            if (m.ne.4) then
              m3=m+1
            else
              m3=1
            endif
            case = castab(sh(m1),sh(m2),sh(m3))
            if (case.ne.0) then
              goto (31,32,33,34,35,36,37,38,39),case
!
!     Case 1 - Line between vertices 1 and 2
!
31            x1=xh(m1)
              y1=yh(m1)
              x2=xh(m2)
              y2=yh(m2)
              goto 40
!
!     Case 2 - Line between vertices 2 and 3
!
32            x1=xh(m2)
              y1=yh(m2)
              x2=xh(m3)
              y2=yh(m3)
              goto 40
!
!     Case 3 - Line between vertices 3 and 1
!
33            x1=xh(m3)
              y1=yh(m3)
              x2=xh(m1)
              y2=yh(m1)
              goto 40
!
!     Case 4 - Line between vertex 1 and side 2-3
!
34            x1=xh(m1)
              y1=yh(m1)
              call xsectsub(x2,h(m2),h(m3),xh(m2),xh(m3))
              call ysectsub(y2,h(m2),h(m3),yh(m2),yh(m3))
              goto 40
!
!     Case 5 - Line between vertex 2 and side 3-1
!
35            x1=xh(m2)
              y1=yh(m2)
              call xsectsub(x2,h(m3),h(m1),xh(m3),xh(m1))
              call ysectsub(y2,h(m3),h(m1),yh(m3),yh(m1))
              goto 40
!
!     Case 6 - Line between vertex 3 and side 1-2
!
36            x1=xh(m3)
              y1=yh(m3)
              call xsectsub(x2,h(m1),h(m2),xh(m1),xh(m2))
              call ysectsub(y2,h(m1),h(m2),yh(m1),yh(m2))
              goto 40
!
!     Case 7 - Line between sides 1-2 and 2-3
!
37            call xsectsub(x1,h(m1),h(m2),xh(m1),xh(m2))
              call ysectsub(y1,h(m1),h(m2),yh(m1),yh(m2))
              call xsectsub(x2,h(m2),h(m3),xh(m2),xh(m3))
              call ysectsub(y2,h(m2),h(m3),yh(m2),yh(m3))
              goto 40
!
!     Case 8 - Line between sides 2-3 and 3-1
!
38            call xsectsub(x1,h(m2),h(m3),xh(m2),xh(m3))
              call ysectsub(y1,h(m2),h(m3),yh(m2),yh(m3))
              call xsectsub(x2,h(m3),h(m1),xh(m3),xh(m1))
              call ysectsub(y2,h(m3),h(m1),yh(m3),yh(m1))
              goto 40
!
!     Case 9 - Line between sides 3-1 and 1-2
!
39            call xsectsub(x1,h(m3),h(m1),xh(m3),xh(m1))
              call ysectsub(y1,h(m3),h(m1),yh(m3),yh(m1))
              call xsectsub(x2,h(m1),h(m2),xh(m1),xh(m2))
              call ysectsub(y2,h(m1),h(m2),yh(m1),yh(m2))
40            continue
              XVAL(ii+1)=x1
              YVAL(ii+1)=y1
              XVAL(ii+2)=x2
              YVAL(ii+2)=y2
              ii=ii+2
            endif
60        continue
        endif
        nbc=ii

        end subroutine

!=====================================================================
! Function used in the conrec subroutine
 subroutine xsectsub(xsect,h1,h2,xh1,xh2)
        !from FUNCTION xsect(h1,h2,xh1,xh2)
        real(kind=nr), intent(out) :: xsect
        real(kind=nr), intent(in) :: h1,h2,xh1,xh2

        if (xh1.eq.xh2) then
          xsect=xh1
        else
          xsect=(h2*xh1-h1*xh2)/(h2-h1)
        endif

        END SUBROUTINE

!=============================================================
! Function used in the conrec subroutine
 subroutine ysectsub(ysect,h1,h2,yh1,yh2)
        !from FUNCTION ysect(h1,h2,yh1,yh2)
        real(kind=nr), intent(out) :: ysect
        real(kind=nr), intent(in) :: h1,h2,yh1,yh2

        if (yh1.eq.yh2) then
          ysect=yh1
        else
          ysect=(h2*yh1-h1*yh2)/(h2-h1)
        endif

        END SUBROUTINE

!===========================================================================
! orgseg is subroutine that orders the points, defined by their longitude and 
! latitude, found by conrec in the west-east direction
        SUBROUTINE orgseg(xseg,yseg,xxf,yyf,nbc,ndim,icompt)

        IMPLICIT NONE
        integer(kind=ni),INTENT(IN) :: ndim,nbc
        integer(kind=ni),INTENT(INOUT) :: icompt
        real(kind=nr),DIMENSION(nbc),INTENT(IN) :: xseg
        real(kind=nr),DIMENSION(nbc),INTENT(IN) :: yseg
        real(kind=nr),DIMENSION(ndim),INTENT(INOUT) :: xxf,yyf
        integer(kind=ni) kk,ll,rr,ptsauv,reste,reste1,llsauv,ptf,rrsuiv,ptsauvb
        integer(kind=ni) faux,icnew,icomptd
 !f2py depend(nbc) xseg,yseg
 !f2py depend(ndim) xxf,yyf

        if (icompt.eq.0) then
          DO ll=1,nbc
            if (xseg(ll).eq.(-180.0)) then
              xxf(1)=xseg(ll)
              yyf(1)=yseg(ll)
              ptsauv=ll
            endif        
          ENDDO
          call parite(ptsauv,ptsauvb)
          xxf(2)=xseg(ptsauvb)
          yyf(2)=yseg(ptsauvb)
          ptf=nbc/2+1
          if (ptf.gt.2) then
            boucleorgpt : &
            DO kk=3,ptf
              faux=0
              DO ll=1,nbc
                if ((xseg(ll).eq.xxf(kk-1)).and.&
               (yseg(ll).eq.YYF(kk-1)).and.(ll.ne.ptsauvb)) then
                  call parite(ll,llsauv)
                  xxf(kk)=xseg(llsauv)
                  yyf(kk)=yseg(llsauv)
                  icnew=kk
                else
                  faux=faux+1
                endif
              ENDDO
              ptsauvb=llsauv
              if (faux.ge.nbc) then
                EXIT boucleorgpt
              endif
            ENDDO boucleorgpt
          endif
          if (faux.ge.nbc) then
            icompt=icnew
          else
            icompt=icompt+ptf
          endif
        else
          ptsauv=0
          DO rr=1,nbc
            if (((xseg(rr).eq.xxf(icompt)).or.&
            (xseg(rr).eq.(-xxf(icompt)))).and.(yseg(rr).eq.yyf(icompt))) then
              call parite(rr,rrsuiv)
              xxf(icompt+1)=xseg(rrsuiv)
              yyf(icompt+1)=yseg(rrsuiv)
              ptsauv=rrsuiv
            endif
          ENDDO
          if (ptsauv.eq.0) then
            goto 111
          endif
          ptf=nbc/2+1
          if (ptf.gt.2) then
            boucleorgpts : &
            DO kk=icompt+2,icompt+ptf-1
              faux=0
              DO ll=1,nbc
                if ((xseg(ll).eq.xxf(kk-1)).and.&
                (yseg(ll).eq.YYF(kk-1)).and.(ll.ne.ptsauv)) then
                  call parite(ll,llsauv)
                  xxf(kk)=xseg(llsauv)
                  yyf(kk)=yseg(llsauv)
                  icnew=kk
                else
                  faux=faux+1
                endif
              ENDDO
              ptsauv=llsauv
              if (faux.ge.nbc) then
                EXIT boucleorgpts
              endif
            ENDDO boucleorgpts
          endif
          if (faux.ge.nbc) then      
            icompt=icnew
          else
            icompt=icompt+ptf-1
          endif
        endif
   
111     CONTINUE

        
        end subroutine

! ===============================================================================================
! PARITE is a subroutine checking if the number is even or odd
! It is used in the orgseg subroutine
        SUBROUTINE parite(nombrea,nombreb)

        IMPLICIT NONE
        integer(kind=ni),INTENT(IN) :: nombrea
        integer(kind=ni),INTENT(OUT) :: nombreb
        integer(kind=ni) :: reste

        reste=0
        reste=mod(nombrea,2)
        if (reste.eq.1) then
          nombreb=nombrea+1
        else
          nombreb=nombrea-1
        endif

        
        end subroutine


! =============================================================================================
! boitesuiv is a subroutine finding the following box where the contour is.
! It uses the last point found of the contour.
        SUBROUTINE boitesuiv(i,j,x,y,NXX,NYY,XVALF,YVALF,ilongl,ilongu,jlatl,jlatu)

        IMPLICIT NONE
        integer(kind=ni),INTENT(IN) :: i,j
        integer(kind=ni) :: inew
        integer(kind=ni),INTENT(OUT) :: ilongl,ilongu,jlatl,jlatu
        integer(kind=ni),INTENT(IN) :: NXX,NYY
        real(kind=nr),DIMENSION(2),INTENT(IN) :: x,y
        real(kind=nr),INTENT(IN) :: XVALF,YVALF
        integer(kind=ni) :: a,b


! 8 cases are possible, 0 is the initial box whose vertices are (i,j),(i+1,j),(i,j+1),(i+1,j+1)
!   |2|8|4|
!   |3|0|6|
!   |1|7|5|
!
        a=1
        b=1
! Case 1
        if ((XVALF.eq.x(a)).and.(YVALF.eq.y(b))) then
          if (i.eq.1) then
            inew=NXX+1
          else
            inew=i
          endif
          ilongl=inew-1
          ilongu=ilongl+1
          jlatl=j-1
          jlatu=jlatl+1
        else
! Case 2
          if ((XVALF.eq.x(a)).and.(YVALF.eq.y(b+1))) then
            if (i.eq.1) then
              inew=NXX+1
            else
              inew=i
            endif
            ilongl=inew-1
            ilongu=ilongl+1
            jlatl=j+1
            jlatu=jlatl+1
          else
! Case 3
            if (XVALF.eq.x(a)) then      
              if (i.eq.1) then
                inew=NXX+1
              else
                inew=i
              endif
              ilongl=inew-1
              ilongu=ilongl+1
              jlatl=j
              jlatu=jlatl+1
            else
! Case 4
              if ((XVALF.eq.x(a+1)).and.(YVALF.eq.y(b+1))) then
                if ((i+1).eq.(NXX+1)) then
                  inew=0
                else
                  inew=i
                endif
                ilongl=inew+1
                ilongu=ilongl+1
                jlatl=j+1
                jlatu=jlatl+1
              else
! Case 5
                if ((XVALF.eq.x(a+1)).and.(YVALF.eq.y(b))) then
                  if ((i+1).eq.(NXX+1)) then
                    inew=0
                  else
                    inew=i
                  endif
                  ilongl=inew+1
                  ilongu=ilongl+1
                  jlatl=j-1
                  jlatu=jlatl+1
                else
! Case 6
                  if (XVALF.eq.x(a+1)) then
                    if ((i+1).eq.(NXX+1)) then
                      inew=0
                    else
                      inew=i
                    endif
                    ilongl=inew+1
                    ilongu=ilongl+1
                    jlatl=j
                    jlatu=jlatl+1
                  else
! Case 7
                    if (YVALF.eq.y(b)) then
                      ilongl=i
                      ilongu=ilongl+1
                      jlatl=j-1
                      jlatu=jlatl+1
                    else
! Case 8
                      if (YVALF.eq.y(b+1)) then
                        ilongl=i
                        ilongu=ilongl+1
                        jlatl=j+1
                        jlatu=jlatl+1
                      endif
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif

        
        end subroutine
! ______________________________________________________________________
! ______________________________________________________________________
! This subroutine determines the points of a contour that belong to an
! AWB or to a CWB
! beta_a_temp is 1 when an AWB is detected at a grid point and 0 otherwise
! beta_c_temp is 1 when a CWB is detected at a grid point and 0 otherwise
        SUBROUTINE eventfunctions(icompt,NX,NY,longmin,latimin,grille, &
        XVAL,YVAL,beta_a_temp,beta_c_temp,nbpta_temp,nbptc_temp)

        integer(kind=ni),INTENT(IN) :: icompt,NX,NY
        real(kind=nr),INTENT(IN) :: longmin,latimin,grille
        real(kind=nr),DIMENSION(icompt),INTENT(IN) :: XVAL
        real(kind=nr),DIMENSION(icompt),INTENT(IN) :: YVAL
        integer(kind=ni),DIMENSION(NX,NY),INTENT(OUT) ::beta_a_temp,beta_c_temp
        integer(kind=ni),DIMENSION(icompt) :: flag1
        integer(kind=ni),INTENT(OUT) :: nbpta_temp,nbptc_temp
        integer(kind=ni) :: aa,bb
        integer(kind=ni) :: kk,ww,ptdpremier,indi,indj
        real(kind=nr) :: diffm1,diffp1
 !f2py depend(icompt) XVAL, YVAL, flag1
 !f2py depend(NX,NY) beta_a_temp,beta_c_temp 

        flag1=0
        aa=0
        bb=0
        beta_a_temp=0
        beta_c_temp=0
        nbpta_temp=0
        nbptc_temp=0
        DO kk=2,icompt-1
          if (flag1(kk).ne.2) then
            diffm1=XVAL(kk)-XVAL(kk-1)
            diffp1=XVAL(kk+1)-XVAL(kk)
            if ((XVAL(kk).eq.-180.0).and.(diffp1.gt.180.0)) then
              diffp1=diffp1-360.0
            endif
            if ((XVAL(kk-1).eq.-180.0).and.(diffm1.gt.180.0)) then
              diffm1=diffm1-360.0
            endif
            if ((diffm1.lt.0.0).and.(diffp1.lt.0.0)) then
              if (XVAL(kk-2).lt.XVAL(kk-1)) then
                ptdpremier=kk
                flag1(kk)=2
                if (YVAL(ptdpremier-1).le.YVAL(ptdpremier)) then
                  indi=int((XVAL(ptdpremier)-longmin)/grille)+1
                  indj=int((YVAL(ptdpremier)-latimin)/grille)+1
                  if ((indi.gt.NX).or.(indj.gt.NY)) then
                    print *,'indi ou indj depassent NX ou NY'
                    STOP
                  endif
                  beta_c_temp(indi,indj)=1
                  aa=aa+1
                  ww=ptdpremier+1
                  diffm1=XVAL(ww)-XVAL(ww-1)
                  diffp1=XVAL(ww+1)-XVAL(ww)
                  if ((XVAL(ww).eq.-180.0).and.(diffp1.gt.180.0)) then
                    diffp1=diffp1-360.0
                  endif
                  if ((XVAL(ww-1).eq.-180.0).and.(diffm1.gt.180.0)) then
                    diffm1=diffm1-360.0
                  endif
                  DO WHILE((diffm1.lt.0.0).and.(diffp1.lt.0.0))
                    aa=aa+1
                    indi=int((XVAL(ww)-longmin)/grille)+1
                    indj=int((YVAL(ww)-latimin)/grille)+1
                    if ((indi.gt.NX).or.(indj.gt.NY)) then
                      STOP
                    endif
                    beta_c_temp(indi,indj)=1
                    flag1(ww)=2
                    ww=ww+1
                    diffm1=XVAL(ww)-XVAL(ww-1)
                    diffp1=XVAL(ww+1)-XVAL(ww)
                    if ((XVAL(ww).eq.-180.0).and.(diffp1.gt.180.0)) then
                      diffp1=diffp1-360.0
                    endif
                    if ((XVAL(ww-1).eq.-180.0).and.(diffm1.gt.180.0)) then
                      diffm1=diffm1-360.0
                    endif
                  ENDDO
                else
                  if (YVAL(ptdpremier-1).gt.YVAL(ptdpremier)) then
                    indi=int((XVAL(ptdpremier)-longmin)/grille)+1
                    indj=int((YVAL(ptdpremier)-latimin)/grille)+1
                    if ((indi.gt.NX).or.(indj.gt.NY)) then
                      STOP
                    endif
                    beta_a_temp(indi,indj)=1
                    bb=bb+1
                    ww=ptdpremier+1
                    diffm1=XVAL(ww)-XVAL(ww-1)
                    diffp1=XVAL(ww+1)-XVAL(ww)
                    if ((XVAL(ww).eq.-180.0).and.(diffp1.gt.180.0)) then
                      diffp1=diffp1-360.0
                    endif
                    if ((XVAL(ww-1).eq.-180.0).and.(diffm1.gt.180.0)) then
                      diffm1=diffm1-360.0
                    endif
                    DO WHILE((diffm1.lt.0.0).and.(diffp1.lt.0.0))
                      bb=bb+1
                      indi=int((XVAL(ww)-longmin)/grille)+1
                      indj=int((YVAL(ww)-latimin)/grille)+1
                      beta_a_temp(indi,indj)=1
                      flag1(ww)=2
                      ww=ww+1
                      diffm1=XVAL(ww)-XVAL(ww-1)
                      diffp1=XVAL(ww+1)-XVAL(ww)
                      if ((XVAL(ww).eq.-180.0).and.(diffp1.gt.180.0)) then
                        diffp1=diffp1-360.0
                      endif
                      if ((XVAL(ww-1).eq.-180.0).and.(diffm1.gt.180.0)) then
                        diffm1=diffm1-360.0
                      endif
                    ENDDO
                  endif
                endif
              endif
            endif
          endif
        ENDDO
        nbptc_temp=aa
        nbpta_temp=bb

        
        END SUBROUTINE
!======================================================================
!======================================================================
!======================================================================
!======================================================================
!======================================================================
! end of this section
!======================================================================


  ! Calculates geostrophic velocity [ug,vg] = v_g(mont,lat)
  subroutine v_g(resx,resy,nx,ny,nz,mont,lat,dx,dy)
    real(kind=nr), intent(in)  :: mont(nz,ny,nx),lat(ny),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: resx(nz,ny,nx), resy(nz,ny,nx)
    real(kind=nr) :: montx(nz,ny,nx), monty(nz,ny,nx),f(nz,ny,nx)
    integer(kind=ni) :: i,j,k, nx,ny,nz
    real(kind=nr), parameter :: Omega = 0.00007292, pi=3.141592

    !f2py depend(nx,ny,nz) resx, resy
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    ! -----------------------------------------------------------------
    !
    call grad(montx,monty,nx,ny,nz,mont,dx,dy)
    forall(i = 1_ni:nx, k = 1_ni:nz)
       f(k,:,i) = 2*Omega*sin(lat*pi/180.)
    end forall
    where (f==0.) f=9E99 !avoid singularity calculating v_g at equator
    resx=-monty/f
    resy=montx/f
  end subroutine

  !Calculates Okubo-Weiss criterion lambda_0=1/4 (sigma^2-omega^2)= 1/4 W
  !This is the square of the eigenvalues in Okubo's paper (assumes small divergence)
  subroutine okuboweiss(res,nx,ny,nz,u,v,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx),v(nz,ny,nx),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: sig(nz,ny,nx), omega(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res,v
    !f2py depend(nx,ny) dx, dy
    call def_total(sig,nx,ny,nz,u,v,dx,dy)
    call vor(omega,nx,ny,nz,u,v,dx,dy)
    res=0.25*(sig**2-omega**2)
    print *,'lam0(1,10,10)=',res(1,10,10)
  end subroutine

  !Calculates Lagrangian acceleration on the isentropic surface
  subroutine laccel(resx,resy,nx,ny,nz,u,v,mont,lat,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx),v(nz,ny,nx),mont(nz,ny,nx),lat(ny),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: resx(nz,ny,nx),resy(nz,ny,nx)
    real(kind=nr) :: a_pressx(nz,ny,nx), a_pressy(nz,ny,nx), f(nz,ny,nx)
    integer(kind=ni) :: i,k,nx,ny,nz
    real(kind=nr), parameter :: Omega = 0.00007292, pi=3.141592
    !f2py depend(nx,ny,nz) resx,resy,v,mont
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    forall(i = 1_ni:nx, k = 1_ni:nz)
       f(k,:,i) = 2.*Omega*sin(lat*pi/180.)
    end forall
    call grad(a_pressx,a_pressy,nx,ny,nz,mont,dx,dy)
    resx=-a_pressx+f*v ! pressure force + coriolis force
    resy=-a_pressy-f*u ! pressure force + coriolis force
  end subroutine


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
    !real(kind=nr) :: uax(nz,ny,nx), uay(nz,ny,nx),vax(nz,ny,nx),vay(nz,ny,nx)
    real(kind=nr) :: accelx(nz,ny,nx),accely(nz,ny,nx)
    real(kind=nr) :: dxaccelx(nz,ny,nx),dyaccelx(nz,ny,nx)
    real(kind=nr) :: dxaccely(nz,ny,nx),dyaccely(nz,ny,nx)
    real(kind=nr) :: tem(2,2)
    real(kind=nr) :: eigensr(2),eigensi(2)
    real(kind=nr) :: dummy1(2,2),dummy2(2,2),dummy3(6)
    integer(kind=ni) :: i,j,k, nx,ny,nz, dummy4
    !real(kind=nr), parameter :: Omega = 0.00007292, pi=3.141592
    !f2py depend(nx,ny,nz) respr,resmr,respi,resmi,v,mont
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat
    ! -----------------------------------------------------------------
    !
    call laccel(accelx,accely,nx,ny,nz,u,v,mont,lat,dx,dy)
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
            !print *,i,j,k    
            !print *,gamma(k,j,i)%t(1:2,1:2)      
            !tem(1,1:2)=(/1.,2./)
            !tem(2,1:2)=(/3.,4./)
            tem=gamma(k,j,i)%t
            call dgeev('N','N',2,tem,2,eigensr,eigensi,dummy1,2,dummy2,2,dummy3,6,dummy4)
            !like eigens(1:2) = eig(gamma(k,j,i)%t)
            respr(k,j,i)=eigensr(1) 
            respi(k,j,i)=eigensi(1) 
            resmr(k,j,i)=eigensr(2) 
            resmi(k,j,i)=eigensi(2) 
          enddo
       enddo
    end do

    !print *,gamma(1,10,10)%t(1:2,1:2)  
    !print *,gamma(1,10,10)%t
    !print *,'tem=',tem 
    !print *,'info=',dummy4
    !print *,'lampr(1,10,10)=',respr(1,10,10)
    !print *,'lampi(1,10,10)=',respi(1,10,10),'i'
    !print *,'lammr(1,10,10)=',resmr(1,10,10)
    !print *,'lammi(1,10,10)=',resmi(1,10,10),'i'
    respr(:,(/1,ny/),:)=0.
    !respr(:,:,(/1,nx/))=0.
    resmr(:,(/1,ny/),:)=0.
    !resmr(:,:,(/1,nx/))=0.
    respi(:,(/1,ny/),:)=0.
    !respi(:,:,(/1,nx/))=0.
    resmi(:,(/1,ny/),:)=0.
    !resmi(:,:,(/1,nx/))=0.

  end subroutine

  !Calculates Lagrangian time derivative of compression axis angle
  ! d(phi)/dt (ref Lapeyre et al 1999)
  subroutine dphidt(res,nx,ny,nz,u,v,mont,lat,dx,dy)
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
    call def_shear(sig_sh,nx,ny,nz,u,v,dx,dy)
    call def_stretch(sig_st,nx,ny,nz,u,v,dx,dy)
    call laccel(accelx,accely,nx,ny,nz,u,v,mont,lat,dx,dy)
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

  !Calculates the r disgnostic of Lapeyre et al (1999)
  ! r is the ratio of effective rotation to strain rate
  ! where effective rotation comprises both vorticity and strain-axes rotation
  subroutine r(res,nx,ny,nz,u,v,mont,lat,dx,dy)
    real(kind=nr), intent(in)  :: u(nz,ny,nx),v(nz,ny,nx),mont(nz,ny,nx),lat(ny),dx(ny,nx),dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    real(kind=nr) :: ddtphi(nz,ny,nx), omega(nz,ny,nx), sig(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res,v,mont
    !f2py depend(nx,ny) dx, dy
    !f2py depend(ny) lat

    call dphidt(ddtphi,nx,ny,nz,u,v,mont,lat,dx,dy)
    call vor(omega,nx,ny,nz,u,v,dx,dy)
    call def_total(sig,nx,ny,nz,u,v,dx,dy)
    
    res=(omega+2.*ddtphi)/sig

  end subroutine

end module
