! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- humidity conversions 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module humidity
  use kind
  use config
  use consts
  !
  implicit none
contains
  !
  ! Calculate saturation water vapour pressure over a liquid surface from temperature
  subroutine e_sat_liq(res,nx,ny,nz,temp,dx,dy)
    real(kind=nr), intent(in)  :: temp(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    res(:,:,:) = 611.2_nr * exp(17.67_nr * (temp(:,:,:) - 273.15)/((temp(:,:,:) - 273.15) + 243.5_nr) )
    !
    return
  end subroutine
  !
  ! Calculate water vapour pressure from specific humidity and pressure
  subroutine e_from_q(res,nx,ny,nz,q,p,dx,dy)
    real(kind=nr), intent(in)  :: q(nz,ny,nx), p(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) p, res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    res(:,:,:) = q(:,:,:)*p(:,:,:)/(0.622_nr + (1.0_nr - 0.622_nr)*q(:,:,:))
    !
    return
  end subroutine
  !
  ! Calculate water vapour pressure from water vapour mixing ratio and pressure
  subroutine e_from_mv(res,nx,ny,nz,mv,p,dx,dy)
    real(kind=nr), intent(in)  :: mv(nz,ny,nx), p(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) p, res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    res(:,:,:) = mv(:,:,:)*p(:,:,:)/(0.622_nr + mv(:,:,:))
    !
    return
  end subroutine
  !
  ! Calculate water vapour pressure from relative humidity (liquid surface) and temperature
  subroutine e_from_rh(res,nx,ny,nz,rh,temp,dx,dy)
    real(kind=nr), intent(in)  :: temp(nz,ny,nx), rh(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, temp
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: e_sat(nz,ny,nx)
    ! -----------------------------------------------------------------
    call e_sat_liq(e_sat, nx,ny,nz, temp,dx,dy)
    res(:,:,:) = rh(:,:,:)/100.0_nr * e_sat(:,:,:)
    !
    return
  end subroutine
  !
  ! Calculate dew point temperature from water vapour pressure
  subroutine tdew_from_e(res,nx,ny,nz,e,dx,dy)
    real(kind=nr), intent(in)  :: e(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    res(:,:,:) = 243.5_nr*log(e(:,:,:)/6.112_nr)  / (17.67_nr - log(e(:,:,:)/6.112_nr)) + 273.15_nr
    !
    return
  end subroutine
  !
  ! Calculate specific humidity from water vapour pressure and pressure
  subroutine q_from_e(res,nx,ny,nz,e,p,dx,dy)
    real(kind=nr), intent(in)  :: e(nz,ny,nx), p(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) p, res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    res(:,:,:) = 0.622_nr * e(:,:,:)/(p(:,:,:) - (1.0_nr - 0.622_nr)*e(:,:,:))
    !
    return
  end subroutine
  !
  ! Calculate water vapour mixing ratio from water vapour pressure and pressure
  subroutine mv_from_e(res,nx,ny,nz,e,p,dx,dy)
    real(kind=nr), intent(in)  :: e(nz,ny,nx), p(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) p, res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    res(:,:,:) = 0.622_nr * e(:,:,:)/(p(:,:,:) - e(:,:,:))
    !
    return
  end subroutine
  !
  ! Calculate relative humidity (liquid surface) from water vapour pressure and temperature
  subroutine rh_from_e(res,nx,ny,nz,e,temp,dx,dy)
    real(kind=nr), intent(in)  :: temp(nz,ny,nx), e(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, temp
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: e_sat(nz,ny,nx)
    ! -----------------------------------------------------------------
    call e_sat_liq(e_sat, nx,ny,nz, temp,dx,dy)
    res(:,:,:) = e(:,:,:)/e_sat(:,:,:) * 100.0_nr
    !
    return
  end subroutine
  !
  ! Calculate mixing ratio (directly) from specific humidity
  subroutine mv_from_q(res,nx,ny,nz,q,dx,dy)
    real(kind=nr), intent(in)  :: q(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res
    !f2py depend(nx,ny) dx, dy
    ! -----------------------------------------------------------------
    res(:,:,:) = q(:,:,:)/(1._nr-q(:,:,:))
    !
    return
  end subroutine
  !
  ! Calculate the (isobaric) equivalent temperature from specific humidity, temperature
  subroutine te_from_q(res,nx,ny,nz,q,temp,dx,dy)
    real(kind=nr), intent(in)  :: temp(nz,ny,nx), q(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, temp
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: mv(nz,ny,nx)
    ! -----------------------------------------------------------------
    call mv_from_q(mv,nx,ny,nz,q,dx,dy)
    res(:,:,:) = temp(:,:,:) + Lv/cp * mv(:,:,:)
    !
    return
  end subroutine
  !
  ! Calculate the equivalent temperature from specific humidity, temperature, pressure
  ! [formula from AMS Glossary: http://glossary.ametsoc.org/wiki/Equivalent_potential_temperature
  subroutine thetae_from_q(res,nx,ny,nz,q,temp,pres,dx,dy)
    real(kind=nr), intent(in)  :: temp(nz,ny,nx), q(nz,ny,nx), pres(nz,ny,nx), dx(ny,nx), dy(ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, temp, pres
    !f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: mv(nz,ny,nx), e(nz,ny,nx), rh(nz,ny,nx)
    ! -----------------------------------------------------------------
    call mv_from_q(mv,nx,ny,nz,q,dx,dy)
    call e_from_q(e,nx,ny,nz,q,pres,dx,dy)
    call rh_from_e(rh,nx,ny,nz,e,temp,dx,dy)
    res(:,:,:) = temp(:,:,:) * (p0/(pres(:,:,:)-e(:,:,:)))**(Rl/cp) &
            &  *   rh(:,:,:) ** (-mv(:,:,:)*(Rv/cp)) * exp(Lv*mv(:,:,:)/(cp*temp(:,:,:))) 
    !
    return
  end subroutine
  !

!--------------------------------------------------------------------------------
  subroutine sat_vapor_press_2D(tt_,es_)
        !calculates saturation vapor pressure given temperature
        double precision, intent(in)    ::  tt_(ny,nz)
        double precision, intent(out)   ::  es_(ny,nz)
	
	!same equations as in WRF
    do j=1,ny
        do k=1,nz
            if(tt_(j,k).gt.T0)then
                 es_(j,k)=e0*exp(17.67*(tt_(j,k)-T0)/(tt_(j,k)-29.65))
            else
                 es_(j,k)=e0*exp(21.8745584*(tt_(j,k)-T0)/(tt_(j,k)-7.66))
            endif
        enddo
    enddo

    end subroutine sat_vapor_press_2D
!--------------------------------------------------------------------------------




end module
