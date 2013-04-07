! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- humidity conversions 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module humidity
  use kind
  use config
  !
  implicit none
contains
  !
  ! Calculate saturation water vapour pressure over a liquid surface from temperature
  subroutine e_sat_liq(res,nx,ny,nz,temp)
    real(kind=nr), intent(in)  :: temp(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res
    ! -----------------------------------------------------------------
    res(:,:,:) = 611.2_nr * exp(17.67_nr * (temp(:,:,:) - 273.15)/((temp(:,:,:) - 273.15) + 243.5_nr) )
    !
    return
  end subroutine
  !
  ! Calculate water vapour pressure from specific humidity and pressure
  subroutine e_from_q(res,nx,ny,nz,q,p)
    real(kind=nr), intent(in)  :: q(nz,ny,nx), p(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) p, res
    ! -----------------------------------------------------------------
    res(:,:,:) = q(:,:,:)*p(:,:,:)/(0.622_nr + (1.0_nr - 0.622_nr)*q(:,:,:))
    !
    return
  end subroutine
  !
  ! Calculate water vapour pressure from water vapour mixing ratio and pressure
  subroutine e_from_mv(res,nx,ny,nz,mv,p)
    real(kind=nr), intent(in)  :: mv(nz,ny,nx), p(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) p, res
    ! -----------------------------------------------------------------
    res(:,:,:) = mv(:,:,:)*p(:,:,:)/(0.622_nr + mv(:,:,:))
    !
    return
  end subroutine
  !
  ! Calculate water vapour pressure from relative humidity (liquid surface) and temperature
  subroutine e_from_rh(res,nx,ny,nz,rh,temp)
    real(kind=nr), intent(in)  :: temp(nz,ny,nx), rh(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, temp
    !
    real(kind=nr) :: e_sat(nz,ny,nx)
    ! -----------------------------------------------------------------
    call e_sat_liq(e_sat, nx,ny,nz, temp)
    res(:,:,:) = rh(:,:,:)/100.0_nr * e_sat(:,:,:)
    !
    return
  end subroutine
  !
  ! Calculate dew point temperature from water vapour pressure
  subroutine tdew_from_e(res,nx,ny,nz,e)
    real(kind=nr), intent(in)  :: e(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res
    ! -----------------------------------------------------------------
    res(:,:,:) = 243.5_nr*log(e(:,:,:)/6.112_nr)  / (17.67_nr - log(e(:,:,:)/6.112_nr)) + 273.15_nr
    !
    return
  end subroutine
  !
  ! Calculate specific humidity from water vapour pressure and pressure
  subroutine q_from_e(res,nx,ny,nz,e,p)
    real(kind=nr), intent(in)  :: e(nz,ny,nx), p(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) p, res
    ! -----------------------------------------------------------------
    res(:,:,:) = 0.622_nr * e(:,:,:)/(p(:,:,:) - (1.0_nr - 0.622_nr)*e(:,:,:))
    !
    return
  end subroutine
  !
  ! Calculate water vapour mixing ratio from water vapour pressure and pressure
  subroutine mv_from_e(res,nx,ny,nz,e,p)
    real(kind=nr), intent(in)  :: e(nz,ny,nx), p(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) p, res
    ! -----------------------------------------------------------------
    res(:,:,:) = 0.622_nr * e(:,:,:)/(p(:,:,:) - e(:,:,:))
    !
    return
  end subroutine
  !
  ! Calculate relative humidity (liquid surface) from water vapour pressure and temperature
  subroutine rh_from_e(res,nx,ny,nz,e,temp)
    real(kind=nr), intent(in)  :: temp(nz,ny,nx), e(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, temp
    !
    real(kind=nr) :: e_sat(nz,ny,nx)
    ! -----------------------------------------------------------------
    call e_sat_liq(e_sat, nx,ny,nz, temp)
    res(:,:,:) = e(:,:,:)/e_sat(:,:,:) * 100.0_nr
    !
    return
  end subroutine
  !
end module
