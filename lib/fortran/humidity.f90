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
  !@ Calculate saturation water vapour pressure over a liquid surface from temperature
  !@
  !@ TODO: Reference for the formula/constants in the formula?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ temp : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Temperature in K.
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
  !@     Calculated saturation water vapour pressure in Pa.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`tdew_from_e`
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
  !@ Calculate water vapour pressure from specific humidity and pressure
  !@
  !@ TODO: Reference for the formula/constants in the formula?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ q : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Specific humidity in kg/kg.
  !@ p : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Pressure in Pa.
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
  !@     Calculated water vapour pressure in Pa.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`q_from_e`, :meth:`e_from_mv`
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
  !@ Calculate water vapour pressure from water vapour mixing ratio and pressure
  !@
  !@ TODO: Reference for the formula/constants in the formula?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ mv : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Water vapour mixing ratio in kg/kg.
  !@ p : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Pressure in Pa.
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
  !@     Calculated water vapour pressure in Pa.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`mv_from_e`, :meth:`e_from_q`
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
  !@ Calculate water vapour pressure from relative humidity (liquid surface) and temperature
  !@
  !@ TODO: Reference for the formula/constants in the formula?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ rh : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Relative humidity in percent.
  !@ temp : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Temperature in K.
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
  !@     Calculated water vapour pressure in Pa.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`q_from_e`, :meth:`e_from_mv`
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
  !@ Calculate dew point temperature from water vapour pressure
  !@
  !@ TODO: Reference for the formula/constants in the formula?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ e : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Water vapour pressure in Pa.
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
  !@     Calculated dew point temperature in K.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`e_sat_liq`
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
  !@ Calculate specific humidity from water vapour pressure and pressure
  !@
  !@ TODO: Reference for the formula/constants in the formula?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ e : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Water vapour pressure in Pa.
  !@ p : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Pressure in Pa.
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
  !@     Calculated specific humidity in kg/kg.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`e_from_q`, :meth:`mv_from_e`
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
  !@ Calculate water vapour mixing ratio from water vapour pressure and pressure
  !@
  !@ TODO: Reference for the formula/constants in the formula?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ e : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Water vapour pressure in Pa.
  !@ p : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Pressure in Pa.
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
  !@     Calculated water vapour mixing ratio in kg/kg.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`e_from_mv`, :meth:`q_from_e`
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
  !@ Calculate relative humidity (liquid surface) from water vapour pressure and temperature
  !@
  !@ TODO: Reference for the formula/constants in the formula?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ e : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Water vapour pressure in Pa.
  !@ temp : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Temperature in K.
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
  !@     Calculated relative humidity in percent.
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
  !@ Calculate mixing ratio (directly) from specific humidity
  !@
  !@ TODO: Reference for the formula/constants in the formula?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ q : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Specific humidity in kg/kg.
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
  !@     Calculated water vapour mixing ratio in kg/kg.
  subroutine mv_from_q(res,nx,ny,nz,q)
    real(kind=nr), intent(in)  :: q(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res
    ! -----------------------------------------------------------------
    res(:,:,:) = q(:,:,:)/(1._nr-q(:,:,:))
    !
    return
  end subroutine
  !
  !@ Calculate the (isobaric) equivalent temperature from specific humidity, temperature
  !@
  !@ TODO: Reference for the formula/constants in the formula?
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ q : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Specific humidity in kg/kg.
  !@ temp : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Temperature in K.
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
  !@     Calculated equivalent temperature in K.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`thetae_from_q`
  subroutine te_from_q(res,nx,ny,nz,q,temp)
    real(kind=nr), intent(in)  :: temp(nz,ny,nx), q(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, temp
    !
    real(kind=nr) :: mv(nz,ny,nx)
    ! -----------------------------------------------------------------
    call mv_from_q(mv,nx,ny,nz,q)
    res(:,:,:) = temp(:,:,:) + Lv/cp * mv(:,:,:)
    !
    return
  end subroutine
  !
  !@ Calculate the equivalent potential temperature from specific humidity, temperature, pressure
  !@
  !@ Formula from AMS Glossary: http://glossary.ametsoc.org/wiki/Equivalent_potential_temperature
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ q : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Specific humidity in kg/kg.
  !@ temp : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Temperature in K.
  !@ p : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Pressure in Pa.
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
  !@     Calculated equivalent potential temperature in K.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`te_from_q`
  subroutine thetae_from_q(res,nx,ny,nz,q,temp,pres)
    real(kind=nr), intent(in)  :: temp(nz,ny,nx), q(nz,ny,nx), pres(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, temp, pres
    !
    real(kind=nr) :: mv(nz,ny,nx), e(nz,ny,nx), rh(nz,ny,nx)
    ! -----------------------------------------------------------------
    call mv_from_q(mv,nx,ny,nz,q)
    call e_from_q(e,nx,ny,nz,q,pres)
    call rh_from_e(rh,nx,ny,nz,e,temp)
    res(:,:,:) = temp(:,:,:) * (p0/(pres(:,:,:)-e(:,:,:)))**(Rl/cp) &
            &  *   rh(:,:,:) ** (-mv(:,:,:)*(Rv/cp)) * exp(Lv*mv(:,:,:)/(cp*temp(:,:,:))) 
    !
    return
  end subroutine
  !
 !
  !@ Calculate the potential temperature from temperature, pressure
  !@
  !@ Formula from AMS Glossary: http://glossary.ametsoc.org/wiki/Potential_temperature
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ temp : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Temperature in K.
  !@ p : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Pressure in Pa.
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
  !@     Calculated potential temperature in K.
  !@
  !@ See Also
  !@ --------
  !@ :meth:`theta`, must compile everything again for it to find the theta function?
  subroutine theta(res,nx,ny,nz,temp,pres)
    real(kind=nr), intent(in)  :: temp(nz,ny,nx), pres(nz,ny,nx)
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) res, pres
    !
    ! -----------------------------------------------------------------
    res(:,:,:) = temp(:,:,:) * ((p0/pres(:,:,:))**(Rl/cp))
    !
    return
  end subroutine
end module
