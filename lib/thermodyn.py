#!/usr/bin/env python

from __future__ import absolute_import, unicode_literals

import sys

from . import dynfor
from . import docutil

# Take over the contents of dynfor.diag to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.thermodyn, 'thermodyn', sys.modules[__name__])


def theta_from_temp(t, pres):
	''' Calculate potential temperature from temperature and pressure

	Uses the Exner function subroutine of dynfor.thermodyn for the calculations.

	Parameters
	----------
	
	t : np.ndarray with shape (nz,ny,nx) and dtype float64
	    Temperature in Kelvin
	pres : np.ndarray with shape (nz,ny,nx) and dtype float64
	    Pressure in Pa
	
	Returns
	-------

	np.ndarray with shape (nz,ny,nx) and dtype float64
	    Potential temperature in Kelvin
	'''

	return t/exner(pres)


def temp_from_theta(pt, pres):
	''' Calculate temperature from potential temperature and pressure

	Uses the Exner function subroutine of dynfor.thermodyn for the calculations.

	Parameters
	----------
	
	pt : np.ndarray with shape (nz,ny,nx) and dtype float64
	    Potential temperature in Kelvin
	pres : np.ndarray with shape (nz,ny,nx) and dtype float64
	    Pressure in Pa
	
	Returns
	-------

	np.ndarray with shape (nz,ny,nx) and dtype float64
	    Temperature in Kelvin
	'''

	return pt*exner(pres)


def density_from_temp(t, pres):
	''' Use ideal gas law to calculate density from temperature and pressure

	Parameters
	----------
	
	t : np.ndarray with shape (nz,ny,nx) and dtype float64
	    Temperature in Kelvin
	pres : np.ndarray with shape (nz,ny,nx) and dtype float64
	    Pressure in Pa
	
	Returns
	-------

	np.ndarray with shape (nz,ny,nx) and dtype float64
	    Density in kg/m3
	'''

	return dynfor.consts.rl*t/pres



#
