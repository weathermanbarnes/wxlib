#!/usr/bin/env python

from __future__ import absolute_import, unicode_literals, division, print_function

import sys
import pickle
import numpy as np
from datetime import datetime as dt, timedelta as td

from . import dynfor
from . import docutil
from . import settings as s
from .metio import metopen, metsave

# Take over the contents of dynfor.diag to this module and inject documentation from the Fortran sources
docutil.takeover(dynfor.detect, 'detect', sys.modules[__name__])


def block_by_grad_rev(q='pv', plev='pt330', lat_band=(30, 70), 
		local_move=(27,36), total_move=(40,54), min_duration=21):
	''' Detect blocking following the Masato et al. (2012) procedure based on persistent gradient reversals

	The procedure involves three main steps:
	 
	 1. Find instantaneous regions with reversed gradients (done by 
	    :meth:`block_indicator_grad_rev`)
	 2. Connext these regions in time and apply stationarity criteria
	 3. Mask those regions with reversed gradients that fullfil the criteria in step 2
	
	The function does not return any data, but saves intermediate and final results in
	``conf.opath``.

	Parameters
	----------

	q : str
	    Which variable to base the detection on. Default: ``'pv'``.
	plev: str or int
	    Which vertical level to base the detection on. Default ``'pv2000'``.
	lat_band: 2-tuple of float or int
	    Which latitudes in degrees to consider for the blocking detection. Default:
	    ``(30, 70)``.
	local_move: 2-tuple of float or int
	    Maximum movement in degrees latitude/longitude of the block center from 
	    one time step to the next. Default: ``(27, 36)``.
	total_move: 2-tuple of float or int
	    Maximum movement in degrees latitude/longitude of the block center during 
	    the entire duration. Default: ``(40, 54)``.
	min_duration: int
	    Minimum number of time steps the block as to persist. Default: ``21``.
	'''

	import scipy.ndimage.filters as ndf

	Q_BI = 'blockint'
	Q_BLOCK = 'block'
	
	print('Stage 1+2: Finding gradient reversals and conenct them in time')
	blocks = []
	prevblocks = {}
	ny, nx = s.conf.gridsize
	dtd = s.conf.timestep
	tlen = 0
	for year in s.conf.years:
		# Stage 1: Load base data, calculate and save blocking indicator
		f, dat, grid = metopen(s.conf.file_std % {'time': year, 'plev': plev, 'qf': s.conf.qf[q]}, q)
		dat = dat.squeeze()
		if year == s.conf.years[0]:
			dynfor.config.grid_cyclic_ew = grid.cyclic_ew
			print(grid.cyclic_ew)
		bi = dynfor.detect.block_indicator_grad_rev(dat, grid.dx, grid.dy)
		metsave(bi[:,np.newaxis,:,:], grid, q=Q_BI, plev=plev)
		tlen += bi.shape[0]
		
		# Translate thresholds from degrees to grid point indexes
		if year == s.conf.years[0]:
			j0 = np.argwhere(grid.y[:,0] == lat_band[0])[0,0]
			j1 = np.argwhere(grid.y[:,0] == lat_band[1])[0,0]
			if j0 > j1: 
				j1, j0 = j0, j1

			dx = abs(sorted(grid.x[0,1:]-grid.x[0,:-1])[nx//2])
			dy = abs(sorted(grid.y[1:,0]-grid.y[:-1,0])[ny//2])

			djl_max, dil_max = local_move[0]/dy, local_move[1]/dx
			djt_max, dit_max = total_move[0]/dy, total_move[1]/dx

		# Stage 2: Connect in time
		bi = bi[:,j0:j1,:]
	
		# Identify local extrema with the appropriate sign, special treatment for pv as it increases with latitude
		if q == 'pv':
			extrema = (ndf.minimum_filter(bi, size=(1,3,3), mode='wrap') == bi)
			extrema[bi >= 0] = False
		else:
			extrema = (ndf.maximum_filter(bi, size=(1,3,3), mode='wrap') == bi)
			extrema[bi <= 0] = False

		# Remove local extrema directly at the border of the considered domain
		extrema[:,0,:] = False
		extrema[:,-1,:] = False

		# Joining into block objects
		for tidx in range(bi.shape[0]):
			done = {}
			curblocks = {}
			# From mask to grid point indexes
			expos = np.argwhere(extrema[tidx,:,:] > 0)
			for prevpos, blocknr in prevblocks.items():
				initpos = blocks[blocknr]['pos'][0]
				dij2 = 9999999
				for pos in expos:
					pos = tuple(pos)
					# Movement since last
					djl = abs(pos[0]-prevpos[0])
					dil = abs(pos[1]-prevpos[1])
					if dil > nx/2:
						dil = nx - dil
					# Total movement
					djt = abs(pos[0]-initpos[0])
					dit = abs(pos[1]-initpos[1])
					if dit > nx/2:
						dit = nx - dit
					# Check if connected to known block
					if djl <= djl_max and djt <= djt_max and \
					   dil <= dil_max and dit <= dit_max:
						# New closest continuation of the block
						if djl**2 + dil**2 < dij2:
							dij2 = djl**2 + dil**2
							savepos = pos
						# Belongs to same block structure, but another maximum is closer
						else:
							pass

						done[pos] = True

				if dij2 < 9999999:
					curblocks[savepos] = blocknr
					blocks[blocknr]['pos'].append(savepos)
					blocks[blocknr]['blockidx'].append(bi[tidx,savepos[0],savepos[1]])

			# Save new blocks
			for pos in expos:
				pos = tuple(pos)
				if not pos in done:
					blocks.append({
						'pos': [pos, ],
						'blockidx': [bi[tidx,pos[0],pos[1]], ],
						'onset': dt(year,1,1) + dtd*tidx,
					})
					curblocks[pos] = len(blocks) -1

			prevblocks = curblocks

	# Save information on time connection
	filename = s.conf.opath+'/'+(s.conf.file_timeless % 
			{'time': str(s.conf.years[0])+'-'+str(s.conf.years[-1]), 
			 'name': plev+'.'+Q_BLOCK})+'.pickle'
	fsave = file(filename, 'w')
	pickle.dump(blocks, fsave)
	fsave.close()

	# Create an list of seeds.
	print('Stage 3: Applying minimum duration criterion and write out block masks')
	seeds = []
	dates = [dt(s.conf.years[0],1,1) + dtd*tidx for tidx in range(tlen)]
	for tidx in range(tlen):
		seeds.append([])
	for block in blocks:
		if len(block['blockidx']) >= min_duration:
			tidx0 = dates.index(block['onset'])
			for dtidx in range(len(block['blockidx'])):
				pos = block['pos'][dtidx]
				# convert to Fortran indexes
				seeds[tidx0+dtidx].append([pos[0]+1, pos[1]+1])
	
	offset = 0
	# Make sure to use the data saved previously, not anything else!
	datapath_ = s.conf.datapath
	s.conf.datapath = [s.conf.opath, ]
	for year in s.conf.years:
		# Open nc file
		fbi, bi, grid  = metopen(s.conf.file_std % {'time': year, 'plev': plev, 'qf': s.conf.qf[Q_BI]}, Q_BI)
		bi = bi.squeeze()
		if q == 'pv':		# PV does climatologically increase polewards
			bi = -bi
		blockmask = np.zeros(bi.shape, dtype='bool')
		bi = bi[:,j0:j1,:]

		for tidx in range(bi.shape[0]):
			if len(seeds[tidx+offset]) > 0:
				blockmask[tidx,j0:j1,:] = dynfor.utils.mask_minimum_connect(
						bi[tidx,:,:], 
						np.array(seeds[tidx+offset], dtype='i4'), 
						0.0
				)

		offset += bi.shape[0]
		
		blockmask = blockmask.astype('i1')
		metsave(blockmask[:,np.newaxis,:,:], grid, q=Q_BLOCK, plev=plev, compress_to_short=False)
	
	# Restore previous data path
	s.conf.datapath = datapath_


def frontalvolume_largescale(tfp, dx, dy):
	''' Detect frontal zones as coherent areas with strong TFP gradients
	
	The large-scale version of this function (applicable for example to ERA-Interim data)
	uses only the large-scale gradient threshold. 

	The gradient threshold can be set by ``dynfor.config.tfp_mingrad_largescale``,
	the minimum size of the frontal zones by ``dynfor.config.tfp_grad_minsize``.

	Parameters
	----------

	tfp : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
	    Thermal front parameter (TFP) field, typically potential or equivalent potential 
	    temperature.
	dx : np.ndarray with shape (ny,nx) and dtype float64
	    The double grid spacing in x-direction to be directly for centered differences.
	    ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``. 
	dy : np.ndarray with shape (ny,nx) and dtype float64
	    The double grid spacing in y-direction to be directly for centered differences.
	    ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
	
	Returns
	-------
	np.ndarray with shape (nz,ny,nx) and dtype int32
	    Detected frontal zones labeled with object IDs > 0, zero everywhere else.

	See also
	--------
	:meth:`frontalvolume_smallscale`, :meth:`frontline`
	'''

	thres = dynfor.config.tfp_mingrad_largescale		# 4.5e-5
	min_size = dynfor.config.minsize			# 75000.0e+6
	cellsize = np.abs(dx)*np.abs(dy)/4.0

	labels = np.empty(tfp.shape, dtype='i4')
	for tidx in range(tfp.shape[0]):
		sys.stdout.write('%d/%d\r' % (tidx+1, tfp.shape[0]))
		sys.stdout.flush()
		ddx, ddy = dynfor.derivatives.grad(tfp[tidx,:,:,:], dx, dy)
		tfp_grad = np.sqrt(ddx**2 + ddy**2)
		mask = tfp_grad > thres

		labels_, sizes = dynfor.utils.label_connected_3d(mask, cellsize, 2500)
		for n, size in zip(range(1,len(sizes)+1), sizes):
			if size == 0:
				break

			if size < (mask.shape[0]*min_size):
				labels_[labels_ == n] = 0

		labels[tidx,:,:,:] = labels_
	
	print('')

	return labels


def frontalvolume_smallscale(tfp, dx, dy):
	''' Detect frontal zones as coherent areas with strong TFP gradients
	
	The small-scale version of this function (applicable for example to NORA10 data)
	uses a combination of a large-scale gradient threshold for a smoothed TFP field
	and a small-scale gradient threshold for the unfiltered TFP field.

	The gradient thresholds can be set by ``dynfor.config.tfp_mingrad_largescale``,
	and ``dynfor.config.tfp_mingrad_smallscale``, the minimum size of the frontal 
	zones by ``dynfor.config.tfp_grad_minsize`` and the amount of smoothing by
	``dynfor.config.nsmooth``.

	Parameters
	----------

	tfp : np.ndarray with shape (nt,nz,ny,nx) and dtype float64
	    Thermal front parameter (TFP) field, typically potential or equivalent potential 
	    temperature.
	dx : np.ndarray with shape (ny,nx) and dtype float64
	    The double grid spacing in x-direction to be directly for centered differences.
	    ``dx(j,i)`` is expected to contain the x-distance between ``(j,i+1)`` and ``(j,i-1)``. 
	dy : np.ndarray with shape (ny,nx) and dtype float64
	    The double grid spacing in y-direction to be directly for centered differences.
	    ``dy(j,i)`` is expected to contain the y-distance between ``(j+1,i)`` and ``(j-1,i)``.
	
	Returns
	-------
	np.ndarray with shape (nz,ny,nx) and dtype int32
	    Detected frontal zones labeled with object IDs > 0, zero everywhere else.

	See also
	--------
	:meth:`frontalvolume_largescale`, :meth:`frontline`
	'''

	nsmooth = dynfor.config.nsmooth
	thres_ls = dynfor.config.tfp_mingrad_largescale		# 4.5e-5
	thres_ss = dynfor.config.tfp_mingrad_smallscale		# 7.5e-5
	min_size = dynfor.config.minsize			# 75000.0e+6
	cellsize = np.abs(dx)*np.abs(dy)/4.0

	labels = np.empty(tfp.shape, dtype='i4')
	for tidx in range(tfp.shape[0]):
		ddx, ddy = dynfor.derivatives.grad(tfp[tidx,:,:,:], dx, dy)
		tfp_grad = np.sqrt(ddx**2 + ddy**2)
	
		tfps = dynfor.utils.smooth_xy(tfp[tidx,:,:,:], nsmooth)
		ddx, ddy = dynfor.derivatives.grad(tfps, dx, dy)
		tfps_grad = np.sqrt(ddx**2 + ddy**2)

		mask = (tfp_grad > thres_ss) & (tfps_grad > thres_ls)

		labels_, sizes = dynfor.utils.label_connected_3d(mask, cellsize, 1000)
		for n, size in zip(range(1,len(sizes)+1), sizes):
			if size == 0:
				break

			if size < (mask.shape[0]*min_size):
				labels_[labels_ == n] = 0

		labels[tidx,:,:,:] = labels_

	return labels



#
