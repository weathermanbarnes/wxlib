#!/usr/bin/env python
# -*- encoding: utf-8


import os
import netCDF4 as nc
import scipy.io.matlab as mat

import numpy as np
import scipy.interpolate as interp
import cftime
import pytz
import copy
from datetime import datetime as dtv, timedelta as td

from ..gridlib import grid_by_static, grid_by_nc
from .. import utils

from .. import version as dynlib_version


# 16 billion values -> 128G of memory
# (Corresponds to a bit more than 40 years of 6 hourly data for a single variable/level at 0.5deg resolution)
MAX_REQUEST_SIZE = 16.0e9 
# 2 billion values -> 16G of memory
WARN_REQUEST_SIZE = 2.0e9


class variable(object):
    def __init__(self, q, long_name, units, qf=None, q_std=None, bins=None, lines=None, objmask=None):
        self.q = q
        self.long_name = long_name
        self.units = units
        
        if not qf:
            self.qf = q
        else:
            self.qf = qf

        if not q_std:
            self.q_std = q
        else:
            self.q_std = q_std
        
        self.bins = bins
        self.lines = lines
        self.objmask = objmask

        if lines:
            self.q_avg = f'{q}_freq'
        elif objmask:
            self.q_avg = f'{objmask}_freq'
        else:
            self.q_avg = q

        return


def plev_to_tuple(plev, with_unit=False):
    ''' Return meta-information about a given plev string

    Returns the type of vertical level and a numeric value

    Parameters
    ----------
    plev : str
        The level specifying string.
    with_unit: bool
        If ``'True'``, return string giving the unit as a third return variable.

    Returns
    -------
    plevtype : str
        The type of vertical level, e.g. ``'isobaric'``, ``'isentropic'``, 
        ``'potential vorticity'``, ``'height above ground'``.
    plevhgt : float
        The height of the vertical level in the appropriate unit for the ``'plevtype'``.
    plevunit : str
        If ``'with_unit'``, return unit for the height of the vertical level.
    '''

    if plev[:2] == 'pt':
        plevtype = 'isentropic'
        plevunit = 'K'
        plevhgt = float(plev[2:])
    elif plev[:2] == 'pv':
        plevtype = 'potential vorticity'
        plevunit = 'PVU'
        plevhgt = float(plev[2:])/1000.0
    elif plev[-1:] == 'm':
        plevtype = 'height above ground'
        plevunit = 'm'
        plevhgt = float(plev[:-1])
    else:
        try: 
            plevhgt = float(plev) * 100.0
        except ValueError:
            raise Exception(f'Cannot interpret vertical level string `{plev}`')
        plevtype = 'isobaric'
        plevunit = 'Pa'

    if with_unit:
        return plevtype, plevhgt, plevunit
    else:
        return plevtype, plevhgt


def tuple_to_plev(plevtype, plevhgt):
    ''' Return brief plev descriptor string based on type and height

    Parameters
    ----------
    plevtype : str
        The type of vertical level, e.g. ``'isobaric'``, ``'isentropic'``, 
        ``'potential vorticity'``, ``'height above ground'``.
    plevhgt : float
        The height of the vertical level in the appropriate unit for the ``'plevtype'``.

    Returns
    -------
    plev : str
        The level specifying string.
    '''

    if plevtype == 'isobaric':
        plev = str(int(plevhgt/100))
    elif plevtype == 'isentropic':
        plev = f'pt{int(plevhgt)}'
    elif plevtype == 'potential vorticity':
        plev = f'pv{int(plevhgt*1000)}'
    elif plevtype == 'height above ground':
        plev = f'{int(plevhgt)}m'
    else:
        raise NotImplementedError(f'Unknown vertical level type `{plevtype}`')

    return plev


class files_by_plevq(object):
    def __init__(self, plevq, start, end):
        plev, q = plevq
        self.plev = plev
        self.q = q
        self.start = start
        self.cur = start
        self.end = end

        return 

    def __iter__(self, ):
        return self
    
    def __next__(self):
        raise StopIteration


def get_from_file(filename, plev, q, **kwargs):
    ''' Get data for level plev and variable q from filename
    
    The file is located by metopen. Optional arguments are passed on to metopen.

    The located files might either be numpy-files, netCDF-files or matlab mat-files. 
    For each of the file types, metopen returns the requested variable and some meta-
    information about the variable, if not suppressed.

    Parameters
    ----------
    filename : str
        The name of the file, excluding the file ending.
    plev : str
        The requested vertical level for the variable. Might be ``'__all__'`` for all
        vertical levels in the file if the data source supports that.
    q : str
        The requested variable within the file. Might be ``'__all__'`` for all
        variables in the file if the data source supports that.

    Keyword arguments
    -----------------
    metopen arguments : all optional
        Optional arguments passed on to calls of metopen within this function.

    Returns
    -------
    (dict of) np.ndarray
        Data of the requested variable for the requested vertical level. If several 
        variables are requested, a dictionary is returned with the variable shorthand q 
        as keys.
    grid.gridlib
        If ``no_static=False`` meta-information about the requested data.
    '''

    raise NotImplementedError('Needs to be implemented per data source')


def get_normalized_from_file_factory(get_from_file, conf):

    def get_normalized_from_file(filename, plev, q, **kwargs):
        ''' A variant of get_from_file which in addition normalises feature detections
        
        Feature detections are normalized such that they appear on a regular grid and can
        directly be averaged. Thus, line features are normalized by ``dynlib.utils.normalize_lines``,
        areal detections are converted to a boolean mask field.

        Parameters
        ----------
        filename : str
            The name of the file, excluding the file ending.
        plev : str
            The requested vertical level for the variable. Might be ``'__all__'`` for all
            vertical levels in the file if the data source supports that.
        q : str
            The requested variable within the file. Might be ``'__all__'`` for all
            variables in the file if the data source supports that.

        Keyword arguments
        -----------------
        metopen arguments : all optional
            Optional arguments passed on to calls of metopen within this function.

        Returns
        -------
        (dict of) np.ndarray
            Data of the requested variable for the requested vertical level. If several 
            variables are requested, a dictionary is returned with the variable shorthand q 
            as keys.
        grid.gridlib
            If ``no_static=False`` meta-information about the requested data.
        '''

        ret = get_from_file(filename, plev, q, **kwargs)

        if type(ret) == tuple:
            dat = ret[0]
            stuff = ret[1:]
        else:
            dat = ret

        # Treat lines
        if q in conf.q_lines:
            datoff, grid = get_from_file(filename, plev, conf.q_lines[q])
            dat = utils.normalize_lines(dat, datoff, grid.dx, grid.dy)

        # Treat objmasks
        if q in conf.q_obj:
            dat = dat > 0.5

        # Binned variables kept as they are
        
        if type(ret) == tuple:
            return (dat, *stuff)
        else:
            return dat

    return get_normalized_from_file
    

# TODO: Return xarray by default, reflect other docu changes (if any?)
def metopen_factory(get_static, conf):
    ''' Create the metopen function based on data source specific helpers '''
    
    # Interpret files without file type designator as netCDF
    NO_ENDING = 'nc'

    def handle_npy(filepath, q=None, cut=slice(None), verbose=False, no_dtype_conversion=False,
            no_static=False, quiet=False, mode='r', no_xarray=False):
        ''' Helper function to read npy files '''

        if not mode == 'r':
            print('WARNING: Can only open npy files in read mode!')
        if q:
            dat = np.load(filepath, mmap_mode='r', encoding='bytes')
            dat = dat[cut]
        else:
            dat = None
        if not quiet:
            print('Found '+filepath)
        f = None

        return f, dat

    def handle_npz(filepath, q=None, cut=slice(None), verbose=False, no_dtype_conversion=False,
            no_static=False, quiet=False, mode='r', no_xarray=False):
        ''' Helper function to read npz files '''

        if not mode == 'r' and not quiet:
            print('WARNING: Can only open npz files in read mode!')
        f = np.load(filepath, encoding='bytes', allow_pickle=True)
        if q:
            if q not in f.files:
                tried.append(filepath)
            dat = f[q][cut]
        else:
            dat = None
        if not quiet:
            print('Found '+filepath)

        return f, dat
    
    def handle_mat(filepath, q=None, cut=slice(None), verbose=False, no_dtype_conversion=False,
            no_static=False, quiet=False, mode='r', no_xarray=False):
        ''' Helper function to read mat files '''

        if not mode == 'r' and not quiet:
            print('WARNING: Can only open mat files in read mode!')
        f = mat.loadmat(filepath)
        if q:
            if q not in f:
                tried.append(filepath)
            dat = f[q][cut]
        else:
            dat = None
        if not quiet:
            print('Found '+filepath)

        return f, dat

    def handle_nc(filepath, q=None, cut=slice(None), verbose=False, no_dtype_conversion=False,
            no_static=False, quiet=False, mode='r', no_xarray=False):
        ''' Helper function to read netCDF files '''

        f = nc.Dataset(filepath, mode)
        f.set_auto_scale(False)
        f.set_auto_mask(False)
        if q:
            var = f.variables[q]
            if q not in f.variables:
                tried.append(filepath)
            if not no_dtype_conversion:
                dat = utils.scale(var, cut=cut)
            else:
                dat = var[cut]
        else:
            dat = None

        if not no_static:
            if q:
                static = grid_by_nc(f, var)
            else:
                static = grid_by_nc(f)
            # TODO: Where to search for topography in nc files?
            static.oro = np.zeros((static.ny, static.nx))
        else:
            static = None

        if not quiet:
            print('Found '+filepath)

        return f, dat, static


    def metopen(filename, q=None, **kwargs):
        ''' Find and open files by name
        
        Uses the conf.datapath list to locate files in a variety of different locations. 

        The located files might either be numpy-files, netCDF-files or matlab mat-files. 
        For each of the file types, metopen returns the requested variable and some meta-
        information about the variable, if not suppressed.

        Parameters
        ----------
        filename : str
            The name of the file, excluding the file ending.
        q : str
            *Optional*. The requested variable within the file.
        cut : slice
            *Optional*, default ``slice(None)``. Limit the request to a given time slice. 
            With the default data layout, only relevant data needs to be read when only 
            a time slice of the entire data is requested. Hence, using cut to limit your 
            data request can make reading the data largely more efficient.
        verbose : bool
            *Optional*, default ``False``. Print debug information on which files are 
            being looked for.
        no_dtype_conversion : bool
            *Optional*, default ``False``. By default, ``metopen`` uncompresses data in the 
            file and converts all data to float64. This behaviour can be suppressed by 
            setting this option to ``True``. This implies, however, that scaling and offset
            cannot be applied automatically.
        no_static : bool
            *Optional*, default ``False``. By default, ``metopen`` does its best to 
            provide meta-information about the requested file, using 
            :module:`grid.gridlib`, and returns the meta-information as a third value. 
            This behaviour can be suppressed by setting this parameter to ``True``.
        mode : str
            *Optional*, default ``'r'``. Only effective for netCDF files. The read/write mode
            with which to open the file. Valid values are ``'r'`` for read-only access, `'a'``
            or ``'r+'` for read-write access and and ``'w'``for replacing the given file. 
        no_xarray : bool
            *Optional*, default ``False``. If mode='r', netCDF files are by default returned 
            as xarray data sets. If ``True`` metopen instead returns an netCDF4 file object.

        Returns
        -------
        data file object
            xarray dataset or python data, netCDF or Matlab file object.
        np.ndarray
            If q given, data of the requested variable.
        grid.gridlib
            If ``no_static=False`` meta-information about the requested data.
        '''

        kwargs['no_static'] = kwargs.get('no_static', False)
        kwargs['verbose'] = kwargs.get('verbose', False)
        kwargs['quiet'] = kwargs.get('quiet', False)
        kwargs['no_dtype_conversion'] = kwargs.get('no_dtype_conversion', False)
        kwargs['cut'] = kwargs.get('cut', slice(None))
        
        if not type(kwargs['cut']) == slice:
            raise ValueError('cut must be a 1-dimensional slice object')

        tried = []
        for path in conf.datapath:
            static = None

            if kwargs['verbose']:
                print('Trying: '+path+'/'+filename+'.*')

            if os.path.exists(path+'/'+filename+'.npy'):
                f, dat = handle_npy(path+'/'+filename+'.npy', q, **kwargs)
            elif os.path.exists(path+'/'+filename) and NO_ENDING == 'npy':
                f, dat = handle_npy(path+'/'+filename, q, **kwargs)

            elif os.path.exists(path+'/'+filename+'.npz'):
                f, dat = handle_npz(path+'/'+filename+'.npz', q, **kwargs)
            elif os.path.exists(path+'/'+filename) and NO_ENDING == 'npz':
                f, dat = handle_npz(path+'/'+filename, q, **kwargs)

            elif os.path.exists(path+'/'+filename+'.mat'):
                f, dat = handle_mat(path+'/'+filename+'.mat', q, **kwargs)
            elif os.path.exists(path+'/'+filename) and NO_ENDING == 'mat':
                f, dat = handle_mat(path+'/'+filename, q, **kwargs)

            elif os.path.exists(path+'/'+filename+'.nc'):
                f, dat, static = handle_nc(path+'/'+filename+'.nc', q, **kwargs)
            elif os.path.exists(path+'/'+filename) and NO_ENDING == 'nc':
                f, dat, static = handle_nc(path+'/'+filename, q, **kwargs)

            else:
                tried.append(path)
                continue
            
            if q and not kwargs['no_dtype_conversion'] and not dat.dtype == 'f8':
                dat = dat.astype('f8')

            if not kwargs['no_static']:
                if not static:
                    static = get_static(kwargs['verbose'], kwargs['no_dtype_conversion'], kwargs['quiet'])
            else:
                if q:
                    return f, dat
                else:
                    return f
            if q:
                return f, dat, static
            else:
                return f, static
        
        raise ValueError('%s.* not found in any data location. \n'
                'Tried the following (in order):\n\t%s' % (filename, '\n\t'.join(tried)) )

    return metopen



def metsave_factory(metopen, conf):
    ''' Create the metsave function based on data source specific helpers '''

    known_vertical_level_units = {
            'Pa': ('pressure', 'down'),
            'K': ('isentropic', 'up'),
            'PVU': ('potential_vorticity', 'up'),
    }

    def _metsave_vargroup(f, all_dims, dat, static, compress_to_short, prefix=''):
        ''' Helper function for metsave: Saving a group of variables given a file with defined dimensions '''

        for key, values in dat.items():
            # Recursively save subgroups of data
            if type(values) == dict:
                _metsave_vargroup(f, all_dims, values, static, compress_to_short, prefix=f'{prefix}/{key}')

                continue
            
            # Expand values if it contains more than the actual numerical data
            elif type(values) == tuple:
                dims, values, varattrs = values
                
            # Otherwise assume that variable covers all defined dimensions and 
            # fill in meta data only from dynlib defaults
            else:
                dims = all_dims
                varattrs = {}
            
            # Variable name with group path, delete old variables if present
            if prefix:
                q = f'{prefix}/{key}'
                if prefix in f.groups and key in f.groups[prefix].variables:
                    i = 1
                    while f'{key}_{i}' in f.groups[prefix].variables:
                        i += 1

                    qold = q
                    q = f'{prefix}/{key}_{i}'
                    print(f'Warning: {qold} already exists, saving as {q} instead.')
            else: 
                q = key
                if key in f.variables:
                    i = 1
                    while f'{key}_{i}' in f.variables:
                        i += 1

                    qold = q
                    q = f'{key}_{i}'
                    print(f'Warning: {qold} already exists, saving as {q} instead.')

            # Optional: conversion from f8 to i2 to save space
            if compress_to_short:
                values, scale, off, fill = utils.unscale(values)
                if fill: 
                    ovar = f.createVariable(q, 'i2', dims, fill_value=fill)
                    ovar.set_auto_scale(False)

                    ovar.missing_value = fill
                else:
                    ovar = f.createVariable(q, 'i2', dims)

                varattrs['add_offset'] = off
                varattrs['scale_factor'] = scale

            else:
                if 'missing_value' in varattrs:
                    fill = varattrs.pop('missing_value')
                    ovar = f.createVariable(q, values.dtype, dims, fill_value=fill)
                    ovar.set_auto_scale(False)

                    ovar.missing_value = fill
                else:
                    ovar = f.createVariable(q, values.dtype, dims)
            
            # Fill in some dynlib default meta information if required
            if 'long_name' not in varattrs and key in conf.q_long:
                varattrs['long_name'] = conf.q_long[key]
            if 'units' not in varattrs and key in conf.q_units:
                varattrs['units'] = conf.q_units[key]

            # Add some attributes to make ncview display the rotated grid correctly
            if getattr(static, 'rotated', False) and static.x_name in dims and static.y_name in dims:
                varattrs['grid_mapping'] = 'rotated_pole'
                varattrs['coordinates'] = f'{static.x_name} {static.y_name}'

            ovar.setncatts(varattrs)
            if not ovar.shape == values.shape:
                raise ValueError(f'Input variable for {q} does not have the required shape. '
                                 f'Expected {ovar.shape}, got {values.shape}.')
            ovar[::] = values

        return 

    def metsave(dat, static, filename, compress_to_short=True, add_dims=[], global_attrs={}):
        ''' Save data in a netCDF file

        The data is saved either to an existing file with a matching name in conf.datapath, or, if such 
        a file does not exist, to a new file in conf.opath. Existing variables will not be overwritten,
        but the new variable will be renamed to resolve the conflict.

        Errors will be raised if, and only if, 
         (a) the shape of the given data does not match the given dimensions or the predefined 
             dimensions, respectively, OR
         (b) the dimensions in the static file and in add_dims do not match the dimensions in a
             prexisting file.

        Parameters
        ----------
        dat : dict, mapping q => np.ndarray or q => (dims, np.ndarray, attributes) or groupname => dict 
            Keys either map the variable name to data or netCDF group name to the contents of a group. 
            If no meta information is passed along with the variable data, the data is expected to cover 
            all given dimensions and attributes are filled as far as possible from dynlib's variable 
            definitions.
        static : gridlib.grid
            Some meta information about the data, like the grid information.
        filename : str
            File name to be saved to.
        compress_to_short : bool
            *Optional*, default ``True``. By default, ``metsave`` compresses the data by converting
            the data field into int16, using the float64 ``add_offset`` and ``scale_factor`` attributes 
            to represent the data.
        add_dims : list of dict
            *Optional*, default empty. List of additional dimensions, in order, to be placed between time 
            and space. The dict must have the keys ``'name'`` and ``'data'``, defining the axis name and 
            its data values. All additional keys will be taken over as attributes for the axis variable.
        global_attrs: 
            *Optional*, default empty. List of additional global attributes to be included in the netCDF file.
        '''
        
        for dim in add_dims:
            if 'name' not in dim or 'data' not in dim:
                raise ValueError('Additional dimensions must have a name and data')

        now = dtv.now(pytz.timezone(conf.local_timezone))
        new_history = '%s by %s' % (now.strftime('%Y-%m-%d %H:%M:%S %Z'), dynlib_version)

        try:
            # Try to locate a matching file
            f = metopen(filename, no_static=True, mode='a')
            
            # If dimensions do not exist, they can and will be created; but if they exist, they should match!
            dims = []
            if static.t_name:
                dims.append((static.t_name, static.t) )
            for dim in add_dims:
                dims.append((dim['name'], dim['data']) )
            if static.z_name:
                dims.append((static.z_name, static.z) )
            if static.y_name:
                dims.append((static.y_name, static.y[:,0]) )
            if static.x_name:
                dims.append((static.x_name, static.x[0,:]) )

            for dim, dimvals in dims:
                if dim and dim in f.dimensions:
                    if not np.all(dimvals == f[dim][::]):
                        print(f'Skipping existing because of mismatching {dim} dimension.')
                        raise ValueError
            
            print('Saving to existing %s' % filename)

        except ValueError: 
            try:
                f.close()
            except:
                pass
            
            # TODO: Check whether exists and if so adapt filename to avoid overwriting data
            f = nc.Dataset(conf.opath+'/'+filename+'.nc', 'w', format='NETCDF4')
            print('Saving to %s/%s.nc' % (conf.opath, filename))


        # Checks whether file dimensions are consistent
        f.setncatts(global_attrs)
        f.setncatts({'Conventions': 'CF-1.0'})

        if hasattr(f, 'history'):
            f.setncatts({'history': f'{f.history}\n{new_history}'})
        else:
            f.setncatts({'history': new_history})

        all_dims = ()
        
        # Time dimension
        if static.t_name:
            if static.t_name not in f.dimensions:
                f.createDimension(static.t_name, len(static.t))
                ot = f.createVariable(static.t_name, 'i', (static.t_name,))
                tattrs = {'long_name': 'time', 'units': static.t_unit}
                if hasattr(static.t_parsed[0], 'calendar'):
                    tattrs['calendar'] = static.t_parsed[0].calendar
                ot.setncatts(tattrs)
                ot[::] = static.t
            all_dims += (static.t_name, )
        
        # Additional dimensions
        for dim in add_dims:
            all_dims += (dim['name'], )

            if dim['name'] not in f.dimensions:
                f.createDimension(dim['name'], len(dim['data']))
                odim = f.createVariable(dim['name'], dim['data'].dtype, (dim['name'], ))
                odim[:] = dim['data']

                del dim['name']
                del dim['data']
                odim.setncatts(dim)
        
        # Vertical dimension
        if static.z_name:
            if static.z_name not in f.dimensions:
                f.createDimension(static.z_name, len(static.z))
                z_longname, z_positive = known_vertical_level_units[static.z_unit]
                olev = f.createVariable(static.z_name, 'f', (static.z_name,))
                olev.setncatts({'long_name': z_longname, 'units': static.z_unit, 'axis': 'Z', 'positive': z_positive})
                olev[::] = static.z[:]
            all_dims += (static.z_name, )
        
        # Horizontal dimensions
        if getattr(static, 'rotated', False):
            raise NotImplementedError('This seems broken, double-check!')

            olat = f.createVariable(static.rot_y_name, 'f', (y_name,))
            olat.setncatts({'long_name': static.rot_y_longname, 'units': static.y_unit, 'axis': 'Y'})
            olat[::] = static.rot_y[:,0]
            olon = f.createVariable(static.rot_x_name, 'f', (x_name,))
            olon.setncatts({'long_name': static.rot_x_longname, 'units': static.x_unit, 'axis': 'X'})
            olon[::] = static.rot_x[0,:]

            olat = f.createVariable(static.y_name, 'f', (y_name, x_name,))
            olat.setncatts({'long_name': static.y_name, 'units': static.y_unit})
            olat[::] = static.y
            olon = f.createVariable(static.x_name, 'f', (y_name, x_name,))
            olon.setncatts({'long_name': static.x_name, 'units': static.x_unit})
            olon[::] = static.x

            orot = f.createVariable('rotated_pole', 'i', ())
            orot.setncatts({
                'grid_north_pole_longitude': static.rot_np[1], 
                'grid_north_pole_latitude': static.rot_np[0],
            })
            all_dims += (static.y_name, static.x_name)

        else:
            if static.y_name:
                if static.y_name not in f.dimensions:
                    f.createDimension(static.y_name, static.y.shape[0])
                    olat = f.createVariable(static.y_name, 'f', (static.y_name,))
                    olat.setncatts({'long_name': static.y_name, 'units': static.y_unit, 'axis': 'Y'})
                    olat[::] = static.y[:,0]
                all_dims += (static.y_name, )
            if static.x_name:
                if static.x_name not in f.dimensions:
                    f.createDimension(static.x_name, static.x.shape[1])
                    olon = f.createVariable(static.x_name, 'f', (static.x_name,))
                    olon.setncatts({'long_name': static.x_name, 'units': static.x_unit, 'axis': 'X'})
                    olon[::] = static.x[0,:]
                all_dims += (static.x_name, )
        
        try:
            _metsave_vargroup(f, all_dims, dat, static, compress_to_short)
        finally:
            f.close()

        return

    def metsave_composite(dat, composites, static, filename, compress_to_short=True, add_dims=[], global_attrs={}):
        ''' Save composite data in a netCDF file

        Convenience interface to metsave
         (1) transforming the data array as returned from get_composite into a data structure as 
             required by metsave, and then
         (2) saving the data using the standard metsave function.

        In the transformation, ``plevs`` are mapped to netCDF groups, and all specified composites
        are concatenated into a first (prepended) data dimension. This allows skipping through the 
        composites in ncview in the order in which they are specified here.

        NOTE: Composites that are not specified are not saved!

        Parameters
        ----------
        dat : dict, mapping (composite_name, plev, q) => np.ndarray
            Composite data as returned from get_composite. 
        composites : list of decider
            Composites, in order, that should be saved.
        static : gridlib.grid
            Some meta information about the data, like the grid information.
        filename : str
            File name to be saved to.
        compress_to_short : bool
            *Optional*, default ``True``. By default, ``metsave`` compresses the data by converting
            the data field into int16, using the float64 ``add_offset`` and ``scale_factor`` attributes 
            to represent the data.
        add_dims : list of dict
            *Optional*, default empty. List of additional dimensions, in order, to be placed between time 
            and space. The dict must have the keys ``'name'`` and ``'data'``, defining the axis name and 
            its data values. All additional keys will be taken over as attributes for the axis variable.
        global_attrs: 
            *Optional*, default empty. List of additional global attributes to be included in the netCDF file.
        '''
        
        # Transform dat array into tosave
        tosave = {}
        names = [composite.name for composite in composites]
        unsaved = 0
        for (name, plev, q), dat_ in dat.items():
            if name not in names:
                unsaved += 1
                continue
            
            cidx = names.index(name)

            if not plev in tosave:
                tosave[plev] = {}

            if not q in tosave[plev]:
                if not q in conf.q_bins:
                    s = dat_['mean'].shape

                    # Remove vertical dimension if necessary
                    if len(s) == (3+len(add_dims)):
                        if s[-3] > 1:
                            raise NotImplementedError('Cannot yet save data with multiple vertical levels in one array')
                        s = s[:-3]+s[-2:]

                    tosave[plev][q] = np.empty((len(names),)+s) * np.nan
                    tosave[plev][f'{q}_stddev'] = np.empty((len(names),)+s) * np.nan
                    tosave[plev][f'{q}_valid'] = np.empty((len(names),)+s, dtype='i4') * np.nan
                else:
                    raise NotImplementedError('Cannot save binned composites yet.')
            
            if not q in conf.q_bins:
                tosave[plev][q][cidx,::] = dat_['mean']
                tosave[plev][f'{q}_stddev'][cidx,::] = dat_['std']
                tosave[plev][f'{q}_valid'][cidx,::] = dat_['valid_cnt']
            else:
                raise NotImplementedError('Cannot save binned composites yet.')

        
        if unsaved > 0:
            print(f'Warning: {unsaved} composite fields left unsaved!')

        # Remove t timension from static if present
        if static.t_name:
            static_ = static
            static = copy.copy(static_)
            static.t_name = ''

        # Prepend composite name dimension to add_dims
        add_dims = (
            [{'name': 'composite', 'data': np.array(names), 'long_name': 'Composite name sequence'},]
            + add_dims
        )

        # Do the actual saving and return whatever metsave might return (None currently)
        return metsave(tosave, static, filename, compress_to_short, add_dims, global_attrs)

    return metsave, metsave_composite


def get_instantaneous_factory(files_by_plevq, metopen, get_from_file, get_static, conf):
    ''' Create the get_instantaneous function based on data source-specific helpers '''

    def get_instantaneous(plevqs, dates, q_special={}, force=False, **kwargs):
        ''' Get instantaneous fields

        Allows general data requests in the configured data base, e.g. ERA-Interim. The request
        can span several files, e.g. by including several vertical levels or by covering several
        years. The returned data can be up to 4-dimensional, with the dimensions (t,z,y,x). 
        
        The method internally uses metopen to locate data files. Hence, it will find data in the 
        locations given conf.datapath.

        Parameters
        ----------
        plevqs : 2-tuple or list of 2-tuples
            Each 2-tuple consists of (1) a string representations of the requested vertical level(s), 
            e.g. ``'700'`` for 700 hPa or ``'pv2000'`` for the PV2-surface, and (2) a variable name 
            identifier, following the ECMWF conventions as far as applicable, e.g. ``'u'`` or ``'msl'``.
            Some data sets might allow to supply ``'__all__'`` instead of either the vertical level
            and/or the variable name, to request all vertical levels/variables available.
        dates : list of datetime
            The minimum and maxmimum dates in this list define the requested time interval. The i
            final date will not be included in the result, i.e. for all time steps in 2016 request
            dates from 2016-01-01 00:00 to 2017-01-01 00:00.
        q_special : dict of callable
            *Optional*, default ``{}``. If given for a specific variable q, this function is used instead 
            of get_from_file to load the data. This allows, for example, to run simple transformations 
            before aggregating the data.
        force : bool
            *Optional*, default ``False``. Turn off the error, if large amounts of data are
            requested at once. **Be sure you know what you are doing, when setting this to 
            ``True``! Your request might make your script occupy a large fraction of the 
            system memory**.
        
        Keyword arguments
        -----------------
        metopen arguments : all optional
            Optional arguments passed on to calls of metopen within this function.

        Returns
        -------
        dict of np.ndarray
            Data for the requested variable.
        grid.gridlib
            If ``no_static=False`` meta-information about the requested data, otherwise ``None``.
        '''
        
        start, end = min(dates), max(dates)

        # Make sure you get something when only requesting one time step
        if start == end:
            end += td(0,1)

        if type(plevqs) == tuple:
            plevqs = [plevqs, ]

        # Dry run to test whether input parameters are valid and to determine the total amount of requested data
        req = {}            # Request info
        datshape = {}       # Shape of the resulting data arrays
        request_size = 0    # Total size of all requested data 
        for plevq in plevqs:
            # req[plevq] contains a list of 4-tuples: 
            # (filename, list of tidx, list of dates, request_size in number of values)
            req[plevq] = list(files_by_plevq(plevq, start=start, end=end))
            datshape[plevq] = req[plevq][0][3]
            for entry in req[plevq][1:]:
                shape = entry[3]
                if not shape[1:] == datshape[plevq][1:]:
                    raise ValueError(f'''Discovered inconsistent data shape across time:
                            plevq: {plevq}
                            file {entry[0]} with shape {shape[1:]}, 
                            preceeding files with shape {datshape[plevq][1:]}.''')
                datshape[plevq] = (datshape[plevq][0] + shape[0],) + shape[1:]
            request_size += sum([np.prod(reqinfo[3]) for reqinfo in req[plevq]])

        # Checking max length
        if request_size > MAX_REQUEST_SIZE:
            if force:
                print(f'Warning: you requested {request_size*8 / 1024**3:.2f}G of data.')
            else:
                raise ValueError(f'''Cowardly refusing to fetch {request_size*8 / 1024**3:.2f}G of data at once.
                    If you are absolutely certain what you're doing, you can use force=True to override.''')

        if not force and request_size > WARN_REQUEST_SIZE:
            print(f'Warning: you requested {request_size*8 / 1024**3:.2f}G of data.')
        
        # Remove no_static if present, this should only effect the output of this function rather than the 
        # metopen calls herein; also static info is sometimes required internally herein.
        no_static = kwargs.pop('no_static', False)
        grid = get_static()

        dat = {}
        dates = {}
        for plev, q in plevqs:
            # Allocate memory for the result
            dat[plev,q] = np.empty(datshape[plev,q])
            dates[plev,q] = []
            
            toff = 0
            for filename, tidxs, dates_, shape in req[plev,q]:
                cut = slice(tidxs[0], tidxs[-1]+1)
                tlen = len(tidxs)
                
                if q in q_special:
                    dat_ = q_special[q](get_from_file, filename, plev, q, cut=cut, no_static=True, **kwargs)
                else:
                    dat_ = get_from_file(filename, plev, q, cut=cut, no_static=True, **kwargs)

                    # Treat special data
                    #
                    #  1. Lines
                    if q in conf.q_lines:
                        ql = conf.q_lines[q]
                        datoff_ = get_from_file(filename, plev, ql, cut=cut, no_static=True, **kwargs)
                        dat_ = utils.normalize_lines(dat_, datoff_, grid.dx, grid.dy)[:,np.newaxis,:,:]

                dat[plev,q][toff:toff+tlen,...] = dat_[...]
                dates[plev,q].extend(dates_)
                
                toff += tlen
        
        if len(dates) > 1:
            # Iterate through all time axes at in parallel, see if all indexes refer to the respective same date
            for dates_tstep in zip(*dates.values()):
                for date in dates_tstep[1:]:
                    if not date == dates_tstep[0]:
                        raise ValueError('Time axes not consistent across the requested variables')
        
        # Prepare grid information, if so requested
        if no_static:
            grid = None
        else:
            grid = grid.new_time(dates[plev,q])
            
        return dat, grid
    
    return get_instantaneous



def get_at_position_factory(files_by_plevq, get_normalized_from_file):
    ''' Create the get_at_position funciton based on data source-specific helpers '''

    def get_at_position(dates, plevs, ys, xs, q):
        ''' Get values of the variable q at given dates, interpolated to given plevs, ys, xs

        The given plevs, ys and xs must all have the same 4-dimensional shape. Dates must be
        one-dimensional, corresponding in length to the first dimension of plevs, ys and xs.
        The returned array will have the same shape as plevs, ys and xs.
        
        Parameters
        ----------
        dates : list/np.ndarray of datetime with dimensions (t,)
            The dates for which data is to be interpolated.
        plevs : np.ndarray with dimensions (z,)
            Pressure levels to be interpolated to.
        ys : np.ndarray with dimensions (t,y,x)
            Meridional distances/latitudes to be interpolated to.
        xs : np.ndarray with dimensions (t,y,x)
            Zonal distances or longitudes to be interpolated to.
        q : str
            A variable name identifier, following the ECMWF conventions as far as applicable,
            e.g. ``'u'`` or ``'msl'``.

        Returns
        -------
        np.ndarray with dimensions (t,z,y,x)
            Interpolated values of variable q at the givens dates/positions.
        '''

        # Consistency checks
        if not xs.shape == ys.shape:
            raise ValueError('xs and ys must have the same shape.')
            
        if not xs.shape[0] == len(dates):
            raise ValueError('Number of dates must be equal length of first dimension of xs and ys.')
        
        # Assuming the chunking of data into files is consistent across plevs, this is a natural 
        # rhythm to iterate through both test and composite data
        #
        # (the below code does not require data to be chunked consistently, 
        #  but the code will much more efficient if it is)
        if type(dates) == list:
            dates = np.array(dates)

        start, end = dates.min(), dates.max() + td(0,1)
        req = list(files_by_plevq((plevs[0], q), start=start, end=end))
        
        # Allocate array to hold the results
        dat = np.empty(xs.shape[:1] + (len(plevs),) + (xs.shape[1:]))
        
        for pidx,plev in enumerate(plevs):

            # Iterate through the data set in chunks natural to the data source
            for filename, tidxs, dates_, shape in req:
                # Skip chunks from which no time step is required
                load_chunk = False
                for date in dates_:
                    if date in dates:
                        load_chunk = True

                if not load_chunk:
                    continue
                
                dat_, grid = get_normalized_from_file(filename, plev, q)
            
                for tidx_in, date in enumerate(dates_):
                    tidxs_out = np.argwhere(dates == date)
                    if tidxs_out.size == 0:
                        continue

                    dat__ = dat_[tidx_in,:,:].squeeze()
                    if np.isnan(dat__).sum() > 0:
                        zonalavg = np.nanmean(dat__, axis=1)
                        dat__ -= zonalavg[:,np.newaxis]
                        dat__, conv = utils.fill_nan(dat__[np.newaxis,:,:])
                        dat__ = dat__[0,:,:]
                        dat__ += zonalavg[:,np.newaxis]

                    ifunc = interp.RectBivariateSpline(grid.y[::-1,0], grid.x[0,:], dat__[::-1,:], kx=1, ky=1)
                    for tidx_out in tidxs_out:
                        tidx_out = tidx_out[0]
                        dat[tidx_out,pidx,...] = ifunc(ys[tidx_out,...], xs[tidx_out,...], grid=False)


        return dat

    
    def get_hor_interpolation_functions(dates, plevs, q):
        ''' Get interpolation functions of for all given plevs of given q

        The given plevs, ys and xs must all have the same 4-dimensional shape. Dates must be
        one-dimensional, corresponding in length to the first dimension of plevs, ys and xs.
        The returned array will have the same shape as plevs, ys and xs.
        
        Parameters
        ----------
        dates : list/np.ndarray of datetime with dimensions (t,)
            The dates for which data is to be interpolated.
        plevs : np.ndarray with dimensions (z,)
            Pressure levels to be interpolated to.
        q : str
            A variable name identifier, following the ECMWF conventions as far as applicable,
            e.g. ``'u'`` or ``'msl'``.

        Returns
        -------
        dict (date, plev) => interpolation function
            The requested interpolation accessible by date and plev.
        '''

        # TODO: Reduce code duplicaiton with get_at_position!

        # Assuming the chunking of data into files is consistent across plevs, this is a natural 
        # rhythm to iterate through both test and composite data
        #
        # (the below code does not require data to be chunked consistently, 
        #  but the code will much more efficient if it is)
        if type(dates) == list:
            dates = np.array(dates)

        # Allocate structure to hold the results
        ifuncs = {}

        start, end = dates.min(), dates.max() + td(0,1)

        for pidx,plev in enumerate(plevs):
            req = list(files_by_plevq((plev, q), start=start, end=end))
            
            # Iterate through the data set in chunks natural to the data source
            for filename, tidxs, dates_, shape in req:
                # Skip chunks from which no time step is required
                load_chunk = False
                for date in dates_:
                    if date in dates:
                        load_chunk = True

                if not load_chunk:
                    continue
                
                dat_, grid = get_normalized_from_file(filename, plev, q)
            
                for tidx_in, date in enumerate(dates_):
                    tidxs_out = np.argwhere(dates == date)
                    if tidxs_out.size == 0:
                        continue

                    dat__ = dat_[tidx_in,:,:].squeeze()
                    if np.isnan(dat__).sum() > 0:
                        zonalavg = np.nanmean(dat__, axis=1)
                        dat__ -= zonalavg[:,np.newaxis]
                        dat__, conv = utils.fill_nan(dat__[np.newaxis,:,:])
                        dat__ = dat__[0,:,:]
                        dat__ += zonalavg[:,np.newaxis]
                    
                    if not date in ifuncs:
                        ifuncs[date] = {}
                    ifuncs[date][plev] = interp.RectBivariateSpline(grid.y[::-1,0], grid.x[0,:], 
                            dat__[::-1,:], kx=1, ky=1)

        return ifuncs
    
    return get_at_position, get_hor_interpolation_functions


def get_time_average_factory(files_by_plevq, get_normalized_from_file, get_static, conf):
    ''' Create the get_time_average function based on data source-specific helpers '''

    def get_time_average(plevqs, dates, q_special={}, **kwargs):
        ''' Get time-average fields

        Allows general data requests in the configured data base, e.g. ERA-Interim. The request
        can span several files, e.g. by including several vertical levels or by covering several
        years. The returned data can be up to 4-dimensional, with the dimensions (t,z,y,x). 
        
        The method internally uses metopen to locate data files. Hence, it will find data in the 
        locations given conf.datapath.

        Parameters
        ----------
        plevqs : 2-tuple or list of 2-tuples
            Each 2-tuple consists of (1) a string representations of the requested vertical level(s), 
            e.g. ``'700'`` for 700 hPa or ``'pv2000'`` for the PV2-surface, and (2) a variable name 
            identifier, following the ECMWF conventions as far as applicable, e.g. ``'u'`` or ``'msl'``.
            Some data sets might allow to supply ``'__all__'`` instead of either the vertical level
            and/or the variable name, to request all vertical levels/variables available.
        dates : list of datetime
            The minimum and maxmimum dates in this list define the requested time interval. The i
            final date will not be included in the result, i.e. for an average of 2016 request
            dates from 2016-01-01 00:00 to 2017-01-01 00:00.
        q_special : dict of callable
            *Optional*, default ``{}``. If given for a specific variable q, this function is used instead 
            of get_normalized_from_file to load the data. This allows, for example, to run simple transformations 
            before aggregating the data.
        
        Keyword arguments
        -----------------
        metopen arguments : all optional
            Optional arguments passed on to calls of metopen within this function.

        Returns
        -------
        dict of np.ndarray
            Time average for the requested variable and vertical level(s).
        grid.gridlib
            If ``no_static=False`` meta-information about the requested data, otherwise ``None``.
        '''

        start, end = min(dates), max(dates)

        if type(plevqs) == tuple:
            plevqs = [plevqs, ]

        # Dry run to test whether input parameters are valid and to determine the total amount of requested data
        req = {}            # Request info
        datshape = {}       # Shape of the resulting data arrays
        for plevq in plevqs:
            # req[plevq] contains a list of 4-tuples: 
            # (filename, list of tidx, list of dates, shape)
            req[plevq] = list(files_by_plevq(plevq, start=start, end=end))
            datshape[plevq] = req[plevq][0][3][1:]
            for entry in req[plevq][1:]:
                shape = entry[3]
                if not shape[1:] == datshape[plevq]:
                    raise ValueError(f'''Discovered inconsistent data shape across time:
                            plevq: {plevq}
                            file {entry[0]} with shape {shape[1:]}, 
                            preceeding files with shape {datshape[plevq]}.''')

        # Remove no_static if present, this should only effect the output of this function rather than the metopen calls herein
        no_static = kwargs.pop('no_static', False)

        dat = {}
        dates = {}
        grid = get_static()
        for plev, q in plevqs:
            qout = conf.q_avg.get(q, q)

            # Allocate memory for the result
            dat[plev,qout] = np.zeros(datshape[plev,q])
            if q in conf.q_bins:
                dat[plev,qout,'hist'] = np.zeros((len(conf.q_bins[q]),) + datshape[plev,q], dtype='i4')
            else:
                dat[plev,qout,'valid'] = np.zeros(datshape[plev,q], dtype='i4')
            
            toff = 0
            for filename, tidxs, dates_, shape in req[plev,q]:
                cut = slice(tidxs[0], tidxs[-1]+1)
                tlen = len(tidxs)
                
                # Allow special function to retrieve data
                if q in q_special:
                    dat_ = q_special[q](get_normalized_from_file, filename, plev, q, cut=cut, no_static=True, **kwargs)

                # Else read data in one of the standard ways
                else:
                    dat_ = get_normalized_from_file(filename, plev, q, cut=cut, no_static=True, **kwargs)

                    # Treat special data (lines and object masks are already taken care of)
                    #
                    #  Binned data
                    if q in conf.q_bins:
                        for bi in range(len(conf.q_bins[q])-1):
                            upper, lower = conf.q_bins[q][bi+1], conf.q_bins[q][bi]
                            if upper > lower:
                                dat[plev,qout,'hist'][bi,:,:] += np.logical_and(dat_ >= lower, dat_ <  upper).sum(axis=0)
                            else:
                                dat[plev,qout,'hist'][bi,:,:] += np.logical_or(dat_ <  upper, dat_ >= lower).sum(axis=0)
                
                    # Summing up non-binned data, taking care of NaNs
                    else:
                        mask = np.isnan(dat_)
                        dat[plev,qout] += dat_.sum(axis=0, where=~mask)
                        dat[plev,qout,'valid'] += tlen - mask.sum(axis=0)
                
                toff += tlen

            dat[plev,qout,'cnt'] = toff

            if q in conf.q_bins:
                dat[plev,qout] = utils.cal_mfv(dat[plev,qout,'hist'], conf.q_bins[q])
            else:
                dat[plev,qout] /= dat[plev,qout,'valid']
        
        # Prepare grid information, if so requested
        if no_static:
            grid = None
        else:
            grid = grid.new_time([start, ])
            grid.t_periods = [(start, end), ]
            
        return dat, grid
    
    return get_time_average


def get_aggregate_factory(files_by_plevq, get_normalized_from_file, get_static, conf):
    ''' Create the get_time_average function based on data source-specific helpers '''

    def get_aggregate(plevqs, dates, agg, q_special={}, **kwargs):
        ''' Get time-aggregates, such as monthly means

        Allows general data requests in the configured data base, e.g. ERA-Interim. The request
        can span several files, e.g. by including several vertical levels or by covering several
        years. The returned data can be up to 4-dimensional, with the dimensions (t,z,y,x). 
        
        The method internally uses metopen to locate data files. Hence, it will find data in the 
        locations given conf.datapath.

        Parameters
        ----------
        plevqs : 2-tuple or list of 2-tuples
            Each 2-tuple consists of (1) a string representations of the requested vertical level(s), 
            e.g. ``'700'`` for 700 hPa or ``'pv2000'`` for the PV2-surface, and (2) a variable name 
            identifier, following the ECMWF conventions as far as applicable, e.g. ``'u'`` or ``'msl'``.
            Some data sets might allow to supply ``'__all__'`` instead of either the vertical level
            and/or the variable name, to request all vertical levels/variables available.
        dates : list of datetime
            The minimum and maxmimum dates in this list define the requested time interval. The i
            final date will not be included in the result, i.e. for all time steps in 2016 request
            dates from 2016-01-01 00:00 to 2017-01-01 00:00.
        agg : str or dynlib.tagg.agg object
            String representation (e.g. ``'cal_month'``, ``'pentad'``, or ``'met_season'``), or 
            time aggregator object representation the aggregation interval.
        q_special : dict of callable
            *Optional*, default ``{}``. If given for a specific variable q, this function is used instead 
            of get_normalized_from_file to load the data. This allows, for example, to run simple transformations 
            before aggregating the data.
        
        Keyword arguments
        -----------------
        metopen arguments : all optional
            Optional arguments passed on to calls of metopen within this function.

        Returns
        -------
        dict of np.ndarray
            Time-aggregated data for the requested variable(s) and vertical level(s).
        grid.gridlib
            If ``no_static=False`` meta-information about the requested data, otherwise ``None``.
        '''

        start, end = min(dates), max(dates)

        if type(plevqs) == tuple:
            plevqs = [plevqs, ]

        # Dry run to test whether input parameters are valid and to determine the total amount of requested data
        req = {}            # Request info
        datshape = {}       # Shape of the resulting data arrays
        for plevq in plevqs:
            # req[plevq] contains a list of 4-tuples: 
            # (filename, list of tidx, list of dates, shape)
            req[plevq] = list(files_by_plevq(plevq, start=start, end=end))
            datshape[plevq] = req[plevq][0][3][1:]
            for entry in req[plevq][1:]:
                shape = entry[3]
                if not shape[1:] == datshape[plevq]:
                    raise ValueError(f'''Discovered inconsistent data shape across time:
                            plevq: {plevq}
                            file {entry[0]} with shape {shape[1:]}, 
                            preceeding files with shape {datshape[plevq]}.''')

        # Remove no_static if present, this should only effect the output of this function rather than the metopen calls herein
        no_static = kwargs.pop('no_static', False)
        grid = get_static() # For conf.q_lines and potentially for output

        dat = {}
        start_dates = {}
        end_dates = {}
        for plev, q in plevqs:
            qout = conf.q_avg.get(q, q)

            # Allocate memory for the result
            dat[plev,qout] = []
            if q in conf.q_bins:
                dat[plev,qout,'hist'] = []
            else:
                dat[plev,qout,'valid'] = []
            dat[plev,qout,'cnt'] = []

            start_dates[plev,q] = []
            end_dates[plev,q] = []

            leftover = None
            for filename, tidxs, dates_, shape in req[plev,q]:
                cut = slice(tidxs[0], tidxs[-1]+1)
                tlen = len(tidxs)
                
                # Allow special function to retrieve data
                if q in q_special:
                    dat_ = q_special[q](get_normalized_from_file, filename, plev, q, cut=cut, no_static=True, **kwargs)

                # Else read data in one of the standard ways
                else:
                    dat_ = get_normalized_from_file(filename, plev, q, cut=cut, no_static=True, **kwargs)

                # Prepend leftover data from the previous file period, if available
                if leftover:
                    dates_ = leftover[0] + dates_
                    dat_ = np.concatenate((leftover[1], dat_), axis=0)
                    leftover = None
                
                # Aggregate data
                start_dates_agg, end_dates_agg, mean_or_mfv_agg, valid_or_hist_agg, cnt_agg = \
                        utils.aggregate(dates_, dat_, agg, bins=conf.q_bins.get(q, None))
                
                # Save results to output arrays (lists for now, concatenated below)
                dat[plev,qout].append(mean_or_mfv_agg)
                if q in conf.q_bins:
                    dat[plev,qout,'hist'].append(valid_or_hist_agg)
                else:
                    dat[plev,qout,'valid'].append(valid_or_hist_agg)
                dat[plev,qout,'cnt'].append(cnt_agg)

                start_dates[plev,q].extend(start_dates_agg)
                end_dates[plev,q].extend(end_dates_agg)
                
                # Some time steps were not (yet) incorporated in an aggregation period. Save them for the next one.
                if end_dates_agg[-1] in dates_:
                    next_tidx = dates_.index(end_dates_agg[-1])
                    leftover = dates_[next_tidx:], dat_[next_tidx:]
            
            # Concatenate the list of arrays to one array
            dat[plev,qout] = np.concatenate(dat[plev,qout], axis=0)
            if q in conf.q_bins:
                dat[plev,qout,'hist'] = np.concatenate(dat[plev,qout,'hist'], axis=0)
            else:
                dat[plev,qout,'valid'] = np.concatenate(dat[plev,qout,'valid'], axis=0)
            dat[plev,qout,'cnt'] = np.concatenate(dat[plev,qout,'cnt'], axis=0)
        
        if len(start_dates) > 1:
            # Iterate through all time axes at in parallel, see if all indexes refer to the respective same date
            for dates_tstep in zip(*start_dates.values()):
                for date in dates_tstep[1:]:
                    if not date == dates_tstep[0]:
                        raise ValueError('Time axes not consistent across the requested variables')
        
        # Prepare grid information, if so requested
        if not no_static:
            grid = grid.new_time(start_dates[plev,q])
            grid.t_periods = list(zip(start_dates[plev,q], end_dates[plev,q]))
            
        return dat, grid
 
    return get_aggregate


def get_composite_factory(files_by_plevq, get_normalized_from_file, get_static, conf):
    ''' Create the get_time_average function based on data source-specific helpers '''

    def _add_chunk(plevq, dat, comp_ts, dat_to_add):
        ''' Internal helper function, adding data for a chunk of time to the relevant composites 
        
        The dat structure of nested dict is modified in-place.
        '''
        
        plev, q = plevq
        qs = conf.q_avg.get(q, q)
        s = dat_to_add.shape[1:]

        if not q in conf.q_bins:
            for name, ts in comp_ts.items():
                if 'mean' not in dat[name,plev,qs]:
                    dat[name,plev,qs]['mean'] = np.zeros(s)
                    dat[name,plev,qs]['valid_cnt'] = np.zeros(s)
                    dat[name,plev,qs]['std'] = np.zeros(s)
                for tidx in np.argwhere(ts)[:,0]:
                    dat_ = dat_to_add[tidx,::]
                    nanmask = ~np.isnan(dat_)
                    dat[name,plev,qs]['mean'][nanmask] += dat_[nanmask]
                    dat[name,plev,qs]['std'][nanmask] += dat_[nanmask]**2
                    dat[name,plev,qs]['valid_cnt'] += nanmask
        
        else:
            for name, ts in comp_ts.items():
                if 'hist' not in dat[name,plev,qs]:
                    dat[name,plev,qs]['hist'] = np.zeros((len(conf.q_bins[q]),)+s)
                    dat[name,plev,qs]['mfv'] = np.zeros(s)
                for tidx in np.argwhere(ts)[:,0]:
                    for bi in range(len(conf.q_bins[q])-1):
                        upper = conf.q_bins[q][bi+1]
                        lower = conf.q_bins[q][bi]
                        if upper > lower:
                            dat[name,plev,qs]['hist'][bi,(dat >= lower).__and__(dat <  upper)] += 1
                        else:
                            dat[name,plev,qs]['hist'][bi,(dat <  upper).__or__(dat >= lower)] += 1

        return
    

    def get_composite(plevqs, dates, composites, **kwargs):
        ''' Get composites, such as multi-year seasonal averages or for NAO+/-

        Allows general data requests in the configured data base, e.g. ERA-Interim. The request
        can span several files, e.g. by including several vertical levels or by covering several
        years. For each composite, the returned data can be up to 3-dimensional, with the 
        dimensions (z,y,x). 
        
        The method internally uses metopen to locate data files. Hence, it will find data in the 
        locations given conf.datapath.

        Parameters
        ----------
        plevqs : 2-tuple or list of 2-tuples
            Each 2-tuple consists of (1) a string representations of the requested vertical level(s), 
            e.g. ``'700'`` for 700 hPa or ``'pv2000'`` for the PV2-surface, and (2) a variable name 
            identifier, following the ECMWF conventions as far as applicable, e.g. ``'u'`` or ``'msl'``.
            Some data sets might allow to supply ``'__all__'`` instead of either the vertical level
            and/or the variable name, to request all vertical levels/variables available.
        dates : list of datetime
            The minimum and maxmimum dates in this list define the requested time interval. The
            final date will not be included in the result, i.e. for all time steps in 2016 request
            dates from 2016-01-01 00:00 to 2017-01-01 00:00.
        composites : composite_test or list of composite_test
            The compositing criteria to be used.
        
        Keyword arguments
        -----------------
        metopen arguments : all optional
            Optional arguments passed on to calls of metopen within this function.

        Returns
        -------
        dict of np.ndarray
            Composite data for the requested variable(s), vertical level(s), and composites.
        '''

        # TODO: Generalise both this function and the composite decider infrastructure to allow returning
        #       Grid information per time step to be taken into the composite

        # Find out if we need any test data and compile an overall set
        #test_plevqs = set([])
        #for composite in composites:
        #    if type(composite.requires) == tuple:
        #        test_plevqs.add(composite.requires)
        
        # Assuming the chunking of data into files is consistent across plevqs, this is a natural rhythm to
        # iterate through both test and composite data
        #
        # (the below code does not require data to be chunked consistently, 
        #  but the code will much more efficient if it is)
        start, end = min(dates), max(dates)
        req = list(files_by_plevq(plevqs[0], start=start, end=end))
        
        # Nested structure to hold the results
        dat = {}
        for composite in composites:
            for plev, q in plevqs:
                qs = conf.q_avg.get(q, q)
                if (composite.name,plev,qs) in dat:
                    raise ValueError(f'Duplicate composite name `{composite.name}`.')
                dat[composite.name,plev,qs] = {}
                    

        for filename, tidxs, dates_, shape in req:
            # 1a. For each composite construct time series
            to_include = {}
            load_chunk = False
            for composite in composites:
                if type(composite.requires) == tuple:
                    # Inject relevant functions to get test data from this data source
                    ts = composite.get_time_series(dates_, files_by_plevq, get_normalized_from_file)
                else:
                    ts = composite.get_time_series(dates_)

                to_include[composite.name] = ts

                # 1b. Evaluate whether any composite requires data from the current chunk
                load_chunk = load_chunk or np.any(ts) 

            if not load_chunk:
                continue

            # Reminder: We need to construct mean, validcnt, std?, hist, mfv, tscnt

            # 2. If yes, request data from the current chunk for all variables, and
            # 3. Add relevant time steps to each composite
            dat_ = get_normalized_from_file(filename, *plevqs[0], no_static=True)
            _add_chunk(plevqs[0], dat, to_include, dat_)
            
            # ... and the same for all the other variables/levels
            for plevq in plevqs[1:]:
                req_ = list(files_by_plevq(plevq, start=min(dates_), end=max(dates_)+td(0,1) )) # the end date should here be included in the requests for this chunk
                for filename_, tidxs, dates_, shape in req_:
                    dat_ = get_normalized_from_file(filename_, *plevq, no_static=True)
                    _add_chunk(plevq, dat, to_include, dat_)
        
        # A bit of post-processing, caculating mean+standard deviation / most-frequent-value
        for composite in composites:
            for plev, q in plevqs:
                qs = conf.q_avg.get(q, q)
                dat_ = dat[composite.name,plev,qs]
                if not q in conf.q_bins:
                    dat_['mean'] /= dat_['valid_cnt']
                    dat_['std'] = np.sqrt((dat_['std'] - dat_['valid_cnt']*dat_['mean']**2)/(dat_['valid_cnt']-1))
                else:
                    # the hist itself does not require any post-processing
                    dat_['mfv'] = dynlib.utils.cal_mfv(dat_['hist'], conf.q_bins[q])

        return dat
    
    return get_composite


# C'est le fin
