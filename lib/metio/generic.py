#!/usr/bin/env python
# -*- encoding: utf-8


from .datasource import *
from ..settings import settings_obj, default_conf


dt = cftime.DatetimeGregorian
conf = settings_obj({
    'q': {},
    'qf': {}, 
    'q_std': {},
    'q_avg': {},
    'q_units': {},
    'q_long': {},
    'q_bins': {},
    'q_lines': {},
    'q_feature_dists': {},
    'q_obj': {},
    'datapath': ['.',], 
    'opath': '.',
    'oformat': 'nc',
    'staticfile': None,
    'epoch': dt(1900,1,1,0),
    'calendar': None,
    'timestep': None,
    'gridsize': (-1,-1),
    'local_timezone': default_conf.local_timezone,
}, [])


# Derive data source-specific version of metopen
metopen = metopen_factory(lambda: None, conf)
metsave, metsave_composite = metsave_factory(metopen, conf)



# C'est le fin
