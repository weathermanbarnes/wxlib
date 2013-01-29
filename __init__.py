#!/usr/bin/env python
# -*- encoding: utf-8

from copy import deepcopy

import dynlib

import dynpie.figures as fig
from dynpie.metopen import metopen
from dynpie.settings import conf 


default_conf = deepcopy(conf)
del conf, deepcopy


# that's it
