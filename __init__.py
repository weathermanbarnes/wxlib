#!/usr/bin/env python
# -*- encoding: utf-8

from copy import deepcopy

import dynlibdoc as dynlib

import dynpie.figures as fig
import dynpie.gridlib as gridlib
import dynpie.utils as utils
import dynpie.settings as settings
from dynpie.metopen import *

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from datetime import datetime as dt, timedelta as td


dynlib_version = (''.join(dynlib.consts.version)).strip()

# Clean up a bit
default_conf = deepcopy(settings.conf)
del deepcopy


# that's it
