#!/usr/bin/env python
# -*- encoding: utf-8

from __future__ import absolute_import

# Importing some widely used functions in a common name space

from . import dynfor

from . import figures as fig
from . import gridlib
from . import utils
from . import settings
from .metio import metopen, metsave, metsave_lines, metsave_timeless, get_aggregate, get_instantaneous, get_static

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from datetime import datetime as dt, timedelta as td


# Saving the default configuration
from copy import deepcopy
default_conf = deepcopy(settings.conf)
del deepcopy


# that's it
