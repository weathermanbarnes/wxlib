#!/usr/bin/env python
# -*- encoding: utf-8

from __future__ import absolute_import

# Importing some widely used functions in a common name space

from . import dynfor

from .metio.generic import metopen, metsave, conf
from . import figures as fig
from . import cm
from . import proj

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from datetime import timedelta as td
from cftime import DatetimeGregorian as dt



# that's it
