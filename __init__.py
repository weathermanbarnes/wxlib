#!/usr/bin/env python
# -*- encoding: utf-8

from copy import deepcopy

import dynlib

import dynpie.figures as fig
import dynpie.gridlib as gridlib
import dynpie.utils as utils
import dynpie.settings as settings
from dynpie.metopen import metopen

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Clean up a bit
default_conf = deepcopy(settings.conf)
del deepcopy


# that's it
