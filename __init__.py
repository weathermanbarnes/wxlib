#!/usr/bin/env python
# -*- encoding: utf-8

from copy import deepcopy

import dynlib

import dynpie.figures as fig
import dynpie.gridlib as gridlib
import dynpie.utils as utils
from dynpie.metopen import metopen
from dynpie.settings import conf 

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Clean up a bit
default_conf = deepcopy(conf)
del conf, deepcopy


# that's it
