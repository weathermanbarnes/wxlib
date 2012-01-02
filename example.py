#!/usr/bin/python

from dynlib import *
import numpy as np

u   = np.ones((3,50,50))
v   = np.zeros((3,50,50))
dx  = np.ones((50,50))
dy  = np.ones((50,50))

div = diag.div(u,v,dx,dx)
print div.shape

# /the/end
