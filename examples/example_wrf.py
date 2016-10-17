#!/usr/bin/env python
# -*- encoding: utf8

from lib.shorthands import *
from lib.settings import *


conf.datapath.insert(1, '/Data/gfi/spengler/ate093/ideal_SFX/data')


f, static = metopen('cyclonePLEV_SLP_D6_SFX')

print static.t
print static.t_parsed
print static.x
print static.y

print static.dx

print f.variables.keys()

#
