#!/usr/bin/env python
# -*- encoding: utf-8


from . import *
from ..erainterim import files_by_plevq, get_static, get_from_file

# Create ERA-Interim specific version of decide_by_data
decide_by_data = decide_by_data_factory(files_by_plevq, get_static, get_from_file)


# C'est le fin
