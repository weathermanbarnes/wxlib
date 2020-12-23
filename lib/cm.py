#!/usr/bin/env python
# -*- encoding: utf-8


''' A collection of useful color maps

These color maps have been devised for some special meteorological applications, 
where they fit better than the matplotlib standard color maps. Some of these
color maps were originally devised for iveret.gfi.uib.no.

The standard matplotlib color maps are also available through this module.
'''

from __future__ import absolute_import, unicode_literals


import matplotlib as mpl
mpl.use('Agg')

# Make the standard matplotlib colorbars available in the same name space
from matplotlib.cm import *



def greys():
    ''' Colorbar from white to 30% grey 

    Returns
    -------
    mpl.colors.LinearSegmentedColormap
        Colormap
    
    See also
    --------
    :meth:`greys_r`
    '''
    cdict = {'red':   ((0.0, 1.0, 1.0), (1.0, 0.3, 0.3)),
         'green': ((0.0, 1.0, 1.0), (1.0, 0.3, 0.3)),
         'blue':  ((0.0, 1.0, 1.0), (1.0, 0.3, 0.3))  }

    return mpl.colors.LinearSegmentedColormap('greys',cdict,256)

def greys_r():
    ''' Colorbar from 30% grey to white

    Returns
    -------
    mpl.colors.LinearSegmentedColormap
        Colormap
    
    See also
    --------
    :meth:`greys`
    '''
    cdict = {'red':   ((0.0, 0.3, 0.3), (1.0, 1.0, 1.0)),
         'green': ((0.0, 0.3, 0.3), (1.0, 1.0, 1.0)),
         'blue':  ((0.0, 0.3, 0.3), (1.0, 1.0, 1.0))  }

    return mpl.colors.LinearSegmentedColormap('my_grey_r',cdict,256)

def greys2():
    ''' Colorbar from white to 10% grey 

    Returns
    -------
    mpl.colors.LinearSegmentedColormap
        Colormap
    '''    
    cdict = {'red':   ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1)),
         'green': ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1)),
         'blue':  ((0.0, 1.0, 1.0), (1.0, 0.1, 0.1))  }

    return mpl.colors.LinearSegmentedColormap('greys2',cdict,256)


def div_bw():
    ''' Colorbar from 60% grey to white to 20% grey

    Divergent colorbar for black-white figures.

    Returns
    -------
    mpl.colors.LinearSegmentedColormap
        Colormap
    
    See also
    --------
    :meth:`div_bw_inv`
    '''    
    cdict = {'red':   ((0.0, 0.6, 0.6), (0.5, 1.0, 1.0), (1.0, 0.2, 0.2)),
         'green': ((0.0, 0.6, 0.6), (0.5, 1.0, 1.0), (1.0, 0.2, 0.2)),
         'blue':  ((0.0, 0.6, 0.6), (0.5, 1.0, 1.0), (1.0, 0.2, 0.2))  }

    return mpl.colors.LinearSegmentedColormap('div_grey',cdict,256)

def div_greys():
    ''' Symmetric diverging colorbar from 50% grey to white to 50% grey

    Divergent colorbar to be used with additional indicator of sign.

    Returns
    -------
    mpl.colors.LinearSegmentedColormap
        Colormap
    
    See also
    --------
    :meth:`div_bw`
    '''    
    cdict = {'red':   ((0.0, 0.5, 0.5), (0.5, 1.0, 1.0), (1.0, 0.5, 0.5)),
         'green': ((0.0, 0.5, 0.5), (0.5, 1.0, 1.0), (1.0, 0.5, 0.5)),
         'blue':  ((0.0, 0.5, 0.5), (0.5, 1.0, 1.0), (1.0, 0.5, 0.5))  }

    return mpl.colors.LinearSegmentedColormap('div_grey',cdict,256)

def div_bw_r():
    ''' Colorbar from 60% grey to 20% grey to white

    Divergent colorbar for black-white figures.

    Returns
    -------
    mpl.colors.LinearSegmentedColormap
        Colormap
    
    See also
    --------
    :meth:`div_bw`
    '''    
    cdict = {'red':   ((0.0, 0.6, 0.6), (0.5, 0.2, 0.2), (1.0, 1.0, 1.0)),
         'green': ((0.0, 0.6, 0.6), (0.5, 0.2, 0.2), (1.0, 1.0, 1.0)),
         'blue':  ((0.0, 0.6, 0.6), (0.5, 0.2, 0.2), (1.0, 1.0, 1.0))  }

    return mpl.colors.LinearSegmentedColormap('div_grey_r',cdict,256)


def defabs():
    ''' Colorbar from white to 15% grey to red

    Often used for total deformation and jet axis frequencies

    Returns
    -------
    mpl.colors.LinearSegmentedColormap
        Colormap
    '''
    cdict = {'red':   ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 1.0, 1.0)),
         'green': ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 0.0, 0.0)),
         'blue':  ((0.0, 1.0, 1.0), (0.75, 0.15, 0.15), (1.0, 0.0, 0.0))  }

    return mpl.colors.LinearSegmentedColormap('defabs',cdict,256)


def q():
    ''' Colorbar from white to dark blue to violet

    Often used for humidity, rain measures, etc.

    Returns
    -------
    mpl.colors.LinearSegmentedColormap
        Colormap
    '''
    cdict = {'red':   ((0.0, 1.0, 1.0), (0.33, 0.30, 0.30),  (0.867, 0.1, 0.1), (1.0, 0.5, 0.5)),
         'green': ((0.0, 1.0, 1.0), (0.33, 0.65, 0.65),  (0.867, 0.2, 0.2), (1.0, 0.2, 0.2)),
         'blue':  ((0.0, 1.0, 1.0), (0.33, 0.80, 0.80),  (0.867, 0.6, 0.6), (1.0, 0.8, 0.8))  }

    return mpl.colors.LinearSegmentedColormap('q',cdict,256)

def periodic():
    ''' Colorbar from white to blue to 50% grey to yellow back to white

    Useful for periodic variablies like the deformation angle or the wind direction.

    Returns
    -------
    mpl.colors.LinearSegmentedColormap
        Colormap
    '''
    cdict = {'red':   ((0.0, 1.0, 1.0), (0.25, 0.3, 0.2), (0.5, 0.5, 0.5), (0.75, 0.65, 0.7), (1.0, 1.0, 1.0)),
         'green': ((0.0, 1.0, 1.0), (0.25, 0.3, 0.2), (0.5, 0.5, 0.5), (0.75, 0.65, 0.7), (1.0, 1.0, 1.0)),
         'blue':  ((0.0, 1.0, 1.0), (0.25, 0.7, 0.7), (0.5, 0.5, 0.5), (0.75, 0.20, 0.3), (1.0, 1.0, 1.0)) }

    return mpl.colors.LinearSegmentedColormap('periodic',cdict,256)


from .cm_dat.KPiRdYl1 import *
from .cm_dat.KPiRdYl2 import *
from .cm_dat.KOr import *
from .cm_dat.KBu import *

from .cm_dat.BuPiYl import *

from .cm_dat.PiBuRdYl import *
from .cm_dat.PiGn2 import *

from .cm_dat.RdBu2 import *
from .cm_dat.RdYlCyBu import *
from .cm_dat.RdYlCyBu2 import *
from .cm_dat.RdYlCyBu3 import *

from .cm_dat.WBuKRd import *
from .cm_dat.WBuPi import *
from .cm_dat.WBuPi2 import *
from .cm_dat.WCyBuPi import *
from .cm_dat.WCyGnRd import *
from .cm_dat.WCyGnYlRd import *
from .cm_dat.WYlBuKPi import *
from .cm_dat.WYlGn import *
from .cm_dat.WYlRd import *

from .cm_dat.parula import *     # New matlab standard color map


dynlib_cms = [
    # Linear BW
    greys, greys2, 

    # Linear color, starting white
    defabs, WBuKRd, 
    q,
    WBuPi, WBuPi2, WYlBuKPi, WCyBuPi, 
    WCyGnRd, WCyGnYlRd, 
    WYlGn, WYlRd, 

    # Linear color, starting black (dark)
    KBu, KOr,
    KPiRdYl1, KPiRdYl2,
    BuPiYl, PiBuRdYl, 
    parula,

    # Diverging BW
    div_bw, div_greys,
    
    # Diverging color
    RdBu2, RdYlCyBu, RdYlCyBu2, RdYlCyBu3, 
    PiGn2, 

    # Misc
    periodic, 
]

if __name__ == '__main__':
    ncms = len(dynlib_cms)

    import matplotlib.pyplot as plt
    import numpy as np
    
    x, y = np.meshgrid(np.arange(-1,1.01,0.1), np.arange(-1,1.01,1))

    plt.figure(figsize=(6, ncms/2.0))
    for idx, cm in enumerate(dynlib_cms):
        if not isinstance(cm, mpl.colors.LinearSegmentedColormap):
            cm = cm()
        ax = plt.subplot(ncms, 1, idx+1)
        ax.set_title(cm.name)
        ax.contourf(x[0,:], y[:,0], x, 256, cmap=cm)
        ax.set_xticks([])
        ax.set_yticks([])
    
    plt.tight_layout()
    plt.savefig('dynlib_cms.pdf')

# that's it
