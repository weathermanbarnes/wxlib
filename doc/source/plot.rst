Using dynlib's plotting facilities
==================================

Motivation and scope
--------------------

With matplotlib and Basemap provide a comprehensive library for general plotting and plotting 
maps with python. Why not use those tools directly, but rather build another library on top
of those? 

The reason is that there are a few sharp edges when plotting maps with matplotlib/Basemap,
that will make it harder than necessary to produce a nice-looking map plot. While matplotlib
and Basemap provide almost every possible feature that you can think of, some of them are 
pretty well hidden, but potentially important for plotting meteorological data. And last but
not least, some plot defaults of matplotlib and Basemap are suboptimal for our applications.

All these issues are relatively easy to circumvent or solve (once you know the solution!). 
These issues nevertheless result in quite a lot of python code for a standard map plot, that
one should not be required to copy&paste every time one needs to plot a map.

The aim of the plotting functions in dynlib is to provide easy-to-use plot commands, that
incoprorate all the workarounds, and that come with meaningful default settings. The resulting
plot should look fine by default, and should easily be customised to produce a publication-
quality plot.

Dynlib currently provides two fundamental plot types

:func:`dynlib.figures.map`
   Maps with arbitrary projection
:func:`dynlib.figures.section`
   Vertical cross sections, including the necessary interpolations

and a new axes type for plotting wind roses, :func:`dynlib.windrose.WindroseAxes`, `written by
Lionel Roubeyrie <http://youarealegend.blogspot.no/2008/09/windrose.html>`_, but adapted to 
work better with dynlib.

Dynlib might be extended with further plot types in the future, if either a correspoding plot
type does not exist in matplotlib by default, or if it typically requires a considerable 
amount of code to make it applicable to meteorological data. A good candidate for future
inclusion would for example be a tephigram.


Using the plot facilities
-------------------------

The basic plot commands :func:`dynlib.figures.map` and :func:`dynlib.figures.section` each 
produce a filled-contours plot of the given data. Additional data can be displayed using 
*overlays*. Both the basic plot commands and the overlays take a extensive list of optional
arguments that can be used to customise the plot. These plot configuration arguments are 
listed in :ref:`plot_configuration`.

The following overlays are available for maps:

:func:`dynlib.figures.map_overlay_barbs`
   Overlay wind barbs.
:func:`dynlib.figures.map_overlay_contours`
   Overlay contours.
:func:`dynlib.figures.map_overlay_dilatation`
   Overlay axes of dilatation, double arrows for orientation and strength.
:func:`dynlib.figures.map_overlay_dots`
   Overlay dots.
:func:`dynlib.figures.map_overlay_fronts`
   Overlay front lines.
:func:`dynlib.figures.map_overlay_latlonbox`
   Overlay a latitude-longitude box.
:func:`dynlib.figures.map_overlay_lines`
   Overlay generic lines, for example jet axes.
:func:`dynlib.figures.map_overlay_quiver`
   Overlay a vector field by arrows.
:func:`dynlib.figures.map_overlay_shading`
   Overlay further filled contours. Currently, no color bar for the overlay can be added.

And the following overlays are available for sections:

:func:`dynlib.figures.section_overlay_contours`
   Overlay contours.


.. _plot configuration:

Plot configuration
------------------

Dynlib plots can be customised by a large number of keyword arguments. Here is a comprehensive list of accepted arguments:

=============================== ======= ======= ======================= ======================= =======================
Name                            Line    Fill    Type                    Default                 Description
=============================== ======= ======= ======================= ======================= =======================
alpha		         	✓       ✓       float                   1.0			Opacity of the layer, 1.0 is fully opaque, 0.0 fully transparent.
cb_disable		        ✓  	✓ 	bool		        ``False``		Do not add a colorbar to the plot. 
cb_expand_fig_fraction		✓ 	✓	float			0.1			Fraction of the figure height or width to extend the figure by, to make room for the colorbar.
cb_label                        ✓ 	✓	list of string                                  Labels for the color bar.
cb_orientation  		✓ 	✓ 	string			``'vertical'``		Orientation of the color bar. 
cb_tickspacing  		✓ 	✓ 	string			``'proportional'``	Spacing policy for ticks at the colorbar. Either ``'proportional'`` or ``'uniform'``.
cb_shrink			✓ 	✓	float			0.8			Fraction of the figure height or width covered by the colorbar.
cmap		         	✓ 	      	mpl-colormap 	 	                        Color map for the contours. 
cmap		        		✓ 	mpl-colormap 	        ``cm.gist_ncar``        Color map. 
coastcolor	        	✓ 	✓ 	mpl-color 	        black 	                Color of the coastlines. 
colors		        	✓ 	      	list of mpl-color       black 	                Presribe color of the contours, instead of using color map. 
colors		        	      	✓ 	list of mpl-color 	 	                Prescribe colors directly instead of using color map. 
contour_labels	        	✓ 	      	bool 		        ``False``	        Label the contours within the plot. 
contour_labels_fontsize         ✓ 	  	int/string 	        ``12`` 	                Font size in point, or a string like ``'smaller'``. 
contour_labels_inline 	        ✓ 	    	bool 	                ``True``                Remove the contour under the label? 
contour_labels_inline_spacing   ✓ 	        int 	                ``2`` 	                Space in pixels around the label where the contour is removed as well. 
contour_labels_format 	        ✓ 	   	string/list of string   ``'%1.1f'``             How to format the numbers, or list of strings used as labels. 
extend				✓	 	string		        ``'neither'``		Extend the color bar on the upper and lower boundaries? Possible values are: ``'neither'``, ``'upper'``, ``'lower'`` and ``'both'``.  
extend			          	✓ 	string		        ``'both'``		Extend the color bar on the upper and lower boundaries? Possible values are: ``'neither'``, ``'upper'``, ``'lower'`` and ``'both'``.
fig_size			✓	✓	float/tuple of float	32.0			If tuple: figure width and height in inches; if float: figure size in in². In that case the appropriate aspect ratio will be determined by the map projection.
fig_dpi				✓	✓	int			150			Pixels per inch for rasterising the figure.
gridcolor		        ✓ 	✓ 	mpl-color 	        black 	                Color lat/lon grid lines. 
hatches				 	✓ 	tuple or list of string	()			List of hatch patterns to fill the contours with.
hook 			        ✓ 	✓ 	callable 			                Function to apply to the data before plotting. 
linestyles		        ✓ 	      	string 	 		                        Line styles for contours. 
linewidths		        ✓ 	      	int/list of int	 		                Line widths for contours. 
m		 	        ✓ 	✓ 	Basemap		        worldmap 	        Map projection. 
maskcolor		        ✓ 	✓ 	mpl-color 	        light grey 	        Color of the parts below orography. 
mark			        ✓ 	✓ 	tuple of list 	 	                        Tuple containing list of x-coordinates and y-coordinates to be marked on the map with circles. 
name				✓	✓	string			''			Name of the plot, used for automatically determining file name and title.
name_prefix			✓	✓	string			''			Prefix for the plot file name.
norm                                            mpl.colors.Normalize    LinearNorm              Defines the mapping of values onto colors. By default a linear mapping is used, such that double the distance in value corresponds to double the distance in the color map. Alternatives are ``mpl.colors.BoundaryNorm`` to achieve an equal color spacing between colour bar ticks, and ``mpl.colors.LogNorm`` for a logarithmic mapping.
oroalpha		        ✓ 	✓ 	float 	                ``0.4``    	        Transparency of the orography isolines, which $0$ being entirely translucent. 
orocolor		        ✓ 	✓ 	mpl-color 	        black 	                Color of the orography isolines. 
oroscale		        ✓ 	✓ 	int/list 	        Δ=1000m                 Anything matplotlib accepts as a scale. 
overlays		        ✓ 	✓ 	list of overlay 	 	                List of overlays (being fronts or countours) to plot on top. 
save			        ✓ 	✓ 	string 		 		                If not ``None``, save plot as the given file name with full path. 
scale			        ✓ 	✓ 	list/int/string         ``'auto'``              If ``'auto'``, use the configurable autoscaling (see scale*-properties), otherwise anything matplotlib accepts. 
scale_exceed_percentiles        ✓       ✓ 	tuple of float          ``(0.01, 0.99)``        Percentiles giving the minimal displayed data coverage of the colorbar. For the default values, the lowermost and uppermost inverval cuts of maximally $2\%$ of the data values. 
scale_intervals 	        ✓ 	✓ 	list of int             ``[1,2,3,5,10]``        Which intervals are considered to be "round"? By default it is this list for any order of magnitude.
scale_intervals_periodic        ✓       ✓ 	bool 	                ``True``                Which intervals are considered to be "round"? If set to ``False`` the above list applies only for the given order of magnitude. 
scale_target_steps 	        ✓ 	✓ 	int 		        ``7`` 		        If scale = ``'auto'``, how many intervals are desired? The actual intervals might be a few more/less, to allow for a "round" interval value. 
scale_symmetric_zero 	        ✓ 	✓ 	bool 		        ``False`` 	        If scale = ``'auto'``, force the scale to be symmetric around zero. 
show			        ✓ 	✓ 	bool 		        ``True`` 	        Show the plot in a window. 
ticks			        ✓ 	✓ 	list 		  		                Where to put the ticks on the colorbar. 
ticklabels		        ✓ 	✓ 	list 		  		                How to label the ticks. 
tile                             	✓       bool                    ``False``               Use pcolormesh instead of contourf to plot the shading, having only one solid color per grid cell.
title			        ✓ 	✓ 	string 		  	``'auto'``		Title for the plot. If ``'auto'`` the plot title is determined automatically from q, plev and name.
Zdata			        ✓ 	✓ 	numpy.ndarray 	 	                        2D geopotential height field, used for masking intersections with orography. 
=============================== ======= ======= ======================= ======================= =======================


In addition the plot functions take two special key word arguments. They can be used to make
dynlib use the default values for a certain variable and vertical level.

=============================== ======= ======= ======================= =======================
Name                            Line    Fill    Type                    Description
=============================== ======= ======= ======================= =======================
plev			        ✓ 	✓ 	string/int 	        Vertical level of the plot. Used also to auto-mask intersections with orography using the ERA-Interim average height of the level. 
q				✓ 	✓ 	string			Variable identifier.
=============================== ======= ======= ======================= =======================


