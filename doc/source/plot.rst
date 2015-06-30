Using dynlib's plotting facilities
==================================

Basic concepts
--------------

To be written.


.. _plot configuration:

Plot configuration
------------------

Dynlib plots can be customised by a large number of keyword arguments. Here is a comprehensive list of accepted arguments:

=============================== ======= ======= ======================= ======================= =======================
Name                            Line    Fill    Type                    Default                 Description
=============================== ======= ======= ======================= ======================= =======================
cb_orientation  			✓ 	string                  ``'vertical'`` 	        Orientation of the color bar. 
cmap		         	✓ 	      	mpl-colormap 	 	                        Color map for the contours. 
cmap		        		✓ 	mpl-colormap 	        ``cm.gist_ncar``        Color map. 
coastcolor	        	✓ 	✓ 	mpl-color 	        black 	                Color of the coastlines. 
colors		        	✓ 	      	list of mpl-color       black 	                Presribe color of the contours, instead of using color map. 
colors		        	      	✓ 	list of mpl-color 	 	                Prescribe colors directly instead of using color map. 
contour_labels	        	✓ 	      	bool 		        ``False``	        Label the contours within the plot. 
contour_labels_fontsize         ✓ 	  	int/string 	        ``12`` 	                Font size in point, or a string like {\tt 'smaller'}. 
contour_labels_inline 	        ✓ 	    	bool 	                ``True``                Remove the contour under the label? 
contour_labels_inline_spacing   ✓ 	        int 	                ``2`` 	                Space in pixels around the label where the contour is removed as well. 
contour_labels_format 	        ✓ 	   	string/list of string   ``'%1.1f'``             How to format the numbers, or list of strings used as labels. 
disable_cb		        ✓ 	✓ 	bool		        ``False`` 	        Do not add a colorbar to the plot. 
extend			          	✓ 	string		        ``'both'``              Extend the color bar on the upper and lower boundaries? Possible values are: {\tt 'neither'}, {\tt 'upper'}, {\tt 'lower'} and {\tt 'both'}.  
gridcolor		        ✓ 	✓ 	mpl-color 	        black 	                Color lat/lon grid lines. 
hook 			        ✓ 	✓ 	callable 			                Function to apply to the data before plotting. 
linestyles		        ✓ 	      	string 	 		                        Line styles for contours. 
linewidths		        ✓ 	      	int/list of int	 		                Line widths for contours. 
m		 	        ✓ 	✓ 	Basemap		        worldmap 	        Map projection. 
maskcolor		        ✓ 	✓ 	mpl-color 	        light grey 	        Color of the parts below orography. 
mark			        ✓ 	✓ 	tupleof list 	 	                        Tuple containing list of x-coordinates and y-coordinates to be marked on the map with circles. 
oroalpha		        ✓ 	✓ 	float 	                ``0.4``    	        Transparency of the orography isolines, which $0$ being entirely translucent. 
orocolor		        ✓ 	✓ 	mpl-color 	        black 	                Color of the orography isolines. 
oroscale		        ✓ 	✓ 	int/list 	        Δ=1000m                 Anything matplotlib accepts as a scale. 
overlays		        ✓ 	✓ 	list of overlay 	 	                List of overlays (being fronts or countours) to plot on top. 
plev			        ✓ 	✓ 	string/int 	                                Vertical level of the plot. Used to auto-mask intersections with orography using the ERA-Interim average height of the level. 
scale			        ✓ 	✓ 	list/int/string         ``'auto'``              If {\tt 'auto'}, use the configurable autoscaling (see {\tt scale*}-properties), otherwise anything {\tt matplotlib} accepts. 
scale_exceed_percentiles        ✓       ✓ 	tuple of float          ``(0.01, 0.99)``        Percentiles giving the minimal displayed data coverage of the colorbar. For the default values, the lowermost and uppermost inverval cuts of maximally $2\%$ of the data values. 
scale_intervals 	        ✓ 	✓ 	list of int             ``[1,2,3,5,10]``        Which intervals are considered to be "round"? By default it is this list for any order of magnitude.
scale_intervals_periodic        ✓       ✓ 	bool 	                ``True``                Which intervals are considered to be "round"? If set to {\tt False} the above list applies only for the given order of magnitude. 
scale_target_steps 	        ✓ 	✓ 	int 		        ``7`` 		        If {\tt scale}$=${\tt 'auto'}, how many intervals are desired? The actual intervals might be a few more/less, to allow for a "round" interval value. 
scale_symmetric_zero 	        ✓ 	✓ 	bool 		        ``False`` 	        If {\tt scale}$=${\tt 'auto'}, force the scale to be symmetric around zero. 
save			        ✓ 	✓ 	string 		 		                If not {\tt None}, save plot as the given file name with full path. 
show			        ✓ 	✓ 	bool 		        ``True`` 	        Show the plot in a window. 
ticks			        ✓ 	✓ 	list 		  		                Where to put the ticks on the colorbar. 
ticklabels		        ✓ 	✓ 	list 		  		                How to label the ticks. 
title			        ✓ 	✓ 	string 		  		                Title for the plot. 
Zdata			        ✓ 	✓ 	numpy.ndarray 	 	                        2D geopotential height field, used for masking intersections with orography. 
=============================== ======= ======= ======================= ======================= =======================


