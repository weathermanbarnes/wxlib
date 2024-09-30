"""
This file is part of WaveBreaking.

WaveBreaking provides indices to detect, classify
and track Rossby Wave Breaking (RWB) in climate and weather data.
The tool was developed during my master thesis at the University of Bern.
Link to thesis: https://occrdata.unibe.ch/students/theses/msc/406.pdf

---

Utility functions for the plotting routines
"""

__author__ = "Severin Kaderli"
__license__ = "MIT"
__email__ = "severin.kaderli@unibe.ch"

# import modules
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.path as mpath


def calculate_periodic_field(da, **kwargs):
    """
    Add first longitude at the end to ensure
    that there is no gap in a periodic field
    """
    return xr.concat(
        [
            da,
            da.isel({kwargs["lon_name"]: 0}).assign_coords(
                {kwargs["lon_name"]: da[kwargs["lon_name"]].max() + 1}
            ),
        ],
        dim=kwargs["lon_name"],
    )


def get_levels(min_freq, max_freq):
    """
    Define default levels
    """

    max_level = np.round(max_freq)
    min_level = np.round(min_freq)
    level_range = max_level - min_level

    if min_level < 0:
        max_both = np.abs([min_level, max_level]).max()
        if level_range / 16 > 0.5:
            return np.round(np.linspace(-max_both, max_both, num=18), 1)[1:-1]
        else:
            return np.round(np.linspace(-max_both, max_both, num=10), 1)[1:-1]
    else:
        return np.round(np.linspace(min_level, max_level, num=8), 1)


def get_new_cmap(color_palette):
    """
    Define default cmaps for climatological plot
    """

    cmap = plt.get_cmap(color_palette)

    return colors.LinearSegmentedColormap.from_list(
        "trunc({n},{a:.2f},{b:.2f})".format(n=cmap.name, a=0.3, b=1),
        cmap(np.linspace(0.3, 1, 100)),
    )


def add_colorbar(plot, caxes, levels, label):
    """
    Define colorbar
    """

    cbar = plt.colorbar(plot, cax=caxes, drawedges=True)
    cbar.ax.set_yticklabels(levels, fontsize=12, weight="bold")
    cbar.set_label(label=label, size=12, fontweight="bold", labelpad=15)
    cbar.outline.set_color("black")
    cbar.outline.set_linewidth(2)
    cbar.dividers.set_color("black")
    cbar.dividers.set_linewidth(2)
    return cbar


def add_grid_lines(axes):
    """
    Define grid lines
    """

    gr = axes.gridlines(
        draw_labels=True, color="black", linestyle="dotted", linewidth=1.1
    )
    gr.xlabel_style = {"size": 12, "color": "black", "rotation": 0}
    gr.ylabel_style = {"size": 12, "color": "black"}
    return gr


def add_circular_boundary(axes):
    """
    Define circular boundary for NorthPolarStereo projection
    """

    return axes.add_patch(
        mpatches.Circle(
            (0.5, 0.5),
            radius=0.5,
            color="k",
            linewidth=5,
            fill=False,
            transform=axes.transAxes,
        )
    )


def add_circular_patch(axes):
    """
    Define circular patch for NorthPolarStereo projection
    """

    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    return axes.set_boundary(circle, transform=axes.transAxes)
