"""This code generates radiative forcing in toa, atm and bot"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors
from cartopy import crs
from cartopy.feature import NaturalEarthFeature, COLORS
from netCDF4 import Dataset
from wrf import (getvar,smooth2d, to_np, get_cartopy, latlon_coords, vertcross,ALL_TIMES,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair)
import cmaps
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import StrMethodFormatter
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from matplotlib.patches import Patch
from pylab import *
from matplotlib.pyplot import figure
import matplotlib as mpl
import scipy as sp
import scipy.ndimage

#========= Load WRF data ================================
jan   = Dataset("wrfout_d01_2017-01-01_00_00_00") # Control

swupt  = getvar(jan,"SWUPT",timeidx=1)
swdnt  = getvar(jan,"SWDNT",timeidx=1)
swuptc  = getvar(jan,"SWUPTC",timeidx=1)
swdntc  = getvar(jan,"SWDNTC",timeidx=1)

lwupt  = getvar(jan,"LWUPT",timeidx=1)
lwdnt  = getvar(jan,"LWDNT",timeidx=1)
lwuptc  = getvar(jan,"LWUPTC",timeidx=1)
lwdntc  = getvar(jan,"LWDNTC",timeidx=1)

swdown = getvar(jan,"SWDOWN",timeidx=1)
glw    =getvar(jan,"GLW",timeidx=1)

net = swdown+glw 



fig, ax =plt.subplots()

#=========== WRF-AOD-Winter ================================
lats, lons = latlon_coords(swdnt)
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax = m.pcolormesh(x,y,net,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,0.5)
plt.title('Jan')
plt.ylabel('SW-RF', labelpad=40,fontsize=12)

plt.show()

