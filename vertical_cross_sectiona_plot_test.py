#===================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors
from cartopy import crs
from cartopy.feature import NaturalEarthFeature, COLORS
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy,get_basemap, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair,ALL_TIMES)
import cmaps
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import StrMethodFormatter
import matplotlib.gridspec as gridspec
import matplotlib as mpl
#===========================================
base =     [Dataset("wrfout_d01_2017-01-01_00_00_00"),
           Dataset("wrfout_d01_2017-02-01_00_00_00"),
           Dataset("wrfout_d01_2017-03-01_00_00_00.nc"),
           Dataset("wrfout_d01_2017-04-01_00_00_00.nc"),
           Dataset("wrfout_d01_2017-05-01_00_00_00.nc"),
           Dataset("wrfout_d01_2017-06-01_00_00_00"),
           Dataset("wrfout_d01_2017-07-01_00_00_00"),
           Dataset("wrfout_d01_2017-08-01_00_00_00"),
           Dataset("wrfout_d01_2017-09-01_00_00_00"),
           Dataset("wrfout_d01_2017-10-01_00_00_00"),
           Dataset("wrfout_d01_2017-11-01_00_00_00"),
           Dataset("wrfout_d01_2017-12-01_00_00_00")]
sen =     [Dataset("wrfout_d01_2017-01-01_nobc"),
           Dataset("wrfout_d01_2017-02-01_nobc"),
           Dataset("wrfout_d01_2017-03-01_nobc"),
           Dataset("wrfout_d01_2017-04-01_nobc"),
           Dataset("wrfout_d01_2017-05-01_nobc"),
           Dataset("wrfout_d01_2017-06-01_nobc"),
           Dataset("wrfout_d01_2017-07-01_nobc"),
           Dataset("wrfout_d01_2017-08-01_nobc"),
           Dataset("wrfout_d01_2017-09-01_nobc"),
           Dataset("wrfout_d01_2017-10-01_nobc"),
           Dataset("wrfout_d01_2017-11-01_nobc"),
           Dataset("wrfout_d01_2017-12-01_nobc")]

#========================================================
nc = Dataset("wrfout_d01_2017-01-01_00_00_00")
#====== Getting PBLH ===========
pblh = getvar(nc,"PBLH", timeidx=-1)
pblh1 = np.mean(pblh,axis=1)

#======= Get data data ==========
t2  = getvar(nc,"tc",timeidx=-1) #Temp data
ter = getvar(nc, "ter", units="km")  # Model Terrain Height
z1  = getvar(nc,"z",units="km")  #Model Height for Mass Grid
dbz = getvar(nc, "dbz")
#========= Get values upto 10KM  from the surface ===============
t   =  t2[0:20,:,:]
z   =  z1[0:20,:,:]
Z   =  10**(t/10.)

#==== Set the start point and end point for the cross section
start_point = CoordPair(lat=30, lon=70)
end_point   = CoordPair(lat=30, lon=118)
# Compute the vertical cross-section interpolation
#z_cross  = vertcross(Z, z,    wrfin=nc, start_point=start_point,end_point=end_point, latlon=True, meta=True)
tc_cross = vertcross(t, z,   wrfin=nc, start_point=start_point,end_point=end_point, latlon=True, meta=True)
ter_line = interpline(ter,    wrfin=nc, start_point=start_point,end_point=end_point,latlon=True, meta=True)

tc_cross   = 10.0 * np.log10(tc_cross)

tc_cross.attrs.update(tc_cross.attrs)
tc_cross.attrs["description"] = "temperature"
tc_cross.attrs["units"] = "degC"

# Make a copy of the z cross data. Let's use regular numpy arrays for this.
tc_cross_filled = np.ma.copy(to_np(tc_cross))

#for i in range(tc_cross_filled.shape[-1]):
#    column_vals = tc_cross_filled[:,i]
#    first_idx = int(np.transpose((column_vals > 0).nonzero())[0])
#    tc_cross_filled[0:first_idx, i] = tc_cross_filled[first_idx, i]

#===Figure
fig, ax = plt.subplots(figsize=(18,9))
#===== Setting contour layers ======
#tc_levels = np.arange(-30, 20, 5)
xs = np.arange(0, tc_cross.shape[-1], 1)
ys = to_np(tc_cross.coords["vertical"])
#=========  Fill in the mountain area =========================================
ht_fill = ax.fill_between(xs, 0, to_np(ter_line),facecolor="grey")
#dbz_cross = 10.0 * np.log10(z_cross)
#pbl_cross = 10.0 * np.log10(pblh1)
# Get the latitude and longitude point
lats, lons = latlon_coords(pblh)
#=============  Create the figure=============
#fig = plt.figure(figsize=(12,6))
#ax = plt.axes()
# Make the contour plot
tc_contours = ax.contourf(to_np(tc_cross), cmap=cmaps.WhiteBlueGreenYellowRed) #WhiteBlueGreenYellowRed) #BlueWhiteOrangeRed)
#pbl = ax.plot(to_np(pblh1),'r')
# Set the x-ticks to use latitude and longitude labels.
coord_pairs = to_np(tc_cross.coords["xy_loc"])
x_ticks     = np.arange(coord_pairs.shape[0])
x_labels    = [pair.latlon_str(fmt="{:.2f}, {:.2f}")
              for pair in to_np(coord_pairs)]

#=== Seeting Number of Xtick Label
num_ticks = 7
thin = int((len(x_ticks) / num_ticks) + .7)
ax.set_xticks(x_ticks[::thin])
ax.set_xticklabels(x_labels[::thin], rotation=45, fontsize=8)

# Set the y-ticks to be height.
vert_vals = to_np(tc_cross.coords["vertical"])
v_ticks   = np.arange(vert_vals.shape[0])

ax.set_yticks(v_ticks[::10])
ax.set_yticklabels(vert_vals[::10], fontsize=8)

# Set the x-axis and  y-axis labels
ax.set_xlabel("Latitude, Longitude", fontsize=12)
ax.set_ylabel("Height (Km)", fontsize=12)

#plt.clim(-2,2)
plt.colorbar(tc_contours, ax=ax)
plt.title("Vertical Cross Section of Temperature (degC)")


plt.show()
