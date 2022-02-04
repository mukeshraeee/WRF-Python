
""" This code generates transport flux for different wind components i.e. c,w and wspd
	=>> written by Mukesh Rai,, 2021/11/24 """

#=== load library
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np
import wrf
from wrf import (getvar, to_np, get_cartopy,get_basemap, latlon_coords, vertcross,ALL_TIMES,cartopy_xlim, cartopy_ylim, interpline, CoordPair)
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
import cmaps
from pylab import *
import scipy as sp
import scipy.ndimage
import matplotlib.gridspec as gridspec
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams

"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"

#===== load files =========================================
data = [Dataset("wrfout_d01_2017-01-01_00_00_00"),
       Dataset("wrfout_d01_2017-02-01_00_00_00"),
       Dataset("wrfout_d01_2017-12-01_00_00_00"),
       Dataset("wrfout_d01_2017-03-01_00_00_00.nc"),
       Dataset("wrfout_d01_2017-04-01_00_00_00.nc"),
       Dataset("wrfout_d01_2017-05-01_00_00_00.nc"),
       Dataset("wrfout_d01_2017-06-01_00_00_00"),
       Dataset("wrfout_d01_2017-07-01_00_00_00"),
       Dataset("wrfout_d01_2017-08-01_00_00_00"),
       Dataset("wrfout_d01_2017-09-01_00_00_00"),
       Dataset("wrfout_d01_2017-10-01_00_00_00"),
       Dataset("wrfout_d01_2017-11-01_00_00_00")]

#====== get variables =========================================
b1  = wrf.getvar(data, 'BC1', timeidx=ALL_TIMES)[:,:,:]    #bc1
b2  = wrf.getvar(data, 'BC2', timeidx=ALL_TIMES)[:,:,:]    #bc2
v   = wrf.getvar(data, 'va',  timeidx=150)[:,:,:]  # meridional wind
w   = wrf.getvar(data, 'wa',  timeidx=200)[:,:,:]  # vertical wind
ws  = wrf.getvar(data, 'wspd',timeidx=150)[:,:,:]  # wind speed

#===== time dimesion average out ====================
bc1 = np.average(b1,axis=0)
bc2 = np.average(b2,axis=0)
bc  = bc1+bc2 

"""=== BC Transport Flux ==>> (tF)======"""
tf1 = bc * v  # TF toward TP 
tf2 = bc * w  # TF  vertical
tf3 = bc * ws # TG due to wind speed

"""======== Now interpolation """
z1 = 5**(tf1/5.)
z2 = 5**(tf2/5.)
z3 = 5**(tf3/5.)

"""========  Define the cross section start and end points========="""

cross_start = CoordPair(lat=22, lon=70)
cross_end   = CoordPair(lat=36, lon=119)

"""===== Get terrain height ============================"""
ter      = getvar(data, "ter", units="km")
"""============= Get model height for mass grid =========================================="""
ht  = getvar(data, "height_agl", units="km")[:,:,:]
"""====  Get the terrain heights along the cross section line ============================="""
ter_line = interpline(ter, wrfin=data, start_point=cross_start,end_point=cross_end)
"""Compute the vertical cross-section interpolation.
   Also, include the lat/lon points along the cross-section in the metadata by setting latlon to True"""
z1_cross = vertcross(z1, ht, wrfin=data,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z2_cross = vertcross(z2, ht, wrfin=data,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z3_cross = vertcross(z3, ht, wrfin=data,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

#== for wind component =====================================================================================
v1_cross  = vertcross(v, ht, wrfin=data,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
w1_cross  = vertcross(w, ht, wrfin=data,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
ws1_cross = vertcross(ws, ht, wrfin=data,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

"""===== Convert back to dBz after interpolation =============================="""
z1_cross  = 5 * np.log10(z1_cross)
z2_cross  = 5 * np.log10(z2_cross)
z3_cross  = 5 * np.log10(z3_cross)

#=== wind compoent ====================
v1_cross   = 5 * np.log10(v1_cross)
w1_cross   = 5 * np.log10(w1_cross)
ws1_cross  = 5 * np.log10(ws1_cross)


"""Make a copy of the z cross data. Let's use regular numpy arrays for this"""
xs1 = np.arange(0, z1_cross.shape[-1], 1)
ys1 = to_np(z1_cross.coords["vertical"])
xs2 = np.arange(0, z2_cross.shape[-1], 1)
ys2 = to_np(z2_cross.coords["vertical"])
xs3 = np.arange(0, z3_cross.shape[-1], 1)
ys3 = to_np(z3_cross.coords["vertical"])

#== For wind =================================
xs4 = np.arange(0, v1_cross.shape[-1], 1)
ys4 = to_np(v1_cross.coords["vertical"])
xs5 = np.arange(0, w1_cross.shape[-1], 1)
ys5 = to_np(w1_cross.coords["vertical"])
xs6 = np.arange(0, ws1_cross.shape[-1], 1)
ys6 = to_np(ws1_cross.coords["vertical"])





"""Make a copy of the z cross data. Let's use regular numpy arrays for this"""
z1_cross_filled = np.ma.copy(to_np(z1_cross))
z2_cross_filled = np.ma.copy(to_np(z2_cross))
z3_cross_filled = np.ma.copy(to_np(z3_cross))

#====wind ====================================
v1_cross_filled  = np.ma.copy(to_np(v1_cross))
w1_cross_filled  = np.ma.copy(to_np(w1_cross))
ws1_cross_filled = np.ma.copy(to_np(ws1_cross))


"""=========================================================================
   For each cross section column, find the first index with non-missing
    values and copy these to the missing elements below
   ========================================================================="""
for i in range(z1_cross_filled.shape[-1]):
    column_vals = z1_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z1_cross_filled[0:first_idx, i] = z1_cross_filled[first_idx, i]

for i in range(z2_cross_filled.shape[-1]):
    column_vals = z2_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z2_cross_filled[0:first_idx, i] = z2_cross_filled[first_idx, i]

for i in range(z3_cross_filled.shape[-1]):
    column_vals = z3_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z3_cross_filled[0:first_idx, i] = z3_cross_filled[first_idx, i]



"""======Plot data =================================================="""

fig, (ax1,ax2,ax3) = plt.subplots(1, 3,figsize=(9,3),dpi=200)

"""===== BC =========================================================="""
p1 = ax1.contourf(xs1,ys1,to_np(z1_cross_filled),50,cmap=cmaps.MPL_RdBu_r)      #v
p2 = ax2.contourf(xs2,ys2,to_np(z2_cross_filled),50,cmap=cmaps.BlueWhiteOrangeRed)    #w
p3 = ax3.contourf(xs3,ys3,to_np(z3_cross_filled),50,cmap=cmaps.WhiteBlueGreenYellowRed)  #ws

#=== wind contour plot ========
#===== v component ==================
p4 = ax1.contour(xs4,ys4,to_np(v1_cross_filled), cmap='hsv',levels=15,linestyles="dashed",linewidths=0.3)
#p4 = ax1.contour(xs4,ys4,to_np(v1_cross_filled), cmap='tab20c',levels=5,linestyles="dashed",linewidths=0.6)
ax1.clabel(p4, inline=True, fontsize=5,inline_spacing=5,fmt='%1.0f')
#===== w component ==================
#p5 = ax2.contour(xs5,ys5,to_np(w1_cross_filled), cmap=cmaps.cmocean_curl,levels=5,linestyles="dashed",linewidths=0.6)
p5 = ax2.contour(xs5,ys5,to_np(w1_cross_filled), cmap='hsv',levels=6, linestyles="dashed",linewidths=0.3)
ax1.clabel(p5, inline=True, fontsize=5,inline_spacing=5,fmt='%1.0f')
#===== v component ==================
#p6 = ax3.contour(xs6,ys6,to_np(ws1_cross_filled), cmap=cmaps.cmocean_curl,levels=5,linestyles="dashed",linewidths=0.5)
p6 = ax3.contour(xs6,ys6,to_np(ws1_cross_filled), cmap='gnuplot',levels=13,linestyles="dashed",linewidths=0.3)
ax1.clabel(p6, inline=True, fontsize=5,inline_spacing=5,fmt='%1.0f')




"""===== Add colorbar=================================================="""
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
#b1 = plt.colorbar(p1, cax=cax,ticks=[-9,-7.5,-6,-4.5,-3,-1.5,0,1.5,3,4.5,6,7.5])
b1 = plt.colorbar(p1, cax=cax,ticks=[-7,-6.5,-6,-5.5,-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4])

b1.ax.tick_params(labelsize=5)
#================================================
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
b2 = plt.colorbar(p2, cax=cax,ticks=[-0.055,-0.05,-0.045,-0.04,-0.035,-0.03,-0.025,-0.02,-0.015,-0.010,-0.005,0,0.005,0.010,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.060])
b2.ax.tick_params(labelsize=5)
#===============================================
divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="5%", pad=0.05)
#b3 = plt.colorbar(p3, cax=cax,ticks=[0,2,4,6,8,10,12,14,16,18])
b3 = plt.colorbar(p3, cax=cax,ticks=[0,1,2,3,4,5,6,7,8,9,10,11])
b3.ax.tick_params(labelsize=5)

"""=========  Fill in the mountain area ========================================="""
ax1.fill_between(xs1, 0, to_np(ter_line),facecolor="#4D4D4D",zorder=5)
ax2.fill_between(xs2, 0, to_np(ter_line),facecolor="#4D4D4D",zorder=5)
ax3.fill_between(xs3, 0, to_np(ter_line),facecolor="#4D4D4D",zorder=5)


"""=====  Set the x-ticks to use latitude and longitude labels==="""
coord_pairs = to_np(z1_cross.coords["xy_loc"])
x_ticks     = np.arange(coord_pairs.shape[0])
x_labels    = [pair.latlon_str(fmt="{:.1f}, {:.1f}")
              for pair in to_np(coord_pairs)]
"""====== Set the desired number of x ticks below===== """
num_ticks = 5
thin = int((len(x_ticks) / num_ticks) + .8)

"""===== Set xtick label ============================================="""
ax1.set_xticklabels(x_labels[::thin],rotation=25, fontsize=7)
ax2.set_xticklabels(x_labels[::thin],rotation=25, fontsize=7)
ax3.set_xticklabels(x_labels[::thin],rotation=25, fontsize=7)
"""====  Set ytick label ============================================="""
ax1.yaxis.set_tick_params(labelsize=7)
ax2.yaxis.set_tick_params(labelsize=7)
ax3.yaxis.set_tick_params(labelsize=7)

""" Add  text """
ax1.text(-35,10,'Elevation [Km]', va='center', rotation='vertical',fontsize=7)
ax1.text(70,18,"[a]",fontsize=7)
ax2.text(70,18,"[b]",fontsize=7)
ax3.text(70,18,"[c]",fontsize=7)

ax1.text(145,-0.8,"$μg m ^{-2} s ^{-1}$",fontsize=5)
ax2.text(145,-0.8,"$μg m ^{-2} s ^{-1}$",fontsize=5)
ax3.text(145,-0.8,"$μg m ^{-2} s ^{-1}$",fontsize=5)



""" Adjust"""
fig.subplots_adjust(top=0.988,
                        bottom=0.124,
                        left=0.06,
                        right=0.96,
                        hspace=0.2,
                        wspace=0.4)

#plt.savefig("/mnt/g/2nd_Paper/tf_2000dpi.png",dpi=2000)
plt.show()
