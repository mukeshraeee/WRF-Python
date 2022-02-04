import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors
from cartopy import crs
from cartopy.feature import NaturalEarthFeature, COLORS
from netCDF4 import Dataset
import wrf
from wrf import (getvar, to_np, get_cartopy,get_basemap, latlon_coords, vertcross,ALL_TIMES,cartopy_xlim, cartopy_ylim, interpline, CoordPair)
import cmaps
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import StrMethodFormatter
import matplotlib.gridspec as gridspec
import cartopy
import cartopy.crs as ccrs
import cartopy.mpl.geoaxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams
"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"



cont = [Dataset("wrfout_d01_2017-01-01_00_00_00"),
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


#===== add data from experiment => no BC  =============
sens = [Dataset("wrfout_d01_2017-01-01_nobc"),
        Dataset("wrfout_d01_2017-02-01_nobc"),
        Dataset("wrfout_d01_2017-12-01_nobc"),
        Dataset("wrfout_d01_2017-03-01_nobc"),
        Dataset("wrfout_d01_2017-04-01_nobc"),
        Dataset("wrfout_d01_2017-05-01_nobc"),
        Dataset("wrfout_d01_2017-06-01_nobc"),
        Dataset("wrfout_d01_2017-07-01_nobc"),
        Dataset("wrfout_d01_2017-08-01_nobc"),
        Dataset("wrfout_d01_2017-09-01_nobc"),
        Dataset("wrfout_d01_2017-10-01_nobc"),
        Dataset("wrfout_d01_2017-11-01_nobc")]

wf = Dataset("wrfout_d01_2017-01-01_00_00_00")
#============= Get model height for mass grid ==========================================
ht  = getvar(cont, "height_agl", units="km",timeidx=-1,method='cat')[0:17,:,:]

#============== Get PBLH from control experiment=====================================================================
pblh_cont1  = getvar(cont, "PBLH",timeidx=ALL_TIMES,method='cat')
pblh_sens1  = getvar(sens, "PBLH",timeidx=ALL_TIMES,method='cat')

#======== tIME aVERAGE===========
pblh_cont2 = np.mean(pblh_cont1,axis=0) # all aerosol
pblh_sens2 = np.mean(pblh_sens1,axis=0) # without BC

#======= Convert meter to KM  ==============
#====== Control============
pbl_cont = (pblh_cont2)*0.001
pbl_sens = (pblh_sens2)*0.001

pblh_bc1 = pbl_cont - pbl_sens
#=========== Get Terrain Height ===========================
ter1 = getvar(cont, "ter", units="km") #control
ter2 = getvar(sens, "ter", units="km") #sensitivity

#===== PBLH height for control ========================================
pblh_con = pbl_cont+ter1
pblh_sen = pbl_sens+ter2

pblh_bc = pblh_bc1+ter2
#========= Get ALT =========================================
alt  = getvar(cont, "ALT",timeidx=-1,method='cat')[0]

#======= Get dbz==========================================
dbz      = getvar(cont, "dbz", timeidx=-1)
#======= Seasonal ter height =============================
con_ter      = getvar(cont, "ter", units="km")
sen_ter      = getvar(sens, "ter", units="km")

#========= get BC ============================================
cont_bc1 = getvar(cont, "BC1", timeidx=ALL_TIMES,method='cat')[:,0:17,:,:] # all time,surface,all lat and lon
cont_bc2 = getvar(cont, "BC2", timeidx=ALL_TIMES,method='cat')[:,0:17,:,:]


#== TIME average ================
con_bc1 = np.mean(cont_bc1,axis=0)
con_bc2 = np.mean(cont_bc2,axis=0)
#==========  ==============
cont_bc3 = con_bc1+con_bc2 
bc = cont_bc3/alt
#======= Get ua and wa = CONTROL ===============================================================
cont_u = getvar(cont, "ua", timeidx=10,method='cat')[0:17,:,:] # all time,surface,all lat and lon
cont_w = getvar(cont, "wa", timeidx=10,method='cat')[0:17,:,:]



"""#=======Now Plot ========================="""

fig, axs = plt.subplots(1, 1,figsize=(8,3),dpi=200)

z1 = 5**(bc/10.) # Use linear Z for interpolation
w1 = 25**(cont_w/25.)
u1 = 25**(cont_u/25.)

#========  Define the cross section start and end points=========
cross_start = CoordPair(lat=22, lon=70)
cross_end   = CoordPair(lat=35, lon=119)
"""Compute the vertical cross-section interpolation.
   Also, include the lat/lon points along the cross-section in the metadata by setting latlon to True"""

#=======  ==============================================================================================
z1_cross = vertcross(z1, ht, wrfin=cont,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
w1_cross = vertcross(w1, ht, wrfin=cont,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
u1_cross = vertcross(u1, ht, wrfin=cont,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

#
z1_cross  = 10.0 * np.log10(z1_cross)
w1_cross  = 10.0 * np.log10(w1_cross)
u1_cross  = 10.0 * np.log10(u1_cross)

"""Make a copy of the z cross data. Let's use regular numpy arrays for this"""
#====== ============================================
z1_cross_filled = np.ma.copy(to_np(z1_cross))
w1_cross_filled = np.ma.copy(to_np(w1_cross))
u1_cross_filled = np.ma.copy(to_np(u1_cross))

#======= SA ==================================================================
for i in range(z1_cross_filled.shape[-1]):
    column_vals = z1_cross_filled[:,i]
   #Let's find the lowest index that isn't filled. The nonzero function
   #finds all unmasked values greater than 0. Since 0 is a valid value
   #for dBZ, let's change that threshold to be -200 dBZ instead.
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z1_cross_filled[0:first_idx, i] = z1_cross_filled[first_idx, i]


#===============  Get the lat/lon points =============
lats, lons = latlon_coords(dbz)
cart_proj = get_cartopy(dbz)

bc_levels = np.arange(0,2,0.01)

"""#====  Get the terrain heights along the cross section line"""
pblh_line_cont  = interpline(pblh_con, wrfin=cont, start_point=cross_start,end_point=cross_end)
pblh_line_sens  = interpline(pblh_sen, wrfin=sens, start_point=cross_start,end_point=cross_end)
pblh_line_bc    = interpline(pblh_bc, wrfin=sens, start_point=cross_start,end_point=cross_end)



xs1 = np.arange(0, z1_cross.shape[-1], 1)
ys1 = to_np(z1_cross.coords["vertical"])

"""======Plot BC data =================================================="""
#=== SA===============
bc_contours = axs.contourf(xs1,ys1,to_np(z1_cross_filled),levels=bc_levels,cmap=cmaps.MPL_YlGnBu)

#===== SET COLOR BAR ======================
cax = fig.add_axes([0.934,0.123,0.010,0.86]) # left, bottom, width, and height
#cbar= fig.colorbar(bc_contours, ax=axs,cax=cax)
#cbar.set_label(r'BC [$μg m ^{-3}$]', rotation=-270, fontsize=7,labelpad=0.8)
#cbar.set_ticks([0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.1])
#cbar.set_ticklabels([0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2])
#cbar.ax.tick_params(labelsize=6)

#===============================================
divider = make_axes_locatable(axs)
#cax = divider.append_axes("right", size="5%", pad=0.05)
b3 = plt.colorbar(bc_contours, cax=cax,ticks=[0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,1.99])
b3.ax.tick_params(labelsize=7)


""" ====== Plot Wind data ============================================"""
#======== SA  =============================================================
axs.quiver(xs1[::5], ys1[::5],to_np(u1_cross_filled[::5, ::5]), to_np(w1_cross_filled[::5, ::5]*300))

#====  Get the terrain heights along the cross section line
ter_line_con = interpline(ter1, wrfin=cont, start_point=cross_start,end_point=cross_end)

"""#=========  Fill in the mountain area ========================================="""
#==== SA ===============================
axs.fill_between(xs1, 0, to_np(ter_line_con),facecolor="#5C5C5C")

"""===== Plot PBLH =============================================================="""
#===== SA ==================================
axs.plot(pblh_line_cont,'magenta',linewidth=0.5)
axs.plot(pblh_line_sens,'black',linewidth=0.5)
axs.plot(pblh_line_bc,'#EE7621',linewidth=0.5)


""" Set the x-ticks to use latitude and longitude labels"""
coord_pairs = to_np(z1_cross.coords["xy_loc"])
x_ticks     = np.arange(coord_pairs.shape[0])
x_labels    = [pair.latlon_str(fmt="{:.1f}, {:.1f}")
              for pair in to_np(coord_pairs)]
#=======  Set the desired number of x ticks below
num_ticks = 8
thin = int((len(x_ticks) / num_ticks) + .8)

axs.set_xticks(x_ticks[::thin])
axs.set_xticklabels(x_labels[::thin],rotation=25, fontsize=8)
axs.yaxis.set_tick_params(labelsize=8)

#===== add text==============
fig.text(0.007,  0.53, 'Elevation [Km]',                 va='center', rotation='vertical',  fontsize=10)
fig.text(0.15,   0.26, 'South Asia [+IGP]',           va='center', rotation='horizontal',fontsize=8)
fig.text(0.46,   0.26, 'Tibetan Plateau \n   Himalayas',va='center', rotation='horizontal',fontsize=8)
fig.text(0.78,   0.26, 'East Asia',                  va='center', rotation='horizontal',fontsize=8)
fig.text(0.08,   0.93, '[c]',                  va='center', rotation='horizontal',fontsize=8)

"""==== Inset axes ===>> ZOOM IN ==>> SA region=============="""
axins1 = zoomed_inset_axes(axs, 1.3, loc='upper left')
axins1.set_xlim(1,45)
axins1.set_ylim(0.02, 1.5)
mark_inset(axs, axins1, loc1=1, loc2=2, fc="none", ec="#3D3D3D",lw=0.8,zorder=3)
mark_inset(axs, axins1, loc1=3, loc2=4, fc="none", ec="#3D3D3D",lw=0.8,zorder=3)


plt.xticks(visible=False)
plt.yticks(visible=False)

"""===Plot 1st inset zoom area  ====="""
bc_contours = axins1.contourf(xs1,ys1,to_np(z1_cross_filled),levels=bc_levels,cmap=cmaps.MPL_YlGnBu)
ter_line_con = interpline(ter1, wrfin=cont, start_point=cross_start,end_point=cross_end)
axins1.fill_between(xs1, 0, to_np(ter_line_con),facecolor="#5C5C5C")
axins1.plot(pblh_line_cont,'magenta',linewidth=0.5)
axins1.plot(pblh_line_sens,'black',linewidth=0.5)
axins1.plot(pblh_line_bc,'#EE7621',linewidth=0.5)


"""==== Inset axes ===>> ZOOM IN ==>> EC region=============="""
axins2 = zoomed_inset_axes(axs, 1.3, loc='upper right')
axins2.set_xlim(100,141)
axins2.set_ylim(0.02, 1.5)
mark_inset(axs, axins2, loc1=1, loc2=2, fc="none", ec="#3D3D3D",lw=0.8,zorder=3)
mark_inset(axs, axins2, loc1=3, loc2=4, fc="none", ec="#3D3D3D",lw=0.8,zorder=3)

plt.xticks(visible=False)
plt.yticks(visible=False)

"""===Plot 2nd inset zoom area  ====="""
bc_contours1 = axins2.contourf(xs1,ys1,to_np(z1_cross_filled),levels=bc_levels,cmap=cmaps.MPL_YlGnBu)
ter_line_con = interpline(ter1, wrfin=cont, start_point=cross_start,end_point=cross_end)
axins2.fill_between(xs1, 0, to_np(ter_line_con),facecolor="#5C5C5C")
axins2.plot(pblh_line_cont,'magenta',linewidth=0.5)
axins2.plot(pblh_line_sens,'black',linewidth=0.5)
axins2.plot(pblh_line_bc,'#EE7621',linewidth=0.5)


""" Insert Legend"""
l1 = mlines.Line2D([], [], color='magenta',   label='PBLH_cont.',linewidth=0.8) #control=>all aerosols
l2 = mlines.Line2D([], [], color='black',     label='PBLH_sens.',linewidth=0.8) #sensiti=>no BC
l3 = mlines.Line2D([], [], color='#EE7621',   label='PBLH_BC',linewidth=0.8) # only====>> BC

axins2.legend(handles=[l1,l2,l3],fontsize = 5)

"""#======= Add cross-sectional map ====="""
sub_axes = plt.axes([0.27, 0.80, 0.45, 0.18])  #left, bottom, width, height
m = Basemap(projection='ortho', lat_0=30, lon_0=86,resolution='h')
m.drawmapboundary(fill_color='#A6CAE0', linewidth=0) #A6CAE0
m.fillcontinents(color='grey', alpha=0.7, lake_color='grey')
m.drawcoastlines(linewidth=0.1, color="white")
m.drawcountries(linewidth=0.1, color="white")
# Cross-section to be plotted ==========================
startlat = 22; startlon = 70
arrlat   = 35; arrlon   = 119
m.drawgreatcircle(startlon,startlat,arrlon,arrlat, linewidth=0.6, color='orange')

#======= adjust 
fig.subplots_adjust(top=0.986,
                        bottom=0.124,
                        left=0.048,
                        right=0.929,
                        hspace=0.105,
                        wspace=0.05)

plt.savefig("/mnt/g/2nd_Paper/vertical_PBL2500.png",dpi=2500)
#plt.show()
