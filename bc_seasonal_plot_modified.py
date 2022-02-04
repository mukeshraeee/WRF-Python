"""code for plotting BC seasonal data"""

from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np
import wrf
from wrf import (to_np, getvar,interplevel, smooth2d, get_basemap, latlon_coords, ALL_TIMES)
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
from matplotlib.legend_handler import HandlerLine2D

#======= Load WRF output file ==========================

win = [Dataset("wrfout_d01_2017-01-01_00_00_00"),
      Dataset("wrfout_d01_2017-02-01_00_00_00"),
      Dataset("wrfout_d01_2017-12-01_00_00_00")]
spr = [Dataset("wrfout_d01_2017-03-01_00_00_00.nc"),
       Dataset("wrfout_d01_2017-04-01_00_00_00.nc"),
       Dataset("wrfout_d01_2017-05-01_00_00_00.nc")]
mon = [Dataset("wrfout_d01_2017-06-01_00_00_00"),
       Dataset("wrfout_d01_2017-07-01_00_00_00"),
       Dataset("wrfout_d01_2017-08-01_00_00_00")]
aut = [Dataset("wrfout_d01_2017-09-01_00_00_00"),
       Dataset("wrfout_d01_2017-10-01_00_00_00"),
       Dataset("wrfout_d01_2017-11-01_00_00_00")]
#================ Load DATA =========================================
nc1 = Dataset("wrfout_d01_2017-01-01_00_00_00")
nc2 = Dataset("wrfout_d01_2017-02-01_00_00_00")
nc3 = Dataset("wrfout_d01_2017-03-01_00_00_00.nc")
nc4 = Dataset("wrfout_d01_2017-04-01_00_00_00.nc")
nc5 = Dataset("wrfout_d01_2017-05-01_00_00_00.nc")
nc6 = Dataset("wrfout_d01_2017-06-01_00_00_00")
nc7 = Dataset("wrfout_d01_2017-07-01_00_00_00")
nc8 = Dataset("wrfout_d01_2017-08-01_00_00_00")
nc9 = Dataset("wrfout_d01_2017-09-01_00_00_00")
nc10 = Dataset("wrfout_d01_2017-10-01_00_00_00")
nc11 = Dataset("wrfout_d01_2017-11-01_00_00_00")
nc12 = Dataset("wrfout_d01_2017-12-01_00_00_00")

#====== Get data =================================
win_bc1 = getvar(win, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
win_bc2 = getvar(win, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
spr_bc1 = getvar(spr, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
spr_bc2 = getvar(spr, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
mon_bc1 = getvar(mon, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
mon_bc2 = getvar(mon, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
aut_bc1 = getvar(aut, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
aut_bc2 = getvar(aut, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]

win_bc3 = win_bc1+win_bc2
spr_bc3 = spr_bc1+spr_bc2
mon_bc3 = mon_bc1+mon_bc2
aut_bc3 = aut_bc1+aut_bc2

win_bc = np.average(win_bc3,axis=0)
spr_bc = np.average(spr_bc3,axis=0)
sum_bc = np.average(mon_bc3,axis=0)
aut_bc = np.average(aut_bc3,axis=0)
#======= Get pressure data ============
p  = wrf.getvar(win, 'pressure', timeidx=-1)
#========== Get WRF lat lon  ===================
lat, lon = latlon_coords(p)
#=== get ALT variables for converting ug/kg to ug/m3 ===
alt1  = getvar(nc1, "ALT", timeidx=-1)[0]
alt2  = getvar(nc2, "ALT", timeidx=-1)[0]
alt3  = getvar(nc3, "ALT", timeidx=-1)[0]
alt4  = getvar(nc4, "ALT", timeidx=-1)[0]
alt5  = getvar(nc5, "ALT", timeidx=-1)[0]
alt6  = getvar(nc6, "ALT", timeidx=-1)[0]
alt8  = getvar(nc8, "ALT", timeidx=-1)[0]
alt9  = getvar(nc9, "ALT", timeidx=-1)[0]
alt10  = getvar(nc10, "ALT", timeidx=-1)[0]
alt11  = getvar(nc11, "ALT", timeidx=-1)[0]
alt12  = getvar(nc12, "ALT", timeidx=-1)[0]
#====== Seasonal alt ===============================
winter_alt = (alt1+alt2+alt12)/3
spring_alt = (alt3+alt4+alt5)/3
summer_alt = (alt6+alt8)/2
autumn_alt = (alt9+alt10+alt11)/3
#====== ug/kg to ug/m3 ===========================
winter_bc = win_bc/winter_alt
spring_bc = spr_bc/spring_alt
summer_bc = sum_bc/summer_alt
autumn_bc = aut_bc/autumn_alt
#====== Smoothing ============
sigma_y = 2
sigma_x = 2
sigma = [sigma_y, sigma_x]
bc_win = sp.ndimage.filters.gaussian_filter(winter_bc, sigma, mode='constant')
bc_spr = sp.ndimage.filters.gaussian_filter(spring_bc, sigma, mode='constant')
bc_sum = sp.ndimage.filters.gaussian_filter(summer_bc, sigma, mode='constant')
bc_aut = sp.ndimage.filters.gaussian_filter(autumn_bc, sigma, mode='constant')
#=============== Get surface wind component ====================
u1 = getvar(win, "ua", units="m s-1",timeidx=ALL_TIMES) #[:,0,:,:]
v1 = getvar(win, "va", units="m s-1",timeidx=ALL_TIMES) #[:,0,:,:]
u2 = getvar(spr, "ua", units="m s-1",timeidx=ALL_TIMES) #[:,0,:,:]
v2 = getvar(spr, "va", units="m s-1",timeidx=ALL_TIMES) #[:,0,:,:]
u3 = getvar(mon, "ua", units="m s-1",timeidx=ALL_TIMES) #[:,0,:,:]
v3 = getvar(mon, "va", units="m s-1",timeidx=ALL_TIMES) #[:,0,:,:]
u4 = getvar(aut, "ua", units="m s-1",timeidx=ALL_TIMES) #[:,0,:,:]
v4 = getvar(aut, "va", units="m s-1",timeidx=ALL_TIMES) #[:,0,:,:]


#=== Now interpolate at 500 hpa  ==============================
winu1 = interplevel(u1, p, 500)
winv1 = interplevel(v1, p, 500)
spru1 = interplevel(u2, p, 500)
sprv1 = interplevel(v2, p, 500)
sumu1 = interplevel(u3, p, 500)
sumv1 = interplevel(v3, p, 500)
autu1 = interplevel(u4, p, 500)
autv1 = interplevel(v4, p, 500)


#===== Take time average ============================
winu = np.average(winu1,axis=0)
winv = np.average(winv1,axis=0)
spru = np.average(spru1,axis=0)
sprv = np.average(sprv1,axis=0)
sumu = np.average(sumu1,axis=0)
sumv = np.average(sumv1,axis=0)
autu = np.average(autu1,axis=0)
autv = np.average(autv1,axis=0)

#===== calculate seasonal wind speed =============
wi_ws1 = np.sqrt(winu**2+winv**2)
sp_ws1 = np.sqrt(spru**2+sprv**2)
su_ws1 = np.sqrt(sumu**2+sumv**2)
au_ws1 = np.sqrt(autu**2+autv**2)

#==== Normalize the data in uniform arrows ============
win_u_wrf = winu/wi_ws1
win_v_wrf = winv/wi_ws1
spr_u_wrf = spru/sp_ws1
spr_v_wrf = sprv/sp_ws1
sum_u_wrf = sumu/su_ws1
sum_v_wrf = sumv/su_ws1
aut_u_wrf = autu/au_ws1
aut_v_wrf = autv/au_ws1


"""
#===============================================================#
#++++++++++++++++++++++ Now Plot +++++++++++++++++++++++++++++++#
===============================================================#"""
fig = plt.figure(figsize=(4.5,3),dpi=300)

#=============== Winter ======================================================
ax = plt.subplot(2, 2,1)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
p1 = m.pcolormesh(x,y,bc_win,cmap=cmaps.WhiteBlueGreenYellowRed) 
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2], labels=[1,0,0,0], fontsize=7, color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2],  fontsize=7, color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.title('Winter',fontsize=6)
plt.clim(0,8)
#==== Plot guiver ===========================================
q1 = plt.quiver(x[::5,::5], y[::5,::5], win_u_wrf[::5,::5], win_v_wrf[::5,::5],
              pivot='middle',scale_units='inches',headwidth=5,headlength=8)

#=============== spring ======================================================
ax = plt.subplot(2, 2,2)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
p2 = m.pcolormesh(x,y,bc_spr,cmap=cmaps.WhiteBlueGreenYellowRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2], fontsize=7, color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2],  fontsize=7, color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.title('Spring',fontsize=6)
plt.clim(0,8)
#==== Plot guiver ===========================================
q2 = plt.quiver(x[::5,::5], y[::5,::5], spr_u_wrf[::5,::5], win_v_wrf[::5,::5],
              pivot='middle',scale_units='inches',headwidth=5,headlength=8)


#=============== summer ======================================================
ax = plt.subplot(2, 2,3)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
p2 = m.pcolormesh(x,y,bc_sum,cmap=cmaps.WhiteBlueGreenYellowRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2], labels=[1,0,0,0],fontsize=7, color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2],  labels=[0,0,0,1],fontsize=7, color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.title('Summer',fontsize=6)
plt.clim(0,8)
#==== Plot guiver ===========================================
q3 = plt.quiver(x[::5,::5], y[::5,::5], sum_u_wrf[::5,::5], win_v_wrf[::5,::5],
              pivot='middle',scale_units='inches',headwidth=5,headlength=8)


#=============== spring ======================================================
ax = plt.subplot(2, 2,4)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
ax1 = m.pcolormesh(x,y,bc_aut,cmap=cmaps.WhiteBlueGreenYellowRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2], fontsize=7, color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2],  labels=[0,0,0,1],fontsize=7, color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.title('Autumn',fontsize=6)
plt.clim(0,8)
#==== Plot guiver ===========================================
q4 = plt.quiver(x[::5,::5], y[::5,::5], aut_u_wrf[::5,::5], win_v_wrf[::5,::5],
              pivot='middle',scale_units='inches',headwidth=5,headlength=8)


#====== Colorbar ===============
cax = fig.add_axes([0.89,0.05,0.02,0.90])   # left, bottom, width, and height
cbar= fig.colorbar(ax1, cax=cax,ticks=[0,1,2,3,4,5,6,7,8,9,10]) #,7,8,9,10])
cbar.set_label(r'BC  [$μg  m^{-3}$]',fontsize=6,labelpad=0)
cbar.ax.tick_params(labelsize=7)

fig.subplots_adjust(top=0.935,
                        bottom=0.085,
                        left=0.045,
                        right=0.895,
                        hspace=0.18,
                        wspace=0.02)
#plt.savefig("/mnt/e/LaTex-PPT/bc_modified_500hpa.png",dpi=300)
plt.show()
