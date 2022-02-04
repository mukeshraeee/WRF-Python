""" This code produce RF and temperatur change due to BC """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors
from cartopy import crs
from cartopy.feature import NaturalEarthFeature, COLORS
from netCDF4 import Dataset
from wrf import (getvar,smooth2d, to_np, get_cartopy, latlon_coords, vertcross,ALL_TIMES,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair)
import wrf
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

#===== Control =======================================
winc = [Dataset("wrfout_d01_2017-01-01_00_00_00"),
      Dataset("wrfout_d01_2017-02-01_00_00_00"),
      Dataset("wrfout_d01_2017-12-01_00_00_00")]
sprc = [Dataset("wrfout_d01_2017-03-01_00_00_00.nc"),
       Dataset("wrfout_d01_2017-04-01_00_00_00.nc"),
       Dataset("wrfout_d01_2017-05-01_00_00_00.nc")]
monc = [Dataset("wrfout_d01_2017-06-01_00_00_00"),
       Dataset("wrfout_d01_2017-07-01_00_00_00"),
       Dataset("wrfout_d01_2017-08-01_00_00_00")]
autc = [Dataset("wrfout_d01_2017-09-01_00_00_00"),
       Dataset("wrfout_d01_2017-10-01_00_00_00"),
       Dataset("wrfout_d01_2017-11-01_00_00_00")]

#====== Sensitivity =====================================
wins = [Dataset("wrfout_d01_2017-01-01_nobc"),
        Dataset("wrfout_d01_2017-02-01_nobc"),
        Dataset("wrfout_d01_2017-12-01_nobc")]
sprs = [Dataset("wrfout_d01_2017-03-01_nobc"),
        Dataset("wrfout_d01_2017-04-01_nobc"),
        Dataset("wrfout_d01_2017-05-01_nobc")]
mons = [Dataset("wrfout_d01_2017-06-01_nobc"),
        Dataset("wrfout_d01_2017-07-01_nobc"),
        Dataset("wrfout_d01_2017-08-01_nobc")]
auts = [Dataset("wrfout_d01_2017-09-01_nobc"),
        Dataset("wrfout_d01_2017-10-01_nobc"),
        Dataset("wrfout_d01_2017-11-01_nobc")]

#======== Winter ===========================================
w_swdwn_c1  = getvar(winc,   "SWDOWN", timeidx=ALL_TIMES)
w_swdwn_s1  = getvar(wins,   "SWDOWN", timeidx=ALL_TIMES)
w_glw_c1    = getvar(winc,   "GLW",    timeidx=ALL_TIMES)
w_glw_s1    = getvar(wins,   "GLW",    timeidx=ALL_TIMES)
w_alb_c1    = getvar(winc,   "ALBEDO", timeidx=ALL_TIMES)
w_alb_s1    = getvar(wins,   "ALBEDO", timeidx=ALL_TIMES)
w_emis_c1   = getvar(winc,   "EMISS",  timeidx=ALL_TIMES)
w_emis_s1   = getvar(wins,   "EMISS",  timeidx=ALL_TIMES)
w_tsk_c1    = getvar(winc,   "TSK",    timeidx=ALL_TIMES)
w_tsk_s1    = getvar(wins,   "TSK",    timeidx=ALL_TIMES)
w_grdflx_c1 = getvar(winc,   "GRDFLX", timeidx=ALL_TIMES)
w_grdflx_s1 = getvar(wins,   "GRDFLX", timeidx=ALL_TIMES)
w_lh_c1     = getvar(winc,   "LH",     timeidx=ALL_TIMES)
w_lh_s1     = getvar(wins,   "LH",     timeidx=ALL_TIMES)
w_sh_c1     = getvar(winc,   "HFX",    timeidx=ALL_TIMES)
w_sh_s1     = getvar(wins,   "HFX",    timeidx=ALL_TIMES)

#===== Time average =======================================
w_swdwn_c  = np.average(w_swdwn_c1, axis=0)
w_swdwn_s  = np.average(w_swdwn_s1, axis=0)
w_glw_c    = np.average(w_glw_c1 ,  axis=0)
w_glw_s    = np.average(w_glw_s1 ,  axis=0)
w_alb_c    = np.average(w_alb_c1,   axis=0)
w_alb_s    = np.average(w_alb_s1,   axis=0)
w_emis_c   = np.average(w_emis_c1,  axis=0)
w_emis_s   = np.average(w_emis_s1,  axis=0)
w_tsk_c    = np.average(w_tsk_c1,   axis=0)
w_tsk_s    = np.average(w_tsk_s1,   axis=0)
w_grdflx_c = np.average(w_grdflx_c1,axis=0)
w_grdflx_s = np.average(w_grdflx_s1,axis=0)
w_lh_c     = np.average(w_lh_c1,    axis=0)
w_lh_s     = np.average(w_lh_s1,    axis=0)
w_sh_c     = np.average(w_sh_c1,    axis=0)
w_sh_s     = np.average(w_sh_s1,    axis=0)


#======== Spring ===========================================
s_swdwn_c1  = getvar(sprc,   "SWDOWN", timeidx=ALL_TIMES)
s_swdwn_s1  = getvar(sprs,   "SWDOWN", timeidx=ALL_TIMES)
s_glw_c1    = getvar(sprc,   "GLW",    timeidx=ALL_TIMES)
s_glw_s1    = getvar(sprs,   "GLW",    timeidx=ALL_TIMES)
s_alb_c1    = getvar(sprc,   "ALBEDO", timeidx=ALL_TIMES)
s_alb_s1    = getvar(sprs,   "ALBEDO", timeidx=ALL_TIMES)
s_emis_c1   = getvar(sprc,   "EMISS",  timeidx=ALL_TIMES)
s_emis_s1   = getvar(sprs,   "EMISS",  timeidx=ALL_TIMES)
s_tsk_c1    = getvar(sprc,   "TSK",    timeidx=ALL_TIMES)
s_tsk_s1    = getvar(sprs,   "TSK",    timeidx=ALL_TIMES)
s_grdflx_c1 = getvar(sprc,   "GRDFLX", timeidx=ALL_TIMES)
s_grdflx_s1 = getvar(sprs,   "GRDFLX", timeidx=ALL_TIMES)
s_lh_c1     = getvar(sprc,   "LH",     timeidx=ALL_TIMES)
s_lh_s1     = getvar(sprs,   "LH",     timeidx=ALL_TIMES)
s_sh_c1     = getvar(sprc,   "HFX",    timeidx=ALL_TIMES)
s_sh_s1     = getvar(sprs,   "HFX",    timeidx=ALL_TIMES)

#===== Time average =======================================
s_swdwn_c  = np.average(s_swdwn_c1, axis=0)
s_swdwn_s  = np.average(s_swdwn_s1, axis=0)
s_glw_c    = np.average(s_glw_c1 ,  axis=0)
s_glw_s    = np.average(s_glw_s1 ,  axis=0)
s_alb_c    = np.average(s_alb_c1,   axis=0)
s_alb_s    = np.average(s_alb_s1,   axis=0)
s_emis_c   = np.average(s_emis_c1,  axis=0)
s_emis_s   = np.average(s_emis_s1,  axis=0)
s_tsk_c    = np.average(s_tsk_c1,   axis=0)
s_tsk_s    = np.average(s_tsk_s1,   axis=0)
s_grdflx_c = np.average(s_grdflx_c1,axis=0)
s_grdflx_s = np.average(s_grdflx_s1,axis=0)
s_lh_c     = np.average(s_lh_c1,    axis=0)
s_lh_s     = np.average(s_lh_s1,    axis=0)
s_sh_c     = np.average(s_sh_c1,    axis=0)
s_sh_s     = np.average(s_sh_s1,    axis=0)

#======== Summer ===========================================
su_swdwn_c1  = getvar(monc,   "SWDOWN", timeidx=ALL_TIMES)
su_swdwn_s1  = getvar(mons,   "SWDOWN", timeidx=ALL_TIMES)
su_glw_c1    = getvar(monc,   "GLW",    timeidx=ALL_TIMES)
su_glw_s1    = getvar(mons,   "GLW",    timeidx=ALL_TIMES)
su_alb_c1    = getvar(monc,   "ALBEDO", timeidx=ALL_TIMES)
su_alb_s1    = getvar(mons,   "ALBEDO", timeidx=ALL_TIMES)
su_emis_c1   = getvar(monc,   "EMISS",  timeidx=ALL_TIMES)
su_emis_s1   = getvar(mons,   "EMISS",  timeidx=ALL_TIMES)
su_tsk_c1    = getvar(monc,   "TSK",    timeidx=ALL_TIMES)
su_tsk_s1    = getvar(mons,   "TSK",    timeidx=ALL_TIMES)
su_grdflx_c1 = getvar(monc,   "GRDFLX", timeidx=ALL_TIMES)
su_grdflx_s1 = getvar(mons,   "GRDFLX", timeidx=ALL_TIMES)
su_lh_c1     = getvar(monc,   "LH",     timeidx=ALL_TIMES)
su_lh_s1     = getvar(mons,   "LH",     timeidx=ALL_TIMES)
su_sh_c1     = getvar(monc,   "HFX",    timeidx=ALL_TIMES)
su_sh_s1     = getvar(mons,   "HFX",    timeidx=ALL_TIMES)
#===== Time average =======================================
su_swdwn_c  = np.average(su_swdwn_c1, axis=0)
su_swdwn_s  = np.average(su_swdwn_s1, axis=0)
su_glw_c    = np.average(su_glw_c1 ,  axis=0)
su_glw_s    = np.average(su_glw_s1 ,  axis=0)
su_alb_c    = np.average(su_alb_c1,   axis=0)
su_alb_s    = np.average(su_alb_s1,   axis=0)
su_emis_c   = np.average(su_emis_c1,  axis=0)
su_emis_s   = np.average(su_emis_s1,  axis=0)
su_tsk_c    = np.average(su_tsk_c1,   axis=0)
su_tsk_s    = np.average(su_tsk_s1,   axis=0)
su_grdflx_c = np.average(su_grdflx_c1,axis=0)
su_grdflx_s = np.average(su_grdflx_s1,axis=0)
su_lh_c     = np.average(su_lh_c1,    axis=0)
su_lh_s     = np.average(su_lh_s1,    axis=0)
su_sh_c     = np.average(su_sh_c1,    axis=0)
su_sh_s     = np.average(su_sh_s1,    axis=0)

#======== Autumn ===========================================
a_swdwn_c1  = getvar(autc,   "SWDOWN", timeidx=ALL_TIMES)
a_swdwn_s1  = getvar(auts,   "SWDOWN", timeidx=ALL_TIMES)
a_glw_c1    = getvar(autc,   "GLW",    timeidx=ALL_TIMES)
a_glw_s1    = getvar(auts,   "GLW",    timeidx=ALL_TIMES)
a_alb_c1    = getvar(autc,   "ALBEDO", timeidx=ALL_TIMES)
a_alb_s1    = getvar(auts,   "ALBEDO", timeidx=ALL_TIMES)
a_emis_c1   = getvar(autc,   "EMISS",  timeidx=ALL_TIMES)
a_emis_s1   = getvar(auts,   "EMISS",  timeidx=ALL_TIMES)
a_tsk_c1    = getvar(autc,   "TSK",    timeidx=ALL_TIMES)
a_tsk_s1    = getvar(auts,   "TSK",    timeidx=ALL_TIMES)
a_grdflx_c1 = getvar(autc,   "GRDFLX", timeidx=ALL_TIMES)
a_grdflx_s1 = getvar(auts,   "GRDFLX", timeidx=ALL_TIMES)
a_lh_c1     = getvar(autc,   "LH",     timeidx=ALL_TIMES)
a_lh_s1     = getvar(auts,   "LH",     timeidx=ALL_TIMES)
a_sh_c1     = getvar(autc,   "HFX",    timeidx=ALL_TIMES)
a_sh_s1     = getvar(auts,   "HFX",    timeidx=ALL_TIMES)


#===== Time average =======================================
a_swdwn_c  = np.average(a_swdwn_c1, axis=0)
a_swdwn_s  = np.average(a_swdwn_s1, axis=0)
a_glw_c    = np.average(a_glw_c1 ,  axis=0)
a_glw_s    = np.average(a_glw_s1 ,  axis=0)
a_alb_c    = np.average(a_alb_c1,   axis=0)
a_alb_s    = np.average(a_alb_s1,   axis=0)
a_emis_c   = np.average(a_emis_c1,  axis=0)
a_emis_s   = np.average(a_emis_s1,  axis=0)
a_tsk_c    = np.average(a_tsk_c1,   axis=0)
a_tsk_s    = np.average(a_tsk_s1,   axis=0)
a_grdflx_c = np.average(a_grdflx_c1,axis=0)
a_grdflx_s = np.average(a_grdflx_s1,axis=0)
a_lh_c     = np.average(a_lh_c1,    axis=0)
a_lh_s     = np.average(a_lh_s1,    axis=0)
a_sh_c     = np.average(a_sh_c1,    axis=0)
a_sh_s     = np.average(a_sh_s1,    axis=0)

#====== Calculate RF =============================================================================

win_rf1 = (w_swdwn_c-w_swdwn_c*w_alb_c+w_emis_c*w_glw_c-w_lh_c-w_grdflx_c-w_sh_c)-\
          (w_swdwn_s-w_swdwn_s*w_alb_s+w_emis_s*w_glw_s-w_lh_s-w_grdflx_s-w_sh_s)
spr_rf1 = (s_swdwn_c-s_swdwn_c*s_alb_c+s_emis_c*s_glw_c-s_lh_c-s_grdflx_c-s_sh_c)-\
          (s_swdwn_s-s_swdwn_s*s_alb_s+s_emis_s*s_glw_s-s_lh_s-s_grdflx_s-s_sh_s)
mon_rf1 = (su_swdwn_c-su_swdwn_c*su_alb_c+su_emis_c*su_glw_c-su_lh_c-su_grdflx_c-su_sh_c)-\
          (su_swdwn_s-su_swdwn_s*su_alb_s+su_emis_s*su_glw_s-su_lh_s-su_grdflx_s-su_sh_s)
aut_rf1 = (a_swdwn_c-a_swdwn_c*a_alb_c+a_emis_c*a_glw_c-a_lh_c-a_grdflx_c-a_sh_c)-\
          (a_swdwn_s-a_swdwn_s*a_alb_s+a_emis_s*a_glw_s-a_lh_s-a_grdflx_s-a_sh_s)

#=======NOW TEMPERATURE ======================================
""" get data from control """
winc_k1 = getvar(winc, 'T2',timeidx=ALL_TIMES)-273
sprc_k1 = getvar(sprc, 'T2',timeidx=ALL_TIMES)-273
monc_k1 = getvar(monc, 'T2',timeidx=ALL_TIMES)-273
autc_k1 = getvar(autc, 'T2',timeidx=ALL_TIMES)-273

""" get data from sensitivity """
wins_k1 = getvar(wins, 'T2',timeidx=ALL_TIMES)-273
sprs_k1 = getvar(sprs, 'T2',timeidx=ALL_TIMES)-273
mons_k1 = getvar(mons, 'T2',timeidx=ALL_TIMES)-273
auts_k1 = getvar(auts, 'T2',timeidx=ALL_TIMES)-273


#===== Take time average =======================
winc_k = np.average(winc_k1,axis=0)
sprc_k = np.average(sprc_k1,axis=0)
monc_k = np.average(monc_k1,axis=0)
autc_k = np.average(autc_k1,axis=0)

wins_k = np.average(wins_k1,axis=0)
sprs_k = np.average(sprs_k1,axis=0)
mons_k = np.average(mons_k1,axis=0)
auts_k = np.average(auts_k1,axis=0)

#===== Temperature due to BC ==================
win_t1 = (winc_k-wins_k)
spr_t1 = (sprc_k-sprs_k)
mon_t1 = (monc_k-mons_k)
aut_t1 = (autc_k-auts_k)



#======== Smoothing RF using Gaussion filter =============
sigma_y   =  2
sigma_x   =  2
sigma     =  [sigma_y, sigma_x]

win_rf   = sp.ndimage.filters.gaussian_filter(win_rf1, sigma, mode='constant')
spr_rf   = sp.ndimage.filters.gaussian_filter(spr_rf1, sigma, mode='constant')
mon_rf   = sp.ndimage.filters.gaussian_filter(mon_rf1, sigma, mode='constant')
aut_rf   = sp.ndimage.filters.gaussian_filter(aut_rf1, sigma, mode='constant')
#======== Smoothing Temp using Gaussion filter =============
win_t    = sp.ndimage.filters.gaussian_filter(win_t1, sigma, mode='constant')
spr_t    = sp.ndimage.filters.gaussian_filter(spr_t1, sigma, mode='constant')
mon_t    = sp.ndimage.filters.gaussian_filter(mon_t1, sigma, mode='constant')
aut_t    = sp.ndimage.filters.gaussian_filter(aut_t1, sigma, mode='constant')


#============= Subplot - Winter  =================================================================
fig, axs = plt.subplots(4, 2,figsize=(6,8),dpi=150)
gridspec.GridSpec(4,2)
plt.subplot2grid((4,2), (0,0))
#============================================================
lats, lons = latlon_coords(w_swdwn_c1)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,win_rf,cmap=cmaps.temp_diff_18lev) #cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_R>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1],  color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.ylabel('Winter', labelpad=5,fontsize=7)
#plt.clim(-8,8)
#============ Spring ====================
plt.subplot2grid((4,2), (1,0))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,spr_rf,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5,dashes=[4, 1],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1],  color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.ylabel('Spring', labelpad=5,fontsize=7)
#plt.clim(-8,8)

#============ Summer ====================
plt.subplot2grid((4,2), (2,0))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,mon_rf,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1],  color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.ylabel('Summer', labelpad=5,fontsize=7)
#plt.clim(-8,8)
#============ Autumn ====================
plt.subplot2grid((4,2), (3,0))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,aut_rf,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=7, color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.ylabel('Autumn', labelpad=5,fontsize=7)
#plt.clim(-8,8)


#==============================================================================================================
#=== Now plot for Temperature ========
plt.subplot2grid((4,2), (0,1))

#========= Winter===================================================
lats, lons = latlon_coords(w_swdwn_c1)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

x, y = m(to_np(lons), to_np(lats))
s1 = m.pcolormesh(x,y,win_t,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1],labels=[0,1,0,0], fontsize=7,  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1],  color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
#plt.clim(-0.2,0.2)
#plt.title('Temperature [°C]',fontsize=8)

#========= spring===================================================
plt.subplot2grid((4,2), (1,1))
lats, lons = latlon_coords(w_swdwn_c1)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

x, y = m(to_np(lons), to_np(lats))
s1 = m.pcolormesh(x,y,spr_t,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1],labels=[0,1,0,0], fontsize=7,  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1],  color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
#plt.clim(-0.2,0.2)

#========== Summer =========================================
plt.subplot2grid((4,2), (2,1))
ats, lons = latlon_coords(w_swdwn_c1)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

x, y = m(to_np(lons), to_np(lats))
s1 = m.pcolormesh(x,y,mon_t,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1],labels=[0,1,0,0], fontsize=7,  color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1],  color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
#plt.clim(-0.2,0.2)

#========== Autumn  =========================================
plt.subplot2grid((4,2), (3,1))
ats, lons = latlon_coords(w_swdwn_c1)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

x, y = m(to_np(lons), to_np(lats))
s1 = m.pcolormesh(x,y,aut_t,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], labels=[0,1,0,0], fontsize=7, color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=7, color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
#plt.clim(-0.2,0.2)

#=== Setting Colorbar for RF=========
#cax1 = fig.add_axes([0.05,0.06,0.4,0.02])   # left, bottom, width, and height
#cb1 = fig.colorbar(s,ax=axs[:,:], cax=cax1,ticks=[-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8], orientation='horizontal')
#cb1.set_label(r'Radiative Forcing [$Wm^{-2}$]',fontsize=7,labelpad=5)
#cb1.ax.tick_params(labelsize=5)

#=== Setting Colorbar for Temperarture=========
#cax2 = fig.add_axes([0.5,0.06,0.42,0.02])   # left, bottom, width, and height
#cb2 = fig.colorbar(s1, ax=axs[:,:], cax=cax2,ticks=[-0.20,-0.15,-0.10,-0.05,0,0.05,0.1,0.15,0.20], orientation='horizontal')
#cb2.set_label('Temperature [°C]',fontsize=7,labelpad=5)
#cb2.ax.tick_params(labelsize=5)

fig.subplots_adjust(top=0.995,bottom=0.1,left=0.04,right=0.93,hspace=0.03,wspace=0.01)

#plt.savefig("/mnt/h/Paper_work/RF_Temp_300dpi.jpeg",dpi=300)
plt.show()









