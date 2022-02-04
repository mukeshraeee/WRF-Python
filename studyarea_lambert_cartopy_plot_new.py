
""" This code generates study domain using Cartopy and etopo1 """

from copy import copy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry as sgeom
import cartopy
import cartopy.feature as cfe
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.offsetbox import AnchoredText
from cartopy.io import shapereader as shpreader
import cartopy.io.img_tiles as cimgt
from shapely.geometry.polygon import LinearRing
import matplotlib.patches as mpatches
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib
import pyproj
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature, OCEAN, LAKES
from mpl_toolkits.basemap import Basemap
from copy import copy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.font_manager import FontProperties
import netCDF4 as nc
import re  # regular expression
import cmaps
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms

from mpl_toolkits.basemap import Basemap
import string
import matplotlib.cm as cm

from matplotlib import rcParams

"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"


#============
def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.
    
    """
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
              'right': [(maxx, miny), (maxx, maxy)],
              'bottom': [(minx, miny), (maxx, miny)],
              'top': [(minx, maxy), (maxx, maxy)],}
    return sgeom.LineString(points[side])


def lambert_xticks(ax, ticks):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
    te = lambda xy: xy[0]
    lc = lambda t, n, b: np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T
    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels])
    

def lambert_yticks(ax, ticks):
    """Draw ricks on the left y-axis of a Lamber Conformal projection."""
    te = lambda xy: xy[1]
    lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels])

def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """Get the tick locations and labels for an axis of a Lambert Conformal projection."""
    outline_patch = sgeom.LineString(ax.outline_patch.get_path().vertices.tolist())
    axis = find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(ccrs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(ccrs.Geodetic(), xy[:, 0], xy[:, 1])
        xyt = proj_xyz[..., :2]
        ls = sgeom.LineString(xyt.tolist())
        locs = axis.intersection(ls)
        if not locs:
            tick = [None]
        else:
            tick = tick_extractor(locs.xy)
        _ticks.append(tick[0])
    # Remove ticks that aren't visible:    
    ticklabels = copy(ticks)
    while True:
        try:
            index = _ticks.index(None)
        except ValueError:
            break
        _ticks.pop(index)
        ticklabels.pop(index)
    return _ticks, ticklabels
#=========================================

et    = 'etopo1_PTPE.nc'
etopo = nc.Dataset(et, 'r', format = 'NETCDF4')

dem     = etopo.variables['z'][:]
dem_lat = etopo.variables['y'][:]
dem_lon = etopo.variables['x'][:]

etopo.close()
dem_lons, dem_lats = np.meshgrid(dem_lon, dem_lat)

for i in np.arange(dem.shape[0]):
    for j in np.arange(dem.shape[1]):
        if dem[i,j]<0:
            dem[i,j]=0


proj = ccrs.LambertConformal(central_longitude=86.3333, central_latitude=27.5,standard_parallels=(27, 30))

# Draw a set of axes with coastlines:
fig = plt.figure(figsize=(7, 7), dpi=260)
fig.tight_layout()
#fig.subplots_adjust(left=0.02, bottom=0.06, right=0.98, top=0.94, wspace=0.05)
ax = fig.add_axes([0.09, 0.09, 0.9, 0.9], projection=proj,frameon=False)
#ax = fig.add_axes([0.1,0.1,0.9,0.9], projection=proj) #,frameon=False)   #l,b,w,h
cs = ax.pcolormesh(dem_lons, dem_lats, dem, cmap=plt.cm.terrain, vmin=0, vmax=8848, alpha=0.3, transform=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.COASTLINE,linewidth=0.3)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-',linewidth=0.3)
ax.add_feature(cfeature.COASTLINE,linewidth=0.3)
img = plt.imread('marble_earth.png')
#ax.add_feature(cfeature.LAND,facecolor='lightgrey',alpha=0.5)
#ax.add_feature(cfeature.OCEAN, color="skyblue", alpha=0.4)
ax.add_feature(OCEAN, edgecolor='k') #,facecolor='deepskyblue')
ax.add_feature(LAKES, edgecolor='k') #,facecolor='deepskyblue')
states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='grey',name='admin_1_states_provinces_shp')
ax.add_feature(states, linewidth=0.4)

#ax.set_extent([30, 135, 0, 60]
#===== Adding bouding box- subregiona rectangular box ========================

##==== Add TP============
##================ 1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,  14,  15,  16,  17,  18,  19,  20,  21, 22,  23, 24,25,26,27,28,29,30,31,32,33])
lat_tp = np.array([35, 40, 40, 37, 37, 36, 36, 37, 37, 38, 38, 40, 40,  38,  38,  37,  37,  32,  32,  25,  25, 28,  28, 29,29,30,30,31,31,32,32,35,35]) 
lon_tp = np.array([70, 70, 75, 75, 80, 80, 84, 84, 86, 86, 90, 90 ,100, 100, 102, 102, 105, 105, 102, 102, 100,100, 86, 86,83,83,81,81,79,79,75,75,70])

poly_tp = np.zeros((len(lat_tp), 2), np.float64)
poly_tp[:,0] = lon_tp
poly_tp[:,1] = lat_tp

poly_tp = mpatches.Polygon(poly_tp, closed=True, ec='orange', fill=False, lw=0.8, ls = 'dashed',fc=None, transform=ccrs.PlateCarree())
ax.add_patch(poly_tp)

#======== Sub-region for East China ===============
lat_ec = np.array([25,  12,  12,  40, 40,40,  38,  38,  37,  37,  32,  32,  25,  25])
lon_ec = np.array([100, 100, 120, 120,100,100, 100, 102, 102, 105, 105, 102, 102, 100]) 
poly_ec = np.zeros((len(lat_ec), 2), np.float64)
poly_ec[:,0] = lon_ec
poly_ec[:,1] = lat_ec

poly_ec = mpatches.Polygon(poly_ec, closed=True, ec='orange', fill=False, lw=0.8, ls = 'dashed',fc=None, transform=ccrs.PlateCarree())
ax.add_patch(poly_ec)

#======== Sub-region for North China ===============
lat_nc = np.array([40, 40, 37, 37, 36, 36, 37, 37, 38, 38, 40, 40, 40, 55, 55,40])
lon_nc = np.array([70, 75, 75, 80, 80, 84, 84, 86, 86, 90, 90 ,100,120,120,70,70])
poly_nc = np.zeros((len(lat_nc), 2), np.float64)
poly_nc[:,0] = lon_nc
poly_nc[:,1] = lat_nc

poly_nc = mpatches.Polygon(poly_nc, closed=True, ec='orange', fill=False, lw=0.8, ls = 'dashed',fc=None, transform=ccrs.PlateCarree())
ax.add_patch(poly_nc)

#======== Sub-region for Central Asia ===============
lat_ca = np.array([12,55,55,12])
lon_ca = np.array([45,45,70,70])
poly_ca = np.zeros((len(lat_ca), 2), np.float64)
poly_ca[:,0] = lon_ca
poly_ca[:,1] = lat_ca

poly_ca = mpatches.Polygon(poly_ca, closed=True, ec='orange', fill=False, lw=0.8, ls = 'dashed',fc=None, transform=ccrs.PlateCarree())
ax.add_patch(poly_ca)

#======== Sub-region for South Asia ===============
lat_sa = np.array([25, 28,  28, 29,29,30,30,31,31,32,32,35,35,12,12,25])
lon_sa = np.array([100,100, 86, 86,83,83,81,81,79,79,75,75,70,70,100,100])
poly_sa = np.zeros((len(lat_sa), 2), np.float64)
poly_sa[:,0] = lon_sa
poly_sa[:,1] = lat_sa

poly_sa = mpatches.Polygon(poly_sa, closed=True, ec='orange', fill=False, lw=0.8, ls = 'dashed',fc=None, transform=ccrs.PlateCarree())
ax.add_patch(poly_sa)

#==== Add locations on the map=================
met     = pd.read_excel('location_metdata_aod_pm25.xlsx',sheet_name='met')
aer     = pd.read_excel('location_metdata_aod_pm25.xlsx',sheet_name='aero')
pm25    = pd.read_excel('location_metdata_aod_pm25.xlsx',sheet_name='pm25')
#bc_data      = pd.read_excel('bc_measured_data.xlsx',sheet_name='Sheet1')

plt.scatter(met.lon,met.lat,color=   'gray',         alpha=0.5,s=70, edgecolors ='white',transform=ccrs.PlateCarree())
#plt.scatter(aer.lon,aer.lat,color=   'lightgrey', alpha=0.9,s=70, edgecolors ="grey",transform=ccrs.PlateCarree())
#plt.scatter(pm25.lon,pm25.lat,color= 'orange',    alpha=0.9,s=60, edgecolors ="orange",transform=ccrs.PlateCarree())

#==== Add text ========
ax.text(110, 30, 'EAS' ,color='silver', fontsize=10, transform=ccrs.PlateCarree())
ax.text(80,  50, 'NAS' ,color='silver', fontsize=10, transform=ccrs.PlateCarree())
ax.text(50,  40, 'CAME' ,color='silver', fontsize=10, transform=ccrs.PlateCarree())
ax.text(80,  20, 'SAS' ,color='silver', fontsize=10,  transform=ccrs.PlateCarree())
ax.text(85,  33, 'TP' ,color='grey',   fontsize=10,  transform=ccrs.PlateCarree())

# *must* call draw in order to get the axis boundary used to add ticks:
fig.canvas.draw()
#ax.stock_img()

# Define gridline locations and draw the lines using cartopy's built-in gridliner:#
xticks = [20,30, 40, 50,60,70,80,90,100,110,120,130,140]
yticks = [0, 10,20,30, 40,50,60]

s = ax.gridlines(xlocs=xticks, ylocs=yticks,linewidth=1, color='grey', alpha=0.4, linestyle='dashed')

s.xlabel_style = {'color': 'red'} #, 'weigh
s.xlabel_style = {'size': 0.5, 'color': 'black'}
s.ylabel_style = {'size': 0.5, 'color': 'black'}
# Label the end-points of the gridlines using the custom tick makers:
ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
lambert_xticks(ax, xticks)
lambert_yticks(ax, yticks)

cax,kw = matplotlib.colorbar.make_axes(ax,location='right',pad=0.01,shrink=0.9)
cbar=fig.colorbar(cs,cax=cax,**kw)
cbar.set_label('Elevation [m]',size=8)
cbar.ax.tick_params(labelsize=7)

#========= Add plot ========
#====== Load Data==============
bc_data = pd.read_excel('bc_measured_data_new.xlsx',sheet_name='Sheet1')
lon = bc_data['lon']
lat = bc_data['lat']
bc =  bc_data['bc']
#===== Make array=========
lons = np.array(lon)
lats = np.array(lat)
bcs  = np.array(bc)

#==== this is an inset axes over the main axes=====================
sub_axes = plt.axes([0.03, .62, 0.3, 0.32]) ##left, bottom, width, height
m = Basemap(projection='ortho', lat_0=30, lon_0=86,resolution='h')
x,y=m(lons,lats)
m.drawmapboundary(fill_color='#4D4D4D') # fill to edge
m.drawcountries(linewidth=0.5)
m.fillcontinents(color='white',lake_color='black',zorder=0)
sub_axes.scatter(x,y,s=bcs**1.3,c='red',marker="o",alpha=0.5)

#=== Legend Assets ========
#======Decorate label box=======================
l1 = plt.scatter([],[], s=1**1.3, edgecolors='none',color='red',alpha=0.5)
l2 = plt.scatter([],[], s=2**1.3, edgecolors='none',color='red',alpha=0.5)
l3 = plt.scatter([],[], s=5**1.3, edgecolors='none',color='red',alpha=0.5)
l4 = plt.scatter([],[], s=7**1.3, edgecolors='none',color='red',alpha=0.5)
l5 = plt.scatter([],[], s=9**1.3, edgecolors='none',color='red',alpha=0.5)
l6 = plt.scatter([],[], s=11**1.3, edgecolors='none',color='red',alpha=0.5)
l7 = plt.scatter([],[], s=14**1.3, edgecolors='none',color='red',alpha=0.5)


labels  = ["1", "2", "5", "7","9","11","14"]
leggend = plt.legend([l1, l2, l3, l4,l5,l6,l7], labels, frameon=True,handlelength=2,handleheight=2,fontsize=4, title_fontsize=4,
          bbox_to_anchor=(0, 0.5), loc = 'center left', facecolor="#8F8F8F",borderpad = 0.3, handletextpad=0.3, title='$BC [μg m ^{-3}]$')



#=============== Draw regional box ================
#fig.subplots_adjust(top=0.995,
#                        bottom=0.1,
#                        left=0.04,
#                        right=0.93,
#                        hspace=0.03,
#                        wspace=0.01)

plt.savefig('/mnt/g/2nd_Paper/studyArea600dpi', dpi=600)
#plt.show()

