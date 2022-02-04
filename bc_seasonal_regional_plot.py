#========== Plotting BC data from WRF output ===================================================
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
from matplotlib import rcParams

"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"

#================ Load control data =========================================
nc1_c = Dataset("wrfout_d01_2017-01-01_00_00_00")
nc2_c = Dataset("wrfout_d01_2017-02-01_00_00_00")
nc3_c = Dataset("wrfout_d01_2017-03-01_00_00_00.nc")
nc4_c = Dataset("wrfout_d01_2017-04-01_00_00_00.nc")
nc5_c = Dataset("wrfout_d01_2017-05-01_00_00_00.nc")
nc6_c = Dataset("wrfout_d01_2017-06-01_00_00_00")
nc7_c = Dataset("wrfout_d01_2017-07-01_00_00_00")
nc8_c = Dataset("wrfout_d01_2017-08-01_00_00_00")
nc9_c = Dataset("wrfout_d01_2017-09-01_00_00_00")
nc10_c = Dataset("wrfout_d01_2017-10-01_00_00_00")
nc11_c = Dataset("wrfout_d01_2017-11-01_00_00_00")
nc12_c = Dataset("wrfout_d01_2017-12-01_00_00_00")
#=============== Load North China (NC) data ==============
nc1_nc = Dataset("wrfout_d01_2017-01-01_00_00_00-NC")
nc2_nc = Dataset("wrfout_d01_2017-02-01_00_00_00-NC")
nc3_nc = Dataset("wrfout_d01_2017-03-01_00_00_00-NC")
nc4_nc = Dataset("wrfout_d01_2017-04-01_00_00_00-NC")
nc5_nc = Dataset("wrfout_d01_2017-05-01_00_00_00-NC")
nc6_nc = Dataset("wrfout_d01_2017-06-01_00_00_00-NC")
nc7_nc = Dataset("wrfout_d01_2017-07-01_00_00_00-NC")
nc8_nc = Dataset("wrfout_d01_2017-08-01_00_00_00-NC")
nc9_nc = Dataset("wrfout_d01_2017-09-01_00_00_00-NC")
nc10_nc = Dataset("wrfout_d01_2017-10-01_00_00_00-NC")
nc11_nc = Dataset("wrfout_d01_2017-11-01_00_00_00-NC")
nc12_nc = Dataset("wrfout_d01_2017-12-01_00_00_00-NC")
#=============== Load Central Asia (CA) data  ==============
nc1_ca = Dataset("wrfout_d01_2017-01-01_00_00_00-CA")
nc2_ca = Dataset("wrfout_d01_2017-02-01_00_00_00-CA")
nc3_ca = Dataset("wrfout_d01_2017-03-01_00_00_00-CA")
nc4_ca = Dataset("wrfout_d01_2017-04-01_00_00_00-CA")
nc5_ca = Dataset("wrfout_d01_2017-05-01_00_00_00-CA")
nc6_ca = Dataset("wrfout_d01_2017-06-01_00_00_00-CA")
nc7_ca = Dataset("wrfout_d01_2017-07-01_00_00_00-CA")
nc8_ca = Dataset("wrfout_d01_2017-08-01_00_00_00-CA")
nc9_ca = Dataset("wrfout_d01_2017-09-01_00_00_00-CA")
nc10_ca = Dataset("wrfout_d01_2017-10-01_00_00_00-CA")
nc11_ca = Dataset("wrfout_d01_2017-11-01_00_00_00-CA")
nc12_ca = Dataset("wrfout_d01_2017-12-01_00_00_00-CA")
#=============== Load South Asia (SA) data  =============
nc1_sa = Dataset("wrfout_d01_2017-01-01_00_00_00-SA")
nc2_sa = Dataset("wrfout_d01_2017-02-01_00_00_00-SA")
nc3_sa = Dataset("wrfout_d01_2017-03-01_00_00_00-SA")
nc4_sa = Dataset("wrfout_d01_2017-04-01_00_00_00-SA")
nc5_sa = Dataset("wrfout_d01_2017-05-01_00_00_00-SA")
nc6_sa = Dataset("wrfout_d01_2017-06-01_00_00_00-SA")
nc7_sa = Dataset("wrfout_d01_2017-07-01_00_00_00-SA")
nc8_sa = Dataset("wrfout_d01_2017-08-01_00_00_00-SA")
nc9_sa = Dataset("wrfout_d01_2017-09-01_00_00_00-SA")
nc10_sa = Dataset("wrfout_d01_2017-10-01_00_00_00-SA")
nc11_sa = Dataset("wrfout_d01_2017-11-01_00_00_00-SA")
nc12_sa = Dataset("wrfout_d01_2017-12-01_00_00_00-SA")
#================ Load East China (EC) data  ================
nc1_ec = Dataset("wrfout_d01_2017-01-01_00_00_00-EC")
nc2_ec = Dataset("wrfout_d01_2017-02-01_00_00_00-EC")
nc3_ec = Dataset("wrfout_d01_2017-03-01_00_00_00-EC")
nc4_ec = Dataset("wrfout_d01_2017-04-01_00_00_00-EC")
nc5_ec = Dataset("wrfout_d01_2017-05-01_00_00_00-EC")
nc6_ec = Dataset("wrfout_d01_2017-06-01_00_00_00-EC")
nc7_ec = Dataset("wrfout_d01_2017-07-01_00_00_00-EC")
nc8_ec = Dataset("wrfout_d01_2017-08-01_00_00_00-EC")
nc9_ec = Dataset("wrfout_d01_2017-09-01_00_00_00-EC")
nc10_ec = Dataset("wrfout_d01_2017-10-01_00_00_00-EC")
nc11_ec = Dataset("wrfout_d01_2017-11-01_00_00_00-EC")
nc12_ec = Dataset("wrfout_d01_2017-12-01_00_00_00-EC")
#========================================================
## =========== Get BC data from Control experiment============
bc1_c = wrf.getvar(nc1_c, 'BC1',timeidx=1)[0]  ## 0 means at surface
bc2_c = wrf.getvar(nc1_c, 'BC2',timeidx=1)[0]
bc3_c = wrf.getvar(nc2_c, 'BC1',timeidx=1)[0]
bc4_c = wrf.getvar(nc2_c, 'BC2',timeidx=1)[0]
bc5_c = wrf.getvar(nc3_c, 'BC1',timeidx=1)[0]
bc6_c = wrf.getvar(nc3_c, 'BC2',timeidx=1)[0]
bc7_c = wrf.getvar(nc4_c, 'BC1',timeidx=1)[0]
bc8_c = wrf.getvar(nc4_c, 'BC2',timeidx=1)[0]
bc9_c = wrf.getvar(nc5_c, 'BC1',timeidx=1)[0]
bc10_c = wrf.getvar(nc5_c, 'BC2',timeidx=1)[0]
bc11_c = wrf.getvar(nc6_c, 'BC1',timeidx=1)[0]
bc12_c = wrf.getvar(nc6_c, 'BC2',timeidx=1)[0]
bc13_c = wrf.getvar(nc7_c, 'BC1',timeidx=1)[0]
bc14_c = wrf.getvar(nc7_c, 'BC2',timeidx=1)[0]
bc15_c = wrf.getvar(nc8_c, 'BC1',timeidx=1)[0]
bc16_c = wrf.getvar(nc8_c, 'BC2',timeidx=1)[0]
bc17_c = wrf.getvar(nc9_c, 'BC1',timeidx=1)[0]
bc18_c = wrf.getvar(nc9_c, 'BC2',timeidx=1)[0]
bc19_c= wrf.getvar(nc10_c, 'BC1',timeidx=1)[0]
bc20_c= wrf.getvar(nc10_c, 'BC2',timeidx=1)[0]
bc21_c= wrf.getvar(nc11_c, 'BC1',timeidx=1)[0]
bc22_c= wrf.getvar(nc11_c, 'BC2',timeidx=1)[0]
bc23_c = wrf.getvar(nc12_c, 'BC1',timeidx=1)[0]
bc24_c = wrf.getvar(nc12_c, 'BC2',timeidx=1)[0]
## =========== Get BC data from NC ============
bc1_nc = wrf.getvar(nc1_nc, 'BC1',timeidx=1)[0]  
bc2_nc = wrf.getvar(nc1_nc, 'BC2',timeidx=1)[0]
bc3_nc = wrf.getvar(nc2_nc, 'BC1',timeidx=1)[0]
bc4_nc = wrf.getvar(nc2_nc, 'BC2',timeidx=1)[0]
bc5_nc = wrf.getvar(nc3_nc, 'BC1',timeidx=1)[0]
bc6_nc = wrf.getvar(nc3_nc, 'BC2',timeidx=1)[0]
bc7_nc = wrf.getvar(nc4_nc, 'BC1',timeidx=1)[0]
bc8_nc = wrf.getvar(nc4_nc, 'BC2',timeidx=1)[0]
bc9_nc = wrf.getvar(nc5_nc, 'BC1',timeidx=1)[0]
bc10_nc = wrf.getvar(nc5_nc, 'BC2',timeidx=1)[0]
bc11_nc = wrf.getvar(nc6_nc, 'BC1',timeidx=1)[0]
bc12_nc = wrf.getvar(nc6_nc, 'BC2',timeidx=1)[0]
bc13_nc = wrf.getvar(nc7_nc, 'BC1',timeidx=1)[0]
bc14_nc = wrf.getvar(nc7_nc, 'BC2',timeidx=1)[0]
bc15_nc = wrf.getvar(nc8_nc, 'BC1',timeidx=1)[0]
bc16_nc = wrf.getvar(nc8_nc, 'BC2',timeidx=1)[0]
bc17_nc = wrf.getvar(nc9_nc, 'BC1',timeidx=1)[0]
bc18_nc = wrf.getvar(nc9_nc, 'BC2',timeidx=1)[0]
bc19_nc= wrf.getvar(nc10_nc, 'BC1',timeidx=1)[0]
bc20_nc= wrf.getvar(nc10_nc, 'BC2',timeidx=1)[0]
bc21_nc= wrf.getvar(nc11_nc, 'BC1',timeidx=1)[0]
bc22_nc= wrf.getvar(nc11_nc, 'BC2',timeidx=1)[0]
bc23_nc = wrf.getvar(nc12_nc, 'BC1',timeidx=1)[0]
bc24_nc = wrf.getvar(nc12_nc, 'BC2',timeidx=1)[0]
## =========== Get BC data from CA ============
bc1_ca = wrf.getvar(nc1_ca, 'BC1',timeidx=1)[0]  
bc2_ca = wrf.getvar(nc1_ca, 'BC2',timeidx=1)[0]
bc3_ca = wrf.getvar(nc2_ca, 'BC1',timeidx=1)[0]
bc4_ca = wrf.getvar(nc2_ca, 'BC2',timeidx=1)[0]
bc5_ca = wrf.getvar(nc3_ca, 'BC1',timeidx=1)[0]
bc6_ca = wrf.getvar(nc3_ca, 'BC2',timeidx=1)[0]
bc7_ca = wrf.getvar(nc4_ca, 'BC1',timeidx=1)[0]
bc8_ca = wrf.getvar(nc4_ca, 'BC2',timeidx=1)[0]
bc9_ca = wrf.getvar(nc5_ca, 'BC1',timeidx=1)[0]
bc10_ca = wrf.getvar(nc5_ca, 'BC2',timeidx=1)[0]
bc11_ca = wrf.getvar(nc6_ca, 'BC1',timeidx=1)[0]
bc12_ca = wrf.getvar(nc6_ca, 'BC2',timeidx=1)[0]
bc13_ca = wrf.getvar(nc7_ca, 'BC1',timeidx=1)[0]
bc14_ca = wrf.getvar(nc7_ca, 'BC2',timeidx=1)[0]
bc15_ca = wrf.getvar(nc8_ca, 'BC1',timeidx=1)[0]
bc16_ca = wrf.getvar(nc8_ca, 'BC2',timeidx=1)[0]
bc17_ca = wrf.getvar(nc9_ca, 'BC1',timeidx=1)[0]
bc18_ca = wrf.getvar(nc9_ca, 'BC2',timeidx=1)[0]
bc19_ca = wrf.getvar(nc10_ca, 'BC1',timeidx=1)[0]
bc20_ca = wrf.getvar(nc10_ca, 'BC2',timeidx=1)[0]
bc21_ca = wrf.getvar(nc11_ca, 'BC1',timeidx=1)[0]
bc22_ca = wrf.getvar(nc11_ca, 'BC2',timeidx=1)[0]
bc23_ca = wrf.getvar(nc12_ca, 'BC1',timeidx=1)[0]
bc24_ca = wrf.getvar(nc12_ca, 'BC2',timeidx=1)[0]
#======== Get BC data from SA ========================
bc1_sa = wrf.getvar(nc1_sa, 'BC1',timeidx=1)[0]  
bc2_sa = wrf.getvar(nc1_sa, 'BC2',timeidx=1)[0]
bc3_sa = wrf.getvar(nc2_sa, 'BC1',timeidx=1)[0]
bc4_sa = wrf.getvar(nc2_sa, 'BC2',timeidx=1)[0]
bc5_sa = wrf.getvar(nc3_sa, 'BC1',timeidx=1)[0]
bc6_sa = wrf.getvar(nc3_sa, 'BC2',timeidx=1)[0]
bc7_sa = wrf.getvar(nc4_sa, 'BC1',timeidx=1)[0]
bc8_sa = wrf.getvar(nc4_sa, 'BC2',timeidx=1)[0]
bc9_sa = wrf.getvar(nc5_sa, 'BC1',timeidx=1)[0]
bc10_sa = wrf.getvar(nc5_sa, 'BC2',timeidx=1)[0]
bc11_sa = wrf.getvar(nc6_sa, 'BC1',timeidx=1)[0]
bc12_sa = wrf.getvar(nc6_sa, 'BC2',timeidx=1)[0]
bc13_sa = wrf.getvar(nc7_sa, 'BC1',timeidx=1)[0]
bc14_sa = wrf.getvar(nc7_sa, 'BC2',timeidx=1)[0]
bc15_sa = wrf.getvar(nc8_sa, 'BC1',timeidx=1)[0]
bc16_sa = wrf.getvar(nc8_sa, 'BC2',timeidx=1)[0]
bc17_sa = wrf.getvar(nc9_sa, 'BC1',timeidx=1)[0]
bc18_sa = wrf.getvar(nc9_sa, 'BC2',timeidx=1)[0]
bc19_sa = wrf.getvar(nc10_sa, 'BC1',timeidx=1)[0]
bc20_sa = wrf.getvar(nc10_sa, 'BC2',timeidx=1)[0]
bc21_sa = wrf.getvar(nc11_sa, 'BC1',timeidx=1)[0]
bc22_sa = wrf.getvar(nc11_sa, 'BC2',timeidx=1)[0]
bc23_sa = wrf.getvar(nc12_sa, 'BC1',timeidx=1)[0]
bc24_sa = wrf.getvar(nc12_sa, 'BC2',timeidx=1)[0]
#======== Get BC data from EC  ========================
bc1_ec = wrf.getvar(nc1_ec, 'BC1',timeidx=2)[0]
bc2_ec = wrf.getvar(nc1_ec, 'BC2',timeidx=2)[0]
bc3_ec = wrf.getvar(nc2_ec, 'BC1',timeidx=2)[0]
bc4_ec = wrf.getvar(nc2_ec, 'BC2',timeidx=2)[0]
bc5_ec = wrf.getvar(nc3_ec, 'BC1',timeidx=2)[0]
bc6_ec = wrf.getvar(nc3_ec, 'BC2',timeidx=2)[0]
bc7_ec = wrf.getvar(nc4_ec, 'BC1',timeidx=2)[0]
bc8_ec = wrf.getvar(nc4_ec, 'BC2',timeidx=2)[0]
bc9_ec = wrf.getvar(nc5_ec, 'BC1',timeidx=2)[0]
bc10_ec = wrf.getvar(nc5_ec, 'BC2',timeidx=2)[0]
bc11_ec = wrf.getvar(nc6_ec, 'BC1',timeidx=2)[0]
bc12_ec = wrf.getvar(nc6_ec, 'BC2',timeidx=2)[0]
bc13_ec = wrf.getvar(nc7_ec, 'BC1',timeidx=2)[0]
bc14_ec = wrf.getvar(nc7_ec, 'BC2',timeidx=2)[0]
bc15_ec = wrf.getvar(nc8_ec, 'BC1',timeidx=2)[0]
bc16_ec = wrf.getvar(nc8_ec, 'BC2',timeidx=2)[0]
bc17_ec = wrf.getvar(nc9_ec, 'BC1',timeidx=2)[0]
bc18_ec = wrf.getvar(nc9_ec, 'BC2',timeidx=2)[0]
bc19_ec = wrf.getvar(nc10_ec, 'BC1',timeidx=2)[0]
bc20_ec = wrf.getvar(nc10_ec, 'BC2',timeidx=2)[0]
bc21_ec = wrf.getvar(nc11_ec, 'BC1',timeidx=2)[0]
bc22_ec = wrf.getvar(nc11_ec, 'BC2',timeidx=2)[0]
bc23_ec = wrf.getvar(nc12_ec, 'BC1',timeidx=2)[0]
bc24_ec = wrf.getvar(nc12_ec, 'BC2',timeidx=2)[0]
#======================================================
alt1  = getvar(nc1_c, "ALT", timeidx=-1)[0]
alt2  = getvar(nc2_c, "ALT", timeidx=-1)[0]
alt3  = getvar(nc3_c, "ALT", timeidx=-1)[0]
alt4  = getvar(nc4_c, "ALT", timeidx=-1)[0]
alt5  = getvar(nc5_c, "ALT", timeidx=-1)[0]
alt6  = getvar(nc6_c, "ALT", timeidx=-1)[0]
#alt7  = getvar(nc7, "ALT", timeidx=1)[0]
alt8  = getvar(nc8_c, "ALT", timeidx=-1)[0]
alt9  = getvar(nc9_c, "ALT", timeidx=-1)[0]
alt10  = getvar(nc10_c, "ALT", timeidx=-1)[0]
alt11  = getvar(nc11_c, "ALT", timeidx=-1)[0]
alt12  = getvar(nc12_c, "ALT", timeidx=-1)[0]

#==
winter_alt = (alt1+alt2+alt12)/3
spring_alt = (alt3+alt4+alt5)/3
summer_alt = (alt6+alt8)/2
autumn_alt = (alt9+alt10+alt11)/3
#============ BC data for all ============================================
winter_c  =  (bc1_c + bc2_c + bc3_c + bc4_c + bc23_c + bc24_c)/3
winter_nc1 = (bc1_nc + bc2_nc + bc3_nc + bc4_nc + bc23_nc + bc24_nc)/3
winter_ca1 = (bc1_ca + bc2_ca + bc3_ca + bc4_ca + bc23_ca + bc24_ca)/3
winter_sa1 = (bc1_sa + bc2_sa + bc3_sa + bc4_sa + bc23_sa + bc24_sa)/3
winter_ec1 = (bc1_ec + bc2_ec + bc3_ec + bc4_ec + bc23_ec + bc24_ec)/3
#=========================================================================
winter_nc2 = (winter_c-winter_nc1)/winter_alt
winter_ca2 = (winter_c-winter_ca1)/winter_alt
winter_sa2 = (winter_c-winter_sa1)/winter_alt
winter_ec2 = (winter_c-winter_ec1)/winter_alt

#====== smoothing winter ============
sigma_y = 2
sigma_x = 2
sigma = [sigma_y, sigma_x]

winter_nc = sp.ndimage.filters.gaussian_filter(winter_nc2, sigma, mode='constant')
winter_ca = sp.ndimage.filters.gaussian_filter(winter_ca2, sigma, mode='constant')
winter_sa = sp.ndimage.filters.gaussian_filter(winter_sa2, sigma, mode='constant')
winter_ec = sp.ndimage.filters.gaussian_filter(winter_ec2, sigma, mode='constant')

#==========================================================================
spring_c =   (bc5_c + bc6_c + bc7_c + bc8_c + bc9_c + bc10_c)/3
spring_nc1 = (bc5_nc + bc6_nc + bc7_nc + bc8_nc + bc9_nc + bc10_nc)/3
spring_ca1 = (bc5_ca + bc6_ca + bc7_ca + bc8_ca + bc9_ca + bc10_ca)/3
spring_sa1 = (bc5_sa + bc6_sa + bc7_sa + bc8_sa + bc9_sa + bc10_sa)/3
spring_ec1 = (bc5_ec + bc6_ec + bc7_ec + bc8_ec + bc9_ec + bc10_ec)/3
#=====================
spring_nc2 = (spring_c - spring_nc1)/spring_alt
spring_ca2 = (spring_c - spring_ca1)/spring_alt
spring_sa2 = (spring_c - spring_sa1)/spring_alt
spring_ec2 = (spring_c - spring_ec1)/spring_alt

#==== Smoothing Spring ==========
spring_nc = sp.ndimage.filters.gaussian_filter(spring_nc2, sigma, mode='constant')
spring_ca = sp.ndimage.filters.gaussian_filter(spring_ca2, sigma, mode='constant')
spring_sa = sp.ndimage.filters.gaussian_filter(spring_sa2, sigma, mode='constant')
spring_ec = sp.ndimage.filters.gaussian_filter(spring_ec2, sigma, mode='constant')
#========================================================================
summer_c =   (bc11_c + bc12_c + bc13_c + bc14_c + bc15_c + bc16_c)/3
summer_nc1 = (bc11_nc + bc12_nc + bc13_nc + bc14_nc + bc15_nc + bc16_nc)/3
summer_ca1 = (bc11_ca + bc12_ca + bc13_ca + bc14_ca + bc15_ca + bc16_ca)/3
summer_sa1 = (bc11_sa + bc12_sa + bc13_sa + bc14_sa + bc15_sa + bc16_sa)/3
summer_ec1 = (bc11_ec + bc12_ec + bc13_ec + bc14_ec + bc15_ec + bc16_ec)/3
#=============================
summer_nc2 = (summer_c - summer_nc1)/summer_alt
summer_ca2 = (summer_c - summer_ca1)/summer_alt
summer_sa2 = (summer_c - summer_sa1)/summer_alt
summer_ec2 = (summer_c - summer_ec1)/summer_alt

summer_nc = sp.ndimage.filters.gaussian_filter(summer_nc2, sigma, mode='constant')
summer_ca = sp.ndimage.filters.gaussian_filter(summer_ca2, sigma, mode='constant')
summer_sa = sp.ndimage.filters.gaussian_filter(summer_sa2, sigma, mode='constant')
summer_ec = sp.ndimage.filters.gaussian_filter(summer_ec2, sigma, mode='constant')
#==========================================================================
autumn_c =  (bc17_c + bc18_c + bc19_c + bc20_c + bc21_c + bc22_c)/3
autumn_nc1 = (bc17_nc + bc18_nc + bc19_nc + bc20_nc + bc21_nc + bc22_nc)/3
autumn_ca1 = (bc17_ca + bc18_ca + bc19_ca + bc20_ca + bc21_ca + bc22_ca)/3
autumn_sa1 = (bc17_sa + bc18_sa + bc19_sa + bc20_sa + bc21_sa + bc22_sa)/3
autumn_ec1 = (bc17_ec + bc18_ec + bc19_ec + bc20_ec + bc21_ec + bc22_ec)/3
#============================
autumn_nc2 = (autumn_c - autumn_nc1)/autumn_alt
autumn_ca2 = (autumn_c - autumn_ca1)/autumn_alt
autumn_sa2 = (autumn_c - autumn_sa1)/autumn_alt
autumn_ec2 = (autumn_c - autumn_ec1)/autumn_alt

#===== Smoothing Autumn =====================
autumn_nc = sp.ndimage.filters.gaussian_filter(autumn_nc2, sigma, mode='constant')
autumn_ca = sp.ndimage.filters.gaussian_filter(autumn_ca2, sigma, mode='constant')
autumn_sa = sp.ndimage.filters.gaussian_filter(autumn_sa2, sigma, mode='constant')
autumn_ec = sp.ndimage.filters.gaussian_filter(autumn_ec2, sigma, mode='constant')
#============ Get BC data based on two season ==================================
## For EC- Monsoon -spring
#mon_ec = (bc9_ec + bc10_ec+bc11_ec + bc12_ec + bc13_ec + bc14_ec + bc15_ec + bc16_ec+bc17_ec + bc18_ec)/5
#nonmon_ec = (bc1_ec + bc2_ec + bc3_ec + bc4_ec + bc23_ec + bc24_ec+bc5_ec + bc6_ec + bc7_ec + bc8_ec+bc19_ec + bc20_ec + bc21_ec + bc22_ec)/7
#mon_c = (bc9_c + bc10_c+bc11_c + bc12_c + bc13_c + bc14_c + bc15_c + bc16_c+bc17_c + bc18_c)/5
#nonmon_c = (bc1_c + bc2_c + bc3_c + bc4_c + bc23_c + bc24_c+bc5_c + bc6_c + bc7_c + bc8_c+bc19_c + bc20_c + bc21_c + bc22_c)/7
#=============================================================================================================================================
"""#====== Plot seasonal BC contribution regionally  ============="""

fig, axs = plt.subplots(4, 4,figsize=(7,4),dpi=300)
gridspec.GridSpec(4,4)

#================ Plot for Nort China ==============================
plt.subplot2grid((4,4), (0,0))
#============= Winter ===========================================
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))

ax1 = m.pcolormesh(x,y,winter_ca,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15),linewidth=0.2, dashes=[2, 2], labels=[1,0,0,0], fontsize=5, color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.2, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.ylabel('CAME', labelpad=15,fontsize=5)
plt.title('Winter',fontsize=5,pad=0)
plt.clim(0,8)
#================Spring  ==============================
plt.subplot2grid((4,4),(0,1))
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,spring_ca ,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2],  fontsize=4,color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2, dashes=[2, 2],  fontsize=4,color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.title('Spring',fontsize=5,pad=0)
plt.clim(0,8)
#================ Summer  ==============================
plt.subplot2grid((4,4), (0,2))
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,summer_ca ,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15),linewidth=0.2, dashes=[2, 2],  color='black',fontsize=4)
m.drawmeridians(np.arange(45, 125,15),linewidth=0.2, dashes=[2, 2],  color='black',fontsize=4)
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.title('Summer',fontsize=5,pad=0)
plt.clim(0,8)
#================ autumn  ==============================
plt.subplot2grid((4,4), (0,3))
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,autumn_ca ,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.title('Autumn',fontsize=5,pad=0)
plt.clim(0,8)
##======== CA ============================================================ 
plt.subplot2grid((4,4), (1,0))
#============= Winter ===========================================
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,winter_nc,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2], labels=[1,0,0,0], fontsize=5, color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.clim(0,8)
plt.ylabel('NA', labelpad=15,fontsize=5)
#================Spring  ==============================
plt.subplot2grid((4,4), (1,1))
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,spring_nc ,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2,  dashes=[2, 2], color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.clim(0,8)
#================ Summer  ==============================
plt.subplot2grid((4,4), (1,2))
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,summer_nc ,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.clim(0,8)
#================ autumn  ==============================
plt.subplot2grid((4,4), (1,3))
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,autumn_nc ,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.clim(0,8)
##======== SA ============================================================
plt.subplot2grid((4,4), (2,0))
#============= Winter ===========================================
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,winter_sa,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2], labels=[1,0,0,0], fontsize=5, color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.clim(0,8)
plt.ylabel('SA', labelpad=15,fontsize=5)
#================Spring  ==============================
plt.subplot2grid((4,4), (2,1))
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,spring_sa ,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.clim(0,8)
#============= summer ===========================================
plt.subplot2grid((4,4), (2,2))
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,summer_sa,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.clim(0,8)
#================Autumn  ==============================
plt.subplot2grid((4,4), (2,3))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,autumn_sa ,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2],  color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2, dashes=[2, 2],  color='black',fontsize=5)
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.clim(0,8)
##======== EC ============================================================
plt.subplot2grid((4,4), (3,0))
#============= Winter ===========================================
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,winter_ec,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2], labels=[1,0,0,0], fontsize=5, color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2, dashes=[2, 2], labels=[0,0,0,1], fontsize=5, color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.clim(0,8)
plt.ylabel('EA', labelpad=15,fontsize=5)
#================Spring  ==============================
plt.subplot2grid((4,4), (3,1))
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,spring_ec ,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2, dashes=[2, 2], labels=[0,0,0,1], fontsize=5, color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.clim(0,8)
#============= summer ===========================================
plt.subplot2grid((4,4), (3,2))
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,summer_ec,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.2, dashes=[2, 2], color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.2, dashes=[2, 2], labels=[0,0,0,1], fontsize=5, color='black')
m.drawcountries(linewidth=0.2)
m.drawcoastlines(linewidth=0.2)
plt.clim(0,8)
#============= autumn ===========================================
plt.subplot2grid((4,4), (3,3))
lats, lons = wrf.latlon_coords(winter_c)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,autumn_ec,cmap=cmaps.WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15),linewidth=0.2, dashes=[2, 2], color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.2, dashes=[2, 2], labels=[0,0,0,1], fontsize=5, color='black')
m.drawcountries(0.2)
m.drawcoastlines(linewidth=0.2)
plt.clim(0,8)


#cax = fig.add_axes([0.88,0.04,0.012,0.9295])   # left, bottom, width, and height
cax = fig.add_axes([0.93,0.04,0.012,0.9295])

cbar= fig.colorbar(ax1, ax=axs[:,:], cax=cax,ticks=[0,1,2,3,4,5,6,7,8,9,10]) #,7,8,9,10])
cbar.set_label(r'BC  [$μg  m^{-3}$]',fontsize=5,labelpad=5)
cbar.ax.tick_params(labelsize=5)

fig.subplots_adjust(top=0.975,
			bottom=0.04,
			left=0.04,
			right=0.935,
			hspace=0.045,
			wspace=0.0)

plt.savefig("/mnt/g/2nd_Paper/bBC_seasonal_regional.png",dpi=300)
#plt.show()
