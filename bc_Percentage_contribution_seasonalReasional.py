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
#================ Load East China (CA) data  ================
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
#================ Load East China (NC) data  ================
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
#================ Load East China (SA) data  ================
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

## =========== Get BC data from Control experiment============
bc1_c = wrf.getvar(nc1_c, 'BC1',timeidx=-1)[0]  ## 0 means at surface
bc2_c = wrf.getvar(nc1_c, 'BC2',timeidx=-1)[0]
bc3_c = wrf.getvar(nc2_c, 'BC1',timeidx=-1)[0]
bc4_c = wrf.getvar(nc2_c, 'BC2',timeidx=-1)[0]
bc5_c = wrf.getvar(nc3_c, 'BC1',timeidx=1)[0]
bc6_c = wrf.getvar(nc3_c, 'BC2',timeidx=1)[0]
bc7_c = wrf.getvar(nc4_c, 'BC1',timeidx=1)[0]
bc8_c = wrf.getvar(nc4_c, 'BC2',timeidx=1)[0]
bc9_c = wrf.getvar(nc5_c, 'BC1',timeidx=1)[0]
bc10_c = wrf.getvar(nc5_c, 'BC2',timeidx=1)[0]
bc11_c = wrf.getvar(nc6_c, 'BC1',timeidx=-1)[0]
bc12_c = wrf.getvar(nc6_c, 'BC2',timeidx=-1)[0]
bc13_c = wrf.getvar(nc7_c, 'BC1',timeidx=-1)[0]
bc14_c = wrf.getvar(nc7_c, 'BC2',timeidx=-1)[0]
bc15_c = wrf.getvar(nc8_c, 'BC1',timeidx=-1)[0]
bc16_c = wrf.getvar(nc8_c, 'BC2',timeidx=-1)[0]
bc17_c = wrf.getvar(nc9_c, 'BC1',timeidx=-1)[0]
bc18_c = wrf.getvar(nc9_c, 'BC2',timeidx=-1)[0]
bc19_c= wrf.getvar(nc10_c, 'BC1',timeidx=-1)[0]
bc20_c= wrf.getvar(nc10_c, 'BC2',timeidx=-1)[0]
bc21_c= wrf.getvar(nc11_c, 'BC1',timeidx=-1)[0]
bc22_c= wrf.getvar(nc11_c, 'BC2',timeidx=-1)[0]
bc23_c = wrf.getvar(nc12_c, 'BC1',timeidx=-1)[0]
bc24_c = wrf.getvar(nc12_c, 'BC2',timeidx=-1)[0]
#======== Get BC data from EC  ========================
bc1_ec = wrf.getvar(nc1_ec, 'BC1',timeidx=-1)[0]
bc2_ec = wrf.getvar(nc1_ec, 'BC2',timeidx=-1)[0]
bc3_ec = wrf.getvar(nc2_ec, 'BC1',timeidx=-1)[0]
bc4_ec = wrf.getvar(nc2_ec, 'BC2',timeidx=-1)[0]
bc5_ec = wrf.getvar(nc3_ec, 'BC1',timeidx=-1)[0]
bc6_ec = wrf.getvar(nc3_ec, 'BC2',timeidx=-1)[0]
bc7_ec = wrf.getvar(nc4_ec, 'BC1',timeidx=-1)[0]
bc8_ec = wrf.getvar(nc4_ec, 'BC2',timeidx=-1)[0]
bc9_ec = wrf.getvar(nc5_ec, 'BC1',timeidx=-1)[0]
bc10_ec = wrf.getvar(nc5_ec, 'BC2',timeidx=-1)[0]
bc11_ec = wrf.getvar(nc6_ec, 'BC1',timeidx=-1)[0]
bc12_ec = wrf.getvar(nc6_ec, 'BC2',timeidx=-1)[0]
bc13_ec = wrf.getvar(nc7_ec, 'BC1',timeidx=-1)[0]
bc14_ec = wrf.getvar(nc7_ec, 'BC2',timeidx=-1)[0]
bc15_ec = wrf.getvar(nc8_ec, 'BC1',timeidx=-1)[0]
bc16_ec = wrf.getvar(nc8_ec, 'BC2',timeidx=-1)[0]
bc17_ec = wrf.getvar(nc9_ec, 'BC1',timeidx=-1)[0]
bc18_ec = wrf.getvar(nc9_ec, 'BC2',timeidx=-1)[0]
bc19_ec = wrf.getvar(nc10_ec, 'BC1',timeidx=-1)[0]
bc20_ec = wrf.getvar(nc10_ec, 'BC2',timeidx=-1)[0]
bc21_ec = wrf.getvar(nc11_ec, 'BC1',timeidx=-1)[0]
bc22_ec = wrf.getvar(nc11_ec, 'BC2',timeidx=-1)[0]
bc23_ec = wrf.getvar(nc12_ec, 'BC1',timeidx=-1)[0]
bc24_ec = wrf.getvar(nc12_ec, 'BC2',timeidx=-1)[0]
#======== Get BC data from CA  ========================
bc1_ca = wrf.getvar(nc1_ca, 'BC1',timeidx=-1)[0]
bc2_ca = wrf.getvar(nc1_ca, 'BC2',timeidx=-1)[0]
bc3_ca = wrf.getvar(nc2_ca, 'BC1',timeidx=-1)[0]
bc4_ca = wrf.getvar(nc2_ca, 'BC2',timeidx=-1)[0]
bc5_ca = wrf.getvar(nc3_ca, 'BC1',timeidx=-1)[0]
bc6_ca = wrf.getvar(nc3_ca, 'BC2',timeidx=-1)[0]
bc7_ca = wrf.getvar(nc4_ca, 'BC1',timeidx=-1)[0]
bc8_ca = wrf.getvar(nc4_ca, 'BC2',timeidx=-1)[0]
bc9_ca = wrf.getvar(nc5_ca, 'BC1',timeidx=-1)[0]
bc10_ca = wrf.getvar(nc5_ca, 'BC2',timeidx=-1)[0]
bc11_ca = wrf.getvar(nc6_ca, 'BC1',timeidx=-1)[0]
bc12_ca = wrf.getvar(nc6_ca, 'BC2',timeidx=-1)[0]
bc13_ca = wrf.getvar(nc7_ca, 'BC1',timeidx=-1)[0]
bc14_ca = wrf.getvar(nc7_ca, 'BC2',timeidx=-1)[0]
bc15_ca = wrf.getvar(nc8_ca, 'BC1',timeidx=-1)[0]
bc16_ca = wrf.getvar(nc8_ca, 'BC2',timeidx=-1)[0]
bc17_ca = wrf.getvar(nc9_ca, 'BC1',timeidx=-1)[0]
bc18_ca = wrf.getvar(nc9_ca, 'BC2',timeidx=-1)[0]
bc19_ca = wrf.getvar(nc10_ca, 'BC1',timeidx=-1)[0]
bc20_ca = wrf.getvar(nc10_ca, 'BC2',timeidx=-1)[0]
bc21_ca = wrf.getvar(nc11_ca, 'BC1',timeidx=-1)[0]
bc22_ca = wrf.getvar(nc11_ca, 'BC2',timeidx=-1)[0]
bc23_ca = wrf.getvar(nc12_ca, 'BC1',timeidx=-1)[0]
bc24_ca = wrf.getvar(nc12_ca, 'BC2',timeidx=-1)[0]
#======== Get BC data from NC  ========================
bc1_nc = wrf.getvar(nc1_nc, 'BC1',timeidx=10)[0]
bc2_nc = wrf.getvar(nc1_nc, 'BC2',timeidx=10)[0]
bc3_nc = wrf.getvar(nc2_nc, 'BC1',timeidx=10)[0]
bc4_nc = wrf.getvar(nc2_nc, 'BC2',timeidx=10)[0]
bc5_nc = wrf.getvar(nc3_nc, 'BC1',timeidx=-1)[0]
bc6_nc = wrf.getvar(nc3_nc, 'BC2',timeidx=-1)[0]
bc7_nc = wrf.getvar(nc4_nc, 'BC1',timeidx=-1)[0]
bc8_nc = wrf.getvar(nc4_nc, 'BC2',timeidx=-1)[0]
bc9_nc = wrf.getvar(nc5_nc, 'BC1',timeidx=-1)[0]
bc10_nc = wrf.getvar(nc5_nc, 'BC2',timeidx=-1)[0]
bc11_nc = wrf.getvar(nc6_nc, 'BC1',timeidx=-1)[0]
bc12_nc = wrf.getvar(nc6_nc, 'BC2',timeidx=-1)[0]
bc13_nc = wrf.getvar(nc7_nc, 'BC1',timeidx=-1)[0]
bc14_nc = wrf.getvar(nc7_nc, 'BC2',timeidx=-1)[0]
bc15_nc = wrf.getvar(nc8_nc, 'BC1',timeidx=-1)[0]
bc16_nc = wrf.getvar(nc8_nc, 'BC2',timeidx=-1)[0]
bc17_nc = wrf.getvar(nc9_nc, 'BC1',timeidx=-1)[0]
bc18_nc = wrf.getvar(nc9_nc, 'BC2',timeidx=-1)[0]
bc19_nc = wrf.getvar(nc10_nc, 'BC1',timeidx=-1)[0]
bc20_nc = wrf.getvar(nc10_nc, 'BC2',timeidx=-1)[0]
bc21_nc = wrf.getvar(nc11_nc, 'BC1',timeidx=-1)[0]
bc22_nc = wrf.getvar(nc11_nc, 'BC2',timeidx=-1)[0]
bc23_nc = wrf.getvar(nc12_nc, 'BC1',timeidx=10)[0]
bc24_nc = wrf.getvar(nc12_nc, 'BC2',timeidx=10)[0]
#======== Get BC data from SA  ========================
bc1_sa = wrf.getvar(nc1_sa, 'BC1',timeidx=-1)[0]
bc2_sa = wrf.getvar(nc1_sa, 'BC2',timeidx=-1)[0]
bc3_sa = wrf.getvar(nc2_sa, 'BC1',timeidx=-1)[0]
bc4_sa = wrf.getvar(nc2_sa, 'BC2',timeidx=-1)[0]
bc5_sa = wrf.getvar(nc3_sa, 'BC1',timeidx=1)[0]
bc6_sa = wrf.getvar(nc3_sa, 'BC2',timeidx=1)[0]
bc7_sa = wrf.getvar(nc4_sa, 'BC1',timeidx=1)[0]
bc8_sa = wrf.getvar(nc4_sa, 'BC2',timeidx=1)[0]
bc9_sa = wrf.getvar(nc5_sa, 'BC1',timeidx=1)[0]
bc10_sa = wrf.getvar(nc5_sa, 'BC2',timeidx=1)[0]
bc11_sa = wrf.getvar(nc6_sa, 'BC1',timeidx=-1)[0]
bc12_sa = wrf.getvar(nc6_sa, 'BC2',timeidx=-1)[0]
bc13_sa = wrf.getvar(nc7_sa, 'BC1',timeidx=-1)[0]
bc14_sa = wrf.getvar(nc7_sa, 'BC2',timeidx=-1)[0]
bc15_sa = wrf.getvar(nc8_sa, 'BC1',timeidx=-1)[0]
bc16_sa = wrf.getvar(nc8_sa, 'BC2',timeidx=-1)[0]
bc17_sa = wrf.getvar(nc9_sa, 'BC1',timeidx=-1)[0]
bc18_sa = wrf.getvar(nc9_sa, 'BC2',timeidx=-1)[0]
bc19_sa = wrf.getvar(nc10_sa, 'BC1',timeidx=-1)[0]
bc20_sa = wrf.getvar(nc10_sa, 'BC2',timeidx=-1)[0]
bc21_sa = wrf.getvar(nc11_sa, 'BC1',timeidx=-1)[0]
bc22_sa = wrf.getvar(nc11_sa, 'BC2',timeidx=-1)[0]
bc23_sa = wrf.getvar(nc12_sa, 'BC1',timeidx=-1)[0]
bc24_sa = wrf.getvar(nc12_sa, 'BC2',timeidx=-1)[0]

##====== Adding up BC1 + BC2 ===========
#= Control===
jan_c = bc1_c + bc2_c
feb_c = bc3_c + bc4_c
mar_c = bc5_c + bc6_c
apr_c = bc7_c + bc8_c
may_c = bc9_c + bc10_c
jun_c = bc11_c + bc12_c
jul_c = bc13_c + bc14_c
aug_c = bc15_c + bc16_c
sep_c = bc17_c + bc18_c
oct_c = bc19_c + bc20_c
nov_c = bc21_c + bc22_c
dec_c = bc23_c + bc24_c
#= EC====================
jan_ec = bc1_ec + bc2_ec
feb_ec = bc3_ec + bc4_ec
mar_ec = bc5_ec + bc6_ec
apr_ec = bc7_ec + bc8_ec
may_ec = bc9_ec + bc10_ec
jun_ec = bc11_ec + bc12_ec
jul_ec = bc13_ec + bc14_ec
aug_ec = bc15_ec + bc16_ec
sep_ec = bc17_ec + bc18_ec
oct_ec = bc19_ec + bc20_ec
nov_ec = bc21_ec + bc22_ec
dec_ec = bc23_ec + bc24_ec

jan_e = (jan_c - jan_ec)/jan_c*100
feb_e = (feb_c - feb_ec)/feb_c*100
mar_e = (mar_c - mar_ec)/mar_c*100
apr_e = (apr_c - apr_ec)/apr_c*100
may_e = (may_c - may_ec)/may_c*100
jun_e = (jun_c - jun_ec)/jun_c*100
jul_e = (jul_c - jul_ec)/jul_c*100
aug_e = (aug_c - aug_ec)/aug_c*100
sep_e = (sep_c - sep_ec)/sep_c*100
oct_e = (oct_c - oct_ec)/oct_c*100
nov_e = (nov_c - nov_ec)/nov_c*100
dec_e = (dec_c - dec_ec)/dec_c*100

winter_ec = (jan_e+feb_e+dec_e)/3
spring_ec = (mar_e+apr_e+may_e)/3
summer_ec = (jun_e+jul_e+aug_e)/3
autumn_ec = (sep_e+oct_e+nov_e)/3

#== NC==
jan_nc = bc1_nc + bc2_nc
feb_nc = bc3_nc + bc4_nc
mar_nc = bc5_nc + bc6_nc
apr_nc = bc7_nc + bc8_nc
may_nc = bc9_nc + bc10_nc
jun_nc = bc11_nc + bc12_nc
jul_nc = bc13_nc + bc14_nc
aug_nc = bc15_nc + bc16_nc
sep_nc = bc17_nc + bc18_nc
oct_nc = bc19_nc + bc20_nc
nov_nc = bc21_nc + bc22_nc
dec_nc = bc23_nc + bc24_nc

jan_n = (jan_c - jan_nc)/jan_c*100
feb_n = (feb_c - feb_nc)/feb_c*100
mar_n = (mar_c - mar_nc)/mar_c*100
apr_n = (apr_c - apr_nc)/apr_c*100
may_n = (may_c - may_nc)/may_c*100
jun_n = (jun_c - jun_nc)/jun_c*100
jul_n = (jul_c - jul_nc)/jul_c*100
aug_n = (aug_c - aug_nc)/aug_c*100
sep_n = (sep_c - sep_nc)/sep_c*100
oct_n = (oct_c - oct_nc)/oct_c*100
nov_n = (nov_c - nov_nc)/nov_c*100
dec_n = (dec_c - dec_nc)/dec_c*100

winter_nc = (jan_n+feb_n+dec_n)/3
spring_nc = (mar_n+apr_n+may_n)/3
summer_nc = (jun_n+jul_n+aug_n)/3
autumn_nc = (sep_n+oct_n+nov_n)/3

#== CA===
jan_ca = bc1_ca + bc2_ca
feb_ca = bc3_ca + bc4_ca
mar_ca = bc5_ca + bc6_ca
apr_ca = bc7_ca + bc8_ca
may_ca = bc9_ca + bc10_ca
jun_ca = bc11_ca + bc12_ca
jul_ca = bc13_ca + bc14_ca
aug_ca = bc15_ca + bc16_ca
sep_ca = bc17_ca + bc18_ca
oct_ca = bc19_ca + bc20_ca
nov_ca = bc21_ca + bc22_ca
dec_ca = bc23_ca + bc24_ca

jan_c1 = (jan_c - jan_ca)/jan_c*100
feb_c1 = (feb_c - feb_ca)/feb_c*100
mar_c1 = (mar_c - mar_ca)/mar_c*100
apr_c1 = (apr_c - apr_ca)/apr_c*100
may_c1 = (may_c - may_ca)/may_c*100
jun_c1 = (jun_c - jun_ca)/jun_c*100
jul_c1 = (jul_c - jul_ca)/jul_c*100
aug_c1 = (aug_c - aug_ca)/aug_c*100
sep_c1 = (sep_c - sep_ca)/sep_c*100
oct_c1 = (oct_c - oct_ca)/oct_c*100
nov_c1 = (nov_c - nov_ca)/nov_c*100
dec_c1 = (dec_c - dec_ca)/dec_c*100

winter_ca = (jan_c1+feb_c1+dec_c1)/3
spring_ca = (mar_c1+apr_c1+may_c1)/3
summer_ca = (jun_c1+jul_c1+aug_c1)/3
autumn_ca = (sep_c1+oct_c1+nov_c1)/3

#== SA ==========
jan_sa = bc1_sa + bc2_sa
feb_sa = bc3_sa + bc4_sa
mar_sa = bc5_sa + bc6_sa
apr_sa = bc7_sa + bc8_sa
may_sa = bc9_sa + bc10_sa
jun_sa = bc11_sa + bc12_sa
jul_sa = bc13_sa + bc14_sa
aug_sa = bc15_sa + bc16_sa
sep_sa = bc17_sa + bc18_sa
oct_sa = bc19_sa + bc20_sa
nov_sa = bc21_sa + bc22_sa
dec_sa = bc23_sa + bc24_sa

jan_s = (jan_c - jan_sa)/jan_c*100
feb_s = (feb_c - feb_sa)/feb_c*100
mar_s = (mar_c - mar_sa)/mar_c*100
apr_s = (apr_c - apr_sa)/apr_c*100
may_s = (may_c - may_sa)/may_c*100
jun_s = (jun_c - jun_sa)/jun_c*100
jul_s = (jul_c - jul_sa)/jul_c*100
aug_s = (aug_c - aug_sa)/aug_c*100
sep_s = (sep_c - sep_sa)/sep_c*100
oct_s = (oct_c - oct_sa)/oct_c*100
nov_s = (nov_c - nov_sa)/nov_c*100
dec_s = (dec_c - dec_sa)/dec_c*100

winter_sa = (jan_s+feb_s+dec_s)/3
spring_sa = (mar_s+apr_s+may_s)/3
summer_sa = (jun_s+jul_s+aug_s)/3
autumn_sa = (sep_s+oct_s+nov_s)/3

#=== Now plot =========================
fig,ax = plt.subplots(4,4, figsize = (12,8))
gridspec.GridSpec(4,4)
lats, lons = wrf.latlon_coords(bc1_c)
#==== Winter==========================
#==== EC========
plt.subplot2grid((4,4), (0,0))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,winter_ec,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed)   #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.title('EA',fontsize=12)
plt.ylabel('Winter', labelpad=40,fontsize=12)
plt.clim(0,100)
#=== NC==========
plt.subplot2grid((4,4), (0,1))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,winter_nc,cmap=cmaps.MPL_GnBu) #WhiteBlueGreenYellowRed)  # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.title('NA',fontsize=12)
plt.clim(0,100)
#== CA=======
plt.subplot2grid((4,4), (0,2))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,winter_ca,cmap=cmaps.MPL_GnBu) #precip3_16lev)  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.title('CAME',fontsize=12)
plt.clim(0,100)
#== SA=======
plt.subplot2grid((4,4), (0,3))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,winter_sa,cmap=cmaps.MPL_GnBu) #cmaps.precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.title('SA',fontsize=12)
plt.clim(0,100)
#== Spring=================
#==== EC========
plt.subplot2grid((4,4), (1,0))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,spring_ec,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.ylabel('Spring', labelpad=40,fontsize=12)
plt.clim(0,100)
#=== NC==========
plt.subplot2grid((4,4), (1,1))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,spring_nc,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,100)
#=== CA==========
plt.subplot2grid((4,4), (1,2))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,spring_ca,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,100)
#=== SA==========
plt.subplot2grid((4,4), (1,3))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,spring_sa,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,100)
#==== Summer ==========
#==== EC========
plt.subplot2grid((4,4), (2,0))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,summer_ec,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.ylabel('Summer', labelpad=40,fontsize=12)
plt.clim(0,100)
#==== NC========
plt.subplot2grid((4,4), (2,1))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,summer_nc,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,100)
#==== CA========
plt.subplot2grid((4,4), (2,2))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,summer_ca,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,100)
#==== SA========
plt.subplot2grid((4,4), (2,3))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,summer_sa,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,100)
#====== Autumn=====================
#=== EC==
plt.subplot2grid((4,4), (3,0))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,autumn_ec,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10,color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.ylabel('Autumn', labelpad=40,fontsize=12)
plt.clim(0,100)
#=== NC==
plt.subplot2grid((4,4), (3,1))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,autumn_nc,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10,color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,100)
#=== CA==
plt.subplot2grid((4,4), (3,2))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,autumn_ca,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10,color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,100)
#=== SA==
plt.subplot2grid((4,4), (3,3))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
           lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
           llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')

#lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax1 = m.pcolormesh(x,y,autumn_sa,cmap=cmaps.MPL_GnBu) #precip3_16lev) #WhiteBlueGreenYellowRed) #amwg) #  ncl_default) # WhiteBlueGreenYellowRed)  #ncl_default)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.5, dashes=[4, 1], color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10,color='black')
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,100)

fig.subplots_adjust(top=0.96,
			bottom=0.145,
			left=0.085,
			right=0.995,
			hspace=0.0,
			wspace=0.02)


#===== Color bar====================
cax = fig.add_axes([0.1, 0.085,0.88,0.03])   # left, bottom, width, and height
cb = fig.colorbar(ax1,  cax=cax, orientation='horizontal',ticks=[0,10,20,30,40,50,60,70,80,90,100])
cb.set_label(label='Contribution Rate [%]',size=12)

plt.savefig("/mnt/g/2nd_Paper/BC_relative_Contribution.png",dpi=1000)

#plt.show()

