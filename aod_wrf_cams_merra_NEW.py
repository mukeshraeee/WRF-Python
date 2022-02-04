#====== Code for AOD plot WRF-CAMS-MERRA-2 ===========================
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
#========= Winter=======================================
nc1 = Dataset("wrfout_d01_2017-01-01_00_00_00")
nc2 = Dataset("wrfout_d01_2017-02-01_00_00_00")
nc3 = Dataset("wrfout_d01_2017-12-01_00_00_00")
#====== Spring =========================================
nc4 = Dataset("wrfout_d01_2017-03-01_00_00_00.nc")
nc5 = Dataset("wrfout_d01_2017-04-01_00_00_00.nc")
nc6 = Dataset("wrfout_d01_2017-05-01_00_00_00.nc")
#=======Summer ========================================
nc7 = Dataset("wrfout_d01_2017-06-01_00_00_00")
nc8 = Dataset("wrfout_d01_2017-07-01_00_00_00")
nc9 = Dataset("wrfout_d01_2017-08-01_00_00_00")
#======= Autumn ========================================
nc10 = Dataset("wrfout_d01_2017-09-01_00_00_00")
nc11 = Dataset("wrfout_d01_2017-10-01_00_00_00")
nc12 = Dataset("wrfout_d01_2017-11-01_00_00_00")
#====== 
jan = getvar(nc1, "EXTCOF55", timeidx=-1)
feb = getvar(nc2, "EXTCOF55", timeidx=-1)
dec = getvar(nc12, "EXTCOF55", timeidx=-1)
mar = getvar(nc3, "EXTCOF55", timeidx=-1)
apr = getvar(nc4, "EXTCOF55", timeidx=-1)
may = getvar(nc5, "EXTCOF55", timeidx=-1)
jun = getvar(nc6, "EXTCOF55", timeidx=-1)
jul = getvar(nc7, "EXTCOF55", timeidx=-1)
aug = getvar(nc8, "EXTCOF55", timeidx=-1)
sep = getvar(nc9, "EXTCOF55", timeidx=-1)
oct = getvar(nc10, "EXTCOF55", timeidx=-1)
nov = getvar(nc11, "EXTCOF55", timeidx=-1)
#==== Getting model height ==============
z1 = getvar(nc1,"z",timeidx=-1)
z2 = getvar(nc2,"z",timeidx=-1)
z3 = getvar(nc12,"z",timeidx=-1)
z4 = getvar(nc3,"z",timeidx=-1)
z5 = getvar(nc4,"z",timeidx=-1)
z6 = getvar(nc5,"z",timeidx=-1)
z7 = getvar(nc6,"z",timeidx=-1)
z8 = getvar(nc7,"z",timeidx=-1)
z9 = getvar(nc8,"z",timeidx=-1)
z10 = getvar(nc9,"z",timeidx=-1)
z11 = getvar(nc10,"z",timeidx=-1)
z12 = getvar(nc11,"z",timeidx=-1)

z = (z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12)/12
ph = getvar(nc1, "PH", timeidx=-1)

#===== Load CAMS AOD Data =====
cams_winter = Dataset('cams_winter_new.nc')
cams_spring = Dataset('cams_spring_new.nc')
cams_summer = Dataset('cams_summer_new.nc')
cams_autumn = Dataset('cams_autumn_new.nc')
## MERRA-2 AOD data ===================
merra_winter = Dataset('merra_winter_new.nc')
merra_spring = Dataset('merra_spring_new.nc')
merra_summer = Dataset('merra_summer_new.nc')
merra_autumn = Dataset('merra_autumn_new.nc')
#=========================================================
aod_cams_winter = cams_winter.variables['aod550'][:]
aod_cams_spring = cams_spring.variables['aod550'][:]
aod_cams_summer = cams_summer.variables['aod550'][:]
aod_cams_autumn = cams_autumn.variables['aod550'][:]
aod_merra_winter = merra_winter.variables['TOTEXTTAU'][:]
aod_merra_spring = merra_spring.variables['TOTEXTTAU'][:]
aod_merra_summer = merra_summer.variables['TOTEXTTAU'][:]
aod_merra_autumn = merra_autumn.variables['TOTEXTTAU'][:]
#======================================================
#========== smooth =========================================
aod_cams_winter_smooth=smooth2d(aod_cams_winter,5, cenweight=5)
aod_cams_spring_smooth=smooth2d(aod_cams_spring,5, cenweight=5)
aod_cams_summer_smooth=smooth2d(aod_cams_summer,5, cenweight=5)
aod_cams_autumn_smooth=smooth2d(aod_cams_autumn,5, cenweight=5)
aod_mera_winter_smooth=smooth2d(aod_merra_winter,5,cenweight=5)
aod_mera_spring_smooth=smooth2d(aod_merra_spring,5,cenweight=5)
aod_mera_summer_smooth=smooth2d(aod_merra_summer,5,cenweight=5)
aod_mera_autumn_smooth=smooth2d(aod_merra_autumn,5,cenweight=5)
#========  Getting lat and long ===============================
lat_cam  = cams_winter.variables['latitude'][:]
lon_cam  = cams_winter.variables['longitude'][:]
lat_merra = merra_winter.variables['lat'][:]
lon_merra = merra_winter.variables['lon'][:]
#=======   average over the first dimension: time ####
aod_cams_winter_av=np.mean(aod_cams_winter_smooth[:,:,:],axis=0)
aod_cams_spring_av=np.mean(aod_cams_spring_smooth[:,:,:],axis=0)
aod_cams_summer_av=np.mean(aod_cams_summer_smooth[:,:,:],axis=0)
aod_cams_autumn_av=np.mean(aod_cams_autumn_smooth[:,:,:],axis=0)
aod_mera_winter_av=np.mean(aod_mera_winter_smooth[:,:,:],axis=0)
aod_mera_spring_av=np.mean(aod_mera_spring_smooth[:,:,:],axis=0)
aod_mera_summer_av=np.mean(aod_mera_summer_smooth[:,:,:],axis=0)
aod_mera_autumn_av=np.mean(aod_mera_autumn_smooth[:,:,:],axis=0)

#====== Get ext for jan From WRF output  ============================
e0_jan = jan[0]
e1_jan = jan[1]
e2_jan = jan[2]
e3_jan = jan[3]
e4_jan = jan[4]
e5_jan = jan[5]
e6_jan = jan[6]
e7_jan = jan[7]
e8_jan = jan[8]
e9_jan = jan[9]
e10_jan = jan[10]
e11_jan = jan[11]
e12_jan = jan[12]
e13_jan = jan[13]
e14_jan=  jan[14]
e15_jan = jan[15]
e16_jan = jan[16]
e17_jan = jan[17]
e18_jan = jan[18]
#====== Get ext for jan ============================
e0_feb = feb[0]
e1_feb = feb[1]
e2_feb = feb[2]
e3_feb = feb[3]
e4_feb = feb[4]
e5_feb = feb[5]
e6_feb = feb[6]
e7_feb = feb[7]
e8_feb = feb[8]
e9_feb = feb[9]
e10_feb = feb[10]
e11_feb = feb[11]
e12_feb = feb[12]
e13_feb = feb[13]
e14_feb=  feb[14]
e15_feb = feb[15]
e16_feb = feb[16]
e17_feb = feb[17]
e18_feb = feb[18]
#====== Get ext for jan ============================
e0_dec = dec[0]
e1_dec = dec[1]
e2_dec = dec[2]
e3_dec = dec[3]
e4_dec = dec[4]
e5_dec = dec[5]
e6_dec = dec[6]
e7_dec = dec[7]
e8_dec = dec[8]
e9_dec = dec[9]
e10_dec = dec[10]
e11_dec = dec[11]
e12_dec = dec[12]
e13_dec = dec[13]
e14_dec=  dec[14]
e15_dec = dec[15]
e16_dec = dec[16]
e17_dec = dec[17]
e18_dec = dec[18]
#====== Get ext for jan ============================
e0_mar = mar[0]
e1_mar = mar[1]
e2_mar = mar[2]
e3_mar = mar[3]
e4_mar = mar[4]
e5_mar = mar[5]
e6_mar = mar[6]
e7_mar = mar[7]
e8_mar = mar[8]
e9_mar = mar[9]
e10_mar = mar[10]
e11_mar = mar[11]
e12_mar = mar[12]
e13_mar = mar[13]
e14_mar=  mar[14]
e15_mar = mar[15]
e16_mar = mar[16]
e17_mar = mar[17]
e18_mar = mar[18]
#====== Get ext for jan ============================
e0_apr = apr[0]
e1_apr = apr[1]
e2_apr = apr[2]
e3_apr = apr[3]
e4_apr = apr[4]
e5_apr = apr[5]
e6_apr = apr[6]
e7_apr = apr[7]
e8_apr = apr[8]
e9_apr = apr[9]
e10_apr = apr[10]
e11_apr = apr[11]
e12_apr = apr[12]
e13_apr = apr[13]
e14_apr=  apr[14]
e15_apr = apr[15]
e16_apr = apr[16]
e17_apr = apr[17]
e18_apr = apr[18]
#====== Get ext for jan ============================
e0_may = may[0]
e1_may = may[1]
e2_may = may[2]
e3_may = may[3]
e4_may = may[4]
e5_may = may[5]
e6_may = may[6]
e7_may = may[7]
e8_may = may[8]
e9_may = may[9]
e10_may = may[10]
e11_may = may[11]
e12_may = may[12]
e13_may = may[13]
e14_may=  may[14]
e15_may = may[15]
e16_may = may[16]
e17_may = may[17]
e18_may = may[18]
#====== Get ext for jan ============================
e0_jun = jun[0]
e1_jun = jun[1]
e2_jun = jun[2]
e3_jun = jun[3]
e4_jun = jun[4]
e5_jun = jun[5]
e6_jun = jun[6]
e7_jun = jun[7]
e8_jun = jun[8]
e9_jun = jun[9]
e10_jun = jun[10]
e11_jun = jun[11]
e12_jun = jun[12]
e13_jun = jun[13]
e14_jun=  jun[14]
e15_jun = jun[15]
e16_jun = jun[16]
e17_jun = jun[17]
e18_jun = jun[18]
#====== Get ext for jan ============================
e0_jul = jul[0]
e1_jul = jul[1]
e2_jul = jul[2]
e3_jul = jul[3]
e4_jul = jul[4]
e5_jul = jul[5]
e6_jul = jul[6]
e7_jul = jul[7]
e8_jul = jul[8]
e9_jul = jul[9]
e10_jul = jul[10]
e11_jul = jul[11]
e12_jul = jul[12]
e13_jul = jul[13]
e14_jul = jul[14]
e15_jul = jul[15]
e16_jul = jul[16]
e17_jul = jul[17]
e18_jul = jul[18]
#====== Get ext for jan ============================
e0_aug = aug[0]
e1_aug = aug[1]
e2_aug = aug[2]
e3_aug = aug[3]
e4_aug = aug[4]
e5_aug = aug[5]
e6_aug = aug[6]
e7_aug = aug[7]
e8_aug = aug[8]
e9_aug = aug[9]
e10_aug = aug[10]
e11_aug = aug[11]
e12_aug = aug[12]
e13_aug = aug[13]
e14_aug=  aug[14]
e15_aug = aug[15]
e16_aug = aug[16]
e17_aug = aug[17]
e18_aug = aug[18]
#====== Get ext for jan ============================
e0_sep = sep[0]
e1_sep = sep[1]
e2_sep = sep[2]
e3_sep = sep[3]
e4_sep = sep[4]
e5_sep = sep[5]
e6_sep = sep[6]
e7_sep = sep[7]
e8_sep = sep[8]
e9_sep = sep[9]
e10_sep = sep[10]
e11_sep = sep[11]
e12_sep = sep[12]
e13_sep = sep[13]
e14_sep=  sep[14]
e15_sep = sep[15]
e16_sep = sep[16]
e17_sep = sep[17]
e18_sep = sep[18]
#====== Get ext for jan ============================
e0_oct = oct[0]
e1_oct = oct[1]
e2_oct = oct[2]
e3_oct = oct[3]
e4_oct = oct[4]
e5_oct = oct[5]
e6_oct = oct[6]
e7_oct = oct[7]
e8_oct = oct[8]
e9_oct = oct[9]
e10_oct = oct[10]
e11_oct = oct[11]
e12_oct = oct[12]
e13_oct = oct[13]
e14_oct = oct[14]
e15_oct = oct[15]
e16_oct = oct[16]
e17_oct = oct[17]
e18_oct = oct[18]
#====== Get ext for jan ============================
e0_nov = nov[0]
e1_nov = nov[1]
e2_nov = nov[2]
e3_nov = nov[3]
e4_nov = nov[4]
e5_nov = nov[5]
e6_nov = nov[6]
e7_nov = nov[7]
e8_nov = nov[8]
e9_nov = nov[9]
e10_nov = nov[10]
e11_nov = nov[11]
e12_nov = nov[12]
e13_nov = nov[13]
e14_nov = nov[14]
e15_nov = nov[15]
e16_nov = nov[16]
e17_nov = nov[17]
e18_nov = nov[18]
#=============================================
#==== Now take every layer thickness =========
z0 = z[0]
z1 = z[1]
z2 = z[2]
z3 = z[3]
z4 = z[4]
z5 = z[5]
z6 = z[6]
z7 = z[7]
z8 = z[8]
z9 = z[9]
z10 = z[10]
z11 = z[11]
z12 = z[12]
z13 = z[13]
z14 = z[14]
z15 = z[15]
z16 = z[16]
z17 = z[17]
z18 = z[18]
#======= AOD -Jan====================
ex0_jan = e0_jan*(z1-z0)*0.001 # since aod is unitless 
ex1_jan = e1_jan*(z2-z1)*0.001
ex2_jan = e2_jan*(z3-z2)*0.001
ex3_jan = e3_jan*(z4-z3)*0.001
ex4_jan = e4_jan*(z5-z4)*0.001
ex5_jan = e5_jan*(z6-z5)*0.001
ex6_jan = e6_jan*(z7-z6)*0.001
ex7_jan = e7_jan*(z8-z7)*0.001
ex8_jan = e8_jan*(z9-z8)*0.001
ex9_jan = e9_jan*(z10-z9)*0.001
ex10_jan = e10_jan*(z11-z10)*0.001
ex11_jan = e11_jan*(z12-z11)*0.001
ex12_jan = e12_jan*(z13-z12)*0.001
ex13_jan = e13_jan*(z14-z13)*0.001
ex14_jan = e14_jan*(z15-z14)*0.001
ex15_jan = e15_jan*(z16-z15)*0.001
ex16_jan = e16_jan*(z17-z0)*0.001
ex17_jan = e17_jan*(z17-z18)*0.001
##===== AOD -Feb =======================
ex0_feb = e0_feb*(z1-z0)*0.001 # since aod is unitless
ex1_feb = e1_feb*(z2-z1)*0.001
ex2_feb = e2_feb*(z3-z2)*0.001
ex3_feb = e3_feb*(z4-z3)*0.001
ex4_feb = e4_feb*(z5-z4)*0.001
ex5_feb = e5_feb*(z6-z5)*0.001
ex6_feb = e6_feb*(z7-z6)*0.001
ex7_feb = e7_feb*(z8-z7)*0.001
ex8_feb = e8_feb*(z9-z8)*0.001
ex9_feb = e9_feb*(z10-z9)*0.001
ex10_feb = e10_feb*(z11-z10)*0.001
ex11_feb = e11_feb*(z12-z11)*0.001
ex12_feb = e12_feb*(z13-z12)*0.001
ex13_feb = e13_feb*(z14-z13)*0.001
ex14_feb = e14_feb*(z15-z14)*0.001
ex15_feb = e15_feb*(z16-z15)*0.001
ex16_feb = e16_feb*(z17-z0)*0.001
ex17_feb = e17_feb*(z17-z18)*0.001
##===== AOD - Dec =======================
ex0_dec = e0_dec*(z1-z0)*0.001 # since aod is unitless
ex1_dec = e1_dec*(z2-z1)*0.001
ex2_dec = e2_dec*(z3-z2)*0.001
ex3_dec = e3_dec*(z4-z3)*0.001
ex4_dec = e4_dec*(z5-z4)*0.001
ex5_dec = e5_dec*(z6-z5)*0.001
ex6_dec = e6_dec*(z7-z6)*0.001
ex7_dec = e7_dec*(z8-z7)*0.001
ex8_dec = e8_dec*(z9-z8)*0.001
ex9_dec = e9_dec*(z10-z9)*0.001
ex10_dec = e10_dec*(z11-z10)*0.001
ex11_dec = e11_dec*(z12-z11)*0.001
ex12_dec = e12_dec*(z13-z12)*0.001
ex13_dec = e13_dec*(z14-z13)*0.001
ex14_dec = e14_dec*(z15-z14)*0.001
ex15_dec = e15_dec*(z16-z15)*0.001
ex16_dec = e16_dec*(z17-z0)*0.001
ex17_dec = e17_dec*(z17-z18)*0.001
#======= AOD - March====================
ex0_mar = e0_mar*(z1-z0)*0.001 # since aod is unitless
ex1_mar = e1_mar*(z2-z1)*0.001
ex2_mar = e2_mar*(z3-z2)*0.001
ex3_mar = e3_mar*(z4-z3)*0.001
ex4_mar = e4_mar*(z5-z4)*0.001
ex5_mar = e5_mar*(z6-z5)*0.001
ex6_mar = e6_mar*(z7-z6)*0.001
ex7_mar = e7_mar*(z8-z7)*0.001
ex8_mar = e8_mar*(z9-z8)*0.001
ex9_mar = e9_mar*(z10-z9)*0.001
ex10_mar = e10_mar*(z11-z10)*0.001
ex11_mar = e11_mar*(z12-z11)*0.001
ex12_mar = e12_mar*(z13-z12)*0.001
ex13_mar = e13_mar*(z14-z13)*0.001
ex14_mar = e14_mar*(z15-z14)*0.001
ex15_mar = e15_mar*(z16-z15)*0.001
ex16_mar = e16_mar*(z17-z0)*0.001
ex17_mar = e17_mar*(z17-z18)*0.001
#======= AOD April -====================
ex0_apr = e0_apr*(z1-z0)*0.001 # since aod is unitless
ex1_apr = e1_apr*(z2-z1)*0.001
ex2_apr = e2_apr*(z3-z2)*0.001
ex3_apr = e3_apr*(z4-z3)*0.001
ex4_apr = e4_apr*(z5-z4)*0.001
ex5_apr = e5_apr*(z6-z5)*0.001
ex6_apr = e6_apr*(z7-z6)*0.001
ex7_apr = e7_apr*(z8-z7)*0.001
ex8_apr = e8_apr*(z9-z8)*0.001
ex9_apr = e9_apr*(z10-z9)*0.001
ex10_apr = e10_apr*(z11-z10)*0.001
ex11_apr = e11_apr*(z12-z11)*0.001
ex12_apr = e12_apr*(z13-z12)*0.001
ex13_apr = e13_apr*(z14-z13)*0.001
ex14_apr = e14_apr*(z15-z14)*0.001
ex15_apr = e15_apr*(z16-z15)*0.001
ex16_apr = e16_apr*(z17-z0)*0.001
ex17_apr = e17_apr*(z17-z18)*0.001
#======= AOD -May====================
ex0_may = e0_may*(z1-z0)*0.001 # since aod is unitless
ex1_may = e1_may*(z2-z1)*0.001
ex2_may = e2_may*(z3-z2)*0.001
ex3_may = e3_may*(z4-z3)*0.001
ex4_may = e4_may*(z5-z4)*0.001
ex5_may = e5_may*(z6-z5)*0.001
ex6_may = e6_may*(z7-z6)*0.001
ex7_may = e7_may*(z8-z7)*0.001
ex8_may = e8_may*(z9-z8)*0.001
ex9_may = e9_may*(z10-z9)*0.001
ex10_may = e10_may*(z11-z10)*0.001
ex11_may = e11_may*(z12-z11)*0.001
ex12_may = e12_may*(z13-z12)*0.001
ex13_may = e13_may*(z14-z13)*0.001
ex14_may = e14_may*(z15-z14)*0.001
ex15_may = e15_may*(z16-z15)*0.001
ex16_may = e16_may*(z17-z0)*0.001
ex17_may = e17_may*(z17-z18)*0.001
#======= AOD -Jun====================
ex0_jun = e0_jun*(z1-z0)*0.001 # since aod is unitless
ex1_jun = e1_jun*(z2-z1)*0.001
ex2_jun = e2_jun*(z3-z2)*0.001
ex3_jun  = e3_jun*(z4-z3)*0.001
ex4_jun = e4_jun*(z5-z4)*0.001
ex5_jun = e5_jun*(z6-z5)*0.001
ex6_jun = e6_jun*(z7-z6)*0.001
ex7_jun = e7_jun*(z8-z7)*0.001
ex8_jun = e8_jun*(z9-z8)*0.001
ex9_jun = e9_jun*(z10-z9)*0.001
ex10_jun = e10_jun*(z11-z10)*0.001
ex11_jun = e11_jun*(z12-z11)*0.001
ex12_jun = e12_jun*(z13-z12)*0.001
ex13_jun = e13_jun*(z14-z13)*0.001
ex14_jun = e14_jun*(z15-z14)*0.001
ex15_jun = e15_jun*(z16-z15)*0.001
ex16_jun = e16_jun*(z17-z0)*0.001
ex17_jun = e17_jun*(z17-z18)*0.001
#======= AOD -Jul====================
ex0_jul = e0_jul*(z1-z0)*0.001 # since aod is unitless
ex1_jul = e1_jul*(z2-z1)*0.001
ex2_jul = e2_jul*(z3-z2)*0.001
ex3_jul = e3_jul*(z4-z3)*0.001
ex4_jul = e4_jul*(z5-z4)*0.001
ex5_jul = e5_jul*(z6-z5)*0.001
ex6_jul = e6_jul*(z7-z6)*0.001
ex7_jul = e7_jul*(z8-z7)*0.001
ex8_jul = e8_jul*(z9-z8)*0.001
ex9_jul = e9_jul*(z10-z9)*0.001
ex10_jul = e10_jul*(z11-z10)*0.001
ex11_jul = e11_jul*(z12-z11)*0.001
ex12_jul = e12_jul*(z13-z12)*0.001
ex13_jul = e13_jul*(z14-z13)*0.001
ex14_jul = e14_jul*(z15-z14)*0.001
ex15_jul = e15_jul*(z16-z15)*0.001
ex16_jul = e16_jul*(z17-z0)*0.001
ex17_jul = e17_jul*(z17-z18)*0.001
#======= AOD -Aug ====================
ex0_aug = e0_aug*(z1-z0)*0.001 # since aod is unitless
ex1_aug = e1_aug*(z2-z1)*0.001
ex2_aug = e2_aug*(z3-z2)*0.001
ex3_aug = e3_aug*(z4-z3)*0.001
ex4_aug = e4_aug*(z5-z4)*0.001
ex5_aug = e5_aug*(z6-z5)*0.001
ex6_aug = e6_aug*(z7-z6)*0.001
ex7_aug = e7_aug*(z8-z7)*0.001
ex8_aug = e8_aug*(z9-z8)*0.001
ex9_aug = e9_aug*(z10-z9)*0.001
ex10_aug = e10_aug*(z11-z10)*0.001
ex11_aug = e11_aug*(z12-z11)*0.001
ex12_aug = e12_aug*(z13-z12)*0.001
ex13_aug = e13_aug*(z14-z13)*0.001
ex14_aug = e14_aug*(z15-z14)*0.001
ex15_aug = e15_aug*(z16-z15)*0.001
ex16_aug = e16_aug*(z17-z0)*0.001
ex17_aug = e17_aug*(z17-z18)*0.001
#======= AOD -Sep ====================
ex0_sep = e0_sep*(z1-z0)*0.001 # since aod is unitless
ex1_sep = e1_sep*(z2-z1)*0.001
ex2_sep = e2_sep*(z3-z2)*0.001
ex3_sep = e3_sep*(z4-z3)*0.001
ex4_sep = e4_sep*(z5-z4)*0.001
ex5_sep = e5_sep*(z6-z5)*0.001
ex6_sep = e6_sep*(z7-z6)*0.001
ex7_sep = e7_sep*(z8-z7)*0.001
ex8_sep = e8_sep*(z9-z8)*0.001
ex9_sep = e9_sep*(z10-z9)*0.001
ex10_sep = e10_sep*(z11-z10)*0.001
ex11_sep = e11_sep*(z12-z11)*0.001
ex12_sep = e12_sep*(z13-z12)*0.001
ex13_sep = e13_sep*(z14-z13)*0.001
ex14_sep = e14_sep*(z15-z14)*0.001
ex15_sep = e15_sep*(z16-z15)*0.001
ex16_sep = e16_sep*(z17-z0)*0.001
ex17_sep = e17_sep*(z17-z18)*0.001
#======= AOD - oct  ====================
ex0_oct = e0_oct*(z1-z0)*0.001 # since aod is unitless
ex1_oct = e1_oct*(z2-z1)*0.001
ex2_oct = e2_oct*(z3-z2)*0.001
ex3_oct = e3_oct*(z4-z3)*0.001
ex4_oct = e4_oct*(z5-z4)*0.001
ex5_oct = e5_oct*(z6-z5)*0.001
ex6_oct = e6_oct*(z7-z6)*0.001
ex7_oct = e7_oct*(z8-z7)*0.001
ex8_oct = e8_oct*(z9-z8)*0.001
ex9_oct = e9_oct*(z10-z9)*0.001
ex10_oct = e10_oct*(z11-z10)*0.001
ex11_oct = e11_oct*(z12-z11)*0.001
ex12_oct = e12_oct*(z13-z12)*0.001
ex13_oct = e13_oct*(z14-z13)*0.001
ex14_oct = e14_oct*(z15-z14)*0.001
ex15_oct = e15_oct*(z16-z15)*0.001
ex16_oct = e16_oct*(z17-z0)*0.001
ex17_oct = e17_oct*(z17-z18)*0.001
#======= AOD - nov  ====================
ex0_nov = e0_nov*(z1-z0)*0.001 # since aod is unitless
ex1_nov = e1_nov*(z2-z1)*0.001
ex2_nov = e2_nov*(z3-z2)*0.001
ex3_nov = e3_nov*(z4-z3)*0.001
ex4_nov = e4_nov*(z5-z4)*0.001
ex5_nov = e5_nov*(z6-z5)*0.001
ex6_nov = e6_nov*(z7-z6)*0.001
ex7_nov = e7_nov*(z8-z7)*0.001
ex8_nov = e8_nov*(z9-z8)*0.001
ex9_nov = e9_nov*(z10-z9)*0.001
ex10_nov = e10_nov*(z11-z10)*0.001
ex11_nov = e11_nov*(z12-z11)*0.001
ex12_nov = e12_nov*(z13-z12)*0.001
ex13_nov= e13_nov*(z14-z13)*0.001
ex14_nov = e14_nov*(z15-z14)*0.001
ex15_nov = e15_nov*(z16-z15)*0.001
ex16_nov = e16_nov*(z17-z0)*0.001
ex17_nov = e17_nov*(z17-z18)*0.001
#============================== Monlthy AOD integreted over column ==========================
aod_jan = (ex0_jan+ex1_jan+ex2_jan+ex3_jan+ex4_jan+ex5_jan+ex6_jan+ex7_jan+ex8_jan+ex9_jan+ex10_jan+ex11_jan+ex12_jan+ex13_jan+ex14_jan+ex15_jan+ex16_jan+ex17_jan)
aod_feb = (ex0_feb+ex1_feb+ex2_feb+ex3_feb+ex4_feb+ex5_feb+ex6_feb+ex7_feb+ex8_feb+ex9_feb+ex10_feb+ex11_feb+ex12_feb+ex13_feb+ex14_feb+ex15_feb+ex16_feb+ex17_feb)
aod_mar = (ex0_mar+ex1_mar+ex2_mar+ex3_mar+ex4_mar+ex5_mar+ex6_mar+ex7_mar+ex8_mar+ex9_mar+ex10_mar+ex11_mar+ex12_mar+ex13_mar+ex14_mar+ex15_mar+ex16_mar+ex17_mar)
aod_apr = (ex0_apr+ex1_apr+ex2_apr+ex3_apr+ex4_apr+ex5_apr+ex6_apr+ex7_apr+ex8_apr+ex9_apr+ex10_apr+ex11_apr+ex12_apr+ex13_apr+ex14_apr+ex15_apr+ex16_apr+ex17_apr)
aod_may = (ex0_may+ex1_may+ex2_may+ex3_may+ex4_may+ex5_may+ex6_may+ex7_may+ex8_may+ex9_may+ex10_may+ex11_may+ex12_may+ex13_may+ex14_may+ex15_may+ex16_may+ex17_may)
aod_jun = (ex0_jun+ex1_jun+ex2_jun+ex3_jun+ex4_jun+ex5_jun+ex6_jun+ex7_jun+ex8_jun+ex9_jun+ex10_jun+ex11_jun+ex12_jun+ex13_jun+ex14_jun+ex15_jun+ex16_jun+ex17_jun)
aod_jul = (ex0_jul+ex1_jul+ex2_jul+ex3_jul+ex4_jul+ex5_jul+ex6_jul+ex7_jul+ex8_jul+ex9_jul+ex10_jul+ex11_jul+ex12_jul+ex13_jul+ex14_jul+ex15_jul+ex16_jul+ex17_jul)
aod_aug = (ex0_aug+ex1_aug+ex2_aug+ex3_aug+ex4_aug+ex5_aug+ex6_aug+ex7_aug+ex8_aug+ex9_aug+ex10_aug+ex11_aug+ex12_aug+ex13_aug+ex14_aug+ex15_aug+ex16_aug+ex17_aug)
aod_sep = (ex0_sep+ex1_sep+ex2_sep+ex3_sep+ex4_sep+ex5_sep+ex6_sep+ex7_sep+ex8_sep+ex9_sep+ex10_sep+ex11_sep+ex12_sep+ex13_sep+ex14_sep+ex15_sep+ex16_sep+ex17_sep)
aod_oct = (ex0_oct+ex1_oct+ex2_oct+ex3_oct+ex4_oct+ex5_oct+ex6_oct+ex7_oct+ex8_oct+ex9_oct+ex10_oct+ex11_oct+ex12_oct+ex13_oct+ex14_oct+ex15_oct+ex16_oct+ex17_oct)
aod_nov = (ex0_nov+ex1_nov+ex2_nov+ex3_nov+ex4_nov+ex5_nov+ex6_nov+ex7_nov+ex8_nov+ex9_nov+ex10_nov+ex11_nov+ex12_nov+ex13_nov+ex14_nov+ex15_nov+ex16_nov+ex17_nov)
aod_dec = (ex0_dec+ex1_dec+ex2_dec+ex3_dec+ex4_dec+ex5_dec+ex6_dec+ex7_dec+ex8_dec+ex9_dec+ex10_dec+ex11_dec+ex12_dec+ex13_dec+ex14_dec+ex15_dec+ex16_dec+ex17_dec)

#================== Seasonal AOD ===============
wrf_winter = (aod_dec+aod_jan+aod_feb)/3
wrf_spring = (aod_mar+aod_apr+aod_may)/3
wrf_summer = (aod_jun+aod_jul+aod_aug)/3
wrf_autumn = (aod_sep+aod_oct+aod_nov)/3
#====== Let's smooth by using Gussian filter  ============
sigma_y = 1.0
sigma_x = 1.0
sigma = [sigma_y, sigma_x]

aod_wrf_winter = sp.ndimage.filters.gaussian_filter(wrf_winter, sigma, mode='constant')
aod_wrf_spring = sp.ndimage.filters.gaussian_filter(wrf_spring, sigma, mode='constant')
aod_wrf_summer = sp.ndimage.filters.gaussian_filter(wrf_summer, sigma, mode='constant')
aod_wrf_autumn = sp.ndimage.filters.gaussian_filter(wrf_autumn, sigma, mode='constant')

#====== Sub - Plot ====================
fig, axs = plt.subplots(3, 4,figsize=(10,6),dpi=180)
gridspec.GridSpec(3,4)
#=========== WRF-AOD-Winter ================================
plt.subplot2grid((3,4), (0,0))
lats, lons = latlon_coords(ph)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
            llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2],   color='black')
x, y = m(to_np(lons), to_np(lats))
ax = m.pcolormesh(x,y,aod_wrf_winter,cmap=cmaps.cmocean_balance) #WhiteBlueGreenYellowRed)
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,0.7)
plt.title('Winter',fontsize=8)
plt.ylabel('WRF-AOD', labelpad=5,fontsize=8)

#=========== WRF-AOD-Spring ================================
plt.subplot2grid((3,4), (0,1))
lats, lons = latlon_coords(ph)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
            llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2],   color='black')

x, y = m(to_np(lons), to_np(lats))
ax = m.pcolormesh(x,y,aod_wrf_spring,cmap=cmaps.cmocean_balance) #WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,0.7)
plt.title('Spring',fontsize=8)

#=========== WRF-AOD-Summer ================================
plt.subplot2grid((3,4), (0,2))
lats, lons = latlon_coords(ph)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
            llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2],   color='black')

x, y = m(to_np(lons), to_np(lats))
ax = m.pcolormesh(x,y,aod_wrf_summer,cmap=cmaps.cmocean_balance) #WhiteBlueGreenYellowRed)
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,0.7)
plt.title('Summer',fontsize=8)
#plt.ylabel('WRF-AOD', labelpad=40,fontsize=12)

#=========== WRF-AOD-Autumn ================================
plt.subplot2grid((3,4), (0,3))
lats, lons = latlon_coords(ph)
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
            llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2], labels=[0,1,0,0],fontsize=6, color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2],    color='black')

x, y = m(to_np(lons), to_np(lats))
ax = m.pcolormesh(x,y,aod_wrf_autumn,cmap=cmaps.cmocean_balance) #WhiteBlueGreenYellowRed)
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,0.7)
plt.title('Autumn',fontsize=8)


#=========  CAMS-winter =========================== 
plt.subplot2grid((3,4), (1,0))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
            llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2] ,color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2],   color='black')
lons,lats = np.meshgrid(lon_cam,lat_cam)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,aod_cams_winter_av,cmap=cmaps.cmocean_balance) #WhiteBlueGreenYellowRed)
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)
plt.ylabel('CAMS-AOD', labelpad=5,fontsize=8)

#=========  CAMS-Spring ===========================
plt.subplot2grid((3,4), (1,1))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
            llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2],  color='black')
lons,lats = np.meshgrid(lon_cam,lat_cam)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,aod_cams_spring_av,cmap=cmaps.cmocean_balance) #WhiteBlueGreenYellowRed)
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)
#=========  CAMS-Summer ===========================
plt.subplot2grid((3,4), (1,2))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
            llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2],  color='black')
lons,lats = np.meshgrid(lon_cam,lat_cam)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,aod_cams_summer_av,cmap=cmaps.cmocean_balance) #WhiteBlueGreenYellowRed)
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)

#=========  CAMS-Autumn ===========================
plt.subplot2grid((3,4), (1,3))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
            llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2],labels=[0,1,0,0], fontsize=6,color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2], color='black')
lons,lats = np.meshgrid(lon_cam,lat_cam)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,aod_cams_autumn_av,cmap=cmaps.cmocean_balance) #WhiteBlueGreenYellowRed)
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)
#=========  MERRA-winter ===========================
plt.subplot2grid((3,4), (2,0))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
            llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2], labels=[0,0,0,1], fontsize=7, color='black')
lons,lats = np.meshgrid(lon_merra,lat_merra)
xi,yi = m(lons, lats)
ax = m.pcolormesh(xi,yi,aod_mera_winter_av,cmap=cmaps.cmocean_balance) #WhiteBlueGreenYellowRed)
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)
plt.ylabel('MERRA-2-AOD', labelpad=5,fontsize=7)
#=========  MERRA-Spring ===========================
plt.subplot2grid((3,4), (2,1))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
            llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2], labels=[0,0,0,1], fontsize=7, color='black')
lons,lats = np.meshgrid(lon_merra,lat_merra)
xi,yi = m(lons, lats)
ax = m.pcolormesh(xi,yi,aod_mera_spring_av,cmap=cmaps.cmocean_balance) #WhiteBlueGreenYellowRed)
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)
#=========  MERRA-Summer ===========================
plt.subplot2grid((3,4), (2,2))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
            llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2], labels=[0,0,0,1], fontsize=7, color='black')
lons,lats = np.meshgrid(lon_merra,lat_merra)
xi,yi = m(lons, lats)
ax = m.pcolormesh(xi,yi,aod_mera_summer_av,cmap=cmaps.cmocean_balance) #WhiteBlueGreenYellowRed)
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)
#=========  MERRA-Autumn ===========================
plt.subplot2grid((3,4), (2,3))
#m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,\
            llcrnrlat=0,urcrnrlon=139,urcrnrlat=55,resolution='h')
m.drawparallels(np.arange(10, 60, 15), linewidth=0.4, dashes=[2, 2], labels=[0,1,0,0], fontsize=6,color='black')
m.drawmeridians(np.arange(45, 125,15),linewidth=0.4, dashes=[2, 2], labels=[0,0,0,1], fontsize=6, color='black')
lons,lats = np.meshgrid(lon_merra,lat_merra)
xi,yi = m(lons, lats)
ax = m.pcolormesh(xi,yi,aod_mera_autumn_av,cmap=cmaps.cmocean_balance) #WhiteBlueGreenYellowRed)
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)

#==== Setting common colorbar ==========
#cbar= fig.colorbar(ax, ax=axs[:,:], location = 'right')
#cbar.set_label('AOD', rotation=-270, fontsize=10)
fig.subplots_adjust(top=0.942,
    			bottom=0.14,
			left=0.015,
			right=.99,
			hspace=0.12,
			wspace=0.0)

cax = fig.add_axes([0.06,0.09,0.89,0.02])   # left, bottom, width, and height
cb = fig.colorbar(ax,  cax=cax, orientation='horizontal',ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
cb.set_label(label='AOD',size=7,labelpad=5)
cb.ax.tick_params(labelsize=7) 


#plt.savefig("/mnt/g/Paper_work/aod_180NEW.png",dpi=180)
plt.show()





