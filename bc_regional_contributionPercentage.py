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
from collections import OrderedDict
import skill_metrics as sm

#====== For Control Case ===========================
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


#===== TP =================

#===== Slice Tibetan Plataeu regiona [lat= 30:40,lon=80:100] ===
bc1_tp = wrf.getvar(nc1_c, 'BC1',timeidx=1)[0,60:93,92:160]  ## 0 means at surface
bc2_tp = wrf.getvar(nc1_c, 'BC2',timeidx=1)[0,60:93,92:160]
bc3_tp = wrf.getvar(nc2_c, 'BC1',timeidx=1)[0,60:93,92:160]
bc4_tp = wrf.getvar(nc2_c, 'BC2',timeidx=1)[0,60:93,92:160]
bc5_tp = wrf.getvar(nc3_c, 'BC1',timeidx=1)[0,60:93,92:160]
bc6_tp = wrf.getvar(nc3_c, 'BC2',timeidx=1)[0,60:93,92:160]
bc7_tp = wrf.getvar(nc4_c, 'BC1',timeidx=1)[0,60:93,92:160]
bc8_tp = wrf.getvar(nc4_c, 'BC2',timeidx=1)[0,60:93,92:160]
bc9_tp= wrf.getvar(nc5_c, 'BC1',timeidx=1)[0,60:93,92:160]
bc10_tp = wrf.getvar(nc5_c, 'BC2',timeidx=1)[0,60:93,92:160]
bc11_tp = wrf.getvar(nc6_c, 'BC1',timeidx=1)[0,60:93,92:160]
bc12_tp = wrf.getvar(nc6_c, 'BC2',timeidx=1)[0,60:93,92:160]
bc13_tp = wrf.getvar(nc7_c, 'BC1',timeidx=1)[0,60:93,92:160]
bc14_tp = wrf.getvar(nc7_c, 'BC2',timeidx=1)[0,60:93,92:160]
bc15_tp = wrf.getvar(nc8_c, 'BC1',timeidx=1)[0,60:93,92:160]
bc16_tp = wrf.getvar(nc8_c, 'BC2',timeidx=1)[0,60:93,92:160]
bc17_tp = wrf.getvar(nc9_c, 'BC1',timeidx=1)[0,60:93,92:160]
bc18_tp = wrf.getvar(nc9_c, 'BC2',timeidx=1)[0,60:93,92:160]
bc19_tp = wrf.getvar(nc10_c, 'BC1',timeidx=1)[0,60:93,92:160]
bc20_tp = wrf.getvar(nc10_c, 'BC2',timeidx=1)[0,60:93,92:160]
bc21_tp = wrf.getvar(nc11_c, 'BC1',timeidx=1)[0,60:93,92:160]
bc22_tp = wrf.getvar(nc11_c, 'BC2',timeidx=1)[0,60:93,92:160]
bc23_tp = wrf.getvar(nc12_c, 'BC1',timeidx=1)[0,60:93,92:160]
bc24_tp = wrf.getvar(nc12_c, 'BC2',timeidx=1)[0,60:93,92:160]

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
bc1_ec = wrf.getvar(nc1_ec, 'BC1',timeidx=1)[0]
bc2_ec = wrf.getvar(nc1_ec, 'BC2',timeidx=1)[0]
bc3_ec = wrf.getvar(nc2_ec, 'BC1',timeidx=1)[0]
bc4_ec = wrf.getvar(nc2_ec, 'BC2',timeidx=1)[0]
bc5_ec = wrf.getvar(nc3_ec, 'BC1',timeidx=1)[0]
bc6_ec = wrf.getvar(nc3_ec, 'BC2',timeidx=1)[0]
bc7_ec = wrf.getvar(nc4_ec, 'BC1',timeidx=1)[0]
bc8_ec = wrf.getvar(nc4_ec, 'BC2',timeidx=1)[0]
bc9_ec = wrf.getvar(nc5_ec, 'BC1',timeidx=1)[0]
bc10_ec = wrf.getvar(nc5_ec, 'BC2',timeidx=1)[0]
bc11_ec = wrf.getvar(nc6_ec, 'BC1',timeidx=1)[0]
bc12_ec = wrf.getvar(nc6_ec, 'BC2',timeidx=1)[0]
bc13_ec = wrf.getvar(nc7_ec, 'BC1',timeidx=1)[0]
bc14_ec = wrf.getvar(nc7_ec, 'BC2',timeidx=1)[0]
bc15_ec = wrf.getvar(nc8_ec, 'BC1',timeidx=1)[0]
bc16_ec = wrf.getvar(nc8_ec, 'BC2',timeidx=1)[0]
bc17_ec = wrf.getvar(nc9_ec, 'BC1',timeidx=1)[0]
bc18_ec = wrf.getvar(nc9_ec, 'BC2',timeidx=1)[0]
bc19_ec = wrf.getvar(nc10_ec, 'BC1',timeidx=1)[0]
bc20_ec = wrf.getvar(nc10_ec, 'BC2',timeidx=1)[0]
bc21_ec = wrf.getvar(nc11_ec, 'BC1',timeidx=1)[0]
bc22_ec = wrf.getvar(nc11_ec, 'BC2',timeidx=1)[0]
bc23_ec = wrf.getvar(nc12_ec, 'BC1',timeidx=1)[0]
bc24_ec = wrf.getvar(nc12_ec, 'BC2',timeidx=1)[0]

#======= Winter =========
winter_tp  =  (bc1_tp + bc2_tp + bc3_tp + bc4_tp + bc23_tp + bc24_tp)/3
winter_c  =  (bc1_c + bc2_c + bc3_c + bc4_c + bc23_c + bc24_c)/3
winter_nc1 = (bc1_nc + bc2_nc + bc3_nc + bc4_nc + bc23_nc + bc24_nc)/3
winter_ca1 = (bc1_ca + bc2_ca + bc3_ca + bc4_ca + bc23_ca + bc24_ca)/3
winter_sa1 = (bc1_sa + bc2_sa + bc3_sa + bc4_sa + bc23_sa + bc24_sa)/3
winter_ec1 = (bc1_ec + bc2_ec + bc3_ec + bc4_ec + bc23_ec + bc24_ec)/3
##==== Spring ===========
spring_tp =   (bc5_tp + bc6_tp + bc7_tp + bc8_tp + bc9_tp + bc10_tp)/3
spring_c =   (bc5_c + bc6_c + bc7_c + bc8_c + bc9_c + bc10_c)/3
spring_nc1 = (bc5_nc + bc6_nc + bc7_nc + bc8_nc + bc9_nc + bc10_nc)/3
spring_ca1 = (bc5_ca + bc6_ca + bc7_ca + bc8_ca + bc9_ca + bc10_ca)/3
spring_sa1 = (bc5_sa + bc6_sa + bc7_sa + bc8_sa + bc9_sa + bc10_sa)/3
spring_ec1 = (bc5_ec + bc6_ec + bc7_ec + bc8_ec + bc9_ec + bc10_ec)/3

##===== Summer===============
summer_tp =   (bc11_tp + bc12_tp + bc13_tp + bc14_tp + bc15_tp + bc16_tp)/3
summer_c =   (bc11_c + bc12_c + bc13_c + bc14_c + bc15_c + bc16_c)/3
summer_nc1 = (bc11_nc + bc12_nc + bc13_nc + bc14_nc + bc15_nc + bc16_nc)/3
summer_ca1 = (bc11_ca + bc12_ca + bc13_ca + bc14_ca + bc15_ca + bc16_ca)/3
summer_sa1 = (bc11_sa + bc12_sa + bc13_sa + bc14_sa + bc15_sa + bc16_sa)/3
summer_ec1 = (bc11_ec + bc12_ec + bc13_ec + bc14_ec + bc15_ec + bc16_ec)/3

##====== Autumn =======================
autumn_tp =  (bc17_tp + bc18_tp + bc19_tp + bc20_tp + bc21_tp + bc22_tp)/3
autumn_c =  (bc17_c + bc18_c + bc19_c + bc20_c + bc21_c + bc22_c)/3
autumn_nc1 = (bc17_nc + bc18_nc + bc19_nc + bc20_nc + bc21_nc + bc22_nc)/3
autumn_ca1 = (bc17_ca + bc18_ca + bc19_ca + bc20_ca + bc21_ca + bc22_ca)/3
autumn_sa1 = (bc17_sa + bc18_sa + bc19_sa + bc20_sa + bc21_sa + bc22_sa)/3
autumn_ec1 = (bc17_ec + bc18_ec + bc19_ec + bc20_ec + bc21_ec + bc22_ec)/3

###========= For each region ===========
winter_nc = (winter_c-winter_nc1)
winter_ca = (winter_c-winter_ca1)
winter_sa = (winter_c-winter_sa1)
winter_ec = (winter_c-winter_ec1)
spring_nc = (spring_c - spring_nc1)
spring_ca = (spring_c - spring_ca1)
spring_sa = (spring_c - spring_sa1)
spring_ec = (spring_c - spring_ec1)
summer_nc = (summer_c - summer_nc1)
summer_ca = (summer_c - summer_ca1)
summer_sa = (summer_c - summer_sa1)
summer_ec = (summer_c - summer_ec1)
autumn_nc = (autumn_c - autumn_nc1)
autumn_ca = (autumn_c - autumn_ca1)
autumn_sa = (autumn_c - autumn_sa1)
autumn_ec = (autumn_c - autumn_ec1)

###======
stats = OrderedDict()
##===== Now take zonal average ===================
##=== For TP ==========
win_tp_av = np.average(winter_tp,axis=1)
spr_tp_av = np.average(spring_tp,axis=1)
sum_tp_av= np.average(summer_tp,axis=1)
aut_tp_av = np.average(autumn_tp,axis=1)

##=== For Winter Season =====
#== For Winter==
win_nc_av      = np.average(winter_nc,axis=1)
win_ca_av      = np.average(winter_ca,axis=1)
win_sa_av      = np.average(winter_sa,axis=1)
win_ec_av      = np.average(winter_ec,axis=1)
#== For Spring==
spr_nc_av      = np.average(spring_nc,axis=1)
spr_ca_av      = np.average(spring_ca,axis=1)
spr_sa_av      = np.average(spring_sa,axis=1)
spr_ec_av      = np.average(spring_ec,axis=1)
#== For Summer==
sum_nc_av      = np.average(summer_nc,axis=1)
sum_ca_av      = np.average(summer_ca,axis=1)
sum_sa_av      = np.average(summer_sa,axis=1)
sum_ec_av      = np.average(summer_ec,axis=1)
#== For autumn==
aut_nc_av      = np.average(autumn_nc,axis=1)
aut_ca_av      = np.average(autumn_ca,axis=1)
aut_sa_av      = np.average(autumn_sa,axis=1)
aut_ec_av      = np.average(autumn_ec,axis=1)
#======  Now take zonal mean from each grid cell====
#== For control===
win_tp = np.average(win_tp_av)
spr_tp = np.average(spr_tp_av)
sum_tp = np.average(sum_tp_av)
aut_tp = np.average(aut_tp_av)
#==winter
nc_bc_win      = np.average(win_nc_av)
ca_bc_win      = np.average(win_ca_av)
sa_bc_win      = np.average(win_sa_av)
ec_bc_win      = np.average(win_ec_av)
#== spring
nc_bc_spr      = np.average(spr_nc_av)
ca_bc_spr      = np.average(spr_ca_av)
sa_bc_spr      = np.average(spr_sa_av)
ec_bc_spr      = np.average(spr_ec_av)
#== summer
nc_bc_sum      = np.average(sum_nc_av)
ca_bc_sum      = np.average(sum_ca_av)
sa_bc_sum      = np.average(sum_sa_av)
ec_bc_sum      = np.average(sum_ec_av)
#=== autumn
nc_bc_aut      = np.average(aut_nc_av)
ca_bc_aut      = np.average(aut_ca_av)
sa_bc_aut      = np.average(aut_sa_av)
ec_bc_aut      = np.average(aut_ec_av)


#=== Now calculate % contribution ======
stats['nc_win'] = (nc_bc_win-win_tp)/nc_bc_win
stats['ca_win'] = (ca_bc_win-win_tp)/ca_bc_win
stats['sa_win'] = (sa_bc_win-win_tp)/sa_bc_win
stats['ec_win'] = (ec_bc_win-win_tp)/ec_bc_win

stats['nc_spr'] = (nc_bc_spr-spr_tp)/nc_bc_spr
stats['ca_spr'] = (ca_bc_spr-spr_tp)/ca_bc_spr
stats['sa_spr'] = (sa_bc_spr-spr_tp)/sa_bc_spr
stats['ec_spr'] = (ec_bc_spr-spr_tp)/ec_bc_spr

stats['nc_sum'] = (nc_bc_sum-sum_tp)/nc_bc_sum
stats['ca_sum'] = (ca_bc_sum-sum_tp)/ca_bc_sum
stats['sa_sum'] = (sa_bc_sum-sum_tp)/sa_bc_sum
stats['ec_sum'] = (ec_bc_sum-sum_tp)/ec_bc_sum

#stats['nc_aut'] = (nc_bc_aut-aut_tp)/nc_bc_aut*100
#stats['ca_aut'] = (ca_bc_aut-aut_tp)/ca_bc_aut*100
#stats['sa_aut'] = (sa_bc_aut-aut_tp)/sa_bc_aut*100
#stats['ec_aut'] = (ec_bc_aut-aut_tp)/ec_bc_aut*100

stats['nc_aut'] = (aut_tp)/nc_bc_aut*100
stats['ca_aut'] = (aut_tp)/ca_bc_aut*100
stats['sa_aut'] = (aut_tp)/sa_bc_aut*100
stats['ec_aut'] = (aut_tp)/ec_bc_aut*100
#===saving===
filename = '/mnt/g/Paper_work/contributionRate_TP14.xlsx'
sm.write_stats(filename,stats,overwrite=True)

