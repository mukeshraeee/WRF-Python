from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np
import wrf
from wrf import (to_np, getvar,interplevel, smooth2d, get_basemap, latlon_coords, ALL_TIMES,ll_to_xy)
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
import cmaps
from pylab import *
import scipy as sp
import scipy.ndimage
import matplotlib.gridspec as gridspec
from collections import OrderedDict
import skill_metrics as sm

import pandas as pd
import numpy as np
from cartopy import crs

#==== Load data =======================================
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

#==== Get BC1 and BC2 data ======= Control case ==================================
win_bc1 = getvar(win, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:] # all time,surface,all lat and lon
win_bc2 = getvar(win, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
spr_bc1 = getvar(spr, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
spr_bc2 = getvar(spr, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
mon_bc1 = getvar(mon, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
mon_bc2 = getvar(mon, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
aut_bc1 = getvar(aut, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
aut_bc2 = getvar(aut, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]

win_bc = (win_bc1+win_bc2)
spr_bc = (spr_bc1+spr_bc2)
mon_bc = (mon_bc1+mon_bc2)
aut_bc = (aut_bc1+aut_bc2)

#==== Take time average ====
win_av = np.average(win_bc,axis=0)
spr_av = np.average(spr_bc,axis=0)
mon_av = np.average(mon_bc,axis=0)
aut_av = np.average(aut_bc,axis=0)
#====== Define TP===============================================================
#lat_tp = np.array([35, 40, 40, 37, 37, 36, 36, 37, 37, 38, 38, 40, 40,  38,  38,  37,  37,  32,  32,  25,  25, 28,  28, 29,29,30,30,31,31,32,32,35,35])
#lon_tp = np.array([70, 70, 75, 75, 80, 80, 84, 84, 86, 86, 90, 90 ,100, 100, 102, 102, 105, 105, 102, 102, 100,100, 86, 86,83,83,81,81,79,79,75,75,70])

#======= TP=================================

tp_win_bc1 = getvar(win, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,60:90,96:155] # all time,surface,all lat and lon
tp_win_bc2 = getvar(win, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,60:90,96:155]
tp_spr_bc1 = getvar(spr, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,60:90,96:155]
tp_spr_bc2 = getvar(spr, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,60:90,96:155]
tp_mon_bc1 = getvar(mon, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,60:90,96:155]
tp_mon_bc2 = getvar(mon, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,60:90,96:155]
tp_aut_bc1 = getvar(aut, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,60:90,96:155]
tp_aut_bc2 = getvar(aut, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,60:90,96:155]

#tp_win_bc = win_bc[lat_tp,lon_tp]
#tp_spr_bc = spr_av[lat_tp,lon_tp] 
#tp_mon_bc = mon_av[lat_tp,lon_tp]
#tp_aut_bc = aut_av[lat_tp,lon_tp]


tp_win_bc = (tp_win_bc1+tp_win_bc2)
tp_spr_bc = (tp_spr_bc1+tp_spr_bc2)
tp_mon_bc = (tp_mon_bc1+tp_mon_bc2)
tp_aut_bc = (tp_aut_bc1+tp_aut_bc2)

#==== Take time average ====
tp_win_av = np.average(tp_win_bc,axis=0)
tp_spr_av = np.average(tp_spr_bc,axis=0)
tp_mon_av = np.average(tp_mon_bc,axis=0)
tp_aut_av = np.average(tp_aut_bc,axis=0)

#==== Load SA data================================
sa_win = [Dataset("wrfout_d01_2017-01-01_00_00_00-SA"),
          Dataset("wrfout_d01_2017-02-01_00_00_00-SA"),
          Dataset("wrfout_d01_2017-12-01_00_00_00-SA")]
sa_spr = [Dataset("wrfout_d01_2017-03-01_00_00_00-SA"),
          Dataset("wrfout_d01_2017-04-01_00_00_00-SA"),
          Dataset("wrfout_d01_2017-05-01_00_00_00-SA")]
sa_mon = [Dataset("wrfout_d01_2017-06-01_00_00_00-SA"),
          Dataset("wrfout_d01_2017-07-01_00_00_00-SA"),
          Dataset("wrfout_d01_2017-08-01_00_00_00-SA")]
sa_aut = [Dataset("wrfout_d01_2017-09-01_00_00_00-SA"),
          Dataset("wrfout_d01_2017-10-01_00_00_00-SA"),
          Dataset("wrfout_d01_2017-11-01_00_00_00-SA")]

sa_win_bc1 = getvar(sa_win, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:] # all time,surface,all lat and lon
sa_win_bc2 = getvar(sa_win, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
sa_spr_bc1 = getvar(sa_spr, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
sa_spr_bc2 = getvar(sa_spr, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
sa_mon_bc1 = getvar(sa_mon, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
sa_mon_bc2 = getvar(sa_mon, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
sa_aut_bc1 = getvar(sa_aut, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
sa_aut_bc2 = getvar(sa_aut, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]

sa_win_bc = (sa_win_bc1+sa_win_bc2)
sa_spr_bc = (sa_spr_bc1+sa_spr_bc2)
sa_mon_bc = (sa_mon_bc1+sa_mon_bc2)
sa_aut_bc = (sa_aut_bc1+sa_aut_bc2)

#==== Take time average ====
sa_win_av = np.average(sa_win_bc,axis=0)
sa_spr_av = np.average(sa_spr_bc,axis=0)
sa_mon_av = np.average(sa_mon_bc,axis=0)
sa_aut_av = np.average(sa_aut_bc,axis=0)

#==== Load EC data================================
ec_win = [Dataset("wrfout_d01_2017-01-01_00_00_00-EC"),
          Dataset("wrfout_d01_2017-02-01_00_00_00-EC"),
          Dataset("wrfout_d01_2017-12-01_00_00_00-EC")]
ec_spr = [Dataset("wrfout_d01_2017-03-01_00_00_00-EC"),
          Dataset("wrfout_d01_2017-04-01_00_00_00-EC"),
          Dataset("wrfout_d01_2017-05-01_00_00_00-EC")]
ec_mon = [Dataset("wrfout_d01_2017-06-01_00_00_00-EC"),
          Dataset("wrfout_d01_2017-07-01_00_00_00-EC"),
          Dataset("wrfout_d01_2017-08-01_00_00_00-EC")]
ec_aut = [Dataset("wrfout_d01_2017-09-01_00_00_00-EC"),
          Dataset("wrfout_d01_2017-10-01_00_00_00-EC"),
          Dataset("wrfout_d01_2017-11-01_00_00_00-EC")]

ec_win_bc1 = getvar(ec_win, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:] # all time,surface,all lat and lon
ec_win_bc2 = getvar(ec_win, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
ec_spr_bc1 = getvar(ec_spr, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
ec_spr_bc2 = getvar(ec_spr, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
ec_mon_bc1 = getvar(ec_mon, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
ec_mon_bc2 = getvar(ec_mon, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
ec_aut_bc1 = getvar(ec_aut, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
ec_aut_bc2 = getvar(ec_aut, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]


ec_win_bc = (ec_win_bc1+ec_win_bc2)
ec_spr_bc = (ec_spr_bc1+ec_spr_bc2)
ec_mon_bc = (ec_mon_bc1+ec_mon_bc2)
ec_aut_bc = (ec_aut_bc1+ec_aut_bc2)

#==== Take time average ====
ec_win_av = np.average(ec_win_bc,axis=0)
ec_spr_av = np.average(ec_spr_bc,axis=0)
ec_mon_av = np.average(ec_mon_bc,axis=0)
ec_aut_av = np.average(ec_aut_bc,axis=0)

#==== Load NC data================================
nc_win = [Dataset("wrfout_d01_2017-01-01_00_00_00-NC"),
          Dataset("wrfout_d01_2017-02-01_00_00_00-NC"),
          Dataset("wrfout_d01_2017-12-01_00_00_00-NC")]
nc_spr = [Dataset("wrfout_d01_2017-03-01_00_00_00-NC"),
          Dataset("wrfout_d01_2017-04-01_00_00_00-NC"),
          Dataset("wrfout_d01_2017-05-01_00_00_00-NC")]
nc_mon = [Dataset("wrfout_d01_2017-06-01_00_00_00-NC"),
          Dataset("wrfout_d01_2017-07-01_00_00_00-NC"),
          Dataset("wrfout_d01_2017-08-01_00_00_00-NC")]
nc_aut = [Dataset("wrfout_d01_2017-09-01_00_00_00-NC"),
          Dataset("wrfout_d01_2017-10-01_00_00_00-NC"),
          Dataset("wrfout_d01_2017-11-01_00_00_00-NC")]

nc_win_bc1 = getvar(nc_win, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:] # all time,surface,all lat and lon
nc_win_bc2 = getvar(nc_win, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
nc_spr_bc1 = getvar(nc_spr, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
nc_spr_bc2 = getvar(nc_spr, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
nc_mon_bc1 = getvar(nc_mon, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
nc_mon_bc2 = getvar(nc_mon, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
nc_aut_bc1 = getvar(nc_aut, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
nc_aut_bc2 = getvar(nc_aut, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]


nc_win_bc = (nc_win_bc1+nc_win_bc2)
nc_spr_bc = (nc_spr_bc1+nc_spr_bc2)
nc_mon_bc = (nc_mon_bc1+nc_mon_bc2)
nc_aut_bc = (nc_aut_bc1+nc_aut_bc2)

#==== Take time average ====
nc_win_av = np.average(nc_win_bc,axis=0)
nc_spr_av = np.average(nc_spr_bc,axis=0)
nc_mon_av = np.average(nc_mon_bc,axis=0)
nc_aut_av = np.average(nc_aut_bc,axis=0)

#==== Load CA data================================
ca_win = [Dataset("wrfout_d01_2017-01-01_00_00_00-CA"),
          Dataset("wrfout_d01_2017-02-01_00_00_00-CA"),
          Dataset("wrfout_d01_2017-12-01_00_00_00-CA")]
ca_spr = [Dataset("wrfout_d01_2017-03-01_00_00_00-CA"),
          Dataset("wrfout_d01_2017-04-01_00_00_00-CA"),
          Dataset("wrfout_d01_2017-05-01_00_00_00-CA")]
ca_mon = [Dataset("wrfout_d01_2017-06-01_00_00_00-CA"),
          Dataset("wrfout_d01_2017-07-01_00_00_00-CA"),
          Dataset("wrfout_d01_2017-08-01_00_00_00-CA")]
ca_aut = [Dataset("wrfout_d01_2017-09-01_00_00_00-CA"),
          Dataset("wrfout_d01_2017-10-01_00_00_00-CA"),
          Dataset("wrfout_d01_2017-11-01_00_00_00-CA")]

ca_win_bc1 = getvar(ca_win, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:] # all time,surface,all lat and lon
ca_win_bc2 = getvar(ca_win, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
ca_spr_bc1 = getvar(ca_spr, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
ca_spr_bc2 = getvar(ca_spr, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
ca_mon_bc1 = getvar(ca_mon, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
ca_mon_bc2 = getvar(ca_mon, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
ca_aut_bc1 = getvar(ca_aut, "BC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
ca_aut_bc2 = getvar(ca_aut, "BC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]

ca_win_bc = (ca_win_bc1+ca_win_bc2)
ca_spr_bc = (ca_spr_bc1+ca_spr_bc2)
ca_mon_bc = (ca_mon_bc1+ca_mon_bc2)
ca_aut_bc = (ca_aut_bc1+ca_aut_bc2)

#==== Take time average ====
ca_win_av = np.average(ca_win_bc,axis=0)
ca_spr_av = np.average(ca_spr_bc,axis=0)
ca_mon_av = np.average(ca_mon_bc,axis=0)
ca_aut_av = np.average(ca_aut_bc,axis=0)

#==== Now take whole area average ==========================

#==== Whole region - Control case  ====
win_ave = np.average(win_av)
spr_ave = np.average(spr_av)
mon_ave = np.average(mon_av)
aut_ave = np.average(aut_av)

#==== TP region ===============================
tp_win_ave = np.average(tp_win_av)
tp_spr_ave = np.average(tp_spr_av)
tp_mon_ave = np.average(tp_mon_av)
tp_aut_ave = np.average(tp_aut_av)

#==== SA region=====================
sa_win_ave = np.average(sa_win_av)
sa_spr_ave = np.average(sa_spr_av)
sa_mon_ave = np.average(sa_mon_av)
sa_aut_ave = np.average(sa_aut_av)

sa_winter  = win_ave-sa_win_ave
sa_spring  = spr_ave-sa_spr_ave
sa_monsoon = mon_ave-sa_mon_ave
sa_autumn  = aut_ave-tp_aut_ave 

#==== EC region================= 
ec_win_ave = np.average(ec_win_av)
ec_spr_ave = np.average(ec_spr_av)
ec_mon_ave = np.average(ec_mon_av)
ec_aut_ave = np.average(ec_aut_av)

ec_winter  = win_ave-ec_win_ave
ec_spring  = spr_ave-ec_spr_ave
ec_monsoon = mon_ave-ec_mon_ave
ec_autumn  = aut_ave-ec_aut_ave

#==== NC region=================
nc_win_ave = np.average(nc_win_av)
nc_spr_ave = np.average(nc_spr_av)
nc_mon_ave = np.average(nc_mon_av)
nc_aut_ave = np.average(nc_aut_av)

nc_winter  = win_ave-nc_win_ave
nc_spring  = spr_ave-nc_spr_ave
nc_monsoon = mon_ave-nc_mon_ave
nc_autumn  = aut_ave-nc_aut_ave

#==== CA region ====================
ca_win_ave = np.average(ca_win_av)
ca_spr_ave = np.average(ca_spr_av)
ca_mon_ave = np.average(ca_mon_av)
ca_aut_ave = np.average(ca_aut_av)

ca_winter  = win_ave-ca_win_ave
ca_spring  = spr_ave-ca_spr_ave
ca_monsoon = mon_ave-ca_mon_ave
ca_autumn  = aut_ave-ca_aut_ave


#=======  Calculate BC Contribution in %  ==========================
stats = OrderedDict()

stats['win_sa'] = (tp_win_ave)/(sa_winter)*100
stats['win_ec'] = (tp_win_ave)/(ec_winter)*100
stats['win_nc'] = (tp_win_ave)/(nc_winter)*100
stats['win_ca'] = (tp_win_ave)/(ca_winter)*100

stats['spr_sa'] = (tp_spr_ave)/(sa_spring)*100
stats['spr_ec'] = (tp_spr_ave)/(ec_spring)*100
stats['spr_nc'] = (tp_spr_ave)/(nc_spring)*100
stats['spr_ca'] = (tp_spr_ave)/(ca_spring)*100

stats['mon_sa'] = (tp_mon_ave)/(sa_monsoon)*100
stats['mon_ec'] = (tp_mon_ave)/(ec_monsoon)*100
stats['mon_nc'] = (tp_mon_ave)/(nc_monsoon)*100
stats['mon_ca'] = (tp_mon_ave)/(ca_monsoon)*100

stats['aut_sa'] = (tp_aut_ave)/(sa_autumn)*100
stats['aut_ec'] = (tp_aut_ave)/(ec_autumn)*100
stats['aut_nc'] = (tp_aut_ave)/(nc_autumn)*100
stats['aut_ca'] = (tp_aut_ave)/(ca_autumn)*100


#===saving===
filename = '/mnt/g/Paper_work/Contribution_BC_TP_newnew.xlsx'
sm.write_stats(filename,stats,overwrite=True)












