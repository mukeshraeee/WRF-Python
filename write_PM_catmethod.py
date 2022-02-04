"""
2020/06/08 == Prepared by RAI. Mukesh
This piece of code helps to extract the
AOD  by integrating the extiction data vertically
in required time period
"""
## load packages =======
import pandas as pd
import numpy as np
from cartopy import crs
from netCDF4 import Dataset
from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,ll_to_xy,cartopy_xlim, cartopy_ylim, ALL_TIMES)
#========= Load WRF data ================================
wrflist = [Dataset("wrfout_d01_2017-01-01_00_00_00"),
           Dataset("wrfout_d01_2017-02-01_00_00_00"),
           Dataset("wrfout_d01_2017-03-01_00_00_00.nc"),
           Dataset("wrfout_d01_2017-04-01_00_00_00.nc"),
           Dataset("wrfout_d01_2017-05-01_00_00_00.nc"),
           Dataset("wrfout_d01_2017-06-01_00_00_00"),
           Dataset("wrfout_d01_2017-07-01_00_00_00"),
           Dataset("wrfout_d01_2017-08-01_00_00_00"),
           Dataset("wrfout_d01_2017-09-01_00_00_00"),
           Dataset("wrfout_d01_2017-10-01_00_00_00"),
           Dataset("wrfout_d01_2017-11-01_00_00_00"),
           Dataset("wrfout_d01_2017-12-01_00_00_00")]



x_y = ll_to_xy(wrflist, 39.8673, 116.366)

#======== get extinction data ==========================
pm25 = getvar(wrflist, "PM2_5_DRY", timeidx=ALL_TIMES,method='cat')
pm10 = getvar(wrflist, "PM10", timeidx=ALL_TIMES,method='cat')

pm_25 = pm25[:,0,x_y[1], x_y[0]]
pm_10 = pm10[:,0,x_y[1], x_y[0]]
df1  = pd.DataFrame(pm_25)
df2  = pd.DataFrame(pm_10)

df1.to_csv('pm2.5_bei.csv')
df2.to_csv('pm10_bei.csv')
