

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

x_y = ll_to_xy(wrflist, 28.36, 86.95) # Change here !! Lat, Lon ##


bc1 = getvar(wrflist, "BC1", timeidx=ALL_TIMES,method='cat')
bc2 = getvar(wrflist, "BC2", timeidx=ALL_TIMES,method='cat')
bc = bc1+bc2

#===== Change into ng/m3==========
#alt = getvar(wrflist, "ALT", timeidx=ALL_TIMES,method='cat')
#bc4 = bc3/alt

bc0 = bc[:,0,x_y[1], x_y[0]]
bc3 = bc[:,3,x_y[1], x_y[0]]
bc4 = bc[:,4,x_y[1], x_y[0]]
bc5 = bc[:,5,x_y[1], x_y[0]]
bc6 = bc[:,6,x_y[1], x_y[0]]
bc7 = bc[:,7,x_y[1], x_y[0]]
bc8 = bc[:,8,x_y[1], x_y[0]]
bc9 = bc[:,9,x_y[1], x_y[0]]
bc10 = bc[:,10,x_y[1], x_y[0]]

df0  = pd.DataFrame(bc0)
df3  = pd.DataFrame(bc3)
df4  = pd.DataFrame(bc4)
df5  = pd.DataFrame(bc5)
df6  = pd.DataFrame(bc6)
df7  = pd.DataFrame(bc7)
df8  = pd.DataFrame(bc8)
df9  = pd.DataFrame(bc9)
df10  = pd.DataFrame(bc10)

df0.to_csv('/mnt/g/Paper_work/BC_extract/qoms/bc_qoms_level0.csv')
#df3.to_csv('/mnt/g/Paper_work/BC_extract/qoms/bc_qoms_level3.csv')
#df4.to_csv('/mnt/g/Paper_work/BC_extract/qoms/bc_qoms_level4.csv')
#df5.to_csv('/mnt/g/Paper_work/BC_extract/qoms/bc_qoms_level5.csv')
#df6.to_csv('/mnt/g/Paper_work/BC_extract/qoms/bc_qoms_level6.csv')
#df7.to_csv('/mnt/g/Paper_work/BC_extract/qoms/bc_qoms_level7.csv')
#df8.to_csv('/mnt/g/Paper_work/BC_extract/qoms/bc_qoms_level8.csv')
#df9.to_csv('/mnt/g/Paper_work/BC_extract/qoms/bc_qoms_level9.csv')
#df10.to_csv('/mnt/g/Paper_work/BC_extract/qoms/bc_qoms_level10.csv')

