## load packages =======
import pandas as pd
import numpy as np
from cartopy import crs
from netCDF4 import Dataset
from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,ll_to_xy,cartopy_xlim, cartopy_ylim, ALL_TIMES)
import wrf
#================ Load control data =========================================
jan = Dataset("wrfout_d01_2017-01-01_00_00_00")
feb = Dataset("wrfout_d01_2017-02-01_00_00_00")
mar = Dataset("wrfout_d01_2017-03-01_00_00_00.nc")
apr = Dataset("wrfout_d01_2017-04-01_00_00_00.nc")
may = Dataset("wrfout_d01_2017-05-01_00_00_00.nc")
jun = Dataset("wrfout_d01_2017-06-01_00_00_00")
jul = Dataset("wrfout_d01_2017-07-01_00_00_00")
aug = Dataset("wrfout_d01_2017-08-01_00_00_00")
sep = Dataset("wrfout_d01_2017-09-01_00_00_00")
oct = Dataset("wrfout_d01_2017-10-01_00_00_00")
nov = Dataset("wrfout_d01_2017-11-01_00_00_00")
dec = Dataset("wrfout_d01_2017-12-01_00_00_00")

#====== Get BC data =======================================
jan_bc1 = getvar(jan, "BC1", timeidx=ALL_TIMES)[:,0,:,:]
jan_bc2 = getvar(jan, "BC2", timeidx=ALL_TIMES)[:,0,:,:]
feb_bc1 = getvar(feb, "BC1", timeidx=ALL_TIMES)[:,0,:,:]
feb_bc2 = getvar(feb, "BC2", timeidx=ALL_TIMES)[:,0,:,:]
mar_bc1 = getvar(mar, "BC1", timeidx=ALL_TIMES)[:,0,:,:]
mar_bc2 = getvar(mar, "BC2", timeidx=ALL_TIMES)[:,0,:,:]
apr_bc1 = getvar(apr, "BC1", timeidx=ALL_TIMES)[:,0,:,:]
apr_bc2 = getvar(apr, "BC2", timeidx=ALL_TIMES)[:,0,:,:]
may_bc1 = getvar(may, "BC1", timeidx=ALL_TIMES)[:,0,:,:]
may_bc2 = getvar(may, "BC2", timeidx=ALL_TIMES)[:,0,:,:]
jun_bc1 = getvar(jun, "BC1", timeidx=ALL_TIMES)[:,0,:,:]
jun_bc2 = getvar(jun, "BC2", timeidx=ALL_TIMES)[:,0,:,:]
jul_bc1 = getvar(jul, "BC1", timeidx=ALL_TIMES)[:,0,:,:]
jul_bc2 = getvar(jul, "BC2", timeidx=ALL_TIMES)[:,0,:,:]
aug_bc1 = getvar(aug, "BC1", timeidx=ALL_TIMES)[:,0,:,:]
aug_bc2 = getvar(aug, "BC2", timeidx=ALL_TIMES)[:,0,:,:]
sep_bc1 = getvar(sep, "BC1", timeidx=ALL_TIMES)[:,0,:,:]
sep_bc2 = getvar(sep, "BC2", timeidx=ALL_TIMES)[:,0,:,:]
oct_bc1 = getvar(oct, "BC1", timeidx=ALL_TIMES)[:,0,:,:]
oct_bc2 = getvar(oct, "BC2", timeidx=ALL_TIMES)[:,0,:,:]
nov_bc1 = getvar(nov, "BC1", timeidx=ALL_TIMES)[:,0,:,:]
nov_bc2 = getvar(nov, "BC2", timeidx=ALL_TIMES)[:,0,:,:]
dec_bc1 = getvar(dec, "BC1", timeidx=ALL_TIMES)[:,0,:,:]
dec_bc2 = getvar(dec, "BC2", timeidx=ALL_TIMES)[:,0,:,:]

#==== BC total ========
jan_bc = jan_bc1+jan_bc2
feb_bc = feb_bc1+feb_bc2
mar_bc = mar_bc1+mar_bc2
apr_bc = apr_bc1+apr_bc2
may_bc = may_bc1+may_bc2
jun_bc = jun_bc1+jun_bc2
jul_bc = jul_bc1+jul_bc2
aug_bc = aug_bc1+aug_bc2
sep_bc = sep_bc1+sep_bc2
oct_bc = oct_bc1+oct_bc2
nov_bc = nov_bc1+nov_bc2
dec_bc = dec_bc1+dec_bc2

#==== BC total ========
#bc_jan = jan_bc1+jan_bc2
#bc_feb = feb_bc1+feb_bc2
#bc_mar = mar_bc1+mar_bc2
#bc_apr = apr_bc1+apr_bc2
#bc_may = may_bc1+may_bc2
#bc_jun = jun_bc1+jun_bc2
#bc_jul = jul_bc1+jul_bc2
#bc_aug = aug_bc1+aug_bc2
#bc_sep = sep_bc1+sep_bc2
#bc_oct = oct_bc1+oct_bc2
#bc_nov = nov_bc1+nov_bc2
#bc_dec = dec_bc1+dec_bc2


"""======= Time extration ============="""
time=wrf.extract_times(jan, timeidx=ALL_TIMES, method='cat', squeeze=True, cache=None, meta=False, do_xtime=False)


#========== Get ALT data ========================================
jan_alt = getvar(jan, "ALT", timeidx=ALL_TIMES)[:,0,:,:]
feb_alt = getvar(feb, "ALT", timeidx=ALL_TIMES)[:,0,:,:]
mar_alt = getvar(mar, "ALT", timeidx=ALL_TIMES)[:,0,:,:]
apr_alt = getvar(apr, "ALT", timeidx=ALL_TIMES)[:,0,:,:]
may_alt = getvar(may, "ALT", timeidx=ALL_TIMES)[:,0,:,:]
jun_alt = getvar(jun, "ALT", timeidx=ALL_TIMES)[:,0,:,:]
#jul_alt = getvar(jun, "ALT", timeidx=ALL_TIMES)[:,0,:,:]
aug_alt = getvar(aug, "ALT", timeidx=ALL_TIMES)[:,0,:,:]
sep_alt = getvar(sep, "ALT", timeidx=ALL_TIMES)[:,0,:,:]
oct_alt = getvar(oct, "ALT", timeidx=ALL_TIMES)[:,0,:,:]
nov_alt = getvar(nov, "ALT", timeidx=ALL_TIMES)[:,0,:,:]
dec_alt = getvar(dec, "ALT", timeidx=ALL_TIMES)[:,0,:,:]

#======== Convert ug/kg to ng/m3 ===================================
bc_jan1 = jan_bc/jan_alt
bc_feb1 = feb_bc/feb_alt
bc_mar1 = mar_bc/mar_alt
bc_apr1 = apr_bc/apr_alt
bc_may1 = may_bc/may_alt
bc_jun1 = jun_bc/jun_alt
bc_jul1 = jul_bc/jun_alt
bc_aug1 = aug_bc/aug_alt
bc_sep1 = sep_bc/sep_alt
bc_oct1 = oct_bc/oct_alt
bc_nov1 = nov_bc/nov_alt
bc_dec1 = dec_bc/dec_alt

#==== Extract BC at QOMS station ============
x_y = ll_to_xy(jan, 28.36, 86.95)

bc_jan = bc_jan1[:,x_y[1], x_y[0]]*1000 #====convert ug to ng 
bc_feb = bc_feb1[:,x_y[1], x_y[0]]*1000
bc_mar = bc_mar1[:,x_y[1], x_y[0]]*1000
bc_apr = bc_apr1[:,x_y[1], x_y[0]]*1000
bc_may = bc_may1[:,x_y[1], x_y[0]]*1000
bc_jun = bc_jun1[:,x_y[1], x_y[0]]*1000
bc_jul = bc_jul1[:,x_y[1], x_y[0]]*1000
bc_aug = bc_aug1[:,x_y[1], x_y[0]]*1000
bc_sep = bc_sep1[:,x_y[1], x_y[0]]*1000
bc_oct = bc_oct1[:,x_y[1], x_y[0]]*1000
bc_nov = bc_nov1[:,x_y[1], x_y[0]]*1000  
bc_dec = bc_dec1[:,x_y[1], x_y[0]]*1000


df1  = pd.DataFrame(bc_jan)
df2  = pd.DataFrame(bc_feb)
df3  = pd.DataFrame(bc_mar)
df4  = pd.DataFrame(bc_apr)
df5  = pd.DataFrame(bc_may)
df6  = pd.DataFrame(bc_jun)
df7  = pd.DataFrame(bc_jul)
df8  = pd.DataFrame(bc_aug)
df9  = pd.DataFrame(bc_sep)
df10 = pd.DataFrame(bc_oct)
df11 = pd.DataFrame(bc_nov)
df12 = pd.DataFrame(bc_dec)

#==== Now saving the data======
df1.to_csv('/mnt/g/Paper_work/BC_extract/qoms/jan.csv')
df2.to_csv('/mnt/g/Paper_work/BC_extract/qoms/feb.csv')
df3.to_csv('/mnt/g/Paper_work/BC_extract/qoms/mar.csv')
df4.to_csv('/mnt/g/Paper_work/BC_extract/qoms/apr.csv')
df5.to_csv('/mnt/g/Paper_work/BC_extract/qoms/may.csv')
df6.to_csv('/mnt/g/Paper_work/BC_extract/qoms/jun.csv')
df7.to_csv('/mnt/g/Paper_work/BC_extract/qoms/jul.csv')
df8.to_csv('/mnt/g/Paper_work/BC_extract/qoms/aug.csv')
df9.to_csv('/mnt/g/Paper_work/BC_extract/qoms/sep.csv')
df10.to_csv('/mnt/g/Paper_work/BC_extract/qoms/oct.csv')
df11.to_csv('/mnt/g/Paper_work/BC_extract/qoms/nov.csv')
df12.to_csv('/mnt/g/Paper_work/BC_extract/qoms/dec.csv')


