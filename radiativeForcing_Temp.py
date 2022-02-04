#====== Code for plotting Radiatie forcing by anthropogenic BC ===========================
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
#========= Load WRF data ================================
jan_c   = Dataset("wrfout_d01_2017-01-01_00_00_00") # Control
jan_s   = Dataset("wrfout_d01_2017-01-01_nobc")   # Sensitivity
feb_c   = Dataset("wrfout_d01_2017-02-01_00_00_00")
feb_s   = Dataset("wrfout_d01_2017-02-01_nobc")
mar_c   = Dataset("wrfout_d01_2017-03-01_00_00_00.nc")
mar_s   = Dataset("wrfout_d01_2017-03-01_nobc")
apr_c   = Dataset("wrfout_d01_2017-04-01_00_00_00.nc")
apr_s   = Dataset("wrfout_d01_2017-04-01_nobc")
may_c   = Dataset("wrfout_d01_2017-05-01_00_00_00.nc")
may_s   = Dataset("wrfout_d01_2017-05-01_nobc")
jun_c   = Dataset("wrfout_d01_2017-06-01_00_00_00")
jun_s   = Dataset("wrfout_d01_2017-06-01_nobc")
jul_c   = Dataset("wrfout_d01_2017-07-01_00_00_00")
jul_s   = Dataset("wrfout_d01_2017-07-01_nobc")
aug_c   = Dataset("wrfout_d01_2017-08-01_00_00_00")
aug_s   = Dataset("wrfout_d01_2017-08-01_nobc")
sep_c   = Dataset("wrfout_d01_2017-09-01_00_00_00")
sep_s   = Dataset("wrfout_d01_2017-09-01_nobc")
oct_c   = Dataset("wrfout_d01_2017-10-01_00_00_00")
oct_s   = Dataset("wrfout_d01_2017-10-01_nobc")
nov_c   = Dataset("wrfout_d01_2017-11-01_00_00_00")
nov_s   = Dataset("wrfout_d01_2017-11-01_nobc")
dec_c   = Dataset("wrfout_d01_2017-12-01_00_00_00")
dec_s   = Dataset("wrfout_d01_2017-12-01_nobc")

#========= Get Radiative Flux data ======================
#====== January ============
jan_swdwn_c  = getvar(jan_c,   "SWDOWN", timeidx=2)
jan_swdwn_s  = getvar(jan_s,   "SWDOWN", timeidx=2)
jan_glw_c    = getvar(jan_c,   "GLW",    timeidx=2)
jan_glw_s    = getvar(jan_s,   "GLW",    timeidx=2)
jan_alb_c    = getvar(jan_c,   "ALBEDO", timeidx=2)
jan_alb_s    = getvar(jan_s,   "ALBEDO", timeidx=2)
jan_emis_c   = getvar(jan_c,   "EMISS",  timeidx=2)
jan_emis_s   = getvar(jan_s,   "EMISS",  timeidx=2)
jan_tsk_c    = getvar(jan_c,   "TSK",    timeidx=2)
jan_tsk_s    = getvar(jan_s,   "TSK",    timeidx=2)
jan_grdflx_c = getvar(jan_c,   "GRDFLX", timeidx=2)
jan_grdflx_s = getvar(jan_s,   "GRDFLX", timeidx=2)
jan_lh_c     = getvar(jan_c,   "LH",     timeidx=2)
jan_lh_s     = getvar(jan_s,   "LH",     timeidx=2)
jan_sh_c     = getvar(jan_c,   "HFX",     timeidx=2)
jan_sh_s     = getvar(jan_s,   "HFX",     timeidx=2)


jan_rf = (jan_swdwn_c-jan_swdwn_c*jan_alb_c+jan_emis_c*jan_glw_c-jan_lh_c-jan_grdflx_c-jan_sh_c)-(jan_swdwn_s-jan_swdwn_s*jan_alb_s+jan_emis_s*jan_glw_s-jan_lh_s-jan_grdflx_s-jan_sh_s)
#====================================================================================================================================================
#jan_rf=jan_swdwn_c+jan_glw_c-jan_swdwn_c*jan_alb_c-0.000000056783*jan_emis_c*jan_tsk_c*jan_tsk_c*jan_tsk_c*jan_tsk_c - (jan_swdwn_s+jan_glw_s-jan_swdwn_s*jan_alb_s-0.000000056783*jan_emis_s*jan_tsk_s*jan_tsk_s*jan_tsk_s*jan_tsk_s)-7
#rf = (swdwn_c-swdwn_c*alb_c+emis_c*glw_c-lh_c-grdflx_c)-(swdwn_s-swdwn_s*alb_s+emis_s*glw_s-lh_s-grdflx_s)
#====== February ============
feb_swdwn_c  = getvar(feb_c,   "SWDOWN", timeidx=2)
feb_swdwn_s  = getvar(feb_s,   "SWDOWN", timeidx=2)
feb_glw_c    = getvar(feb_c,   "GLW",    timeidx=2)
feb_glw_s    = getvar(feb_s,   "GLW",    timeidx=2)
feb_alb_c    = getvar(feb_c,   "ALBEDO", timeidx=2)
feb_alb_s    = getvar(feb_s,   "ALBEDO", timeidx=2)
feb_emis_c   = getvar(feb_c,   "EMISS",  timeidx=2)
feb_emis_s   = getvar(feb_s,   "EMISS",  timeidx=2)
feb_tsk_c    = getvar(feb_c,   "TSK",    timeidx=2)
feb_tsk_s    = getvar(feb_s,   "TSK",    timeidx=2)
feb_grdflx_c = getvar(feb_c,   "GRDFLX", timeidx=2)
feb_grdflx_s = getvar(feb_s,   "GRDFLX", timeidx=2)
feb_lh_c     = getvar(feb_c,   "LH",     timeidx=2)
feb_lh_s     = getvar(feb_s,   "LH",     timeidx=2)
feb_sh_c     = getvar(feb_c,   "HFX",     timeidx=2)
feb_sh_s     = getvar(feb_s,   "HFX",     timeidx=2)


feb_rf = (feb_swdwn_c-feb_swdwn_c*feb_alb_c+feb_emis_c*feb_glw_c-feb_lh_c-feb_grdflx_c-feb_sh_c)-(feb_swdwn_s-feb_swdwn_s*feb_alb_s+feb_emis_s*feb_glw_s-feb_lh_s-feb_grdflx_s-feb_sh_s)
#feb_rf=feb_swdwn_c+feb_glw_c-feb_swdwn_c*feb_alb_c-0.000000056783*feb_emis_c*feb_tsk_c*feb_tsk_c*feb_tsk_c*feb_tsk_c - (feb_swdwn_s+feb_glw_s-feb_swdwn_s*feb_alb_s-0.000000056783*feb_emis_s*feb_tsk_s*feb_tsk_s*feb_tsk_s*feb_tsk_s)-7
#==========================================================================
#====== March ============
mar_swdwn_c  = getvar(mar_c,   "SWDOWN", timeidx=2)
mar_swdwn_s  = getvar(mar_s,   "SWDOWN", timeidx=2)
mar_glw_c    = getvar(mar_c,   "GLW",    timeidx=2)
mar_glw_s    = getvar(mar_s,   "GLW",    timeidx=2)
mar_alb_c    = getvar(mar_c,   "ALBEDO", timeidx=2)
mar_alb_s    = getvar(mar_s,   "ALBEDO", timeidx=2)
mar_emis_c   = getvar(mar_c,   "EMISS",  timeidx=2)
mar_emis_s   = getvar(mar_s,   "EMISS",  timeidx=2)
mar_tsk_c    = getvar(mar_c,   "TSK",    timeidx=2)
mar_tsk_s    = getvar(mar_s,   "TSK",    timeidx=2)
mar_grdflx_c = getvar(mar_c,   "GRDFLX", timeidx=2)
mar_grdflx_s = getvar(mar_s,   "GRDFLX", timeidx=2)
mar_lh_c     = getvar(mar_c,   "LH",     timeidx=2)
mar_lh_s     = getvar(mar_s,   "LH",     timeidx=2)
mar_sh_c     = getvar(mar_c,   "HFX",     timeidx=2)
mar_sh_s     = getvar(mar_s,   "HFX",     timeidx=2)

mar_rf = (mar_swdwn_c-mar_swdwn_c*mar_alb_c+mar_emis_c*mar_glw_c-mar_lh_c-mar_grdflx_c-mar_sh_c)-(mar_swdwn_s-mar_swdwn_s*mar_alb_s+mar_emis_s*mar_glw_s-mar_lh_s-mar_grdflx_s-mar_sh_s)
#mar_rf=mar_swdwn_c+mar_glw_c-mar_swdwn_c*mar_alb_c-0.000000056783*mar_emis_c*mar_tsk_c*mar_tsk_c*mar_tsk_c*mar_tsk_c-(mar_swdwn_s+mar_glw_s-mar_swdwn_s*mar_alb_s-0.000000056783*mar_emis_s*mar_tsk_s*mar_tsk_s*mar_tsk_s*mar_tsk_s)-7
#====== April ============
apr_swdwn_c  = getvar(apr_c,   "SWDOWN", timeidx=2)
apr_swdwn_s  = getvar(apr_s,   "SWDOWN", timeidx=2)
apr_glw_c    = getvar(apr_c,   "GLW",    timeidx=2)
apr_glw_s    = getvar(apr_s,   "GLW",    timeidx=2)
apr_alb_c    = getvar(apr_c,   "ALBEDO", timeidx=2)
apr_alb_s    = getvar(apr_s,   "ALBEDO", timeidx=2)
apr_emis_c   = getvar(apr_c,   "EMISS",  timeidx=2)
apr_emis_s   = getvar(apr_s,   "EMISS",  timeidx=2)
apr_tsk_c    = getvar(apr_c,   "TSK",    timeidx=2)
apr_tsk_s    = getvar(apr_s,   "TSK",    timeidx=2)
apr_grdflx_c = getvar(apr_c,   "GRDFLX", timeidx=2)
apr_grdflx_s = getvar(apr_s,   "GRDFLX", timeidx=2)
apr_lh_c     = getvar(apr_c,   "LH",     timeidx=2)
apr_lh_s     = getvar(apr_s,   "LH",     timeidx=2)
apr_sh_c     = getvar(apr_c,   "HFX",     timeidx=2)
apr_sh_s     = getvar(apr_s,   "HFX",     timeidx=2)


apr_rf = (apr_swdwn_c-apr_swdwn_c*apr_alb_c+apr_emis_c*apr_glw_c-apr_lh_c-apr_grdflx_c-apr_sh_c)-(apr_swdwn_s-apr_swdwn_s*apr_alb_s+apr_emis_s*apr_glw_s-apr_lh_s-apr_grdflx_s-apr_sh_s)
#apr_rf=apr_swdwn_c+apr_glw_c-apr_swdwn_c*apr_alb_c-0.000000056783*apr_emis_c*apr_tsk_c*apr_tsk_c*apr_tsk_c*apr_tsk_c-(apr_swdwn_s+apr_glw_s-apr_swdwn_s*apr_alb_s-0.000000056783*apr_emis_s*apr_tsk_s*apr_tsk_s*apr_tsk_s*apr_tsk_s)-7
#====== May ============
may_swdwn_c  = getvar(may_c,   "SWDOWN", timeidx=2)
may_swdwn_s  = getvar(may_s,   "SWDOWN", timeidx=2)
may_glw_c    = getvar(may_c,   "GLW",    timeidx=2)
may_glw_s    = getvar(may_s,   "GLW",    timeidx=2)
may_alb_c    = getvar(may_c,   "ALBEDO", timeidx=2)
may_alb_s    = getvar(may_s,   "ALBEDO", timeidx=2)
may_emis_c   = getvar(may_c,   "EMISS",  timeidx=2)
may_emis_s   = getvar(may_s,   "EMISS",  timeidx=2)
may_tsk_c    = getvar(may_c,   "TSK",    timeidx=2)
may_tsk_s    = getvar(may_s,   "TSK",    timeidx=2)
may_grdflx_c = getvar(may_c,   "GRDFLX", timeidx=2)
may_grdflx_s = getvar(may_s,   "GRDFLX", timeidx=2)
may_lh_c     = getvar(may_c,   "LH",     timeidx=2)
may_lh_s     = getvar(may_s,   "LH",     timeidx=2)
may_sh_c     = getvar(may_c,   "HFX",     timeidx=2)
may_sh_s     = getvar(may_s,   "HFX",     timeidx=2)


may_rf = (may_swdwn_c-may_swdwn_c*may_alb_c+may_emis_c*may_glw_c-may_lh_c-may_grdflx_c-may_sh_c)-(may_swdwn_s-may_swdwn_s*may_alb_s+may_emis_s*may_glw_s-may_lh_s-may_grdflx_s-may_sh_s)
#may_rf=may_swdwn_c+may_glw_c-may_swdwn_c*may_alb_c-0.000000056783*may_emis_c*may_tsk_c*may_tsk_c*may_tsk_c*may_tsk_c-(may_swdwn_s+may_glw_s-may_swdwn_s*may_alb_s-0.000000056783*may_emis_s*may_tsk_s*may_tsk_s*may_tsk_s*may_tsk_s)-7
#====== June ============
jun_swdwn_c  = getvar(jun_c,   "SWDOWN", timeidx=2)
jun_swdwn_s  = getvar(jun_s,   "SWDOWN", timeidx=2)
jun_glw_c    = getvar(jun_c,   "GLW",    timeidx=2)
jun_glw_s    = getvar(jun_s,   "GLW",    timeidx=2)
jun_alb_c    = getvar(jun_c,   "ALBEDO", timeidx=2)
jun_alb_s    = getvar(jun_s,   "ALBEDO", timeidx=2)
jun_emis_c   = getvar(jun_c,   "EMISS",  timeidx=2)
jun_emis_s   = getvar(jun_s,   "EMISS",  timeidx=2)
jun_tsk_c    = getvar(jun_c,   "TSK",    timeidx=2)
jun_tsk_s    = getvar(jun_s,   "TSK",    timeidx=2)
jun_grdflx_c = getvar(jun_c,   "GRDFLX", timeidx=2)
jun_grdflx_s = getvar(jun_s,   "GRDFLX", timeidx=2)
jun_lh_c     = getvar(jun_c,   "LH",     timeidx=2)
jun_lh_s     = getvar(jun_s,   "LH",     timeidx=2)
jun_sh_c     = getvar(jun_c,   "HFX",     timeidx=2)
jun_sh_s     = getvar(jun_s,   "HFX",     timeidx=2)



jun_rf = (jun_swdwn_c-jun_swdwn_c*jun_alb_c+jun_emis_c*jun_glw_c-jun_lh_c-jun_grdflx_c-jun_sh_c)-(jun_swdwn_s-jun_swdwn_s*jun_alb_s+jun_emis_s*jun_glw_s-jun_lh_s-jun_grdflx_s-jun_sh_s)
#jun_rf=jun_swdwn_c+jun_glw_c-jun_swdwn_c*jun_alb_c-0.000000056783*jun_emis_c*jun_tsk_c*jun_tsk_c*jun_tsk_c*jun_tsk_c-(jun_swdwn_s+jun_glw_s-jun_swdwn_s*jun_alb_s-0.000000056783*jun_emis_s*jun_tsk_s*jun_tsk_s*jun_tsk_s*jun_tsk_s)-7
#====== July ============
jul_swdwn_c   = getvar(jul_c,  "SWDOWN", timeidx=2)
jul_swdwn_s  = getvar(jul_s,   "SWDOWN", timeidx=2)
jul_glw_c    = getvar(jul_c,   "GLW",    timeidx=2)
jul_glw_s    = getvar(jul_s,   "GLW",    timeidx=2)
jul_alb_c    = getvar(jul_c,   "ALBEDO", timeidx=2)
jul_alb_s    = getvar(jul_s,   "ALBEDO", timeidx=2)
jul_emis_c   = getvar(jul_c,   "EMISS",  timeidx=2)
jul_emis_s   = getvar(jul_s,   "EMISS",  timeidx=2)
jul_tsk_c    = getvar(jul_c,   "TSK",    timeidx=2)
jul_tsk_s    = getvar(jul_s,   "TSK",    timeidx=2)
jul_grdflx_c = getvar(jul_c,   "GRDFLX", timeidx=2)
jul_grdflx_s = getvar(jul_s,   "GRDFLX", timeidx=2)
jul_lh_c     = getvar(jul_c,   "LH",     timeidx=2)
jul_lh_s     = getvar(jul_s,   "LH",     timeidx=2)
jul_sh_c     = getvar(jul_c,   "HFX",     timeidx=2)
jul_sh_s     = getvar(jul_s,   "HFX",     timeidx=2)



jul_rf = (jul_swdwn_c-jul_swdwn_c*jul_alb_c+jul_emis_c*jul_glw_c-jul_lh_c-jul_grdflx_c-jul_sh_c)-(jul_swdwn_s-jul_swdwn_s*jul_alb_s+jul_emis_s*jul_glw_s-jul_lh_s-jul_grdflx_s-jul_sh_s)
#jul_rf=jul_swdwn_c+jul_glw_c-jul_swdwn_c*jul_alb_c-0.000000056783*jul_emis_c*jul_tsk_c*jul_tsk_c*jul_tsk_c*jul_tsk_c-(jul_swdwn_s+jul_glw_s-jul_swdwn_s*jul_alb_s-0.000000056783*jul_emis_s*jul_tsk_s*jul_tsk_s*jul_tsk_s*jul_tsk_s)-7
#====== August ============
aug_swdwn_c  = getvar(aug_c,   "SWDOWN", timeidx=2)
aug_swdwn_s  = getvar(aug_s,   "SWDOWN", timeidx=2)
aug_glw_c    = getvar(aug_c,   "GLW",    timeidx=2)
aug_glw_s    = getvar(aug_s,   "GLW",    timeidx=2)
aug_alb_c    = getvar(aug_c,   "ALBEDO", timeidx=2)
aug_alb_s    = getvar(aug_s,   "ALBEDO", timeidx=2)
aug_emis_c   = getvar(aug_c,   "EMISS",  timeidx=2)
aug_emis_s   = getvar(aug_s,   "EMISS",  timeidx=2)
aug_tsk_c    = getvar(aug_c,   "TSK",    timeidx=2)
aug_tsk_s    = getvar(aug_s,   "TSK",    timeidx=2)
aug_grdflx_c = getvar(aug_c,   "GRDFLX", timeidx=2)
aug_grdflx_s = getvar(aug_s,   "GRDFLX", timeidx=2)
aug_lh_c     = getvar(aug_c,   "LH",     timeidx=2)
aug_lh_s     = getvar(aug_s,   "LH",     timeidx=2)
aug_sh_c     = getvar(aug_c,   "HFX",     timeidx=2)
aug_sh_s     = getvar(aug_s,   "HFX",     timeidx=2)



aug_rf = (aug_swdwn_c-aug_swdwn_c*aug_alb_c+aug_emis_c*aug_glw_c-aug_lh_c-aug_grdflx_c-aug_sh_c)-(aug_swdwn_s-aug_swdwn_s*aug_alb_s+aug_emis_s*aug_glw_s-aug_lh_s-aug_grdflx_s-aug_sh_s)
#aug_rf=aug_swdwn_c+aug_glw_c-aug_swdwn_c*aug_alb_c-0.000000056783*aug_emis_c*aug_tsk_c*aug_tsk_c*aug_tsk_c*aug_tsk_c-(aug_swdwn_s+aug_glw_s-aug_swdwn_s*aug_alb_s-0.000000056783*aug_emis_s*aug_tsk_s*aug_tsk_s*aug_tsk_s*aug_tsk_s)-7
#====== September ============
sep_swdwn_c  = getvar(sep_c,   "SWDOWN", timeidx=2)
sep_swdwn_s  = getvar(sep_s,   "SWDOWN", timeidx=2)
sep_glw_c    = getvar(sep_c,   "GLW",    timeidx=2)
sep_glw_s    = getvar(sep_s,   "GLW",    timeidx=2)
sep_alb_c    = getvar(sep_c,   "ALBEDO", timeidx=2)
sep_alb_s    = getvar(sep_s,   "ALBEDO", timeidx=2)
sep_emis_c   = getvar(sep_c,   "EMISS",  timeidx=2)
sep_emis_s   = getvar(sep_s,   "EMISS",  timeidx=2)
sep_tsk_c    = getvar(sep_c,   "TSK",    timeidx=2)
sep_tsk_s    = getvar(sep_s,   "TSK",    timeidx=2)
sep_grdflx_c = getvar(sep_c,   "GRDFLX", timeidx=2)
sep_grdflx_s = getvar(sep_s,   "GRDFLX", timeidx=2)
sep_lh_c     = getvar(sep_c,   "LH",     timeidx=2)
sep_lh_s     = getvar(sep_s,   "LH",     timeidx=2)
sep_sh_c     = getvar(sep_c,   "HFX",     timeidx=2)
sep_sh_s     = getvar(sep_s,   "HFX",     timeidx=2)



sep_rf = (sep_swdwn_c-sep_swdwn_c*sep_alb_c+sep_emis_c*sep_glw_c-sep_lh_c-sep_grdflx_c-sep_sh_c)-(sep_swdwn_s-sep_swdwn_s*sep_alb_s+sep_emis_s*sep_glw_s-sep_lh_s-sep_grdflx_s-sep_sh_s)
#sep_rf=sep_swdwn_c+sep_glw_c-sep_swdwn_c*sep_alb_c-0.000000056783*sep_emis_c*sep_tsk_c*sep_tsk_c*sep_tsk_c*sep_tsk_c-(sep_swdwn_s+sep_glw_s-sep_swdwn_s*sep_alb_s-0.000000056783*sep_emis_s*sep_tsk_s*sep_tsk_s*sep_tsk_s*sep_tsk_s)-7
#====== October ============
oct_swdwn_c  = getvar(oct_c,   "SWDOWN", timeidx=2)
oct_swdwn_s  = getvar(oct_s,   "SWDOWN", timeidx=2)
oct_glw_c    = getvar(oct_c,   "GLW",    timeidx=2)
oct_glw_s    = getvar(oct_s,   "GLW",    timeidx=2)
oct_alb_c    = getvar(oct_c,   "ALBEDO", timeidx=2)
oct_alb_s    = getvar(oct_s,   "ALBEDO", timeidx=2)
oct_emis_c   = getvar(oct_c,   "EMISS",  timeidx=2)
oct_emis_s   = getvar(oct_s,   "EMISS",  timeidx=2)
oct_tsk_c    = getvar(oct_c,   "TSK",    timeidx=2)
oct_tsk_s    = getvar(oct_s,   "TSK",    timeidx=2)
oct_grdflx_c = getvar(oct_c,   "GRDFLX", timeidx=2)
oct_grdflx_s = getvar(oct_s,   "GRDFLX", timeidx=2)
oct_lh_c     = getvar(oct_c,   "LH",     timeidx=2)
oct_lh_s     = getvar(oct_s,   "LH",     timeidx=2)
oct_sh_c     = getvar(oct_c,   "HFX",     timeidx=2)
oct_sh_s     = getvar(oct_s,   "HFX",     timeidx=2)


oct_rf = (oct_swdwn_c-oct_swdwn_c*oct_alb_c+oct_emis_c*oct_glw_c-oct_lh_c-oct_grdflx_c-oct_sh_c)-(oct_swdwn_s-oct_swdwn_s*oct_alb_s+oct_emis_s*oct_glw_s-oct_lh_s-oct_grdflx_s-oct_sh_s)
#oct_rf=oct_swdwn_c+oct_glw_c-oct_swdwn_c*oct_alb_c-0.000000056783*oct_emis_c*oct_tsk_c*oct_tsk_c*oct_tsk_c*oct_tsk_c-(oct_swdwn_s+oct_glw_s-oct_swdwn_s*oct_alb_s-0.000000056783*oct_emis_s*oct_tsk_s*oct_tsk_s*oct_tsk_s*oct_tsk_s)-7
#====== November  ============
nov_swdwn_c  = getvar(nov_c,   "SWDOWN", timeidx=2)
nov_swdwn_s  = getvar(nov_s,   "SWDOWN", timeidx=2)
nov_glw_c    = getvar(nov_c,   "GLW",    timeidx=2)
nov_glw_s    = getvar(nov_s,   "GLW",    timeidx=2)
nov_alb_c    = getvar(nov_c,   "ALBEDO", timeidx=2)
nov_alb_s    = getvar(nov_s,   "ALBEDO", timeidx=2)
nov_emis_c   = getvar(nov_c,   "EMISS",  timeidx=2)
nov_emis_s   = getvar(nov_s,   "EMISS",  timeidx=2)
nov_tsk_c    = getvar(nov_c,   "TSK",    timeidx=2)
nov_tsk_s    = getvar(nov_s,   "TSK",    timeidx=2)
nov_grdflx_c = getvar(nov_c,   "GRDFLX", timeidx=2)
nov_grdflx_s = getvar(nov_s,   "GRDFLX", timeidx=2)
nov_lh_c     = getvar(nov_c,   "LH",     timeidx=2)
nov_lh_s     = getvar(nov_s,   "LH",     timeidx=2)
nov_sh_c     = getvar(nov_c,   "HFX",     timeidx=2)
nov_sh_s     = getvar(nov_s,   "HFX",     timeidx=2)



nov_rf = (nov_swdwn_c-nov_swdwn_c*nov_alb_c+nov_emis_c*nov_glw_c-nov_lh_c-nov_grdflx_c-nov_sh_c)-(nov_swdwn_s-nov_swdwn_s*nov_alb_s+nov_emis_s*nov_glw_s-nov_lh_s-nov_grdflx_s-nov_sh_s)
#nov_rf=nov_swdwn_c+nov_glw_c-nov_swdwn_c*nov_alb_c-0.000000056783*nov_emis_c*nov_tsk_c*nov_tsk_c*nov_tsk_c*nov_tsk_c-(nov_swdwn_s+nov_glw_s-nov_swdwn_s*nov_alb_s-0.000000056783*nov_emis_s*nov_tsk_s*nov_tsk_s*nov_tsk_s*nov_tsk_s)-7
#====== December  ============
dec_swdwn_c  = getvar(dec_c,   "SWDOWN", timeidx=2)
dec_swdwn_s  = getvar(dec_s,   "SWDOWN", timeidx=2)
dec_glw_c    = getvar(dec_c,   "GLW",    timeidx=2)
dec_glw_s    = getvar(dec_s,   "GLW",    timeidx=2)
dec_alb_c    = getvar(dec_c,   "ALBEDO", timeidx=2)
dec_alb_s    = getvar(dec_s,   "ALBEDO", timeidx=2)
dec_emis_c   = getvar(dec_c,   "EMISS",  timeidx=2)
dec_emis_s   = getvar(dec_s,   "EMISS",  timeidx=2)
dec_tsk_c    = getvar(dec_c,   "TSK",    timeidx=2)
dec_tsk_s    = getvar(dec_s,   "TSK",    timeidx=2)
dec_grdflx_c = getvar(dec_c,   "GRDFLX", timeidx=2)
dec_grdflx_s = getvar(dec_s,   "GRDFLX", timeidx=2)
dec_lh_c     = getvar(dec_c,   "LH",     timeidx=2)
dec_lh_s     = getvar(dec_s,   "LH",     timeidx=2)
dec_sh_c     = getvar(dec_c,   "HFX",     timeidx=2)
dec_sh_s     = getvar(dec_s,   "HFX",     timeidx=2)

dec_rf = (dec_swdwn_c-dec_swdwn_c*dec_alb_c+dec_emis_c*dec_glw_c-dec_lh_c-dec_grdflx_c-dec_sh_c)-(dec_swdwn_s-dec_swdwn_s*dec_alb_s+dec_emis_s*dec_glw_s-dec_lh_s-dec_grdflx_s-dec_sh_s)
#dec_rf=dec_swdwn_c+dec_glw_c-dec_swdwn_c*dec_alb_c-0.000000056783*dec_emis_c*dec_tsk_c*dec_tsk_c*dec_tsk_c*dec_tsk_c-(dec_swdwn_s+dec_glw_s-dec_swdwn_s*dec_alb_s-0.000000056783*dec_emis_s*dec_tsk_s*dec_tsk_s*dec_tsk_s*dec_tsk_s)-7
#======= Seasonal radiative forcing ====================
rf_winter  = (jan_rf+feb_rf+dec_rf)/3
rf_spring  = (mar_rf+apr_rf+may_rf)/3
rf_summer  = (jun_rf+jul_rf+aug_rf)/3
rf_autumn  = (sep_rf+oct_rf+nov_rf)/3
#=============================================================
#==== Temperature data - control============
tk1   = wrf.getvar(jan_c, 'T2',timeidx=5)-273
tk2   = wrf.getvar(feb_c, 'T2',timeidx=5)-273
tk3   = wrf.getvar(mar_c, 'T2',timeidx=5)-273
tk4   = wrf.getvar(apr_c, 'T2',timeidx=5)-273
tk5   = wrf.getvar(may_c, 'T2',timeidx=5)-273
tk6   = wrf.getvar(jun_c, 'T2',timeidx=5)-273
tk7   = wrf.getvar(jul_c, 'T2',timeidx=5)-273
tk8   = wrf.getvar(aug_c, 'T2',timeidx=5)-273
tk9   = wrf.getvar(sep_c, 'T2',timeidx=5)-273
tk10  = wrf.getvar(oct_c, 'T2',timeidx=5)-273
tk11  = wrf.getvar(nov_c, 'T2',timeidx=5)-273
tk12  = wrf.getvar(dec_c, 'T2',timeidx=5)-273
#==== Temperature data - control============
tk13   = wrf.getvar(jan_s, 'T2',timeidx=5)-273
tk14   = wrf.getvar(feb_s, 'T2',timeidx=5)-273
tk15   = wrf.getvar(mar_s, 'T2',timeidx=5)-273
tk16   = wrf.getvar(apr_s, 'T2',timeidx=5)-273
tk17   = wrf.getvar(may_s, 'T2',timeidx=5)-273
tk18   = wrf.getvar(jun_s, 'T2',timeidx=5)-273
tk19   = wrf.getvar(jul_s, 'T2',timeidx=5)-273
tk20   = wrf.getvar(aug_s, 'T2',timeidx=5)-273
tk21   = wrf.getvar(sep_s, 'T2',timeidx=5)-273
tk22   = wrf.getvar(oct_s, 'T2',timeidx=5)-273
tk23   = wrf.getvar(nov_s, 'T2',timeidx=5)-273
tk24   = wrf.getvar(dec_s, 'T2',timeidx=5)-273
#================================================
jan_tc = (tk1-tk13)
feb_tc = (tk2-tk14)
mar_tc = (tk3-tk15)
apr_tc = (tk4-tk16)
may_tc = (tk5-tk17)
jun_tc = (tk6-tk18)
jul_tc = (tk7-tk19)
aug_tc = (tk8-tk20)
sep_tc = (tk9-tk21)
oct_tc = (tk10-tk22)
nov_tc = (tk11-tk23)
dec_tc = (tk12-tk24)
#==== Seasonal temperature=================
win_t = (jan_tc+feb_tc+dec_tc)/3
spr_t = (mar_tc+apr_tc+may_tc)/3
sum_t = (jun_tc+jul_tc+aug_tc)/3
aut_t = (sep_tc+oct_tc+nov_tc)/3
#======== Smoothing RF using Gaussion filter =============
sigma_y = 3
sigma_x = 3
sigma = [sigma_y, sigma_x]
winter_rf = sp.ndimage.filters.gaussian_filter(rf_winter, sigma, mode='constant')
spring_rf = sp.ndimage.filters.gaussian_filter(rf_spring, sigma, mode='constant')
summer_rf = sp.ndimage.filters.gaussian_filter(rf_summer, sigma, mode='constant')
autumn_rf = sp.ndimage.filters.gaussian_filter(rf_autumn, sigma, mode='constant')
#======== Smoothing Temp using Gaussion filter =============
winter_t = sp.ndimage.filters.gaussian_filter(win_t, sigma, mode='constant')
spring_t = sp.ndimage.filters.gaussian_filter(spr_t, sigma, mode='constant')
summer_t = sp.ndimage.filters.gaussian_filter(sum_t, sigma, mode='constant')
autumn_t = sp.ndimage.filters.gaussian_filter(aut_t, sigma, mode='constant')
#============= Subplot - Winter  =================================================================
#fig, axs = plt.subplots(2, 2,dpi=150,figsize=(4,4),sharex=True, sharey=True)
fig, axs = plt.subplots(4, 2,figsize=(6,8),dpi=150)
gridspec.GridSpec(4,2)
plt.subplot2grid((4,2), (0,0))
#============================================================
lats, lons = latlon_coords(rf_winter)
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,winter_rf,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1],  color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1],  color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.ylabel('Winter', labelpad=5,fontsize=7)
#plt.clim(-0.3,0.3)
plt.clim(-8,8)
#plt.title(r'Radiative Forcing [$Wm^{-2}$]',fontsize=8)
#============ Spring ====================
plt.subplot2grid((4,2), (1,0))
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,spring_rf,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5,dashes=[4, 1],  color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1],  color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.ylabel('Spring', labelpad=5,fontsize=7)
plt.clim(-8,8)
#plt.clim(-0.3,0.3)
#============ Summer ====================
plt.subplot2grid((4,2), (2,0))
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,summer_rf,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1],  color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1],  color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.ylabel('Summer', labelpad=5,fontsize=7)
plt.clim(-8,8)
#plt.clim(-0.3,0.3)
#============ Autumn ====================
plt.subplot2grid((4,2), (3,0))
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,autumn_rf,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1],  color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=7, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.ylabel('Autumn', labelpad=5,fontsize=7)
plt.clim(-8,8)
#==============================================================================================================
#=== Now plot for Temperature ========
plt.subplot2grid((4,2), (0,1))
#========= Winter===================================================
lats, lons = latlon_coords(rf_winter)
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
s1 = m.pcolormesh(x,y,winter_t,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1],labels=[0,1,0,0], fontsize=7,  color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1],  color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(-0.3,0.3)
#plt.title('Temperature [°C]',fontsize=8)
#========= spring===================================================
plt.subplot2grid((4,2), (1,1))
lats, lons = latlon_coords(rf_winter)
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
s1 = m.pcolormesh(x,y,spring_t,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1],labels=[0,1,0,0], fontsize=7,  color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1],  color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(-0.3,0.3)
#========== Summer =========================================
plt.subplot2grid((4,2), (2,1))
ats, lons = latlon_coords(rf_winter)
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
s1 = m.pcolormesh(x,y,summer_t,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1],labels=[0,1,0,0], fontsize=7,  color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1],  color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(-0.3,0.3)
#========== Autumn  =========================================
plt.subplot2grid((4,2), (3,1))
ats, lons = latlon_coords(rf_winter)
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
s1 = m.pcolormesh(x,y,autumn_t,cmap=cmaps.temp_diff_18lev) #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MPL_RdBu)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[0,1,0,0], fontsize=7, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=7, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(-0.3,0.3)

#=== Setting Colorbar for RF=========
cax1 = fig.add_axes([0.05,0.06,0.4,0.02])   # left, bottom, width, and height
cb1 = fig.colorbar(s,ax=axs[:,:], cax=cax1,ticks=[-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8], orientation='horizontal')
#cb1 = fig.colorbar(s,ax=axs[:,:], cax=cax1,ticks=[-80,-77,-66,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80], orientation='horizontal')
cb1.set_label(r'Radiative Forcing [$Wm^{-2}$]',fontsize=7,labelpad=5)
cb1.ax.tick_params(labelsize=5)
#=== Setting Colorbar for Temperarture=========
cax2 = fig.add_axes([0.5,0.06,0.42,0.02])   # left, bottom, width, and height
cb2 = fig.colorbar(s1, ax=axs[:,:], cax=cax2,ticks=[-0.20,-0.15,-0.10,-0.05,0,0.05,0.1,0.15,0.20], orientation='horizontal')
#cb2 = fig.colorbar(s1, ax=axs[:,:], cax=cax2,ticks=[-0.6,-0.4,-0.2,0,0.5,1,1.5], orientation='horizontal')
cb2.set_label('Temperature [°C]',fontsize=7,labelpad=5)
cb2.ax.tick_params(labelsize=5)

fig.subplots_adjust(top=0.995,
			bottom=0.1,
			left=0.04,
			right=0.93,
			hspace=0.03,
			wspace=0.01)

#plt.savefig("/mnt/h/Paper_work/RF_Temp_300dpi.jpeg",dpi=300)	
plt.show()



