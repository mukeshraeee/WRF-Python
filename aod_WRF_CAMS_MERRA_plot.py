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
nc3 = Dataset("wrfout_d01_2017-03-01_00_00_00.nc")
nc4 = Dataset("wrfout_d01_2017-04-01_00_00_00.nc")
nc5 = Dataset("wrfout_d01_2017-05-01_00_00_00.nc")
nc6 = Dataset("wrfout_d01_2017-06-01_00_00_00")
nc7 = Dataset("wrfout_d01_2017-07-01_00_00_00")
nc8 = Dataset("wrfout_d01_2017-08-01_00_00_00")
nc9 = Dataset("wrfout_d01_2017-09-01_00_00_00")
nc10 = Dataset("wrfout_d01_2017-10-01_00_00_00")
nc11 = Dataset("wrfout_d01_2017-11-01_00_00_00")
nc12 = Dataset("wrfout_d01_2017-12-01_00_00_00")

#====== Get Extinction Profile for all Months ======
jan = getvar(nc1, "EXTCOF55", timeidx=1)
feb = getvar(nc2, "EXTCOF55", timeidx=1)
mar = getvar(nc3, "EXTCOF55", timeidx=1)
apr = getvar(nc4, "EXTCOF55", timeidx=1)
may = getvar(nc5, "EXTCOF55", timeidx=1)
jun = getvar(nc6, "EXTCOF55", timeidx=1)
jul = getvar(nc7, "EXTCOF55", timeidx=1)
aug = getvar(nc8, "EXTCOF55", timeidx=1)
sep = getvar(nc9, "EXTCOF55", timeidx=1)
oct = getvar(nc10, "EXTCOF55", timeidx=1)
nov = getvar(nc11, "EXTCOF55", timeidx=1)
dec = getvar(nc12, "EXTCOF55", timeidx=1)
#======= Get geopotential Height and ======
ph1  = getvar(nc1, "PH", timeidx=1)
ph2  = getvar(nc2, "PH", timeidx=1)
ph3  = getvar(nc3, "PH", timeidx=1)
ph4  = getvar(nc4, "PH", timeidx=1)
ph5  = getvar(nc5, "PH", timeidx=1)
ph6  = getvar(nc6, "PH", timeidx=1)
ph7  = getvar(nc7, "PH", timeidx=1)
ph8  = getvar(nc8, "PH", timeidx=1)
ph9  = getvar(nc9, "PH", timeidx=1)
ph10  = getvar(nc10, "PH", timeidx=1)
ph11  = getvar(nc11, "PH", timeidx=1)
ph12  = getvar(nc12, "PH", timeidx=1)
#===== Get Pertubed geoptential Heitght======= 
phb1  = getvar(nc1, "PHB", timeidx=1)
phb2  = getvar(nc2, "PHB", timeidx=1)
phb3  = getvar(nc3, "PHB", timeidx=1)
phb4  = getvar(nc4, "PHB", timeidx=1)
phb5  = getvar(nc5, "PHB", timeidx=1)
phb6  = getvar(nc6, "PHB", timeidx=1)
phb7  = getvar(nc7, "PHB", timeidx=1)
phb8  = getvar(nc8, "PHB", timeidx=1)
phb9  = getvar(nc9, "PHB", timeidx=1)
phb10  = getvar(nc10, "PHB", timeidx=1)
phb11  = getvar(nc11, "PHB", timeidx=1)
phb12  = getvar(nc12, "PHB", timeidx=1)
#===== Calculate  thikness /depth of the layer  =====
z   = (ph1+phb1)/9.81 #jan
b   = (ph2+phb2)/9.81 #feb
c   = (ph3+phb3)/9.81 #mar
d   = (ph4+phb4)/9.81 #apr
e   = (ph5+phb5)/9.81 #may
f   = (ph6+phb6)/9.81 #jun
g   = (ph7+phb7)/9.81 #jul
h   = (ph8+phb8)/9.81 #aug
i   = (ph9+phb9)/9.81 #sep
j   = (ph10+phb10)/9.81 #oct
k   = (ph11+phb11)/9.81 #nov
l   = (ph12+phb12)/9.81 #dec

#=======  Get in each model layers  =====
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
z19 = z[19]
z20 = z[20]
z21 = z[21]
z22 = z[22]
z23 = z[23]
z24 = z[24]
z25 = z[25]
z26 = z[26]
z27 = z[27]
z28 = z[28]
z29 = z[29]
z30 = z[30]
z31 = z[31]
z32 = z[32]
z33 = z[33]
z34 = z[34]
z35 = z[35]
z36 = z[36]
z37 = z[37]
z38 = z[38]
z39 = z[39]

b0 = b[0]
b1 = b[1]
b2 = b[2]
b3 = b[3]
b4 = b[4]
b5 = b[5]
b6 = b[6]
b7 = b[7]
b8 = b[8]
b9 = b[9]
b10 = b[10]
b11 = b[11]
b12 = b[12]
b13 = b[13]
b14 = b[14]
b15 = b[15]
b16 = b[16]
b17 = b[17]
b18 = b[18]
b19 = b[19]
b20 = b[20]
b21 = b[21]
b22 = b[22]
b23 = b[23]
b24 = b[24]
b25 = b[25]
b26 = b[26]
b27 = b[27]
b28 = b[28]
b29 = b[29]
b30 = b[30]
b31 = b[31]
b32 = b[32]
b33 = b[33]
b34 = b[34]
b35 = b[35]
b36 = b[36]
b37 = b[37]
b38 = b[38]
b39 = b[39]

c0 = c[0]
c1 = c[1]
c2 = c[2]
c3 = c[3]
c4 = c[4]
c5 = c[5]
c6 = c[6]
c7 = c[7]
c8 = c[8]
c9 = c[9]
c10 = c[10]
c11 = c[11]
c12 = c[12]
c13 = c[13]
c14 = c[14]
c15 = c[15]
c16 = c[16]
c17 = c[17]
c18 = c[18]
c19 = c[19]
c20 = c[20]
c21 = c[21]
c22 = c[22]
c23 = c[23]
c24 = c[24]
c25 = c[25]
c26 = c[26]
c27 = c[27]
c28 = c[28]
c29 = c[29]
c30 = c[30]
c31 = c[31]
c32 = c[32]
c33 = c[33]
c34 = c[34]
c35 = c[35]
c36 = c[36]
c37 = c[37]
c38 = c[38]
c39 = c[39]

d0 = d[0]
d1 = d[1]
d2 = d[2]
d3 = d[3]
d4 = d[4]
d5 = d[5]
d6 = d[6]
d7 = d[7]
d8 = d[8]
d9 = d[9]
d10 = d[10]
d11 = d[11]
d12 = d[12]
d13 = d[13]
d14 = d[14]
d15 = d[15]
d16 = d[16]
d17 = d[17]
d18 = d[18]
d19 = d[19]
d20 = d[20]
d21 = d[21]
d22 = d[22]
d23 = d[23]
d24 = d[24]
d25 = d[25]
d26 = d[26]
d27 = d[27]
d28 = d[28]
d29 = d[29]
d30 = d[30]
d31 = d[31]
d32 = d[32]
d33 = d[33]
d34 = d[34]
d35 = d[35]
d36 = d[36]
d37 = d[37]
d38 = d[38]
d39 = d[39]

e0 = e[0]
e1 = e[1]
e2 = e[2]
e3 = e[3]
e4 = e[4]
e5 = e[5]
e6 = e[6]
e7 = e[7]
e8 = e[8]
e9 = e[9]
e10 =e[10]
e11 = e[11]
e12 = e[12]
e13 = e[13]
e14 = e[14]
e15 = e[15]
e16 = e[16]
e17 = e[17]
e18 = e[18]
e19 = e[19]
e20 = e[20]
e21 = e[21]
e22 = e[22]
e23 = e[23]
e24 = e[24]
e25 = e[25]
e26 = e[26]
e27 = e[27]
e28 = e[28]
e29 = e[29]
e30 = e[30]
e31 = e[31]
e32 = e[32]
e33 = e[33]
e34 = e[34]
e35 = e[35]
e36 = e[36]
e37 = e[37]
e38 = e[38]
e39 = e[39]

f0 = f[0]
f1 = f[1]
f2 = f[2]
f3 = f[3]
f4 = f[4]
f5 = f[5]
f6 = f[6]
f7 = f[7]
f8 = f[8]
f9 = f[9]
f10 = f[10]
f11 = f[11]
f12 = f[12]
f13 = f[13]
f14 = f[14]
f15 = f[15]
f16 = f[16]
f17 = f[17]
f18 = f[18]
f19 = f[19]
f20 = f[20]
f21 = f[21]
f22 = f[22]
f23 = f[23]
f24 = f[24]
f25 = f[25]
f26 = f[26]
f27 = f[27]
f28 = f[28]
f29 = f[29]
f30 = f[30]
f31 = f[31]
f32 = f[32]
f33 = f[33]
f34 = f[34]
f35 = f[35]
f36 = f[36]
f37 = f[37]
f38 = f[38]
f39 = f[39]

g0 = g[0]
g1 = g[1]
g2 = g[2]
g3 = g[3]
g4 = g[4]
g5 = g[5]
g6 = g[6]
g7 = g[7]
g8 = g[8]
g9 = g[9]
g10 = g[10]
g11 = g[11]
g12 = g[12]
g13 = g[13]
g14 = g[14]
g15 = g[15]
g16 = g[16]
g17 = g[17]
g18 = g[18]
g19 = g[19]
g20 = g[20]
g21 = g[21]
g22 = g[22]
g23 = g[23]
g24 = g[24]
g25 = g[25]
g26 = g[26]
g27 = g[27]
g28 = g[28]
g29 = g[29]
g30 = g[30]
g31 = g[31]
g32 = g[32]
g33 = g[33]
g34 = g[34]
g35 = g[35]
g36 = g[36]
g37 = g[37]
g38 = g[38]
g39 = g[39]

h0 = h[0]
h1 = h[1]
h2 = h[2]
h3 = h[3]
h4 = h[4]
h5 = h[5]
h6 = h[6]
h7 = h[7]
h8 = h[8]
h9 = h[9]
h10 = h[10]
h11 = h[11]
h12 = h[12]
h13 = h[13]
h14 = h[14]
h15 = h[15]
h16 = h[16]
h17 = h[17]
h18 = h[18]
h19 = h[19]
h20 = h[20]
h21 = h[21]
h22 = h[22]
h23 = h[23]
h24 = h[24]
h25 = h[25]
h26 = h[26]
h27 = h[27]
h28 = h[28]
h29 = h[29]
h30 = h[30]
h31 = h[31]
h32 = h[32]
h33 = h[33]
h34 = h[34]
h35 = h[35]
h36 = h[36]
h37 = h[37]
h38 = h[38]
h39 = h[39]

i0 = i[0]
i1 = i[1]
i2 = i[2]
i3 = i[3]
i4 = i[4]
i5 = i[5]
i6 = i[6]
i7 = i[7]
i8 = i[8]
i9 = i[9]
i10 = i[10]
i11 = i[11]
i12 = i[12]
i13 = i[13]
i14 = i[14]
i15 = i[15]
i16 = i[16]
i17 = i[17]
i18 = i[18]
i19 = i[19]
i20 = i[20]
i21 = i[21]
i22 = i[22]
i23 = i[23]
i24 = i[24]
i25 = i[25]
i26 = i[26]
i27 = i[27]
i28 = i[28]
i29 = i[29]
i30 = i[30]
i31 = i[31]
i32 = i[32]
i33 = i[33]
i34 = i[34]
i35 = i[35]
i36 = i[36]
i37 = i[37]
i38 = i[38]
i39 = i[39]

j0 = j[0]
j1 = j[1]
j2 = j[2]
j3 = j[3]
j4 = j[4]
j5 = j[5]
j6 = j[6]
j7 = j[7]
j8 = j[8]
j9 = j[9]
j10 = j[10]
j11 = j[11]
j12 = j[12]
j13 = j[13]
j14 = j[14]
j15 = j[15]
j16 = j[16]
j17 = j[17]
j18 = j[18]
j19 = j[19]
j20 = j[20]
j21 = j[21]
j22 = j[22]
j23 = j[23]
j24 = j[24]
j25 = j[25]
j26 = j[26]
j27 = j[27]
j28 = j[28]
j29 = j[29]
j30 = j[30]
j31 = j[31]
j32 = j[32]
j33 = j[33]
j34 = j[34]
j35 = j[35]
j36 = j[36]
j37 = j[37]
j38 = j[38]
j39 = j[39]

k0 = k[0]
k1 = k[1]
k2 = k[2]
k3 = k[3]
k4 = k[4]
k5 = k[5]
k6 = k[6]
k7 = k[7]
k8 = k[8]
k9 = k[9]
k10 = k[10]
k11 = k[11]
k12 = k[12]
k13 = k[13]
k14 = k[14]
k15 = k[15]
k16 = k[16]
k17 = k[17]
k18 = k[18]
k19 = k[19]
k20 = k[20]
k21 = k[21]
k22 = k[22]
k23 = k[23]
k24 = k[24]
k25 = k[25]
k26 = k[26]
k27 = k[27]
k28 = k[28]
k29 = k[29]
k30 = k[30]
k31 = k[31]
k32 = k[32]
k33 = k[33]
k34 = k[34]
k35 = k[35]
k36 = k[36]
k37 = k[37]
k38 = k[38]
k39 = k[39]

l0 = l[0]
l1 = l[1]
l2 = l[2]
l3 = l[3]
l4 = l[4]
l5 = l[5]
l6 = l[6]
l7 = l[7]
l8 = l[8]
l9 = l[9]
l10 = l[10]
l11 = l[11]
l12 = l[12]
l13 = l[13]
l14 = l[14]
l15 = l[15]
l16 = l[16]
l17 = l[17]
l18 = l[18]
l19 = l[19]
l20 = l[20]
l21 = l[21]
l22 = l[22]
l23 = l[23]
l24 = l[24]
l25 = l[25]
l26 = l[26]
l27 = l[27]
l28 = l[28]
l29 = l[29]
l30 = l[30]
l31 = l[31]
l32 = l[32]
l33 = l[33]
l34 = l[34]
l35 = l[35]
l36 = l[36]
l37 = l[37]
l38 = l[38]
l39 = l[39]

#==== Get level wise exticntion coef.  ==========
##====== January =================
jan0 = jan[0]
jan1 = jan[1]
jan2 = jan[2]
jan3 = jan[3]
jan4 = jan[4]
jan5 = jan[5]
jan6 = jan[6]
jan7 = jan[7]
jan8 = jan[8]
jan9 = jan[9]
jan10 = jan[10]
jan11 = jan[11]
jan12 = jan[12]
jan13 = jan[13]
jan14 = jan[14]
jan15 = jan[15]
jan16 = jan[16]
jan17 = jan[17]
jan18 = jan[18]
jan19 = jan[19]
jan20 = jan[20]
jan21 = jan[21]
jan22 = jan[22]
jan23 = jan[23]
jan24 = jan[24]
jan25 = jan[25]
jan26 = jan[26]
jan27 = jan[27]
jan28 = jan[28]
jan29 = jan[29]
jan30 = jan[30]
jan31 = jan[31]
jan32 = jan[32]
jan33 = jan[33]
jan34 = jan[34]
jan35 = jan[35]
jan36 = jan[36]
jan37 = jan[37]
jan38 = jan[38]

#===== February =================
feb0 = feb[0]
feb1 = feb[1]
feb2 = feb[2]
feb3 = feb[3]
feb4 = feb[4]
feb5 = feb[5]
feb6 = feb[6]
feb7 = feb[7]
feb8 = feb[8]
feb9 = feb[9]
feb10 = feb[10]
feb11 = feb[11]
feb12 = feb[12]
feb13 = feb[13]
feb14 = feb[14]
feb15 = feb[15]
feb16 = feb[16]
feb17 = feb[17]
feb18 = feb[18]
feb19 = feb[19]
feb20 = feb[20]
feb21 = feb[21]
feb22 = feb[22]
feb23 = feb[23]
feb24 = feb[24]
feb25 = feb[25]
feb26 = feb[26]
feb27 = feb[27]
feb28 = feb[28]
feb29 = feb[29]
feb30 = feb[30]
feb31 = feb[31]
feb32 = feb[32]
feb33 = feb[33]
feb34 = feb[34]
feb35 = feb[35]
feb36 = feb[36]
feb37 = feb[37]
feb38 = feb[38]

#=== March =======================
mar0 = mar[0]
mar1 = mar[1]
mar2 = mar[2]
mar3 = mar[3]
mar4 = mar[4]
mar5 = mar[5]
mar6 = mar[6]
mar7 = mar[7]
mar8 = mar[8]
mar9 = mar[9]
mar10 =mar[10]
mar11 = mar[11]
mar12 = mar[12]
mar13 = mar[13]
mar14 = mar[14]
mar15 = mar[15]
mar16 = mar[16]
mar17 = mar[17]
mar18 = mar[18]
mar19 = mar[19]
mar20 = mar[20]
mar21 = mar[21]
mar22 = mar[22]
mar23 = mar[23]
mar24 = mar[24]
mar25 = mar[25]
mar26 = mar[26]
mar27 = mar[27]
mar28 = mar[28]
mar29 = mar[29]
mar30 = mar[30]
mar31 = mar[31]
mar32 = mar[32]
mar33 = mar[33]
mar34 = mar[34]
mar35 = mar[35]
mar36 = mar[36]
mar37 = mar[37]
mar38 = mar[38]

#=== April ================
apr0 = apr[0]
apr1 = apr[1]
apr2 = apr[2]
apr3 = apr[3]
apr4 = apr[4]
apr5 = apr[5]
apr6 = apr[6]
apr7 = apr[7]
apr8 = apr[8]
apr9 = apr[9]
apr10 = apr[10]
apr11 = apr[11]
apr12 = apr[12]
apr13 = apr[13]
apr14 = apr[14]
apr15 = apr[15]
apr16 = apr[16]
apr17 = apr[17]
apr18 = apr[18]
apr19 = apr[19]
apr20 = apr[20]
apr21 = apr[21]
apr22 = apr[22]
apr23 = apr[23]
apr24 = apr[24]
apr25 = apr[25]
apr26 = apr[26]
apr27 = apr[27]
apr28 = apr[28]
apr29 = apr[29]
apr30 = apr[30]
apr31 = apr[31]
apr32 = apr[32]
apr33 = apr[33]
apr34 = apr[34]
apr35 = apr[35]
apr36 = apr[36]
apr37 = apr[37]
apr38 = apr[38]

#=== May ==============
may0 = may[0]
may1 = may[1]
may2 = may[2]
may3 = may[3]
may4 = may[4]
may5 = may[5]
may6 = may[6]
may7 = may[7]
may8 = may[8]
may9 = may[9]
may10 =may[10]
may11 = may[11]
may12 = may[12]
may13 = may[13]
may14 = may[14]
may15 = may[15]
may16 = may[16]
may17 = may[17]
may18 = may[18]
may19 = may[19]
may20 = may[20]
may21 = may[21]
may22 = may[22]
may23 = may[23]
may24 = may[24]
may25 = may[25]
may26 = may[26]
may27 = may[27]
may28 = may[28]
may29 = may[29]
may30 = may[30]
may31 = may[31]
may32 = may[32]
may33 = may[33]
may34 = may[34]
may35 = may[35]
may36 = may[36]
may37 = may[37]
may38 = may[38]

##===== June ==============
jun0 = jun[0]
jun1 = jun[1]
jun2 = jun[2]
jun3 = jun[3]
jun4 = jun[4]
jun5 = jun[5]
jun6 = jun[6]
jun7 = jun[7]
jun8 = jun[8]
jun9 = jun[9]
jun10 =jun[10]
jun11 = jun[11]
jun12 = jun[12]
jun13 = jun[13]
jun14 = jun[14]
jun15 = jun[15]
jun16 = jun[16]
jun17 = jun[17]
jun18 = jun[18]
jun19 = jun[19]
jun20 = jun[20]
jun21 = jun[21]
jun22 = jun[22]
jun23 = jun[23]
jun24 = jun[24]
jun25 = jun[25]
jun26 = jun[26]
jun27 = jun[27]
jun28 = jun[28]
jun29 = jun[29]
jun30 = jun[30]
jun31 = jun[31]
jun32 = jun[32]
jun33 = jun[33]
jun34 = jun[34]
jun35 = jun[35]
jun36 = jun[36]
jun37 = jun[37]
jun38 = jun[38]

##=== July ==============
jul0 = jul[0]
jul1 = jul[1]
jul2 = jul[2]
jul3 = jul[3]
jul4 = jul[4]
jul5 = jul[5]
jul6 = jul[6]
jul7 = jul[7]
jul8 = jul[8]
jul9 = jul[9]
jul10 = jul[10]
jul11 = jul[11]
jul12 = jul[12]
jul13 = jul[13]
jul14 = jul[14]
jul15 = jul[15]
jul16 = jul[16]
jul17 = jul[17]
jul18 = jul[18]
jul19 = jul[19]
jul20 = jul[20]
jul21 = jul[21]
jul22 = jul[22]
jul23 = jul[23]
jul24 = jul[24]
jul25 = jul[25]
jul26 = jul[26]
jul27 = jul[27]
jul28 = jul[28]
jul29 = jul[29]
jul30 = jul[30]
jul31 = jul[31]
jul32 = jul[32]
jul33 = jul[33]
jul34 = jul[34]
jul35 = jul[35]
jul36 = jul[36]
jul37 = jul[37]
jul38 = jul[38]

##===== August ========
aug0 = aug[0]
aug1 = aug[1]
aug2 = aug[2]
aug3 = aug[3]
aug4 = aug[4]
aug5 = aug[5]
aug6 = aug[6]
aug7 = aug[7]
aug8 = aug[8]
aug9 = aug[9]
aug10 = aug[10]
aug11 = aug[11]
aug12 = aug[12]
aug13 = aug[13]
aug14 = aug[14]
aug15 = aug[15]
aug16 = aug[16]
aug17 = aug[17]
aug18 = aug[18]
aug19 = aug[19]
aug20 = aug[20]
aug21 = aug[21]
aug22 = aug[22]
aug23 = aug[23]
aug24 = aug[24]
aug25 = aug[25]
aug26 = aug[26]
aug27 = aug[27]
aug28 = aug[28]
aug29 = aug[29]
aug30 = aug[30]
aug31 = aug[31]
aug32 = aug[32]
aug33 = aug[33]
aug34 = aug[34]
aug35 = aug[35]
aug36 = aug[36]
aug37 = aug[37]
aug38 = aug[38]

###==== September =========
sep0 = sep[0]
sep1 = sep[1]
sep2 = sep[2]
sep3 = sep[3]
sep4 = sep[4]
sep5 = sep[5]
sep6 = sep[6]
sep7 = sep[7]
sep8 = sep[8]
sep9 = sep[9]
sep10 = sep[10]
sep11 = sep[11]
sep12 = sep[12]
sep13 = sep[13]
sep14 = sep[14]
sep15 = sep[15]
sep16 = sep[16]
sep17 = sep[17]
sep18 = sep[18]
sep19 = sep[19]
sep20 = sep[20]
sep21 = sep[21]
sep22 = sep[22]
sep23 = sep[23]
sep24 = sep[24]
sep25 = sep[25]
sep26 = sep[26]
sep27 = sep[27]
sep28 = sep[28]
sep29 = sep[29]
sep30 = sep[30]
sep31 = sep[31]
sep32 = sep[32]
sep33 = sep[33]
sep34 = sep[34]
sep35 = sep[35]
sep36 = sep[36]
sep37 = sep[37]
sep38 = sep[38]

##===== Ocotber ==============
oct0 = oct[0]
oct1 = oct[1]
oct2 = oct[2]
oct3 = oct[3]
oct4 = oct[4]
oct5 = oct[5]
oct6 = oct[6]
oct7 = oct[7]
oct8 = oct[8]
oct9 = oct[9]
oct10 =oct[10]
oct11 = oct[11]
oct12 = oct[12]
oct13 = oct[13]
oct14 = oct[14]
oct15 = oct[15]
oct16 = oct[16]
oct17 = oct[17]
oct18 = oct[18]
oct19 = oct[19]
oct20 = oct[20]
oct21 = oct[21]
oct22 = oct[22]
oct23 = oct[23]
oct24 = oct[24]
oct25 = oct[25]
oct26 = oct[26]
oct27 = oct[27]
oct28 = oct[28]
oct29 = oct[29]
oct30 = oct[30]
oct31 = oct[31]
oct32 = oct[32]
oct33 = oct[33]
oct34 = oct[34]
oct35 = oct[35]
oct36 = oct[36]
oct37 = oct[37]
oct38 = oct[38]

##=== November ===============
nov0 = nov[0]
nov1 = nov[1]
nov2 = nov[2]
nov3 = nov[3]
nov4 = nov[4]
nov5 = nov[5]
nov6 = nov[6]
nov7 = nov[7]
nov8 = nov[8]
nov9 = nov[9]
nov10 = nov[10]
nov11 = nov[11]
nov12 = nov[12]
nov13 = nov[13]
nov14 = nov[14]
nov15 = nov[15]
nov16 = nov[16]
nov17 = nov[17]
nov18 = nov[18]
nov19 = nov[19]
nov20 = nov[20]
nov21 = nov[21]
nov22 = nov[22]
nov23 = nov[23]
nov24 = nov[24]
nov25 = nov[25]
nov26 = nov[26]
nov27 = nov[27]
nov28 = nov[28]
nov29 = nov[29]
nov30 = nov[30]
nov31 = nov[31]
nov32 = nov[32]
nov33 = nov[33]
nov34 = nov[34]
nov35 = nov[35]
nov36 = nov[36]
nov37 = nov[37]
nov38 = nov[38]

##==== December ============
dec0 = dec[0]
dec1 = dec[1]
dec2 = dec[2]
dec3 = dec[3]
dec4 = dec[4]
dec5 = dec[5]
dec6 = dec[6]
dec7 = dec[7]
dec8 = dec[8]
dec9 = dec[9]
dec10 = dec[10]
dec11 = dec[11]
dec12 = dec[12]
dec13 = dec[13]
dec14 = dec[14]
dec15 = dec[15]
dec16 = dec[16]
dec17 = dec[17]
dec18 = dec[18]
dec19 = dec[19]
dec20 = dec[20]
dec21 = dec[21]
dec22 = dec[22]
dec23 = dec[23]
dec24 = dec[24]
dec25 = dec[25]
dec26 = dec[26]
dec27 = dec[27]
dec28 = dec[28]
dec29 = dec[29]
dec30 = dec[30]
dec31 = dec[31]
dec32 = dec[32]
dec33 = dec[33]
dec34 = dec[34]
dec35 = dec[35]
dec36 = dec[36]
dec37 = dec[37]
dec38 = dec[38]

##=====  Now AOD calculation ======================
jan_aod = jan0*(z1-z0)+jan1*(z2-z1)+jan2*(z3-z2)+jan3*(z4-z3)+jan4*(z5-z4)+jan5*(z6-z5)+jan6*(z7-z6)+jan7*(z8-z7)+jan8*(z9-z8)+jan9*(z10-z9)+jan10*(z11-z10)+jan11*(z12-z11)+jan12*(z13-z12)+jan13*(z14-z13)+jan14*(z15-z14)+jan15*(z16-z15)+jan16*(z17-z16)+jan17*(z18-z17)+jan18*(z19-z18)+jan19*(z20-z19)+jan20*(z21-z20)+jan21*(z22-z21)+jan22*(z23-z22)+jan23*(z24-z23)+jan24*(z25-z24)+jan25*(z26-z25)+jan26*(z27-z26)+jan27*(z28-z27)+jan28*(z29-z28)+jan29*(z30-z29)+jan30*(z31-z30)+jan31*(z32-z31)+jan32*(z33-z32)+jan33*(z34-z33)+jan34*(z35-z34)+jan35*(z36-z35)+jan36*(z37-z36)+jan37*(z38-z37)+jan38*(z39-z38)
feb_aod = feb0*(b1-b0)+feb1*(b2-b1)+feb2*(b3-b2)+feb3*(b4-b3)+feb4*(b5-b4)+feb5*(b6-b5)+feb6*(b7-b6)+feb7*(b8-b7)+feb8*(b9-b8)+feb9*(b10-b9)+feb10*(b11-b10)+feb11*(b12-b11)+feb12*(b13-b12)+feb13*(b14-b13)+feb14*(b15-b14)+feb15*(b16-b15)+feb16*(b17-b16)+feb17*(b18-b17)+feb18*(b19-b18)+feb19*(b20-b19)+feb20*(b21-b20)+feb21*(b22-b21)+feb22*(b23-b22)+feb23*(b24-b23)+feb24*(b25-b24)+feb25*(b26-b25)+feb26*(b27-b26)+feb27*(b28-b27)+feb28*(b29-b28)+feb29*(b30-b29)+feb30*(b31-b30)+feb31*(b32-b31)+feb32*(b33-b32)+feb33*(b34-b33)+feb34*(b35-b34)+feb35*(b36-b35)+feb36*(b37-b36)+feb37*(b38-b37)+feb38*(b39-b38)
mar_aod = mar0*(c1-c0)+mar1*(c2-c1)+mar2*(c3-c2)+mar3*(c4-c3)+mar4*(c5-c4)+mar5*(c6-c5)+mar6*(c7-c6)+mar7*(c8-c7)+mar8*(c9-c8)+mar9*(c10-c9)+mar10*(c11-c10)+mar11*(c12-c11)+mar12*(c13-c12)+mar13*(c14-c13)+mar14*(c15-c14)+mar15*(c16-c15)+mar16*(c17-c16)+mar17*(c18-c17)+mar18*(c19-c18)+mar19*(c20-c19)+mar20*(c21-c20)+mar21*(c22-c21)+mar22*(c23-c22)+mar23*(c24-c23)+mar24*(c25-c24)+mar25*(c26-c25)+mar26*(c27-c26)+mar27*(c28-c27)+mar28*(c29-c28)+mar29*(c30-c29)+mar30*(c31-c30)+mar31*(c32-c31)+mar32*(c33-c32)+mar33*(c34-c33)+mar34*(c35-c34)+mar35*(c36-c35)+mar36*(c37-c36)+mar37*(c38-c37)+mar38*(c39-c38)
apr_aod = apr0*(d1-d0)+apr1*(d2-d1)+apr2*(d3-d2)+apr3*(d4-d3)+apr4*(d5-d4)+apr5*(d6-z5)+apr6*(d7-d6)+apr7*(d8-d7)+apr8*(d9-d8)+apr9*(d10-d9)+apr10*(d11-d10)+apr11*(d12-d11)+apr12*(d13-d12)+apr13*(d14-d13)+apr14*(d15-d14)+apr15*(d16-d15)+apr16*(d17-d16)+apr17*(d18-d17)+apr18*(d19-d18)+apr19*(d20-d19)+apr20*(d21-d20)+apr21*(d22-d21)+apr22*(d23-d22)+apr23*(d24-d23)+apr24*(d25-d24)+apr25*(d26-d25)+apr26*(d27-d26)+apr27*(d28-d27)+apr28*(d29-d28)+apr29*(d30-d29)+apr30*(d31-d30)+apr31*(d32-d31)+apr32*(d33-d32)+apr33*(d34-d33)+apr34*(d35-d34)+apr35*(d36-d35)+apr36*(d37-d36)+apr37*(d38-d37)+apr38*(d39-d38)
may_aod = may0*(e1-e0)+may1*(e2-e1)+may2*(e3-e2)+may3*(e4-e3)+may4*(e5-e4)+may5*(e6-e5)+may6*(e7-e6)+may7*(e8-e7)+may8*(e9-e8)+may9*(e10-e9)+may10*(e11-e10)+may11*(e12-e11)+may12*(e13-e12)+may13*(e14-e13)+may14*(e15-e14)+may15*(e16-e15)+may16*(e17-e16)+may17*(e18-e17)+may18*(e19-e18)+may19*(e20-e19)+may20*(e21-e20)+may21*(e22-e21)+may22*(e23-e22)+may23*(e24-e23)+may24*(e25-e24)+may25*(e26-e25)+may26*(e27-e26)+may27*(e28-e27)+may28*(e29-e28)+may29*(e30-e29)+may30*(e31-e30)+may31*(e32-e31)+may32*(e33-e32)+may33*(e34-e33)+may34*(e35-e34)+may35*(e36-e35)+may36*(e37-e36)+may37*(e38-e37)+may38*(e39-e38)
jun_aod = jun0*(f1-f0)+jun1*(f2-f1)+jun2*(f3-f2)+jun3*(f4-f3)+jun4*(f5-f4)+jun5*(f6-f5)+jun6*(f7-f6)+jun7*(f8-f7)+jun8*(f9-f8)+jun9*(f10-f9)+jun10*(f11-f10)+jun11*(f12-f11)+jun12*(f13-f12)+jun13*(f14-f13)+jun14*(f15-f14)+jun15*(f16-f15)+jun16*(f17-f16)+jun17*(f18-f17)+jun18*(f19-f18)+jun19*(f20-f19)+jun20*(f21-f20)+jun21*(f22-f21)+jun22*(f23-f22)+jun23*(f24-f23)+jun24*(f25-f24)+jun25*(f26-f25)+jun26*(f27-f26)+jun27*(f28-f27)+jun28*(f29-f28)+jun29*(f30-f29)+jun30*(f31-f30)+jun31*(f32-f31)+jun32*(f33-f32)+jun33*(f34-f33)+jun34*(f35-f34)+jun35*(f36-f35)+jun36*(f37-f36)+jun37*(f38-f37)+jun38*(f39-f38)
jul_aod = jul0*(g1-g0)+jul1*(g2-g1)+jul2*(g3-g2)+jul3*(g4-g3)+jul4*(g5-g4)+jul5*(g6-g5)+jul6*(g7-g6)+jul7*(g8-g7)+jul8*(g9-g8)+jul9*(g10-g9)+jul10*(g11-g10)+jul11*(g12-g11)+jul12*(g13-g12)+jul13*(g14-g13)+jul14*(g15-g14)+jul15*(g16-g15)+jul16*(g17-g16)+jul17*(g18-g17)+jul18*(g19-g18)+jul19*(g20-g19)+jul20*(g21-g20)+jul21*(g22-g21)+jul22*(g23-g22)+jul23*(g24-g23)+jul24*(g25-g24)+jul25*(g26-g25)+jul26*(g27-g26)+jul27*(g28-g27)+jul28*(g29-g28)+jul29*(g30-g29)+jul30*(g31-g30)+jul31*(g32-g31)+jul32*(g33-g32)+jul33*(g34-g33)+jul34*(g35-g34)+jul35*(g36-g35)+jul36*(g37-g36)+jul37*(g38-g37)+jul38*(g39-g38)
aug_aod = aug0*(h1-h0)+aug1*(h2-h1)+aug2*(h3-h2)+aug3*(h4-h3)+aug4*(h5-h4)+aug5*(h6-h5)+aug6*(h7-h6)+aug7*(h8-h7)+aug8*(h9-h8)+aug9*(h10-h9)+aug10*(h11-h10)+aug11*(h12-h11)+aug12*(h13-h12)+aug13*(h14-h13)+aug14*(h15-h14)+aug15*(h16-h15)+aug16*(h17-h16)+aug17*(h18-h17)+aug18*(h19-h18)+aug19*(h20-h19)+aug20*(h21-h20)+aug21*(h22-h21)+aug22*(h23-h22)+aug23*(h24-h23)+aug24*(h25-h24)+aug25*(h26-h25)+aug26*(h27-h26)+aug27*(h28-h27)+aug28*(h29-h28)+aug29*(h30-h29)+aug30*(h31-h30)+aug31*(h32-h31)+aug32*(h33-h32)+aug33*(h34-h33)+aug34*(h35-h34)+aug35*(h36-h35)+aug36*(h37-h36)+aug37*(h38-h37)+aug38*(h39-h38)
sep_aod = sep0*(i1-i0)+sep1*(i2-i1)+sep2*(i3-i2)+sep3*(i4-i3)+sep4*(i5-i4)+sep5*(i6-i5)+sep6*(i7-i6)+sep7*(i8-i7)+sep8*(i9-i8)+sep9*(i10-i9)+sep10*(i11-i10)+sep11*(i12-i11)+sep12*(i13-i12)+sep13*(i14-i13)+sep14*(i15-i14)+sep15*(i16-i15)+sep16*(i17-i16)+sep17*(i18-i17)+sep18*(i19-i18)+sep19*(i20-i19)+sep20*(i21-i20)+sep21*(i22-i21)+sep22*(i23-i22)+sep23*(i24-i23)+sep24*(i25-i24)+sep25*(i26-i25)+sep26*(i27-i26)+sep27*(i28-i27)+sep28*(i29-i28)+sep29*(i30-i29)+sep30*(i31-i30)+sep31*(i32-i31)+sep32*(i33-i32)+sep33*(i34-i33)+sep34*(i35-i34)+sep35*(i36-i35)+sep36*(i37-i36)+sep37*(i38-i37)+sep38*(i39-i38)
oct_aod = oct0*(j1-j0)+oct1*(j2-j1)+oct2*(j3-j2)+oct3*(j4-j3)+oct4*(j5-j4)+oct5*(j6-j5)+oct6*(j7-j6)+oct7*(j8-j7)+oct8*(j9-j8)+oct9*(j10-j9)+oct10*(j11-j10)+oct11*(j12-j11)+oct12*(j13-j12)+oct13*(j14-j13)+oct14*(j15-j14)+oct15*(j16-j15)+oct16*(j17-j16)+oct17*(j18-j17)+oct18*(j19-j18)+oct19*(j20-j19)+oct20*(j21-j20)+oct21*(j22-j21)+oct22*(j23-j22)+oct23*(j24-j23)+oct24*(j25-j24)+oct25*(j26-j25)+oct26*(j27-j26)+oct27*(j28-j27)+oct28*(j29-j28)+oct29*(j30-j29)+oct30*(j31-j30)+oct31*(j32-j31)+oct32*(j33-j32)+oct33*(j34-j33)+oct34*(j35-j34)+oct35*(j36-j35)+oct36*(j37-j36)+oct37*(j38-j37)+oct38*(j39-j38)
nov_aod = nov0*(k1-k0)+nov1*(k2-k1)+nov2*(k3-k2)+nov3*(k4-k3)+nov4*(k5-k4)+nov5*(k6-k5)+nov6*(k7-k6)+nov7*(k8-k7)+nov8*(k9-k8)+nov9*(k10-k9)+nov10*(k11-k10)+nov11*(k12-k11)+nov12*(k13-k12)+nov13*(k14-k13)+nov14*(k15-k14)+nov15*(k16-k15)+nov16*(k17-k16)+nov17*(k18-k17)+nov18*(k19-k18)+nov19*(k20-k19)+nov20*(k21-k20)+nov21*(k22-k21)+nov22*(k23-k22)+nov23*(k24-k23)+nov24*(k25-k24)+nov25*(k26-k25)+nov26*(k27-k26)+nov27*(k28-k27)+nov28*(k29-k28)+nov29*(k30-k29)+nov30*(k31-k30)+nov31*(k32-k31)+nov32*(k33-k32)+nov33*(k34-k33)+nov34*(k35-k34)+nov35*(k36-k35)+nov36+(k37-k36)+nov37*(k38-k37)+nov38*(k39-k38)
dec_aod = dec0*(l1-l0)+dec1*(l2-l1)+dec2*(l3-l2)+dec3*(l4-l3)+dec4*(l5-l4)+dec5*(l6-l5)+dec6*(l7-l6)+dec7*(l8-l7)+dec8*(l9-l8)+dec9*(l10-l9)+dec10*(l11-l10)+dec11*(l12-l11)+dec12*(l13-l12)+dec13*(l14-l13)+dec14*(l15-l14)+dec15*(l16-l15)+dec16*(l17-l16)+dec17*(l18-l17)+dec18*(l19-l18)+dec19*(l20-l19)+dec20*(l21-l20)+dec21*(l22-l21)+dec22*(l23-l22)+dec23*(l24-l23)+dec24*(l25-l24)+dec25*(l26-l25)+dec26*(l27-l26)+dec27*(l28-l27)+dec28*(l29-l28)+dec29*(l30-l29)+dec30*(l31-l30)+dec31*(l32-l31)+dec32*(l33-l32)+dec33*(l34-l33)+dec34*(l35-l34)+dec35*(l36-l35)+dec36*(l37-l36)+dec37*(l38-l37)+dec38*(l39-l38)

##===== Seasonal AOD ========================
wrf_winter = (dec_aod+jan_aod+feb_aod)/3 * 0.001
wrf_spring = (mar_aod+apr_aod+may_aod)/3 * 0.001
wrf_summer = (jun_aod+jul_aod+aug_aod)/3 * 0.001 
wrf_autumn = (sep_aod+oct_aod+nov_aod)/3 * 0.001

#====== Let's smooth by using Gussian filter  ============
sigma_y = 2.0
sigma_x = 2.0
sigma = [sigma_y, sigma_x]

aod_wrf_winter = sp.ndimage.filters.gaussian_filter(wrf_winter, sigma, mode='constant')
aod_wrf_spring = sp.ndimage.filters.gaussian_filter(wrf_spring, sigma, mode='constant')
aod_wrf_summer = sp.ndimage.filters.gaussian_filter(wrf_summer, sigma, mode='constant')
aod_wrf_autumn = sp.ndimage.filters.gaussian_filter(wrf_autumn, sigma, mode='constant')

#===== Load CAMS AOD Data =====
cams_winter = Dataset('cams_winter.nc')
cams_spring = Dataset('cams_spring.nc')
cams_summer = Dataset('cams_summer.nc')
cams_autumn = Dataset('cams_autumn.nc')
## MERRA-2 AOD data ===================
merra_winter = Dataset('merra_winter.nc')
merra_spring = Dataset('merra_spring.nc')
merra_summer = Dataset('merra_summer.nc')
merra_autumn = Dataset('merra_autumn.nc')
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

#====== Sub - Plot ====================
fig, axs = plt.subplots(3, 4,figsize=(18,9))
gridspec.GridSpec(3,4)
#=========== WRF-AOD-Winter ================================
plt.subplot2grid((3,4), (0,0))
lats, lons = latlon_coords(ph1)
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax = m.pcolormesh(x,y,aod_wrf_winter,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)
plt.title('Winter')
plt.ylabel('WRF-AOD', labelpad=40,fontsize=12)

#=========== WRF-AOD-Spring ================================
plt.subplot2grid((3,4), (0,1))
lats, lons = latlon_coords(ph1)
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax = m.pcolormesh(x,y,aod_wrf_spring,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)
plt.title('Spring')

#=========== WRF-AOD-Summer ================================
plt.subplot2grid((3,4), (0,2))
lats, lons = latlon_coords(ph1)
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax = m.pcolormesh(x,y,aod_wrf_summer,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)
plt.title('Summer')
#=========== WRF-AOD-Autumn ================================
plt.subplot2grid((3,4), (0,3))
lats, lons = latlon_coords(ph1)
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lon, lat = np.meshgrid(lons, lats)
x, y = m(to_np(lons), to_np(lats))
ax = m.pcolormesh(x,y,aod_wrf_autumn,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)
plt.title('Autumn')

#=========  CAMS-winter =========================== 
plt.subplot2grid((3,4), (1,0))
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lons,lats = np.meshgrid(lon_cam,lat_cam) # set the latitude and longitude variables from the data
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,aod_cams_winter_av,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)
plt.ylabel('CAMS-AOD', labelpad=40,fontsize=12)

#=========  CAMS-Spring ===========================
plt.subplot2grid((3,4), (1,1))
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lons,lats = np.meshgrid(lon_cam,lat_cam)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,aod_cams_spring_av,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)

#=========  CAMS-Summer ===========================
plt.subplot2grid((3,4), (1,2))
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lons,lats = np.meshgrid(lon_cam,lat_cam)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,aod_cams_summer_av,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)


#=========  CAMS-Autumn ===========================
plt.subplot2grid((3,4), (1,3))
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lons,lats = np.meshgrid(lon_cam,lat_cam)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,aod_cams_autumn_av,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)


#=========  MERRA-winter ===========================
plt.subplot2grid((3,4), (2,0))
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lons,lats = np.meshgrid(lon_merra,lat_merra)
xi,yi = m(lons, lats)
ax = m.pcolormesh(xi,yi,aod_mera_winter_av,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)
plt.ylabel('MERRA-2-AOD', labelpad=40,fontsize=12)

#=========  MERRA-Spring ===========================
plt.subplot2grid((3,4), (2,1))
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lons,lats = np.meshgrid(lon_merra,lat_merra)
xi,yi = m(lons, lats)
ax = m.pcolormesh(xi,yi,aod_mera_spring_av,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)


#=========  MERRA-Summer ===========================
plt.subplot2grid((3,4), (2,2))
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lons,lats = np.meshgrid(lon_merra,lat_merra)
xi,yi = m(lons, lats)
ax = m.pcolormesh(xi,yi,aod_mera_summer_av,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)


#=========  MERRA-Autumn ===========================
plt.subplot2grid((3,4), (2,3))
m = Basemap(projection='merc',llcrnrlon=45,llcrnrlat=12,urcrnrlon=118,urcrnrlat=52,resolution='h')
lons,lats = np.meshgrid(lon_merra,lat_merra)
xi,yi = m(lons, lats)
ax = m.pcolormesh(xi,yi,aod_mera_autumn_av,cmap=cmaps.WhiteBlueGreenYellowRed)
#========== Add Grid lines =======================================================================================
m.drawparallels(np.arange(10, 55, 10), linewidth=0.5, dashes=[4, 1], labels=[1,0,0,0], fontsize=10, color='black')
m.drawmeridians(np.arange(45, 120,15),linewidth=0.5, dashes=[4, 1], labels=[0,0,0,1], fontsize=10, color='black')
#========= Add countries boundasy =============================================================================
m.drawcountries()
m.drawcoastlines(linewidth=0.3)
plt.clim(0,1)




cbar= fig.colorbar(ax, ax=axs[:,:], location = 'right')
cbar.set_label('AOD', rotation=-270, fontsize=10)

plt.show()









