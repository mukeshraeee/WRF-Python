#====== Code for vertical cros-seection plot===========================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors
from cartopy import crs
from cartopy.feature import NaturalEarthFeature, COLORS
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy,get_basemap, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair)
import cmaps
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import StrMethodFormatter
import matplotlib.gridspec as gridspec
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
#========  Define the cross section start and end points=========
cross_start = CoordPair(lat=30, lon=60)
cross_end   = CoordPair(lat=30, lon=110)
#======= Get WRF variables ======================================
ht       = getvar(nc1, "z", timeidx=-1)    # Model Height for Mass Grid
ter      = getvar(nc1, "ter", units="km") #timeidx=-1)  # Model Terrain Height
dbz      = getvar(nc1, "dbz", timeidx=-1)  # Reflectivity
max_dbz  = getvar(nc1, "mdbz", timeidx=-1) # Maximum Reflectivity
#====== Winter  ======================================
bc1 = getvar(nc1, "BC1", timeidx=-1)
bc2 = getvar(nc1, "BC2", timeidx=-1)
bc3 = getvar(nc2, "BC1", timeidx=-1)
bc4 = getvar(nc2, "BC2", timeidx=-1)
bc5 = getvar(nc3, "BC1", timeidx=-1)
bc6 = getvar(nc3, "BC2", timeidx=-1)
#====== spring  ======================================
bc7  = getvar(nc4, "BC1", timeidx=-1)
bc8  = getvar(nc4, "BC2", timeidx=-1)
bc9  = getvar(nc5, "BC1", timeidx=-1)
bc10 = getvar(nc5, "BC2", timeidx=-1)
bc11 = getvar(nc6, "BC1", timeidx=-1)
bc12 = getvar(nc6, "BC2", timeidx=-1)
#====== summer  ======================================
bc13  = getvar(nc7, "BC1", timeidx=-1)
bc14  = getvar(nc7, "BC2", timeidx=-1)
bc15  = getvar(nc8, "BC1", timeidx=-1)
bc16  = getvar(nc8, "BC2", timeidx=-1)
bc17  = getvar(nc9, "BC1", timeidx=-1)
bc18  = getvar(nc9, "BC2", timeidx=-1)
#====== autumn  ======================================
bc19  = getvar(nc10, "BC1", timeidx=-1)
bc20  = getvar(nc10, "BC2", timeidx=-1)
bc21  = getvar(nc11, "BC1", timeidx=-1)
bc22  = getvar(nc11, "BC2", timeidx=-1)
bc23  = getvar(nc12, "BC1", timeidx=-1)
bc24  = getvar(nc12, "BC2", timeidx=-1)

jan = bc1 + bc2
feb = bc3 + bc4
dec = bc5 + bc6
mar = bc7 + bc8
apr = bc9 + bc10
may = bc11 + bc12
jun = bc13 + bc14 
jul = bc15 + bc16
aug = bc17 + bc18
sep = bc19+ bc20
oct = bc21 + bc22
nov = bc23 + bc24

winter_bc = (jan+feb+dec)/3
spring_bc = (mar+apr+may)/3
summer_bc = (jun+jul+aug)/3
autumn_bc = (sep+oct+nov)/3

#=== get ALT variables for converting ug/kg to ug/m3 ===
alt1  = getvar(nc1, "ALT", timeidx=1)[0]
alt2  = getvar(nc2, "ALT", timeidx=1)[0]
alt3  = getvar(nc3, "ALT", timeidx=1)[0]
alt4  = getvar(nc4, "ALT", timeidx=1)[0]
alt5  = getvar(nc5, "ALT", timeidx=1)[0]
alt6  = getvar(nc6, "ALT", timeidx=1)[0]
alt7  = getvar(nc7, "ALT", timeidx=1)[0]
#alt8  = getvar(nc8, "ALT", timeidx=1)[0]
alt9  = getvar(nc9, "ALT", timeidx=1)[0]
alt10  = getvar(nc10, "ALT", timeidx=1)[0]
alt11  = getvar(nc11, "ALT", timeidx=1)[0]
alt12  = getvar(nc12, "ALT", timeidx=1)[0]

winter_alt = (alt1+alt2+alt3)/3
spring_alt = (alt4+alt5+alt6)/3
summer_alt = (alt7+alt9)/2
autumn_alt = (alt10+alt11+alt12)/3

winter = winter_bc/winter_alt
spring = spring_bc/spring_alt
summer = summer_bc/summer_alt
autumn = autumn_bc/autumn_alt

#======Get wind data ================================================
w1 = getvar(nc1, "wa", timeidx=-1)
u1 = getvar(nc1, "ua", timeidx=-1)
w2 = getvar(nc2, "wa", timeidx=-1)
u2 = getvar(nc2, "ua", timeidx=-1)
w3 = getvar(nc3, "wa", timeidx=-1)
u3 = getvar(nc3, "ua", timeidx=-1)
w4 = getvar(nc4, "wa", timeidx=-1)
u4 = getvar(nc4, "ua", timeidx=-1)
w5 = getvar(nc5, "wa", timeidx=-1)
u5 = getvar(nc5, "ua", timeidx=-1)
u6 = getvar(nc6, "ua", timeidx=-1)
w6 = getvar(nc6, "wa", timeidx=-1)
u7 = getvar(nc7, "ua", timeidx=-1)
w7 = getvar(nc7, "wa", timeidx=-1)
u8 = getvar(nc8, "ua", timeidx=-1)
w8 = getvar(nc8, "wa", timeidx=-1)
u9 = getvar(nc9, "ua", timeidx=-1)
w9 = getvar(nc9, "wa", timeidx=-1)
u10 = getvar(nc10, "ua", timeidx=-1)
w10 = getvar(nc10, "wa", timeidx=-1)
u11 = getvar(nc11, "ua", timeidx=-1)
w11 = getvar(nc11, "wa", timeidx=-1)
u12 = getvar(nc12, "ua", timeidx=-1)
w12 = getvar(nc12, "wa", timeidx=-1)
##===  Seasonal U and W wind components =======
w_winter  = (w1+w2+w3)/3 
u_winter  = (u1+u2+u3)/3
w_spring  = (w4+w5+w6)/3
u_spring  = (u4+u5+u6)/3
w_summer  = (w7+w8+w9)/3
u_summer  = (u7+u8+u9)/3
w_autumn  = (w10+w11+w12)/3
u_autumn  = (u10+u11+u12)/3

w_win   = w_winter[0:17,:,:]
u_win   = u_winter[0:17,:,:]
w_sp    = w_spring[0:17,:,:]
u_sp    = u_spring[0:17,:,:]
w_su    = w_summer[0:17,:,:]
u_su    = u_summer[0:17,:,:]
w_au    = w_autumn[0:17,:,:]
u_au    = u_autumn[0:17,:,:]

#=== Getting  model height for mass grid ===========
ht1       = getvar(nc1, "height_agl", units="km") 
ht2       = getvar(nc2, "height_agl", units="km") 
ht3       = getvar(nc3, "height_agl", units="km") 
ht4       = getvar(nc4, "height_agl", units="km") 
ht5       = getvar(nc5, "height_agl", units="km")
ht6       = getvar(nc6, "height_agl", units="km")
ht7       = getvar(nc7, "height_agl", units="km")
ht8       = getvar(nc8, "height_agl", units="km")
ht9       = getvar(nc9, "height_agl", units="km")
ht10       = getvar(nc10, "height_agl", units="km")
ht11       = getvar(nc11, "height_agl", units="km")
ht12       = getvar(nc12, "height_agl", units="km")
# ==== Seasonal ht ===============
ht_winter =  (ht1+ht2+ht3)/3 
ht_spring =  (ht4+ht5+ht6)/3 
ht_summer =  (ht7+ht8+ht9)/3 
ht_autumn =  (ht10+ht11+ht12)/3 

ht_win =  ht_winter[0:17,:,:]
ht_spr =  ht_spring[0:17,:,:]
ht_sum =  ht_summer[0:17,:,:]
ht_aut =  ht_autumn[0:17,:,:]

#ht_win =  ht_winter[0:10,:,:]
#ht_spr =  ht_spring[0:10,:,:]
#ht_sum =  ht_summer[0:10,:,:]
#ht_aut =  ht_autumn[0:10,:,:]

#====== Get seasonal BC ==============
winter_bc = winter[0:17,:,:]
spring_bc = spring[0:17,:,:]
summer_bc = summer[0:17,:,:]
autumn_bc = autumn[0:17,:,:]

##========  Creat subplot figures====================
fig, axs = plt.subplots(2, 2,figsize=(18,9)) #,constrained_layout=True) 
#========  interpolation ==================================================
Z1 = 10**(winter_bc/10.) # Use linear Z for interpolation
W1 = 10**(w_win/10.)
U1 = 10**(u_win/10.)

Z2 = 10**(spring_bc/10.) # Use linear Z for interpolation
W2 = 10**(w_sp/10.)
U2 = 10**(u_sp/10.)

Z3 = 10**(summer_bc/10.) # Use linear Z for interpolation
W3 = 10**(w_su/10.)
U3 = 10**(u_su/10.)

Z4 = 10**(autumn_bc/10.) # Use linear Z for interpolation
W4 = 10**(w_au/10.)
U4 = 10**(u_au/10.)
#========= Compute the vertical cross-section interpolation.  Also, include the
#========= lat/lon points along the cross-section in the metadata by setting latlon to True.
z1_cross = vertcross(Z1, ht_win, wrfin=nc1,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
w1_cross = vertcross(W1, ht_win, wrfin=nc1,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
u1_cross = vertcross(U1, ht_win, wrfin=nc1,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

z2_cross = vertcross(Z2, ht_spr, wrfin=nc5,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
w2_cross = vertcross(W2, ht_spr, wrfin=nc5,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
u2_cross = vertcross(U2, ht_spr, wrfin=nc5,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

z3_cross = vertcross(Z3, ht_sum, wrfin=nc12,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
w3_cross = vertcross(W3, ht_sum, wrfin=nc12,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
u3_cross = vertcross(U3, ht_sum, wrfin=nc12,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

z4_cross = vertcross(Z4, ht_aut, wrfin=nc9,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
w4_cross = vertcross(W4, ht_aut, wrfin=nc9,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
u4_cross = vertcross(U4, ht_aut, wrfin=nc9,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

#===== Convert back to dBz after interpolation ==============================
bc1_cross   = 10.0 * np.log10(z1_cross)
w1_cross    = 10.0 * np.log10(w1_cross)
u1_cross    = 10.0 * np.log10(u1_cross)

bc2_cross   = 10.0 * np.log10(z2_cross)
w2_cross    = 10.0 * np.log10(w2_cross)
u2_cross    = 10.0 * np.log10(u2_cross)

bc3_cross   = 10.0 * np.log10(z3_cross)
w3_cross    = 10.0 * np.log10(w3_cross)
u3_cross    = 10.0 * np.log10(u3_cross)

bc4_cross   = 10.0 * np.log10(z4_cross)
w4_cross    = 10.0 * np.log10(w4_cross)
u4_cross    = 10.0 * np.log10(u4_cross)

#====== Add back the attributes that xarray dropped from the operations above ==
bc1_cross.attrs.update(bc1_cross.attrs)
bc1_cross.attrs["description"] = "Hydrophobic Black Carbon"
bc1_cross.attrs["units"] = "ug/kg-dryair"

bc2_cross.attrs.update(bc2_cross.attrs)
bc2_cross.attrs["description"] = "Hydrophobic Black Carbon"
bc2_cross.attrs["units"] = "ug/kg-dryair"

bc3_cross.attrs.update(bc3_cross.attrs)
bc3_cross.attrs["description"] = "Hydrophobic Black Carbon"
bc3_cross.attrs["units"] = "ug/kg-dryair"

bc4_cross.attrs.update(bc4_cross.attrs)
bc4_cross.attrs["description"] = "Hydrophobic Black Carbon"
bc4_cross.attrs["units"] = "ug/kg-dryair"
#==== Add back the attributes that xarray dropped from the operations above ==
w1_cross.attrs.update(w1_cross.attrs)
w1_cross.attrs["description"] = "destaggered w-wind component"
w1_cross.attrs["units"] = "m s-1"

w2_cross.attrs.update(w2_cross.attrs)
w2_cross.attrs["description"] = "destaggered w-wind component"
w2_cross.attrs["units"] = "m s-1"

w3_cross.attrs.update(w3_cross.attrs)
w3_cross.attrs["description"] = "destaggered w-wind component"
w3_cross.attrs["units"] = "m s-1"

w4_cross.attrs.update(w4_cross.attrs)
w4_cross.attrs["description"] = "destaggered w-wind component"
w4_cross.attrs["units"] = "m s-1"
#===  Add back the attributes that xarray dropped from the operations above ===
u1_cross.attrs.update(u1_cross.attrs)
u1_cross.attrs["description"] = "destaggered u-wind component"
u1_cross.attrs["units"] = "m s-1"

u2_cross.attrs.update(u2_cross.attrs)
u2_cross.attrs["description"] = "destaggered u-wind component"
u2_cross.attrs["units"] = "m s-1"

u3_cross.attrs.update(u3_cross.attrs)
u3_cross.attrs["description"] = "destaggered u-wind component"
u3_cross.attrs["units"] = "m s-1"

u4_cross.attrs.update(u4_cross.attrs)
u4_cross.attrs["description"] = "destaggered u-wind component"
u4_cross.attrs["units"] = "m s-1"

# Make a copy of the z cross data. Let's use regular numpy arrays for this.
bc1_cross_filled = np.ma.copy(to_np(bc1_cross))
w1_cross_filled = np.ma.copy(to_np(w1_cross))
u1_cross_filled = np.ma.copy(to_np(u1_cross))

bc2_cross_filled = np.ma.copy(to_np(bc2_cross))
w2_cross_filled = np.ma.copy(to_np(w2_cross))
u2_cross_filled = np.ma.copy(to_np(u2_cross))

bc3_cross_filled = np.ma.copy(to_np(bc3_cross))
w3_cross_filled = np.ma.copy(to_np(w3_cross))
u3_cross_filled = np.ma.copy(to_np(u3_cross))

bc4_cross_filled = np.ma.copy(to_np(bc4_cross))
w4_cross_filled = np.ma.copy(to_np(w4_cross))
u4_cross_filled = np.ma.copy(to_np(u4_cross))

# For each cross section column, find the first index with non-missing
# values and copy these to the missing elements below.
for i in range(bc1_cross_filled.shape[-1]):
    column_vals = bc1_cross_filled[:,i]
   #Let's find the lowest index that isn't filled. The nonzero function
   #finds all unmasked values greater than 0. Since 0 is a valid value
   #for dBZ, let's change that threshold to be -200 dBZ instead.
    first_idx = int(np.transpose((column_vals > 0).nonzero())[0])
    bc1_cross_filled[0:first_idx, i] = bc1_cross_filled[first_idx, i]

for i in range(bc2_cross_filled.shape[-1]):
    column_vals = bc2_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > 0).nonzero())[0])
    bc2_cross_filled[0:first_idx, i] = bc2_cross_filled[first_idx, i]

for i in range(bc3_cross_filled.shape[-1]):
    column_vals = bc3_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > 0).nonzero())[0])
    bc3_cross_filled[0:first_idx, i] = bc3_cross_filled[first_idx, i]

for i in range(bc4_cross_filled.shape[-1]):
    column_vals = bc4_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > 0).nonzero())[0])
    bc4_cross_filled[0:first_idx, i] = bc4_cross_filled[first_idx, i]

#for i in range(w_cross_filled.shape[-1]):
#    column_vals = w_cross_filled[:,i]
#    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
#    w_cross_filled[0:first_idx, i] = w_cross_filled[first_idx, i]

#for i in range(u_cross_filled.shape[-1]):
#    column_vals = u_cross_filled[:,i]
#    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
#    u_cross_filled[0:first_idx, i] = u_cross_filled[first_idx, i]

#====  Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=nc1, start_point=cross_start,end_point=cross_end)
#===============  Get the lat/lon points =============
lats, lons = latlon_coords(dbz)
cart_proj = get_cartopy(dbz)
#============= Make the cross section plot =============
bc_levels = np.arange(0, 3, 0.1)

xs1 = np.arange(0, bc1_cross.shape[-1], 1)
ys1 = to_np(bc1_cross.coords["vertical"])

xs2 = np.arange(0, bc2_cross.shape[-1], 1)
ys2 = to_np(bc2_cross.coords["vertical"])

xs3 = np.arange(0, bc3_cross.shape[-1], 1)
ys3 = to_np(bc3_cross.coords["vertical"])

xs4 = np.arange(0, bc4_cross.shape[-1], 1)
ys4 = to_np(bc4_cross.coords["vertical"])

#===============================================================================================================
bc_contours = axs[0,0].contourf(xs1,ys1,to_np(bc1_cross_filled),levels=bc_levels,cmap=cmaps.WhiteBlueGreenYellowRed)
bc_contours = axs[0,1].contourf(xs2,ys2,to_np(bc2_cross_filled),levels=bc_levels,cmap=cmaps.WhiteBlueGreenYellowRed)
bc_contours = axs[1,0].contourf(xs3,ys3,to_np(bc3_cross_filled),levels=bc_levels,cmap=cmaps.WhiteBlueGreenYellowRed)
bc_contours = axs[1,1].contourf(xs4,ys4,to_np(bc4_cross_filled),levels=bc_levels,cmap=cmaps.WhiteBlueGreenYellowRed)
#=========  Fill in the mountain area =========================================
ht_fill = axs[0,0].fill_between(xs1, 0, to_np(ter_line),facecolor="grey")
ht_fill = axs[0,1].fill_between(xs2, 0, to_np(ter_line),facecolor="grey")
ht_fill = axs[1,0].fill_between(xs3, 0, to_np(ter_line),facecolor="grey")
ht_fill = axs[1,1].fill_between(xs4, 0, to_np(ter_line),facecolor="grey")
#======== plot wind  =============================================================
axs[0,0].quiver(xs1[::4], ys1[::4],to_np(u1_cross_filled[::4, ::4]), to_np(w1_cross_filled[::4, ::4]*500))
axs[0,1].quiver(xs2[::4], ys2[::4],to_np(u2_cross_filled[::4, ::4]), to_np(w2_cross_filled[::4, ::4]*500))
axs[1,0].quiver(xs3[::4], ys3[::4],to_np(u3_cross_filled[::4, ::4]), to_np(w3_cross_filled[::4, ::4]*500))
axs[1,1].quiver(xs4[::4], ys4[::4],to_np(u4_cross_filled[::4, ::4]), to_np(w4_cross_filled[::4, ::4]*500))

# Set the x-ticks to use latitude and longitude labels
coord_pairs = to_np(bc1_cross.coords["xy_loc"])
x_ticks     = np.arange(coord_pairs.shape[0])
x_labels    = [pair.latlon_str(fmt="{:.2f}, {:.2f}") 
              for pair in to_np(coord_pairs)]
#=======  Set the desired number of x ticks below
num_ticks = 7
thin = int((len(x_ticks) / num_ticks) + .10)

#ax[0,0].set_xticks(x_ticks[::thin])
#ax[0,0].set_xticklabels(x_labels[::thin], rotation=35, fontsize=10)

#ax[0,1].set_xticks(x_ticks[::thin])
#ax[0,1].set_xticklabels(x_labels[::thin], rotation=35, fontsize=10)
axs[0,0].set_xticklabels([])
axs[0,1].set_xticklabels([])
axs[0,1].set_yticklabels([])
axs[1,1].set_yticklabels([])

axs[1,0].set_xticks(x_ticks[::thin])
axs[1,0].set_xticklabels(x_labels[::thin],rotation=25, fontsize=10)
#axs[1,0].xaxis.set_major_locator(MaxNLocator(integer=True))
#plt.gca().axs[1,0].set_major_formatter(StrMethodFormatter('{x:,.2f}'))
#plt.gca().ticklabel_format(axis='both', style='plain', useOffset=False)

axs[1,1].set_xticks(x_ticks[::thin])
axs[1,1].set_xticklabels(x_labels[::thin],rotation=25,fontsize=10)
#===== Set the x-axis and  y-axis labels
#axs[1,0].set_xlabel("Latitude", fontsize=12)
#axs[1,0].set_ylabel("Height [m]", fontsize=12)

#fig.text(0.5, 0.06, 'Latitude', ha='center', va='center',fontsize=15)
fig.text(0.08, 0.5, 'Height [Km]', va='center', rotation='vertical',fontsize=12)

#ax.xaxis.set_major_locator(MaxNLocator(integer=False))
#bar = fig.colorbar(bc_contours, ax=ax)
cbar= fig.colorbar(bc_contours, ax=axs[:,:], location = 'right')
#fig.colorbar((bc_contours, ax=ax[0,1])
#cbar.ax[0,1].tick_params(labelsize=12)
#cbar.ax[1,1].tick_params(labelsize=12)
cbar.set_label(r'$BC [μg m ^{-3}]$', rotation=-270, fontsize=10)
#cbar.set_label(r'$Extinction Coef.at 550nm [km^{-1}]$', rotation=-270, fontsize=12)
#cbar.set_label(r'$Extinction Coef.at 550nm [km^{-1}]$', rotation=-270, fontsize=12)

#==== PLot title ===========
axs[0,0].set_title('Winter',fontsize=12)
axs[0,1].set_title('Spring',fontsize=12)
axs[1,0].set_title('Summer',fontsize=12)
axs[1,1].set_title('Autumn',fontsize=12)

plt.show()

