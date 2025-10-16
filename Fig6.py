import xarray as xr

data_path = '../data/'


#reading the temperature dataset

dset=xr.open_dataset(data_path+'SAT_anomaly.nc')
temp=dset.__xarray_dataarray_variable__

years=xr.date_range(start='1979-01-01', end='2023-12-31', freq='YE')

# Select months of interest for NDJFMA and MJJASO
months_ndjfma = [11, 12, 1, 2, 3, 4]
months_mjjaso = [5, 6, 7, 8, 9, 10]

# Select data for NDJFMA and MJJASO
temp_ndjfma1 = temp.sel(time=temp['time.month'].isin(months_ndjfma))
temp_mjjaso1 = temp.sel(time=temp['time.month'].isin(months_mjjaso))

#calculating yearly mean values of temperatures for every NDJFMA and MJJASO
temp_ndjfma=temp_ndjfma1[0]
for i in range(len(years.year)):
    temp = temp_ndjfma1.sel(time=slice(str(years.year[i])+'-11-01', str(years.year[i]+1)+'-04-30')).mean('time')
    temp_ndjfma =xr.concat([temp_ndjfma, temp],dim='years')

temp_ndjfma=temp_ndjfma[1:].assign_coords(years=years)

temp_mjjaso=temp_mjjaso1[0]
for i in range(len(years.year)):
    temp = temp_mjjaso1.sel(time=slice(str(years.year[i])+'-05-01', str(years.year[i])+'-10-30')).mean('time')
    temp_mjjaso =xr.concat([temp_mjjaso, temp],dim='years')

temp_mjjaso=temp_mjjaso[1:].assign_coords(years=years)


#reading the precipitation dataset
dset=xr.open_dataset(data_path+'precip_anomaly.nc')
precip=dset.__xarray_dataarray_variable__

years=xr.date_range(start='1979-01-01', end='2023-12-31', freq='YE')

# Select months of interest for NDJFMA and MJJASO
months_ndjfma = [11, 12, 1, 2, 3, 4]
months_mjjaso = [5, 6, 7, 8, 9, 10]

# Select data for NDJFMA and MJJASO
precip_ndjfma1 = precip.sel(time=precip['time.month'].isin(months_ndjfma))
precip_mjjaso1 = precip.sel(time=precip['time.month'].isin(months_mjjaso))

#calculating yearly mean values of precipitation for every NDJFMA and MJJASO
precip_ndjfma=precip_ndjfma1[0]
for i in range(len(years.year)):
    precip = precip_ndjfma1.sel(time=slice(str(years.year[i])+'-11-01', str(years.year[i]+1)+'-04-30')).mean('time')
    precip_ndjfma =xr.concat([precip_ndjfma, precip],dim='years')

precip_ndjfma=precip_ndjfma[1:].assign_coords(years=years)

precip_mjjaso=precip_mjjaso1[0]
for i in range(len(years.year)):
    precip = precip_mjjaso1.sel(time=slice(str(years.year[i])+'-05-01', str(years.year[i])+'-10-30')).mean('time')
    precip_mjjaso =xr.concat([precip_mjjaso, precip],dim='years')

precip_mjjaso=precip_mjjaso[1:].assign_coords(years=years)

#reading the sea ice dataset

dset_sie=xr.open_dataset(data_path+'seaice_conc.nc').rename({'tdim':'time'})
sea_ice_extent = dset_sie.cdr_seaice_conc_monthly.sel(time=slice('1979-01-01','2023-12-01'))

months_of_interest = [11,12,1,2,3,4] 
sic_ndjfma1 = sea_ice_extent.sel(time=sea_ice_extent['time.month'].isin(months_of_interest))
months_of_interest = [5,6,7,8,9,10] 
sic_mjjaso1 = sea_ice_extent.sel(time=sea_ice_extent['time.month'].isin(months_of_interest))
years=xr.date_range(start='1979-01-01', end='2023-12-31', freq='YE')

sic_ndjfma=sic_ndjfma1[0]
for i in range(len(years.year)):
    sic = sic_ndjfma1.sel(time=slice(str(years.year[i])+'-11-01', str(years.year[i]+1)+'-04-30')).mean('time')
    sic_ndjfma =xr.concat([sic_ndjfma, sic],dim='years')

sic_ndjfma=sic_ndjfma[1:].assign_coords(years=years)

sic_mjjaso=sic_mjjaso1[0]
for i in range(len(years.year)):
    sic = sic_mjjaso1.sel(time=slice(str(years.year[i])+'-05-01', str(years.year[i])+'-10-30')).mean('time')
    sic_mjjaso =xr.concat([sic_mjjaso, sic],dim='years')

sic_mjjaso=sic_mjjaso[1:].assign_coords(years=years)

#reading SAM index
dset=xr.open_dataset(data_path+'SAM_GW_1m_1979-2023.nc')
sam=dset.SAM

months_of_interest = [11,12,1,2,3,4] 
sam_ndjfma1 = sam.sel(time=sam['time.month'].isin(months_of_interest))

months_of_interest = [5,6,7,8,9,10] 
sam_mjjaso1 = sam.sel(time=sam['time.month'].isin(months_of_interest))

sam_ndjfma=sam_ndjfma1[0]
for i in range(len(years.year)):
    temp = sam_ndjfma1.sel(time=slice(str(years.year[i])+'-11-01', str(years.year[i]+1)+'-04-30')).mean('time')
    sam_ndjfma =xr.concat([sam_ndjfma, temp],dim='years')

sam_ndjfma=sam_ndjfma[1:].assign_coords(years=years)

sam_mjjaso=sam_mjjaso1[0]
for i in range(len(years.year)):
    temp = sam_mjjaso1.sel(time=slice(str(years.year[i])+'-05-01', str(years.year[i])+'-10-30')).mean('time')
    sam_mjjaso =xr.concat([sam_mjjaso, temp],dim='years')

sam_mjjaso=sam_mjjaso[1:].assign_coords(years=years)


#reading wind anomalies
dset=xr.open_dataset(data_path+'u850_anomaly.nc')
wind=dset.__xarray_dataarray_variable__
months_of_interest = [11,12,1,2,3,4] 
wind_ndfma = wind.sel(time=wind['time.month'].isin(months_of_interest))
months_of_interest = [5,6,7,8,9,10] 
wind_mjjaso = wind.sel(time=wind['time.month'].isin(months_of_interest))

# Select data for NDJFMA and MJJASO
wind_ndjfma1 = wind.sel(time=wind['time.month'].isin(months_ndjfma))
wind_mjjaso1 = wind.sel(time=wind['time.month'].isin(months_mjjaso))

wind_ndjfma=wind_ndjfma1[0]
for i in range(len(years.year)):
    wind = wind_ndjfma1.sel(time=slice(str(years.year[i])+'-11-01', str(years.year[i]+1)+'-04-30')).mean('time')
    wind_ndjfma =xr.concat([wind_ndjfma, wind],dim='years')

wind_ndjfma=wind_ndjfma[1:].assign_coords(years=years)

wind_mjjaso=wind_mjjaso1[0]
for i in range(len(years.year)):
    wind = wind_mjjaso1.sel(time=slice(str(years.year[i])+'-05-01', str(years.year[i])+'-10-30')).mean('time')
    wind_mjjaso =xr.concat([wind_mjjaso, wind],dim='years')

wind_mjjaso=wind_mjjaso[1:].assign_coords(years=years)

dset=xr.open_dataset(data_path+'v850_anomaly.nc')
vwind=dset.__xarray_dataarray_variable__
months_of_interest = [11,12,1,2,3,4] 
vwind_ndfma = vwind.sel(time=vwind['time.month'].isin(months_of_interest))
months_of_interest = [5,6,7,8,9,10] 
vwind_mjjaso = vwind.sel(time=vwind['time.month'].isin(months_of_interest))

# Select data for NDJFMA and MJJASO
vwind_ndjfma1 = vwind.sel(time=vwind['time.month'].isin(months_ndjfma))
vwind_mjjaso1 = vwind.sel(time=vwind['time.month'].isin(months_mjjaso))

vwind_ndjfma=vwind_ndjfma1[0]
for i in range(len(years.year)):
    vwind = vwind_ndjfma1.sel(time=slice(str(years.year[i])+'-11-01', str(years.year[i]+1)+'-04-30')).mean('time')
    vwind_ndjfma =xr.concat([vwind_ndjfma, vwind],dim='years')

vwind_ndjfma=vwind_ndjfma[1:].assign_coords(years=years)

vwind_mjjaso=vwind_mjjaso1[0]
for i in range(len(years.year)):
    vwind = vwind_mjjaso1.sel(time=slice(str(years.year[i])+'-05-01', str(years.year[i])+'-10-30')).mean('time')
    vwind_mjjaso =xr.concat([vwind_mjjaso, vwind],dim='years')

vwind_mjjaso=vwind_mjjaso[1:].assign_coords(years=years)

import numpy as np
import sacpy as scp
import matplotlib.pyplot as plt
import sacpy.Map # need cartopy or you can just not import
import cartopy.crs as ccrs
from numpy import ma

#Linear regression is performed by regressing the SAM index onto field anomalies using the Sacpy linear regression function

temp_jja = scp.LinReg(sam_mjjaso[0:],temp_mjjaso[0:])
temp_djf = scp.LinReg(sam_ndjfma[0:],temp_ndjfma[0:])

precip_jja = scp.LinReg(sam_mjjaso[0:],precip_mjjaso[0:])
precip_djf = scp.LinReg(sam_ndjfma[0:],precip_ndjfma[0:])

wind_jja = scp.LinReg(sam_mjjaso[0:],wind_mjjaso[0:])
wind_djf = scp.LinReg(sam_ndjfma[0:],wind_ndjfma[0:])

vwind_jja = scp.LinReg(sam_mjjaso[0:],vwind_mjjaso[0:])
vwind_djf = scp.LinReg(sam_ndjfma[0:],vwind_ndjfma[0:])

sic_jja = scp.LinReg(sam_mjjaso[0:],sic_mjjaso[0:])
sic_djf = scp.LinReg(sam_ndjfma[0:],sic_ndjfma[0:])


#mask smaller and insignificant values
def masking(px):
    M =((abs(px)<.1))
    px_mask = ma.array(px,mask = M)
    
    #apply mask
    px_mask[px_mask.mask] = np.nan
    return(px_mask,)
vwind_djf.slope=masking(vwind_djf.slope)
wind_djf.slope=masking(wind_djf.slope)
vwind_jja.slope=masking(vwind_jja.slope)
wind_jja.slope=masking(wind_jja.slope)

def masking(px):
    M =((abs(px)<.01))
    px_mask = ma.array(px,mask = M)
    
    #apply mask
    px_mask[px_mask.mask] = np.nan
    return(px_mask,)

sic_djf.slope=sic_djf.slope.where(sic_djf.p_value<.1)
sic_djf.slope=masking(sic_djf.slope)
sic_jja.slope=sic_jja.slope.where(sic_jja.p_value<.1)
sic_jja.slope=masking(sic_jja.slope)

precip_jja.slope=precip_jja.slope*(10**6)
precip_djf.slope=precip_djf.slope*(10**6)

#importing colorbars using cmap library

import matplotlib as mpl
import cmaps
cmap_temp_ncl = cmaps.precip_diff_12lev
cmap_prec_ncl = cmaps.MPL_BrBG 
cmap_temp_ncl = cmap_temp_ncl(np.arange(cmap_temp_ncl.N))
locs = [5,6,7]

for n in locs:
    cmap_temp_ncl[n] = np.array([1]) # White -> 125, 126, 127, 128, 129
cmap_temp_ncl = mpl.colors.ListedColormap(cmap_temp_ncl)

cmap_temp_ncl1 = cmaps.BlueWhiteOrangeRed
cmap_prec_ncl1 = cmaps.MPL_BrBG 
cmap_temp_ncl1 = cmap_temp_ncl1(np.arange(cmap_temp_ncl1.N))
locs = [120,121,122]

for n in locs:
    cmap_temp_ncl1[n] = np.array([1]) # White -> 125, 126, 127, 128, 129
cmap_temp_ncl1 = mpl.colors.ListedColormap(cmap_temp_ncl1)
cmap_temp_ncl1

cmap_temp_ncl2 = cmaps.GreenMagenta16
cmap_prec_ncl2 = cmaps.MPL_BrBG 
cmap_temp_ncl2 = cmap_temp_ncl2(np.arange(cmap_temp_ncl2.N))
locs = [7,8]

for n in locs:
    cmap_temp_ncl2[n] = np.array([1]) # White -> 125, 126, 127, 128, 129
cmap_temp_ncl2 = mpl.colors.ListedColormap(cmap_temp_ncl2)
cmap_temp_ncl2

import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import numpy as np
import matplotlib.path as mpath
import cartopy.crs as ccrs
from matplotlib.ticker import MaxNLocator, FormatStrFormatter

# Enable LaTeX rendering for all text in plots
#plt.rcParams['text.usetex'] = True

# Set LaTeX preamble to include AMS math package and mathpazo font with old-style figures and small caps
#plt.rcParams['text.latex.preamble'] = "\\usepackage{amsmath} \\usepackage[osf,sc]{mathpazo}"

# Set the default font family to sans-serif
#mpl.rcParams['font.family'] = 'sans-serif'
# Specify Helvetica as the default sans-serif font
#plt.rcParams['font.sans-serif'] = ['Helvetica']


# Load the dataset containing surface air temperature anomalies
dset = xr.open_dataset(data_path+'SAT_anomaly.nc')

# Select the data variable automatically detected in xarray dataset
temp = dset.__xarray_dataarray_variable__.sel(latitude=slice(-30, -90))  # Select southern latitudes

# Get longitude and latitude coordinates
lon, lat = temp.longitude, temp.latitude

# Create a circular boundary path for plotting polar stereographic projection
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

# Initialize figure with specific size and DPI
fig = plt.figure(figsize=[8, 13], dpi=100)

# First subplot: NDJFMA surface air temperature trend on South Polar Stereographic projection
ax = fig.add_subplot(321, projection=ccrs.SouthPolarStereo(central_longitude=0))
m = ax.scontourf(lon, lat, temp_djf.slope.where(temp_djf.p_value < 0.1).sel(latitude=slice(-30, -90)),
                 transform=ccrs.PlateCarree(), levels=np.arange(-.8, .9, .1),
                 cmap=cmaps.BlueWhiteOrangeRed, extend='both')
ax.set_extent([-180, 180, -90, -30], crs=ccrs.PlateCarree())
ax.set_boundary(circle, transform=ax.transAxes)  # Apply circular boundary to plot area
ax.coastlines()
title_australia_ndjfma = 'a  NDJFMA, surface air temperature'  
ax.set_title(title_australia_ndjfma, fontsize=10, fontweight='bold')

# Second subplot: MJJASO surface air temperature trend
ax1 = fig.add_subplot(322, projection=ccrs.SouthPolarStereo(central_longitude=0))
m = ax1.scontourf(lon, lat, temp_jja.slope.where(temp_jja.p_value < 0.1).sel(latitude=slice(-30, -90)),
                  transform=ccrs.PlateCarree(), levels=np.arange(-.8, .9, .1),
                  cmap=cmaps.BlueWhiteOrangeRed, extend='both')
ax1.set_extent([-180, 180, -90, -30], crs=ccrs.PlateCarree())
ax1.set_boundary(circle, transform=ax1.transAxes)
ax1.coastlines()
title_australia_ndjfma = 'b  MJJASO, surface air temperature'
ax1.set_title(title_australia_ndjfma, fontsize=10, fontweight='bold')

# Add colorbar for surface air temperature plots
cb_ax = fig.add_axes([0.93, 0.70, 0.009, 0.15])
cbar1 = fig.colorbar(m, cax=cb_ax, orientation="vertical", shrink=0.5, cmap=cmaps.BlueWhiteOrangeRed)
locator = MaxNLocator(nbins=4, integer=True)  # Locator for tick positions with integer constraint
cbar1.ax.xaxis.set_major_locator(locator)
formatter = FormatStrFormatter('%0.1f')  # Format tick label with one decimal place
cbar1.ax.xaxis.set_major_formatter(formatter)
cbar1.ax.set_ylabel('Surface air temperature ($^\circ$C)', rotation=270, fontsize=12, labelpad=15)  # Label with degree symbol
cbar1.ax.tick_params(labelsize=11)  # Colorbar tick label font size


# Third subplot: NDJFMA precipitation trend with wind vectors overlaid
ax2 = fig.add_subplot(323, projection=ccrs.SouthPolarStereo(central_longitude=0))
m1 = ax2.scontourf(lon, lat, precip_djf.slope.where(precip_djf.p_value < 0.1).sel(latitude=slice(-30, -90)),
                   transform=ccrs.PlateCarree(), levels=np.arange(-4, 5, 1),
                   cmap=cmaps.precip_diff_12lev, extend='both')
ax2.set_extent([-180, 180, -90, -30], crs=ccrs.PlateCarree())

# Create meshgrid for wind data coordinates for plotting arrows
lons, lats = np.meshgrid(wind.longitude, wind.latitude[:])
x, y = (lons, lats)

quiverscale = 20
quiverwidth = 0.003
qstep = 25

# Plot quiver arrows for wind components with scaling and width parameters
q = ax2.quiver(x[::qstep, ::qstep], y[::qstep, ::qstep],
               wind_djf.slope[0][::qstep, ::qstep], vwind_djf.slope[0][::qstep, ::qstep],
               transform=ccrs.PlateCarree(), scale=quiverscale, width=quiverwidth)

# Add quiver key (reference arrow) with label above arrow showing wind speed scale
ax2.quiverkey(q, X=0.93, Y=0.89, U=2,
              label=r'2 ms$^{-1}$', labelpos='N', labelsep=0.05)

ax2.set_boundary(circle, transform=ax2.transAxes)
ax2.coastlines()
title_australia_ndjfma = r'c  NDJFMA, precipitation'
ax2.set_title(title_australia_ndjfma, fontsize=10, fontweight='bold')


# Fourth subplot: MJJASO precipitation with wind vectors
ax3 = fig.add_subplot(324, projection=ccrs.SouthPolarStereo(central_longitude=0))
m1 = ax3.scontourf(lon, lat, precip_jja.slope.where(precip_jja.p_value < 0.1).sel(latitude=slice(-30, -90)),
                   transform=ccrs.PlateCarree(), levels=np.arange(-4, 5, 1),
                   cmap=cmaps.precip_diff_12lev, extend='both')
ax3.set_extent([-180, 180, -90, -30], crs=ccrs.PlateCarree())

lons, lats = np.meshgrid(wind.longitude, wind.latitude[:])
x, y = (lons, lats)

q = ax3.quiver(x[::qstep, ::qstep], y[::qstep, ::qstep],
               wind_jja.slope[0][::qstep, ::qstep], vwind_jja.slope[0][::qstep, ::qstep],
               transform=ccrs.PlateCarree(), scale=quiverscale, width=quiverwidth)

ax3.quiverkey(q, X=.93, Y=0.89, U=2,
              label=r'2 ms$^{-1}$', labelpos='N', labelsep=0.05)

ax3.set_boundary(circle, transform=ax3.transAxes)
ax3.coastlines()
title_australia_ndjfma = r'd  MJJASO, precipitation'
ax3.set_title(title_australia_ndjfma, fontsize=10, fontweight='bold')

# Add colorbar for precipitation plots
cc_ax = fig.add_axes([.93, 0.43, 0.009, 0.15])
cbar = fig.colorbar(m1, cax=cc_ax, orientation="vertical", shrink=0.5, cmap=cmaps.precip_diff_12lev)
locator = MaxNLocator(integer=True)
cbar.ax.xaxis.set_major_locator(locator)
formatter = FormatStrFormatter('%0.1f')
cbar.ax.set_ylabel(r'Precipitation (mm day$^{-1}$ $\times 10^{-6}$)',
                   rotation=270, fontsize=12, labelpad=16)
cbar.ax.yaxis.set_label_coords(8, 0.4)  # Fine-tune label position (x, y)
cbar.ax.tick_params(labelsize=11)


# Load sea ice concentration dataset
#dset = xr.open_dataset('/g/data/v45/SAMworkshop2024/NSIDC/seaice_conc_monthly_sh_197811_202312_v04r00.nc')

#sea_ice_extent = dset.cdr_seaice_conc_monthly
data_proj = ccrs.epsg(dset_sie.projection.attrs["srid"].split("::")[-1])
xgrid = dset_sie.xgrid.values
ygrid = dset_sie.ygrid.values

# Fifth subplot: NDJFMA sea ice concentration
ax4 = fig.add_subplot(325, projection=ccrs.SouthPolarStereo(central_longitude=0))
l = ax4.scontourf(xgrid, ygrid, sic_djf.slope[0][:,:]*100, transform=data_proj,
                  levels=np.arange(-8, 9, 1), cmap=cmap_temp_ncl2)
ax4.set_extent([-180, 180, -50, -90], ccrs.PlateCarree())
ax4.coastlines()
ax4.set_boundary(circle, transform=ax4.transAxes)
ax4.coastlines()
title_australia_ndjfma = r'e NDJFMA, sea ice concentration'
ax4.set_title(title_australia_ndjfma, fontsize=10, fontweight='bold')

# Sixth subplot: MJJASO sea ice concentration
ax5 = fig.add_subplot(326, projection=ccrs.SouthPolarStereo(central_longitude=0))
l = ax5.scontourf(xgrid, ygrid, sic_jja.slope[0][:,:]*100, transform=data_proj,
                  levels=np.arange(-8, 9, 1), cmap=cmap_temp_ncl2)
ax5.set_extent([-180, 180, -50, -90], ccrs.PlateCarree())
ax5.coastlines()
ax5.set_boundary(circle, transform=ax5.transAxes)
ax5.coastlines()
title_australia_ndjfma = r'f  MJJASO, sea ice concentration'
ax5.set_title(title_australia_ndjfma, fontsize=10, fontweight='bold')

# Add colorbar for sea ice concentration
cc_ax = fig.add_axes([.93, 0.16, 0.009, 0.15])
cbar2 = fig.colorbar(l, cax=cc_ax, orientation="vertical", shrink=0.4, cmap=cmap_temp_ncl2)
locator = MaxNLocator(integer=True)
cbar2.ax.yaxis.set_major_locator(locator)
formatter = FormatStrFormatter('%0.1f')  # Format tick labels
cbar2.ax.xaxis.set_major_formatter(formatter)
cbar2.ax.get_xaxis().labelpad = 2
cbar2.ax.set_ylabel('Sea ice concentration (%)', rotation=270, fontsize=12, labelpad=16)
cbar2.ax.tick_params(labelsize=11)

plt.savefig('Fig6.pdf',bbox_inches='tight')
