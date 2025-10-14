import xarray as xr

datadir = '../data/'

dset=xr.open_dataset(datadir+'SAT_anomaly.nc')
temp=dset.__xarray_dataarray_variable__

years=xr.cftime_range(start='1979-01-01', end='2023-12-31', freq='Y')

# Select months of interest for NDJFMA and MJJASO
months_ndjfma = [11, 12, 1, 2, 3, 4]
months_mjjaso = [5, 6, 7, 8, 9, 10]

# Select data for NDJFMA and MJJASO
temp_ndjfma1 = temp.sel(time=temp['time.month'].isin(months_ndjfma))
temp_mjjaso1 = temp.sel(time=temp['time.month'].isin(months_mjjaso))

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

dset=xr.open_dataset(datadir+'precip_anomaly.nc')
precip=dset.__xarray_dataarray_variable__

years=xr.cftime_range(start='1979-01-01', end='2023-12-31', freq='Y')

# Select months of interest for NDJFMA and MJJASO
months_ndjfma = [11, 12, 1, 2, 3, 4]
months_mjjaso = [5, 6, 7, 8, 9, 10]

# Select data for NDJFMA and MJJASO
precip_ndjfma1 = precip.sel(time=precip['time.month'].isin(months_ndjfma))
precip_mjjaso1 = precip.sel(time=precip['time.month'].isin(months_mjjaso))

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

import cartopy.crs as ccrs

dset=xr.open_dataset('/scratch/w97/rr8442/seaice_conc.nc')
sea_ice_extent = dset.cdr_seaice_conc_monthly.sel(tdim=slice('1979-01-01','2023-12-01')).values 
data_proj = ccrs.epsg(dset.projection.attrs["srid"].split("::")[-1])
xgrid = dset.xgrid.values
ygrid = dset.ygrid.values

sea_ice_extent=xr.DataArray(sea_ice_extent)
sea_ice_extent['dim_0']=dset.sel(tdim=slice('1979-01-01','2023-12-01')).tdim.values
months_of_interest = [11,12,1,2,3,4] 
sic_ndjfma1 = sea_ice_extent.sel(dim_0=sea_ice_extent['dim_0.month'].isin(months_of_interest))
months_of_interest = [5,6,7,8,9,10] 
sic_mjjaso1 = sea_ice_extent.sel(dim_0=sea_ice_extent['dim_0.month'].isin(months_of_interest))

sic_ndjfma=sic_ndjfma1[0]
for i in range(len(years.year)):
    sic = sic_ndjfma1.sel(dim_0=slice(str(years.year[i])+'-11-01', str(years.year[i]+1)+'-04-30')).mean('dim_0')
    sic_ndjfma =xr.concat([sic_ndjfma, sic],dim='years')

sic_ndjfma=sic_ndjfma[1:].assign_coords(years=years)

sic_mjjaso=sic_mjjaso1[0]
for i in range(len(years.year)):
    sic = sic_mjjaso1.sel(dim_0=slice(str(years.year[i])+'-05-01', str(years.year[i])+'-10-30')).mean('dim_0')
    sic_mjjaso =xr.concat([sic_mjjaso, sic],dim='years')

sic_mjjaso=sic_mjjaso[1:].assign_coords(years=years)

dset=xr.open_dataset(datadir+'u850_anomaly.nc')
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

dset=xr.open_dataset(datadir+'v850_anomaly.nc')
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

dset=xr.open_dataset(datadir+'SAM_GW_1m_1979-2023.nc')
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

import numpy as np
import sacpy as scp
import matplotlib.pyplot as plt
import sacpy.Map # need cartopy or you can just not import
import cartopy.crs as ccrs

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

import numpy as np
from numpy import ma  #mask smaller values
def masking(px):
    M =((abs(px)<.05))
    px_mask = ma.array(px,mask = M)
    
    #apply mask: very inefficient way
    px_mask[px_mask.mask] = np.nan
    return(px_mask,)

vwind_djf.slope=masking(vwind_djf.slope)
wind_djf.slope=masking(wind_djf.slope)
vwind_jja.slope=masking(vwind_jja.slope)
wind_jja.slope=masking(wind_jja.slope)

import cartopy.crs as ccrs
import sacpy.Map
from matplotlib import rc
rc('text', usetex=True)
plt.rc('axes', linewidth=1.5)
plt.rc('font', weight='bold',size=15)
rc('font',**{'family':'serif','serif':['Times']})


dset=xr.open_dataset(datadir+'SAT_anomaly.nc')
temp=dset.__xarray_dataarray_variable__

lon ,lat = temp.longitude , temp.latitude
fig = plt.figure(figsize=[13,13],dpi=100)
ax = fig.add_subplot(321,projection=ccrs.PlateCarree())
m = ax.scontourf(lon,lat,temp_djf.slope.where(temp_djf.p_value<0.1),cmap=cmaps.BlueWhiteOrangeRed ,extend='both')
ax.init_map(stepx=50,smally=2.5)
ax.set_title('NDJFMA',)

ax1 = fig.add_subplot(322,projection=ccrs.PlateCarree())
m = ax1.scontourf(lon,lat,temp_jja.slope.where(temp_jja.p_value<0.1),cmap=cmaps.BlueWhiteOrangeRed ,extend='both',)
ax1.init_map(stepx=50,smally=2.5)
ax1.set_title('MJJASO',)

ax2 = fig.add_subplot(323,projection=ccrs.PlateCarree())
m1 = ax2.scontourf(lon,lat,precip_djf.slope.where(precip_djf.p_value<0.1),cmap=cmaps.precip_diff_12lev,extend='both')
ax2.squiver(lon,lat,wind_djf.slope.where(wind_djf.p_value<0.1),vwind_djf.slope.where(vwind_djf.p_value<0.1),stepx=30,stepy=15,headwidth=2,headaxislength=8,width=.0055,minlength=2)
ax2.coastlines()
ax2.init_map(smally=2.5)
#stepx=40,stepy=20,headwidth=2,headaxislength=3,width=.006,minlength=2
ax3 = fig.add_subplot(324,projection=ccrs.PlateCarree())
m1 = ax3.scontourf(lon,lat,precip_jja.slope.where(precip_jja.p_value<0.1),cmap=cmaps.precip_diff_12lev,extend='both')
ax3.squiver(lon,lat,wind_jja.slope.where(wind_jja.p_value<0.1),vwind_jja.slope.where(vwind_jja.p_value<0.1),stepx=30,stepy=15,headwidth=2,headaxislength=8,width=.0055,minlength=2)
ax3.coastlines()
ax3.init_map(smally=2.5)

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

levels=np.array([-0.04,  -0.03,  -0.02, -0.01, 
        0.   , 0.01,  0.02,  0.03, 0.04, ])
levels = levels
#MaxNLocator(nbins=15).tick_values(np.nanmin(sic_jja.slope[0]), np.nanmax(sic_jja.slope[0]))
cmap=cmaps.BlueWhiteOrangeRed
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)


ax4 = fig.add_subplot(325,projection=ccrs.Orthographic(180, -90))
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], .5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax4.set_boundary(circle, transform=ax4.transAxes)
ax4.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
ax4.add_feature(cfeature.COASTLINE, edgecolor='black')
m2=ax4.pcolormesh(xgrid, ygrid, sic_djf.slope.where(sic_djf.p_value<0.1),transform=data_proj,norm=norm,cmap=cmaps.BlueWhiteOrangeRed,)
ax4.coastlines()
center_new, radius_new = [0.5, 0.5], 0.60  # New center and radius
circle_new = mpath.Path(verts * radius_new + center_new)
ax4.set_boundary(circle_new, transform=ax4.transAxes)
gl = ax4.gridlines(draw_labels=False, xlocs=None, ylocs=None)


ax5 = fig.add_subplot(326,projection=ccrs.Orthographic(180, -90))
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], .5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax5.set_boundary(circle, transform=ax5.transAxes)
ax5.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
ax5.add_feature(cfeature.COASTLINE, edgecolor='black')
m2=ax5.pcolormesh(xgrid, ygrid, sic_jja.slope.where(sic_jja.p_value<0.1),transform=data_proj,norm=norm,cmap=cmaps.BlueWhiteOrangeRed,)
ax5.coastlines()
center_new, radius_new = [0.5, 0.5], 0.60  # New center and radius
circle_new = mpath.Path(verts * radius_new + center_new)
ax5.set_boundary(circle_new, transform=ax5.transAxes)
gl = ax5.gridlines(draw_labels=False, xlocs=None, ylocs=None)



cc_ax = fig.add_axes([.945, 0.40, 0.01, 0.2])
cbar=fig.colorbar(m1,cax=cc_ax,orientation="vertical")
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('Precip (mm/day)', rotation=90)

cb_ax = fig.add_axes([.945, 0.68, 0.01, 0.2])
cbar1=fig.colorbar(m,cax=cb_ax,orientation="vertical",)
cbar1.ax.get_xaxis().labelpad = 15
cbar1.ax.set_ylabel('TA ($^\circ$C) ', rotation=90)

cb_ax = fig.add_axes([.945, 0.15, 0.01, 0.2])
cbar2=fig.colorbar(m2,cax=cb_ax,orientation="vertical",extend='both')
cbar2.ax.get_xaxis().labelpad = 15
cbar2.ax.set_ylabel('Sea Ice Concentration ', rotation=90)
