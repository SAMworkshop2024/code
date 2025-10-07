import xarray as xr

datadir = '../data/'

lat_lim = -40

# Select months of interest for NDJFMA and MJJASO
months_ndjfma = [11, 12, 1, 2, 3, 4]
months_mjjaso = [5, 6, 7, 8, 9, 10]


years=xr.date_range(start='1979-01-01', end='2023-12-31', freq='YS',use_cftime=True)


dset=xr.open_dataset(datadir+'u850_anomaly.nc')
dset = dset.sel(latitude = slice(lat_lim,-90))
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
dset = dset.sel(latitude = slice(lat_lim,-90))
lon,lat = dset.longitude,dset.latitude
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
import cartopy.crs as ccrs

wind_jja = scp.LinReg(sam_mjjaso[0:],wind_mjjaso[0:])
wind_djf = scp.LinReg(sam_ndjfma[0:],wind_ndjfma[0:])

vwind_jja = scp.LinReg(sam_mjjaso[0:],vwind_mjjaso[0:])
vwind_djf = scp.LinReg(sam_ndjfma[0:],vwind_ndjfma[0:])


#import numpy as np
#from numpy import ma  #mask smaller values
#def masking(px):
#    M =((abs(px)<.05))
#    px_mask = ma.array(px,mask = M)
#    
#    #apply mask: very inefficient way
#    px_mask[px_mask.mask] = np.nan
#    return (px_mask,)

# mj different way - currently doesn't work when trying to apply .where() below
def masking(px):
    M = ((abs(px)>=0.05))
    return px.where(M)

#vwind_djf.slope=masking(vwind_djf.slope)
#wind_djf.slope=masking(wind_djf.slope)
#vwind_jja.slope=masking(vwind_jja.slope)
#wind_jja.slope=masking(wind_jja.slope)

wind_djf.slope.name = 'u'
vwind_djf.slope.name = 'v'
uv_djf = xr.merge([wind_djf.slope,vwind_djf.slope])

wind_jja.slope.name = 'u'
vwind_jja.slope.name = 'v'
uv_jja = xr.merge([wind_jja.slope,vwind_jja.slope])


import cartopy.crs as ccrs

proj = ccrs.SouthPolarStereo()

quivArgs = {
    'x' : 'longitude',
    'y' : 'latitude',
    'u' : 'u',
    'v' : 'v',
    'transform': ccrs.PlateCarree(),
    'add_guide': False,
    'scale':2e1,
    }
    

ii = 30
jj = 15

fig = plt.figure(figsize=[12,6],dpi=100)


ax = fig.add_subplot(121,projection=proj)

sig = ( wind_djf.p_value<0.1 ) + ( vwind_djf.p_value<0.1 )

uv_djf.where(sig).isel(longitude=slice(None,None,ii),latitude=slice(None,None,jj)).plot.quiver(ax=ax,**quivArgs)
ax.coastlines()


ax = fig.add_subplot(122,projection=proj)

sig = ( wind_jja.p_value<0.1 )*( vwind_jja.p_value<0.1 )

uv_jja.where(sig).isel(longitude=slice(None,None,ii),latitude=slice(None,None,jj)).plot.quiver(ax=ax,**quivArgs)
ax.coastlines()

