import numpy as np
import xarray as xr
from aostools import climate as ac
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

data_dir ='../data_minimal/'


# read SH data only
def ReadSH(ds):
    ds = ac.StandardGrid(ds,rename=True)
    return ds.sel(lat=slice(-90,-20))

# read sam index netcdf file
sam = xr.open_dataset(f'{data_dir}SAM_GW_1m_1979-2023.nc', engine='netcdf4')

# get MSLP from ERA5 montlhy
#msl = xr.open_mfdataset(era_base+'single-levels/monthly-averaged/msl/*/*.nc', preprocess=ReadSH)
msl = ReadSH(xr.open_dataset(data_dir+'msl_anomaly.nc'))


# time period to use
period = slice('1979', '2023')
msl = msl.sel(time=period)
sam = sam.sel(time=period)

# change name to gong_wang if name is SAM
if 'SAM' in sam.data_vars:
    sam = sam.rename({'SAM':'gong_wang'})

# choose one index (e.g, Gong and Wang)
# compute correlation between SAM index and SLP
idx_name = 'gong_wang'
corr = xr.corr(sam[idx_name], msl.msl, dim='time').load()
corr.name = 'Total'

# get symmetric/asymmetric components
# asymmetric
corr_asymm = corr - corr.mean('lon')
corr_asymm.name ='Asymmetric'
# symmetric
corr_symm = corr - corr_asymm
corr_symm.name = 'Symmetric'

# Merge all into one Dataset
corrmap = xr.merge([corr, corr_symm, corr_asymm])
corrmap = ac.CloseGlobe(corrmap)

## Plot components

from string import ascii_lowercase
import matplotlib.path as mpath

# to make circular map only (no square box)
theta = np.linspace(0, 2*np.pi, 100)
map_circle = mpath.Path(np.vstack([np.sin(theta), np.cos(theta)]).T * 0.5 + [0.5, 0.5])
vmin1 = -0.6
vmax1 = 0.6
vmin2 = -0.3
vmax2 = 0.3
clevs1 = 16 
clevs2 = 16 
transf = ccrs.PlateCarree()
fs = 18

figArgs = {'transform':transf, 'extend':'both', 'add_colorbar':False}
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(6*3,6), subplot_kw={'projection': ccrs.SouthPolarStereo()})

# titles for panels
titles = ['Total SAM and MSLP corr.', 'Symmetric SAM and MSLP corr.', 'Asymmetric SAM and MSLP corr.']

# total
cf1 = corrmap['Total'].plot.contourf(ax=axs[0], levels=clevs1, vmin=vmin1, **figArgs)
ttl_a = titles[0]
axs[0].set_title(r'$\bf{a\ }$' + f' {ttl_a}', fontsize=fs, loc='center')

# symmetric
cf1 = corrmap['Symmetric'].plot.contourf(ax=axs[1], levels=clevs1, vmin=vmin1, **figArgs)
ttl_b = titles[1]
axs[1].set_title(r'$\bf{b\ }$' + f' {ttl_b}', fontsize=fs, loc='center')

# asymmetric separately
cf3 = corrmap['Asymmetric'].plot.contourf(ax=axs[2], levels=clevs2, vmin=vmin2, **figArgs) 
ttl_c = titles[2]
axs[2].set_title(r'$\bf{c\ }$' + f' {ttl_c}', fontsize=fs, loc='center')

for ax in axs.flatten():
    ax.set_extent([-180, 180, -90, -20], ccrs.PlateCarree())
    ax.set_boundary(map_circle, transform=ax.transAxes)
    ax.coastlines()
    
cbarArgs1 = {'orientation':'horizontal', 'extend':'both', 'fraction':0.06, 'pad':0.05}
cbarArgs2 = {'orientation':'horizontal', 'extend':'both', 'fraction':0.06, 'pad':0.05}

# colorbar 1
cb1 = fig.colorbar(cf1, ax=axs[:2], **cbarArgs1)
#cb1.ax.locator_params(nbins=8)
cb1.ax.tick_params(which='major', labelsize=14, length=7)
ticks1 = [-0.6, -0.36, -0.12, 0.12, 0.36, 0.6]
cb1.set_ticks(ticks1)

# colorbar 2
cb3 = fig.colorbar(cf3, ax=axs[2], **cbarArgs2)
cb3.ax.locator_params(nbins=6)
cb3.ax.tick_params(which='major', labelsize=14, length=7)
ticks3 = [-0.3, -0.18, -0.06, 0.06, 0.18, 0.3]
cb3.set_ticks(ticks3)

# change here path/name
figFile = 'Fig3abc.pdf'
fig.savefig(figFile, bbox_inches='tight', facecolor='white', transparent=False)
print(figFile)
