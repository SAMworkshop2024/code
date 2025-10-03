import numpy as np
import xarray as xr
from aostools import climate as ac
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# read SH data only
def ReadSH(ds):
    ds = ac.StandardGrid(ds,rename=True)
    return ds.sel(lat=slice(-90,-20))

# Workshop data dir
data_dir ='/scratch/v45/SAMworkshop2024/data'
# ERA5 data dir
era_base = '/g/data/rt52/era5/'

# read sam index netcdf file

# mj standardization
#sam = xr.open_dataset(f'{data_dir}/sam_indices_mxj.nc', engine='netcdf4')

# the other one
sam = xr.open_dataset(f'{data_dir}/SAM_GW_1m_1979-2023.nc', engine='netcdf4')

# get MSLP from ERA5 montlhy
msl = xr.open_mfdataset(era_base+'single-levels/monthly-averaged/msl/*/*.nc', preprocess=ReadSH)

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

## Plot components

from string import ascii_lowercase
import matplotlib.path as mpath

# to make circular map only (no square box)
theta = np.linspace(0, 2*np.pi, 100)
map_circle = mpath.Path(np.vstack([np.sin(theta), np.cos(theta)]).T * 0.5 + [0.5, 0.5])

transf = ccrs.PlateCarree()
#clevs = 19
clevs1 = np.linspace(-0.6, 0.6, 19)
clevs2 = np.linspace(-0.3, 0.3, 19)
figArgs = {'transform':transf, 'extend':'both', 'add_colorbar':False}
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(5*3,5), subplot_kw={'projection': ccrs.SouthPolarStereo()})
#fig.suptitle(f'Correlation between MSLP and SAM Index, monthly ERA5 {period.start}-{period.stop}')

# total
cf1 = corrmap['Total'].plot.contourf(ax=axs[0], levels=clevs1, **figArgs)
axs[0].set_title(corrmap['Total'].name)
# symmetric
cf1 = corrmap['Symmetric'].plot.contourf(ax=axs[1], levels=clevs1, **figArgs)
axs[1].set_title(corrmap['Symmetric'].name)
# asymmetric separately
cf3 = corrmap['Asymmetric'].plot.contourf(ax=axs[2], levels=clevs2, **figArgs) 
axs[2].set_title(corrmap['Asymmetric'].name)

for ax in axs.flatten():
    ax.set_extent([-180, 180, -90, -20], ccrs.PlateCarree())
    ax.set_boundary(map_circle, transform=ax.transAxes)
    ax.coastlines()
    ax.coastlines()
    
# Add panel labels    
for a,ax in enumerate(axs):
    ax.text(x=0.01, y=1.06, s=ascii_lowercase[a],weight='bold', fontsize=20,transform=ax.transAxes)
    
cbarArgs = {'orientation':'horizontal', 'extend':'both', 'fraction':0.05, 'pad':0.05}

# colorbar 1
cb1 = fig.colorbar(cf1, ax=axs[:2], **cbarArgs)
cb1.ax.locator_params(nbins=8)

# colorbar 2
cb3 = fig.colorbar(cf3, ax=axs[2], **cbarArgs)
cb3.ax.locator_params(nbins=8)

plt.show()


# # change here path/name
# figFile = 'Fig3.png'
# fig.savefig(figFile, bbox_inches='tight', facecolor='white', transparent=False)
