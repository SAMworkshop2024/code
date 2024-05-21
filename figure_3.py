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
corr.name = 'total'

# get components
# asymmetric
corr_asymm = corr - corr.mean('lon')
corr_asymm.name ='asymmetric'
# symmetric
corr_symm = corr - corr_asymm
corr_symm.name = 'symmetric'

# Merge all into one Dataset
corrmap = xr.merge([corr, corr_symm, corr_asymm])


## Plot components

transf = ccrs.PlateCarree()
#clevs = 19
clevs = np.linspace(-0.64, 0.64, 19)
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(5*3,5), subplot_kw={'projection': ccrs.SouthPolarStereo()})
#fig.suptitle(f'Corr. MSLP with SAM Index ({sam[idx_name].name}), monthly ERA5 {period.start}-{period.stop}')
fig.suptitle(f'Correlation between MSLP and SAM Index, monthly ERA5 {period.start}-{period.stop}')

for var, ax in zip(corrmap.data_vars, axs.flatten()):
    cf = corrmap[var].plot.contourf(ax=ax, transform=transf, levels=clevs, extend='both', add_colorbar=False)
    ax.set_title(corrmap[var].name)
    ax.coastlines()
    ax.gridlines()

cb = fig.colorbar(cf, ax=axs.ravel().tolist(), orientation='horizontal', extend='both', fraction=0.1, 
             shrink=0.7, pad=0.06)
cb.ax.locator_params(nbins=8)

plt.show()

# ## save figure
# # hange here path/name
# figFile = 'Fig3.png'
# fig.savefig(figFile, bbox_inches='tight', facecolor='white', transparent=False)
