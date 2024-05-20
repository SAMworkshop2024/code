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
sam = xr.open_dataset(f'{data_dir}/sam_indices_mxj.nc', engine='netcdf4')
# get MSLP from ERA5 montlhy
msl = xr.open_mfdataset(era_base+'single-levels/monthly-averaged/msl/*/*.nc', preprocess=ReadSH)

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

# Plot components

transf = ccrs.PlateCarree()
clevs = 19
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(5*3,5), subplot_kw={'projection': ccrs.SouthPolarStereo()})
fig.suptitle(f'Corr. MSLP with SAM Index ({sam[idx_name].name}), ERA5 monthly')

for var, ax in zip(corrmap.data_vars, axs.flatten()):
    corrmap[var].plot.contourf(ax=ax, transform=transf, levels=clevs)
    ax.set_title(corrmap[var].name)
    ax.coastlines()
    ax.gridlines()
    
## save figure
# # change here path/name
# figFile = 'figure_5.png'
# fig.savefig(figFile, bbox_inches='tight', facecolor='white', transparent=False)
