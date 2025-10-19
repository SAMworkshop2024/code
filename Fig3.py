import numpy as np
import pandas as pd
import xarray as xr
from aostools import climate as ac
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

data_path ='../data/'

# -----------------------------------------------------------
# 1) CORRELATION OF SAM WITH SLP
# -----------------------------------------------------------

# read SH data only
def ReadSH(ds):
    ds = ac.StandardGrid(ds,rename=True)
    return ds.sel(lat=slice(-90,-20))

# read sam index netcdf file
sam = xr.open_dataset(f'{data_path}/SAM_GW_1m_1979-2023.nc', engine='netcdf4')

# get MSLP from ERA5 montlhy
#msl = xr.open_mfdataset(era_base+'single-levels/monthly-averaged/msl/*/*.nc', preprocess=ReadSH)
msl = ReadSH(xr.open_dataset(data_path+'msl_anomaly.nc'))

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


# -----------------------------------------------------------
# 2) VARIANCE EXPLAINED BY SYMMETRIC/ASYMMETRIC (ELIO'S DATA)
# -----------------------------------------------------------

# sam_variances file created by Elio Campitelli
#dir_elio = '/g/data/v45/ec0044'
sam_variances = pd.read_csv(f'{data_path}sam_variances.csv')

# Pivot the dataframe to have 'month' as index and 'sym' and 'asym' as columns
pivot = sam_variances.pivot(index='month', columns='name', values='value')  # still a df
sam_variances_ds = xr.Dataset.from_dataframe(pivot)

# Prepare data for plotting
months = range(1, 13)
explained_sym = sam_variances_ds.sym * 100
explained_asym = sam_variances_ds.asym * 100


# -----------------------------------------------------------
# 3) CLIMATE DRIVERS TIME SERIES
# -----------------------------------------------------------
# read data of climate drivers time series
df = pd.read_csv(f'{data_path}/drivers_timeseries.csv', index_col=0)

# convert to xarray.DataSet
ds_drivers = df.to_xarray()
# rename dimension
ds_drivers = ds_drivers.rename({'Unnamed: 0': 'time'})


# -----------------------------------------------------------
#   PLOTTING — ALL PANELS IN ONE FIGURE
# -----------------------------------------------------------

import matplotlib.path as mpath

fs = 18
fig = plt.figure(figsize=(18, 18))
gs = fig.add_gridspec(ncols=12, nrows=12, figure=fig)

# titles for top panels
titles = ['Total SAM and MSLP corr.', 'Symmetric SAM and MSLP corr.', 'Asymmetric SAM and MSLP corr.']

# ################ TOP PANEL ################

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
figArgs_abc = {'transform':transf, 'extend':'both', 'add_colorbar':False}

ax1 = fig.add_subplot(gs[:4,:4], projection=ccrs.SouthPolarStereo())
ax2 = fig.add_subplot(gs[:4,4:8], projection=ccrs.SouthPolarStereo())
ax3 = fig.add_subplot(gs[:4,8:], projection=ccrs.SouthPolarStereo())

# total
cf1 = corrmap['Total'].plot.contourf(ax=ax1, levels=clevs1, vmin=vmin1, **figArgs_abc)
ttl_a = titles[0]
ax1.set_title(r'$\bf{a\ }$' + f' {ttl_a}', fontsize=fs, loc='center')

# Symmetric
cf2 = corrmap['Symmetric'].plot.contourf(ax=ax2, levels=clevs1, vmin=vmin1, **figArgs_abc)
ttl_b = titles[1]
ax2.set_title(r'$\bf{b\ }$' + f' {ttl_b}', fontsize=fs, loc='center')

# asymmetric separately
cf3 = corrmap['Asymmetric'].plot.contourf(ax=ax3, levels=clevs2, vmin=vmin2, **figArgs_abc) 
ttl_c = titles[2]
ax3.set_title(r'$\bf{c\ }$' + f' {ttl_c}', fontsize=fs, loc='center')

for ax in [ax1,ax2,ax3]:
    ax.set_extent([-180, 180, -90, -20], ccrs.PlateCarree())
    ax.set_boundary(map_circle, transform=ax.transAxes)
    ax.coastlines()
    
# add colorbars
cbarArgs = {'orientation':'horizontal', 'extend':'both', 'fraction':0.06, 'pad':0.05}
# colorbar 1
cb1 = fig.colorbar(cf1, ax=[ax1, ax2], **cbarArgs)
#cb1.ax.locator_params(nbins=8)
cb1.ax.tick_params(which='major', labelsize=14, length=7)
ticks1 = [-0.6, -0.36, -0.12, 0.12, 0.36, 0.6]
cb1.set_ticks(ticks1)


# colorbar 2
cb3 = fig.colorbar(cf3, ax=ax3, **cbarArgs)
cb3.ax.locator_params(nbins=6)
cb3.ax.tick_params(which='major', labelsize=14, length=7)
ticks3= [-0.3, -0.18, -0.06, 0.06, 0.18, 0.3]
cb3.set_ticks(ticks3)


# # ################ MIDDLE PANEL ################

# Panel d : Barplot
ttl_d = 'Monthly variance explained by the symmetric and asymmetric components of the SAM'
months = range(1, 13)

ax4 = fig.add_subplot(gs[5:8, 1:11])

# Stacked bar plot
ax4.bar(months, explained_sym, label='Symmetric', color='darkorange', alpha=0.5)
ax4.bar(months, explained_asym, bottom=explained_sym, label='Asymmetric', color='teal', alpha=0.5)
# Adding labels and title
ax4.set_xlabel('Month', fontsize=16)
ax4.set_ylabel('Explained Variance (%)', fontsize=16)
ax4.set_xticks(months)
ax4.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax4.legend(loc=(0.96,0.45), fontsize=12)
ax4.set_title(r'$\bf{d\ }$' + f' {ttl_d}', fontsize=18, loc='left')
ax4.tick_params(which='major', labelsize=14, length=7)


# ################ BOTTOM PANEL ################
# Panel e : drivers time series
# title for panel e
ttl_e = 'Correlation and predictability of seasonal-mean SAM with climate drivers'

# legends -> STILL NEED TO BE UPDATED/CORRECTED -> WAIT FOR EUN-PA INPUT
legends = {'SPV':'SH strat. polar vortex',
           'NINO3':'Nino 3 SST index', 
           'EMI':'El Nino-Modoki SST index', 
           'UpperPred': 'Max. pred. from drivers', 
           'ACCESS-S2 LT0': 'ACCESS-S2 forecast skill'}

clrs = ['C0', 'C1', 'C2'] 

# thresholds to add statistical significance
th_drivers = 0.3
th_access = 0.33

FigArgs_e = dict(lw=3, marker='o', markersize=7.)

ax5 = fig.add_subplot(gs[9:, 1:11])

# line 1
var1 = 'SPV'
ax5.plot(ds_drivers.time, ds_drivers[var1], lw=2, label=legends[var1])
ax5.plot(ds_drivers.time, ds_drivers[var1].where(abs(ds_drivers[var1])>th_drivers), c='C0', **FigArgs_e)
# line 2
var2 = 'NINO3'
ax5.plot(ds_drivers.time, ds_drivers[var2], lw=2, label=legends[var2])
ax5.plot(ds_drivers.time, ds_drivers[var2].where(abs(ds_drivers[var2])>th_drivers), c='C1', **FigArgs_e)
# line 3
var3 = 'EMI'
ax5.plot(ds_drivers.time, ds_drivers[var3], lw=2, label=legends[var3])
ax5.plot(ds_drivers.time, ds_drivers[var3].where(abs(ds_drivers[var3])>th_drivers), c='C2', **FigArgs_e)
# line 4
var4 = 'UpperPred'
ax5.plot(ds_drivers.time, ds_drivers[var4], c='0.2', lw=2, label=legends[var4]) 
# line 5
var5 = 'ACCESS-S2 LT0'
ax5.plot(ds_drivers.time, ds_drivers[var5], c='0.2', lw=2,  ls='--', label=legends[var5])
ax5.plot(ds_drivers.time, ds_drivers[var5].where(abs(ds_drivers[var5])>th_access), ls='', marker='o', c='0.2')

# Adding labels, legend, title
ax5.set_xlabel('Season', fontsize=16)
ax5.set_ylabel('Corr. coefficient', fontsize=16)
ax5.legend(loc=(0.96,0.3), fontsize=12)
ax5.set_title(r'$\bf{e\ }$' + f'{ttl_e}', fontsize=18, loc='left')
ax5.axhline(0, color='k', linestyle='--', linewidth=1)
ax5.tick_params(which='major', labelsize=14, length=7)
ax5.grid(which='major', alpha=0.3, linestyle='dashed')

figFile = 'Fig3.pdf'
fig.savefig(figFile, bbox_inches='tight', format='pdf', dpi=400, facecolor='white', transparent=False)

