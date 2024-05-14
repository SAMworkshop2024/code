from aostools import climate as ac
from aostools import inout as ai
import xarray as xr
from matplotlib import pyplot as plt
import numpy as np

era_base = '/g/data/rt52/era5/'

# read SH data only
def ReadSH(ds):
    ds = ac.StandardGrid(ds,rename=True)
    return ds.sel(lat=slice(-90,-20))

sams = []
#####
# first, SAM based on SLP
#####
msl = xr.open_mfdataset(era_base+'single-levels/monthly-averaged/msl/*/*.nc',preprocess=ReadSH)
time = msl.time


# Marshall index
sam_marshall = ac.ClimateIndex(msl.msl,index='sam')
sam_marshall.name = 'marshall'
sam_marshall = sam_marshall.assign_coords(time=time)
sams.append(sam_marshall)

# EOF based
slpa = ac.Anomaly(msl.msl,'time.month')
eof,pc,E,_,_,_ = ac.eof(slpa.values)
# convert np.arrays into xarray.DataArrays
eofx = xr.DataArray(eof,coords=[msl.lat,msl.lon,('mode',np.arange(1,eof.shape[-1]+1))],name='EOF')
pcx = xr.DataArray(pc,coords=[('mode',np.arange(1,eof.shape[-1]+1)),msl.time],name='PC')
# the SAM is standardized PC1
sam_slp = ac.Standardize(pcx.sel(mode=1),'time.month')
sam_slp.name = 'slp_eof'
sam_slp = sam_slp.assign_coords(time=time)
sams.append(sam_slp)

#####
# second, SAM based on Zg
#####

# EOF based
z = xr.open_mfdataset(era_base+'pressure-levels/monthly-averaged/z/*/*.nc',preprocess=ReadSH).z/9.81
am = ac.ComputeAnnularMode(z.lat.values, z.pres.values, z.mean('lon').values)
amx = xr.DataArray(am,coords=[time,z.pres],name='annular_mode')
sams.append(amx)

# Antarctic polar cap
zpc = ac.GlobalAvgXr(z.mean('lon'),lats=[-90,-60])
zpc = ac.Standardize(zpc,'time.month')
zpc.name = 'polar_cap'
zpc = zpc.assign_coords(time=time)
sams.append(-zpc)

#######
# Combine all indices
######
sam = xr.merge(sams).load()

######
# Write them all into a file
######
# compress output variables
enc = ai.DefCompress(sam)
outFile = 'sams.nc'
sam.to_netcdf(outFile,encoding=enc)
print('written file {0}'.format(outFile))

######
# plot all the indices
######
try:
    import seaborn as sns
    sns.set_theme('notebook')
    sns.set_style('whitegrid')
    sns.color_palette('tab10')
except:
    pass
plot_level = 500
fig,ax = plt.subplots()
for var in sam.data_vars:
    tmp = sam[var]
    if 'pres' in tmp.coords:
        tmp = tmp.sel(pres=plot_level)
        vname = var+'_{0}'.format(plot_level)
    else:
        vname = var
    ax.plot(time,tmp,label=vname)
ax.set_ylabel('SAM index')
ax.set_xlabel('time')
ax.legend()
ylims = ax.get_ylim()
ax.set_ylim(min(min(ylims),-max(ylims)),max(max(ylims),-min(ylims)))
try:
    sns.despine(bottom=True)
except:
    ax.grid()
figFile = 'sams.pdf'
fig.savefig(figFile,bbox_inches='tight',transparent=True)
print('saved figure {0}'.format(figFile))



