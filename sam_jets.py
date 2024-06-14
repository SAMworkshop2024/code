import xarray as xr
from aostools import climate as ac
import numpy as np
import scipy.signal as sg
import os

level = 300
proj = 'Robinson'
latslice = {'lat':slice(-90,0)}

threshold = 0.0

plot_seasons = ['DJF','JJA']

detrend = False # makes no difference

def detrend_x(x):
    return ac.Detrend(x,'time','linear',keep_mean=True)


overlap = True

era_base = '/g/data/rt52/era5/'

sam = xr.open_dataarray('/scratch/v45/SAMworkshop2024/data/SAM_GW_1m_1979-2023.nc')
if detrend:
    sam = sam.groupby('time.season').map(detrend_x)

def ReadERA(ds):
    ds = ac.StandardGrid(ds,rename=True)
    return ds.sel(pres=level).sel(latslice)

ufile = 'u_monthly.nc'
if os.path.isfile(ufile):
    print('reading',ufile)
    u = xr.open_dataarray(ufile)
else:
    u = xr.open_mfdataset(era_base+'pressure-levels/monthly-averaged/u/*/*',preprocess=ReadERA)
    u = u.sel(time=slice(sam.time[0],sam.time[-1])).u.load()
    if detrend:
        u = u.groupby('time.season').map(detrend_x)
    u.to_netcdf('u_monthly.nc')
    print(ufile)
    
MAMJJA = (sam['time.month']>=3)*(sam['time.month']<=8)
SONDJF = np.invert(MAMJJA)

DJF = sam['time.season'] == 'DJF'
MAM = sam['time.season'] == 'MAM'
JJA = sam['time.season'] == 'JJA'
SON = sam['time.season'] == 'SON'

SAMneg = sam < -threshold
SAMpos = sam >  threshold

filtr = {'SONDJF, SAM < 0' : SONDJF*SAMneg,
         'SONDJF,_SAM > 0' : SONDJF*SAMpos,
         'MAMJJA,_SAM > 0' : MAMJJA*SAMpos,
         'MAMJJA, SAM < 0' : MAMJJA*SAMneg
        }
nrows = 2

filtr2 = {'DJF, SAM < {0}'.format(-threshold) : DJF*SAMneg,
          'DJF, SAM > {0}'.format(threshold) : DJF*SAMpos,
          'MAM, SAM < {0}'.format(-threshold) : MAM*SAMneg,
          'MAM, SAM > {0}'.format(threshold) : MAM*SAMpos,
          'JJA, SAM < {0}'.format(-threshold) : JJA*SAMneg,
          'JJA, SAM > {0}'.format(threshold) : JJA*SAMpos,
          'SON, SAM < {0}'.format(-threshold) : SON*SAMneg,
          'SON, SAM > {0}'.format(threshold) : SON*SAMpos
          }
nrows = len(plot_seasons)

#if nrows == 4:
fltr = filtr2
#elif nrows == 2:
#    fltr = filtr

#if overlap:
#    nrows = nrows//2
    
ncols = 3
nlevs = 10
fig,axs,transf = ac.Projection(proj,ncols=ncols,nrows=nrows,kw_args={'central_longitude':155})
fig.set_figheight(nrows*1.1)
fig.set_figwidth(ncols*4*1.2)

keys = sorted(list(fltr.keys()))

if overlap:
    seasons = np.unique([k.split(',')[0] for k in keys])
    sams = np.unique([s.split(', ')[-1] for s in keys])
    for a in range(len(plot_seasons)):
        season = seasons[a]
        polarity = 'SAM > {0}'.format(threshold)
        key = '{0}, {1}'.format(season,polarity)
        ax = axs[a][0]
        upos = u.isel(time=fltr[key]).mean('time')
        cl = upos.plot.contourf(levels=nlevs,ax=ax,vmin=10,vmax=40,cmap='Reds',add_colorbar=False,**transf) 
        ax.set_title(key)
        ax.coastlines()  
        key = key.replace('>','<').replace('{0}'.format(threshold),'{0}'.format(-threshold))
        ax = axs[a][1]
        uneg = u.isel(time=fltr[key]).mean('time')
        cl = uneg.plot.contourf(levels=nlevs,ax=ax,vmin=10,vmax=40,cmap='Reds',add_colorbar=False,**transf)
        ax.set_title(key)
        ax.coastlines()  
        ax = axs[a][2]
        cd = (upos-uneg).plot.contourf(levels=nlevs,ax=ax,vmin=-6,cmap='RdBu_r',extend='both',add_colorbar=False,**transf) 
        ax.set_title('{0}, difference'.format(season))
        ax.coastlines()
        # there's a bug where the last axes show the whole globe
        axs[a][2].set_ylim(axs[0][0].get_ylim())
else:
    for a,ax in enumerate(axs.flat):
        key = keys[a]
        season = key.split(',')[0]
        if season not in plot_seasons:
            continue
        cl = u.isel(time=fltr[key]).mean('time').plot.contourf(levels=11,ax=ax,vmin=10,vmax=50,cmap='Reds',add_colorbar=False,**transf)
        ax.set_title(key)
        ax.coastlines()
if overlap:
    ac.AddColorbar(fig,axs,cd,shrink=0.8,cbar_args={'label':'zonal wind difference [ms-1]'})
ac.AddColorbar(fig,axs,cl,shrink=0.8,cbar_args={'label':'zonal wind [ms-1]','location':'left'})
ac.AddPanelLabels(axs,'lower left')
fig.suptitle('Seasonal {0}hPa zonal wind by SAM phase'.format(level))
outFile = 'SAM_u{1}_SAM{0}.pdf'.format(threshold,level)
if detrend:
    outFile = outFile.replace('.pdf','_detrend.pdf')
fig.savefig(outFile,bbox_inches='tight')
print(outFile)
