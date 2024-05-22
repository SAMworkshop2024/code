import xarray as xr
from aostools import climate as ac
import numpy as np

level = 300
proj = 'Robinson'
latslice = {'lat':slice(-90,20)}
overlap = True

plot_seasons = ['DJF','JJA']

era_base = '/g/data/rt52/era5/'

sam = xr.open_dataarray('/scratch/v45/SAMworkshop2024/data/SAM_GW_1m_1979-2023.nc')

def ReadERA(ds):
    ds = ac.StandardGrid(ds,rename=True)
    return ds.sel(pres=level).sel(latslice)

u = xr.open_mfdataset(era_base+'pressure-levels/monthly-averaged/u/*/*',preprocess=ReadERA)
u = u.sel(time=slice(sam.time[0],sam.time[-1])).u.load()

MAMJJA = (sam['time.month']>=3)*(sam['time.month']<=8)
SONDJF = np.invert(MAMJJA)

DJF = sam['time.season'] == 'DJF'
MAM = sam['time.season'] == 'MAM'
JJA = sam['time.season'] == 'JJA'
SON = sam['time.season'] == 'SON'

SAMneg = sam < -0.2
SAMpos = sam >  0.2

filtr = {'SONDJF, SAM < 0' : SONDJF*SAMneg,
         'SONDJF,_SAM > 0' : SONDJF*SAMpos,
         'MAMJJA,_SAM > 0' : MAMJJA*SAMpos,
         'MAMJJA, SAM < 0' : MAMJJA*SAMneg
        }
nrows = 2

filtr2 = {'DJF, SAM < 0' : DJF*SAMneg,
          'DJF, SAM > 0' : DJF*SAMpos,
          'MAM, SAM < 0' : MAM*SAMneg,
          'MAM, SAM > 0' : MAM*SAMpos,
          'JJA, SAM < 0' : JJA*SAMneg,
          'JJA, SAM > 0' : JJA*SAMpos,
          'SON, SAM < 0' : SON*SAMneg,
          'SON, SAM > 0' : SON*SAMpos
          }
nrows = len(plot_seasons)

#if nrows == 4:
fltr = filtr2
#elif nrows == 2:
#    fltr = filtr

#if overlap:
#    nrows = nrows//2
    
ncols = 3
fig,axs,transf = ac.Projection(proj,ncols=ncols,nrows=nrows,kw_args={'central_longitude':155})
fig.set_figheight(nrows*2)
fig.set_figwidth(ncols*4*1.2)

keys = sorted(list(fltr.keys()))

if overlap:
    seasons = np.unique([k.split(',')[0] for k in keys])
    sams = np.unique([s.split(', ')[-1] for s in keys])
    for a in range(len(plot_seasons)):
        season = seasons[a]
        polarity = 'SAM > 0'
        key = '{0}, {1}'.format(season,polarity)
        ax = axs[a][0]
        upos = u.isel(time=fltr[key]).mean('time')
        cl = upos.plot.contourf(levels=11,ax=ax,vmin=10,vmax=40,cmap='Reds',add_colorbar=False,**transf) 
        ax.set_title(key)
        ax.coastlines()  
        key = key.replace('>','<')
        ax = axs[a][1]
        uneg = u.isel(time=fltr[key]).mean('time')
        cl = uneg.plot.contourf(levels=11,ax=ax,vmin=10,vmax=40,cmap='Reds',add_colorbar=False,**transf)
        ax.set_title(key)
        ax.coastlines()  
        ax = axs[a][2]
        cd = (upos-uneg).plot(levels=13,ax=ax,cmap='RdBu_r',vmin=-6,add_colorbar=False,**transf) 
        ax.set_title('difference')
        ax.coastlines()  
else:
    for a,ax in enumerate(axs.flat):
        key = keys[a]
        season = key.split(',')[0]
        if season not in plot_seasons:
            continue
        cl = u.isel(time=fltr[key]).mean('time').plot(ax=ax,vmin=10,vmax=50,cmap='Reds',add_colorbar=False,**transf)
        ax.set_title(key)
        ax.coastlines()
if overlap:
    ac.AddColorbar(fig,axs,cd,shrink=0.8,cbar_args={'label':'u300 [m/s]'})
ac.AddColorbar(fig,axs,cl,shrink=0.8,cbar_args={'label':'u300 [m/s]','location':'left'})
ac.AddPanelLabels(axs,'lower left')
fig.suptitle('300hPa zonal wind, ERA5 1979-2023')
outFile = 'SAM_u300.pdf'
fig.savefig(outFile,bbox_inches='tight')
print(outFile)
