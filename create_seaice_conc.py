import xarray as xr

datadir = '../data/'

dset = xr.open_dataset(datadir+'seaice_conc_monthly_sh_197811_202312_v04r00.nc')

projx = dset.projection
sice = xr.DataArray(dset.cdr_seaice_conc_monthly.values,coords=[dset.time,dset.ygrid,dset.xgrid],name='cdr_seaice_conc_monthly')

data = xr.merge([sice,projx]).rename({'x':'xgrid','y':'ygrid'})

data.to_netcdf(datadir+'seaice_conc.nc')
