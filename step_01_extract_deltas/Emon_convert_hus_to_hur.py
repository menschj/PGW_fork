#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     Convert GCM specific humidity to relative humidity.
author			Christoph Heim
date created    08.02.2021
usage           no args
"""
###############################################################################
import os, argparse
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
###############################################################################

def specific_to_relative_humidity(QV, P, T):
    """
    Compute relative humidity from specific humidity.
    """
    RH = 0.263 * P * QV *(np.exp(17.67*(T - 273.15)/(T-29.65)))**(-1)
    return(RH)


###############################################################################
###############################################################################
###############################################################################


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 
                'Convert GCM specific humidity to relative humidity')
    # input specific humidity file
    parser.add_argument('hus_file', type=str)
    # input temperature file
    parser.add_argument('ta_file', type=str)
    # output relative humidity file
    parser.add_argument('hur_file', type=str)
    # Amon relative humidity file
    parser.add_argument('-a', '--amon_hur_file', type=None)
    args = parser.parse_args()

    #xr.set_options(keep_attrs=True)

    # temperature
    var_name = 'ta'
    ta = xr.open_dataset(args.ta_file, decode_cf=False).ta
    # specific humidity
    var_name = 'hus'
    ds = xr.open_dataset(args.hus_file, decode_cf=False)
    hus = ds.hus
    # pressure
    pa = ds.plev.expand_dims(
            dim={'lon':ds.lon, 'lat':ds.lat, 'time':ds.time})
    pa = pa.transpose('time','plev','lat','lon')

    if hus.shape != ta.shape:
        print(hus.shape)
        print(ta.shape)
        raise ValueError()

    hur = specific_to_relative_humidity(hus, pa, ta)

    ## If hur is given from the coarser dataset Amon
    ## it is here used to interpolate hur to the higher
    ## resolution using information from the high-resolved hur.
    ## The reason is that the high-resolved hur is computed
    ## based on monthly mean values and deviates a lot from 
    ## the true hur in the coarse data set on the coarse levels.
    ## Nevertheless, it contains information about the vertical
    ## variability and this is exploited here for a better
    ## informed vertical interpolation of the coarse hur
    ## to high resolution.
    ## This means that the above computed high-resolved hur
    ## is only indirectly used for the final hur output.
    amon_hur = xr.open_dataset(args.amon_hur_file, decode_cf=False).hur

    hur_interp = hur.copy()
    #print(amon_hur.plev.values)
    #print(hur.plev.values)

    for plev in hur.plev.values:
        if plev not in amon_hur.plev.values:
            print('{}: interpolate'.format(plev))
            plev_below = amon_hur.plev.where((amon_hur.plev-plev) > 0, np.nan)
            plev_below = amon_hur.plev.isel(plev=plev_below.argmin(dim='plev').values).values

            plev_above = amon_hur.plev.where((amon_hur.plev-plev) < 0, np.nan)
            plev_above = amon_hur.plev.isel(plev=plev_above.argmax(dim='plev').values).values

            #print(plev_above)
            #print(plev_below)

            hur_plev = hur.sel(plev=plev)
            hur_above = hur.sel(plev=plev_above)
            hur_below = hur.sel(plev=plev_below)

            #print(hur_above.isel(time=0,lon=10,lat=10).values)
            #print(hur_plev.isel(time=0,lon=10,lat=10).values)
            #print(hur_below.isel(time=0,lon=10,lat=10).values)

            weight_above = 1 - np.abs(hur_plev - hur_above)/(
                                np.abs(hur_plev - hur_above) +
                                np.abs(hur_plev - hur_below))
            weight_below = 1 - np.abs(hur_plev - hur_below)/(
                                np.abs(hur_plev - hur_above) +
                                np.abs(hur_plev - hur_below))

            #print(weight_above.isel(time=0,lon=10,lat=10).values)
            #print(weight_below.isel(time=0,lon=10,lat=10).values)

            interp = (
                amon_hur.sel(plev=plev_above) * weight_above +
                amon_hur.sel(plev=plev_below) * weight_below 
            )
            #print(interp.isel(time=0,lon=10,lat=10).values)
            #print()
            hur_interp.loc[dict(plev=plev)] = interp
            #quit()
        else:
            print('{}: take from Amon'.format(plev))
            hur_interp.loc[dict(plev=plev)] = amon_hur.sel(plev=plev)
    #quit()


    handles = []
    handle, = plt.plot(hur.mean(dim=['time','lon','lat']),
             hur.plev, label='Emon hur=f(hus,ta)')
    handles.append(handle)
    handle, = plt.plot(amon_hur.mean(dim=['time','lon','lat']),
             amon_hur.plev, label='Amon hur')
    handles.append(handle)
    handle, = plt.plot(hur_interp.mean(dim=['time','lon','lat']),
             hur_interp.plev, label='Amon hur interpolated using Emon hur')
    handles.append(handle)
    plt.legend(handles=handles)
    plt.ylim((100000,5000))
    plt.ylabel('p [Pa]')
    plt.xlabel('RH [%]')
    plt.show()
    ##quit()

    ds_out = ds.copy()
    ds_out['hur'] = hur_interp
    del ds_out['hus']
    ds_out.attrs['variable_id'] = 'hur'

    ## make sure to keep time encoding identical
    for key,val in ds.time.attrs.items():
        ds_out.time.attrs[key] = val
    for key,val in ds.lon.attrs.items():
        ds_out.lon.attrs[key] = val
    for key,val in ds.lat.attrs.items():
        ds_out.lat.attrs[key] = val
    for key,val in ds.hus.attrs.items():
        if key == 'standard_name':
            ds_out.hur.attrs[key] = 'relative_humidity'
        if key == 'long_name':
            ds_out.hur.attrs[key] = 'Relative Humidity'
        else:
            ds_out.hur.attrs[key] = val


    ds_out.to_netcdf(args.hur_file)

