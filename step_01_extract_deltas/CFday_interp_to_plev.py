#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     Extract CFday data
author			Christoph Heim
date created    21.02.2021
usage           no args
"""
###############################################################################
import os, argparse
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
from functions import interp_logp_4d
from settings import (
    TIME_GCM, LEV_GCM, PLEV_GCM, LON_GCM, LAT_GCM,
)
###############################################################################


###############################################################################
###############################################################################
###############################################################################


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 
                'Interpolate CFday output to pressure levels')
    # variable name 
    parser.add_argument('var_names', type=str)
    # CMIP6 experiment
    parser.add_argument('experiment', type=str)
    args = parser.parse_args()

    inp_dir = '/net/o3/hymet_nobackup/heimc/data/pgw/download/subdomain'
    out_base_dir ='/net/o3/hymet_nobackup/heimc/data/pgw/download/interp_plev' 

    model_name = 'MPI-ESM1-2-HR'
    #var_name = 'ta'
    #experiment = 'ssp585'
    #experiment = 'historical'

    var_names = args.var_names.split(',')
    experiment = args.experiment

    times = {
        'ssp585': [
            '20700101-20741231',
            '20750101-20791231',
            '20800101-20841231',
            '20850101-20891231',
            '20900101-20941231',
            '20950101-20991231',
        ],
        'historical': [
            '19850101-19891231',
            '19900101-19941231',
            '19950101-19991231',
            '20000101-20041231',
            '20050101-20091231',
            '20100101-20141231',
        ],
    }
    for var_name in var_names:
        for time_ind in range(0,len(times[experiment])):
            print(time_ind)

            inp_file_name = '{}_CFday_{}_{}_r1i1p1f1_gn_{}.nc'.format(
                                    var_name, model_name, experiment, 
                                    times[experiment][time_ind])
            out_dir = os.path.join(out_base_dir, model_name)
            Path(out_dir).mkdir(parents=True, exist_ok=True)
            out_file_name = '{}_CFday_{}_{}_r1i1p1f1_gn_{}.nc'.format(
                                    var_name, model_name, experiment, 
                                    times[experiment][time_ind])

            # create input and output file paths
            out_file_path = os.path.join(out_dir, out_file_name)
            inp_file_path = os.path.join(inp_dir, inp_file_name)

            print('Process input file: \n{}\nto output file: \n{}'.format(
                    inp_file_path, out_file_path))

            ds = xr.open_dataset(inp_file_path)

            # sort pressure ascending
            ds = ds.reindex({LEV_GCM:ds['lev'][::-1]})
            # compute pressure on full levels
            source_P = (ds.ap + ds.b * ds.ps).transpose(
                                TIME_GCM, LEV_GCM, LAT_GCM, LON_GCM)
            var = ds[var_name]



            ### Determine target pressure using tropical ocea-only domain
            #############################################################
            #mean_p = p.mean(dim=['lon','lat','time']).values
            #print(np.around(mean_p[30:], -2))
            #p_integ = np.arange(101000, 100000, -1000)
            #p_integ = np.append(p_integ, np.arange(100000, 80000, -1000))
            #p_integ = np.append(p_integ, np.arange(80000, 30000, -2500))
            #p_integ = np.append(p_integ, np.arange(30000, 20000, -2000))
            #p_integ = np.append(p_integ, np.arange(20000, 10000, -1000))
            #p_integ = np.append(p_integ, mean_p[30:-22])
            #print(p_integ)
            ##plt.scatter(np.arange(len(p.lev.values)),
            ##            p.mean(dim=['lon','lat','time']).values)
            ##plt.show()
            #quit()

            # load target pressure levels
            targ_plev = np.sort(np.loadtxt('CFday_target_p_MPI-ESM1-2-HR.dat'))
            targ_P = xr.DataArray(targ_plev, dims=[PLEV_GCM])

            # create 4d target pressure data array
            targ_P = targ_P.expand_dims(
                            dim={LON_GCM:var[LON_GCM],
                                 LAT_GCM:var[LAT_GCM],
                                 TIME_GCM:var[TIME_GCM]}).transpose(
                                        TIME_GCM, PLEV_GCM, LAT_GCM, LON_GCM)

            # run interpolation from GCM model levels to constant pressure levels
            var_out = interp_logp_4d(var, source_P, targ_P, extrapolate='constant',
                                    time_key=TIME_GCM, lat_key=LAT_GCM,
                                    lon_key=LON_GCM)

            # set pressure levels as coordinate
            var_out = var_out.assign_coords(coords={PLEV_GCM:targ_plev})
            
            # sort for descending pressure
            var_out = var_out.reindex(
                    {PLEV_GCM:list(reversed(var_out[PLEV_GCM]))})
            #print(var_out)
            #var_out.to_netcdf('test.nc')

            # convert to dataset
            ds_out = var_out.to_dataset(name=var_name)

            ## make sure to keep time encoding identical
            ds_out.time.encoding = ds.time.encoding
            ## make sure to keep attributes identical
            for key,val in ds.time.attrs.items():
                ds_out.time.attrs[key] = val
            for key,val in ds.lon.attrs.items():
                ds_out.lon.attrs[key] = val
            for key,val in ds.lat.attrs.items():
                ds_out.lat.attrs[key] = val
            for key,val in ds[var_name].attrs.items():
                ds_out[var_name].attrs[key] = val

            # save output file
            ds_out.to_netcdf(out_file_path)
