# -*- coding: utf-8 -*-
"""
Omon
"""
import numpy as np
import xarray as xr
import time
import os
import re
#%%
omon_data_hist=xr.open_dataset('../data/MPI-ESM1-2-HR_Omon/historical/msftmz/msftmz_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_185001-185412.nc')

#%%
realization_id=np.int32(np.linspace(1,10,10))
initialization_id=np.ones(10)
physics_id=np.ones(10)
forcings_id=np.ones(10)*2

#%%
msftmz_hist_yr=np.zeros([165,41,90,10])
msftmz_ssp370_yr=np.zeros([86,41,90,10])
#%% lat and lev
lat=omon_data_hist.lat.to_numpy()[90:]
lev=omon_data_hist.lev.to_numpy()
#%%
filenames_hist=sorted(os.listdir('../data/MPI-ESM1-2-HR_Omon/historical/msftmz'))
def extract_number(filename):
    match = re.search('historical_r(\d+)', filename)
    return int(match.group(1)) if match else 0

# Creating a new sorted list using the custom key
filenames_hist = sorted(filenames_hist, key=extract_number)

filenames_ssp370=sorted(os.listdir('../data/MPI-ESM1-2-HR_Omon/ssp370/msftmz'))
def extract_number(filename):
    match = re.search('ssp370_r(\d+)', filename)
    return int(match.group(1)) if match else 0

# Creating a new sorted list using the custom key
filenames_ssp370 = sorted(filenames_ssp370, key=extract_number)

path=['../data/MPI-ESM1-2-HR_Omon/historical/msftmz/','../data/MPI-ESM1-2-HR_Omon/ssp370/msftmz/']

#%%
for basin in range(3):
    start_time=time.time()
    for member in range(10):
        #%%historical
        year_init=0
        for file in range(33):
            #import u data from file
            msftmz_temp_mon_hist= xr.open_dataset( path[0] + filenames_hist[np.int32(33*member + file)] ).msftmz.to_numpy()[:,basin,:,90:]
            print('imported from:'+filenames_hist[np.int32(33*member + file)])
            
            yrs_hist=np.int32(msftmz_temp_mon_hist.shape[0]/12)
            
            #initialise array for annual mean
            msftmz_temp_yr_hist=np.zeros([yrs_hist,41,90])
            for yr in range(yrs_hist):
                msftmz_temp_yr_hist[yr,...]=np.mean(msftmz_temp_mon_hist[yr*12:(yr*12)+12,:,:],axis=0)
            
            print('annual mean from '+filenames_hist[np.int32(33*member + file)])
            
            msftmz_hist_yr[year_init:year_init+yrs_hist,:,:,member]=msftmz_temp_yr_hist
            year_init=year_init+yrs_hist
       #%% ssp370
        year_init=0
        for file in range(18):
            #import u data from file          
            msftmz_temp_mon_ssp370= xr.open_dataset( path[1] + filenames_ssp370[np.int32(18*member + file)] ).msftmz.to_numpy()[:,basin,:,90:]
            print('imported from:'+filenames_ssp370[np.int32(18*member + file)])
            yrs_ssp370=np.int32(msftmz_temp_mon_ssp370.shape[0]/12)
            
            #initialise array for annual mean
            msftmz_temp_yr_ssp370=np.zeros([yrs_ssp370,41,90])
            for yr in range(yrs_ssp370):
                msftmz_temp_yr_ssp370[yr,...]=np.mean(msftmz_temp_mon_ssp370[yr*12:(yr*12)+12,:,:],axis=0)
            
            print('annual mean from '+filenames_ssp370[np.int32(18*member + file)])
            
            msftmz_ssp370_yr[year_init:year_init+yrs_ssp370,:,:,member]=msftmz_temp_yr_ssp370
            year_init=year_init+yrs_ssp370
        print('done for member',member)
    #%% combine historical and ssp370
    msftmz_yr=np.concatenate((msftmz_hist_yr,msftmz_ssp370_yr),axis=0)
    #%%
    directory='../data/MPI-ESM1-2-HR_Omon'
    filename_yr='msftmz_Omon_from_basin'+str(basin) +'_yr_MPI-ESM1-2-HR'
    filepath_yr=os.path.join(directory,filename_yr)
    #%%
    np.savez(filepath_yr,msftmz_yr=msftmz_yr,lat=lat,lev=lev)






