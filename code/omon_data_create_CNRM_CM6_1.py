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
omon_data_hist=xr.open_dataset('../data/CNRM-CM6-1_Omon/historical/msftyz/msftyz_Omon_CNRM-CM6-1_historical_r1i1p1f2_gn_185001-201412.nc')

#%%
realization_id=np.int32(np.linspace(1,6,6))
initialization_id=np.ones(6)
physics_id=np.ones(6)
forcings_id=np.ones(6)*2

#%%
msftyz_hist_yr=np.zeros([165,75,180,6])
msftyz_ssp370_yr=np.zeros([86,75,180,6])
#%% lat and lev
lat=omon_data_hist.variables['j-mean'].to_numpy()-180
lev=omon_data_hist.lev.to_numpy()
#%%
filenames_hist=sorted(os.listdir('../data/CNRM-CM6-1_Omon/historical/msftyz'))
def extract_number(filename):
    match = re.search('historical_r(\d+)', filename)
    return int(match.group(1)) if match else 0

# Creating a new sorted list using the custom key
filenames_hist = sorted(filenames_hist, key=extract_number)

filenames_ssp370=sorted(os.listdir('../data/CNRM-CM6-1_Omon/ssp370/msftyz'))
def extract_number(filename):
    match = re.search('ssp370_r(\d+)', filename)
    return int(match.group(1)) if match else 0

# Creating a new sorted list using the custom key
filenames_ssp370 = sorted(filenames_ssp370, key=extract_number)

path=['../data/CNRM-CM6-1_Omon/historical/msftyz/','../data/CNRM-CM6-1_Omon/ssp370/msftyz/']

#%%
for basin in range(3):
    start_time=time.time()
    for member in range(6):
        #%%historical
        year_init=0
        for file in range(33):
            #import u data from file
            msftyz_temp_mon_hist= xr.open_dataset( path[0] + filenames_hist[np.int32(33*member + file)] ).msftyz.to_numpy()[:,basin,:,:]
            print('imported from:'+filenames_hist[np.int32(33*member + file)])
            
            yrs_hist=np.int32(msftyz_temp_mon_hist.shape[0]/12)
            
            #initialise array for annual mean
            msftyz_temp_yr_hist=np.zeros([yrs_hist,41,180])
            for yr in range(yrs_hist):
                msftyz_temp_yr_hist[yr,...]=np.mean(msftyz_temp_mon_hist[yr*12:(yr*12)+12,:,:],axis=0)
            
            print('annual mean from '+filenames_hist[np.int32(33*member + file)])
            
            msftyz_hist_yr[year_init:year_init+yrs_hist,:,:,member]=msftyz_temp_yr_hist
            year_init=year_init+yrs_hist
       #%% ssp370
        year_init=0
        for file in range(18):
            #import u data from file          
            msftyz_temp_mon_ssp370= xr.open_dataset( path[1] + filenames_ssp370[np.int32(18*member + file)] ).msftyz.to_numpy()[:,basin,:,:]
            print('imported from:'+filenames_ssp370[np.int32(18*member + file)])
            yrs_ssp370=np.int32(msftyz_temp_mon_ssp370.shape[0]/12)
            
            #initialise array for annual mean
            msftyz_temp_yr_ssp370=np.zeros([yrs_ssp370,41,180])
            for yr in range(yrs_ssp370):
                msftyz_temp_yr_ssp370[yr,...]=np.mean(msftyz_temp_mon_ssp370[yr*12:(yr*12)+12,:,:],axis=0)
            
            print('annual mean from '+filenames_ssp370[np.int32(18*member + file)])
            
            msftyz_hist_yr[year_init:year_init+yrs_ssp370,:,:,member]=msftyz_temp_yr_ssp370
            year_init=year_init+yrs_ssp370
        print('done for member',member)
    #%% combine historical and ssp370
    msftyz_yr=np.concatenate((msftyz_hist_yr,msftyz_ssp370_yr),axis=0)
    #%%
    directory='../data/CNRM-CM6-1_Omon'
    filename_yr='msftyz_Omon_from_basin'+str(basin) +'_yr_CNRM-CM6-1'
    filepath_yr=os.path.join(directory,filename_yr)
    #%%
    np.savez(filepath_yr,msftyz_yr=msftyz_yr,lat=lat,lev=lev)






