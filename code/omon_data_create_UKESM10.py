# -*- coding: utf-8 -*-
"""
Omon
"""

import numpy as np
import xarray as xr
import time
import matplotlib.pyplot as plt
import os
#%%
omon_data_hist=xr.open_dataset('../data/UKESM1-0-LL_Omon/historical/msftyz/msftyz_Omon_UKESM1-0-LL_historical_r4i1p1f2_gn_185001-194912.nc')


#%%

realization_id=np.int32(np.concatenate((np.linspace(1,4,4),np.linspace(8,12,5),np.array([16,17,18,19]))))
initialization_id=np.ones(13)
physics_id=np.ones(13)
forcings_id=np.ones(13)*2

#%%
msftyz_hist_mon1=np.zeros([1200,75,145,13])
msftyz_hist_mon2=np.zeros([780,75,145,13])
msftyz_ssp370_mon1=np.zeros([420,75,145,13])
msftyz_ssp370_mon2=np.zeros([612,75,145,13])
#%% For a single basin
# basin_no=1
# ref_lat_want=65
# j_lat_want = np.argmin(abs(np.array(omon_data_hist.variables["lat"]) - ref_lat_want))
lat=omon_data_hist.rlat.to_numpy()[185:]
lev=omon_data_hist.lev.to_numpy()
#%%
for basin in range(3):
    start_time=time.time()
    for i in range(len(realization_id)):                 
        member_name='r'+ str(realization_id[i]) + 'i1p1f2'
        #Conversion to xarray
        omon_hist_mon1=xr.open_dataset('../data/UKESM1-0-LL_Omon/historical/msftyz/msftyz_Omon_UKESM1-0-LL_historical_' + member_name + '_gn_185001-194912.nc')
        omon_hist_mon2=xr.open_dataset('../data/UKESM1-0-LL_Omon/historical/msftyz/msftyz_Omon_UKESM1-0-LL_historical_' + member_name + '_gn_195001-201412.nc')
        current_time= time.time()-start_time    
        print(member_name+'_hist completed. Time =', current_time)
        
        omon_ssp370_mon1=xr.open_dataset('../data/UKESM1-0-LL_Omon/ssp370/msftyz/msftyz_Omon_UKESM1-0-LL_ssp370_' + member_name + '_gn_201501-204912.nc')
        omon_ssp370_mon2=xr.open_dataset('../data/UKESM1-0-LL_Omon/ssp370/msftyz/msftyz_Omon_UKESM1-0-LL_ssp370_' + member_name + '_gn_205001-210012.nc')
        current_time= time.time()-start_time
        print(member_name+'_ssp370 completed. Time =', current_time)
        
        #conversion to numpy
        msftyz_hist_mon1[:,:,:,i]=omon_hist_mon1.msftyz.to_numpy()[:,basin,:,185:]
        msftyz_hist_mon2[:,:,:,i]=omon_hist_mon2.msftyz.to_numpy()[:,basin,:,185:]
        current_time= time.time()-start_time
        print(member_name+'_hist numpy converstion completed. Time =', current_time)
        
        msftyz_ssp370_mon1[:,:,:,i]=omon_ssp370_mon1.msftyz.to_numpy()[:,basin,:,185:]
        msftyz_ssp370_mon2[:,:,:,i]=omon_ssp370_mon2.msftyz.to_numpy()[:,basin,:,185:]
        current_time= time.time()-start_time
        print(member_name+'_ssp370 numpy conversion completed completed. Time =', current_time)
        
    msftyz_mon=np.concatenate((msftyz_hist_mon1,msftyz_hist_mon2, msftyz_ssp370_mon1,msftyz_ssp370_mon2),axis=0)
    msftyz_yr=np.zeros([251,75,145,13])
    for i in range(251):
        msftyz_yr[i,:,:,:]= np.mean(msftyz_mon[i*12:(i*12)+12,:,:,:],axis=0)
    #%%
    directory='../data/UKESM1-0-LL_Omon'
    filename_mon='msftyz_Omon_from_basin'+str(basin) +'_mon_UKESM1-0-LL'
    filename_yr='msftyz_Omon_from_basin'+str(basin) +'_yr_UKESM1-0-LL'
    filepath_mon=os.path.join(directory,filename_mon)
    filepath_yr=os.path.join(directory,filename_yr)
    #%%
    np.savez(filepath_mon,msftyz_mon=msftyz_mon,lat=lat,lev=lev,realization_id=realization_id,initialization_id=initialization_id,physics_id=physics_id,forcings_id=forcings_id)
    np.savez(filepath_yr,msftyz_yr=msftyz_yr,lat=lat,lev=lev,realization_id=realization_id,initialization_id=initialization_id,physics_id=physics_id,forcings_id=forcings_id)






