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
omon_data_hist=xr.open_dataset('../data/CanESM5_Omon/historical/msftmz/msftmz_Omon_CanESM5_historical_r1i1p2f1_gn_185001-201412.nc')
# omon_data_ssp370=xr.open_dataset('../data/CanESM5_Omon/ssp370/msftmz/msftmz_Omon_CanESM5_ssp370_r5i1p2f1_gn_201501-210012.nc')

#%%

realization_id=np.int32(np.linspace(1,25,25))
initialization_id=np.ones(25)
physics_id=np.ones(25)*2
forcings_id=np.ones(25)

#%%
msftmz_hist_mon=np.zeros([1980,45,146,25])
msftmz_ssp370_mon=np.zeros([1032,45,146,25])
#%% For a single basin
# basin_no=1
# ref_lat_want=65
# j_lat_want = np.argmin(abs(np.array(omon_data_hist.variables["lat"]) - ref_lat_want))
lat=omon_data_hist.lat.to_numpy()[144:]
lev=omon_data_hist.lev.to_numpy()
#%%
for basin in range(3):
    start_time=time.time()
    for i in range(len(realization_id)):                 
        member_name='r'+ str(realization_id[i]) + 'i1p2f1'
        #Conversion to xarray
        omon_hist_mon=xr.open_dataset('../data/CanESM5_Omon/historical/msftmz/msftmz_Omon_CanESM5_historical_' + member_name + '_gn_185001-201412.nc')
        
        current_time= time.time()-start_time    
        print(member_name+'_hist completed. Time =', current_time)
        
        omon_ssp370_mon=xr.open_dataset('../data/CanESM5_Omon/ssp370/msftmz/msftmz_Omon_CanESM5_ssp370_' + member_name + '_gn_201501-210012.nc')
        
        current_time= time.time()-start_time
        print(member_name+'_ssp370 completed. Time =', current_time)
        
        #conversion to numpy
        msftmz_hist_mon[:,:,:,i]=omon_hist_mon.msftmz.to_numpy()[:,basin,:,144:]
        
        current_time= time.time()-start_time
        print(member_name+'_hist numpy converstion completed. Time =', current_time)
        
        msftmz_ssp370_mon[:,:,:,i]=omon_ssp370_mon.msftmz.to_numpy()[:,basin,:,144:]
        current_time= time.time()-start_time
        print(member_name+'_ssp370 numpy conversion completed completed. Time =', current_time)
        
    msftmz_mon=np.concatenate((msftmz_hist_mon, msftmz_ssp370_mon),axis=0)
    msftmz_yr=np.zeros([251,45,146,25])
    for i in range(251):
        msftmz_yr[i,:,:,:]= np.mean(msftmz_mon[i*12:(i*12)+12,:,:,:],axis=0)
    #%%
    directory='../data/CanESM5_Omon'
    filename_mon='msftmz_Omon_from_basin'+str(basin) +'_mon_CanESM5'
    filename_yr='msftmz_Omon_from_basin'+str(basin) +'_yr_CanESM5'
    filepath_mon=os.path.join(directory,filename_mon)
    filepath_yr=os.path.join(directory,filename_yr)
    #%%
    np.savez(filepath_mon,msftmz_mon=msftmz_mon,lat=lat,lev=lev,realization_id=realization_id,initialization_id=initialization_id,physics_id=physics_id,forcings_id=forcings_id)
    np.savez(filepath_yr,msftmz_yr=msftmz_yr,lat=lat,lev=lev,realization_id=realization_id,initialization_id=initialization_id,physics_id=physics_id,forcings_id=forcings_id)






