# -*- coding: utf-8 -*-
"""
obs data
"""
import numpy as np
from scipy.stats import linregress
import h5py
import matplotlib.pyplot as plt
import xarray as xr
#%% OHT observation data
data=h5py.File('../data/ECCOv4r5_rc2_analysis/MHT/MHT.jld2',"r")

oht_obs_mon=data['single_stored_object'][:]

yr_oht        = np.linspace(1992, 2019,28)
latitude        = np.arange(-89.0, 90.0, 1.0)
nc_var_name_oht = "single_stored_object"
oht_unit_factor = 1000.0

oht_obs_yr=np.zeros([28,179])
for i in range(28):
    oht_obs_yr[i,:]= np.mean(oht_obs_mon[i*12:(i*12)+12,:],axis=0)
#%%
slope, intercept, r_value, p_value, std_err = linregress(yr_oht, oht_obs_yr[:,154])
trendline = (slope * yr_oht) + intercept
#%%
plt.plot(yr_oht,oht_obs_yr[:,154])
plt.plot(yr_oht,trendline,linestyle='--',label='Linear Trendline')
plt.title('OHT observations')
plt.legend()
plt.grid()
plt.show()
#%% iel Observations

iel_obs_nasa=xr.open_dataset('../data/_archives/passive_microwave/iel_zm_yr_05deg_bil/iel_zm_yr_05deg_bil_nh_passive_microwave_NSIDC-0051.nc').iel_zm.to_numpy()
# iel_obs_bootstrap=xr.open_dataset('../data/_archives/passive_microwave/iel_zm_yr_05deg_bil/iel_zm_yr_05deg_bil_nh_passive_microwave_NSIDC-0079.nc')
yr_iel=np.linspace(1979,2022,44)
plt.plot(yr_iel,iel_obs_nasa)