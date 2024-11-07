# -*- coding: utf-8 -*-
"""
single mode time series combined plots
"""

import numpy as np
from Functions import *
from Funcs_importing import *
from Funcs_plots import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import xarray as xr
import os
#%%
modelname=['CanESM5','UKESM1-0-LL','MPI-ESM1-2-HR','CNRM-CM6-1']
ref_lat_want=65
j_lat=np.zeros(4)
years=np.linspace(1850,2100,251)

#%%
colors=['Orange','Blue','Green','Red']
fontsize=22
fig,ax=plt.subplots(3,2,sharex=True,figsize=(15,18))
ax=ax.ravel()
ax[5].axis('off')
plt.suptitle('Time Series Plots',fontsize=fontsize+3)
oht_obs_mon,oht_obs_yr,yr_oht_obs,lat_oht_obs=import_oht_obs()
for m in range(4):
    #%% Importing and unit change
    if modelname[m]=='CanESM5': modelfunc=import_CanESM5()
    elif modelname[m]=='UKESM1-0-LL': modelfunc=import_UKESM1_0()
    elif modelname[m]=='MPI-ESM1-2-HR': modelfunc=import_MPIESM1_2_HR()
    elif modelname[m]=='CNRM-CM6-1': modelfunc=import_CNRM_CM6_1()
    oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tas= modelfunc
    if m==3:
        lat,lev,msft=np.linspace(0,90,10),np.linspace(0,3000,10),np.zeros([251,10,10,6])
    else:
        lat,lev,msft=import_Omon(0,'yr',modelname[m])
    msft_sv=msft*(10**-9) #converting to svedrups
    msft_sv_total=np.sum(msft_sv,axis=1) # depth intergrated
    years=np.linspace(1850,2100,251)
    avg_p=20 # years
    #%% ref lat index
    if modelname[m]=='CNRM-CM6-1':
        j_lat[0] = np.argmin(abs(np.array(oht_hist.variables["ref_lat_n"]) - ref_lat_want))
    else:
        j_lat[0] = np.argmin(abs(np.array(oht_hist.variables["ref_lat"]) - ref_lat_want))
    j_lat[1] = np.argmin(abs(np.array(aht_hist.variables["ref_lat_n"]) - ref_lat_want))
    j_lat[2] = np.argmin(abs(np.array(tas_hist.variables["ref_lat_n"]) - ref_lat_want))
    j_lat[3] = np.argmin(abs(lat - 26.5))
    j_lat=j_lat.astype(int)
    #%% data at ref lat
    amoc_index=np.nanmax(msft_sv[:,:,j_lat[3],:],axis=1)
    amoc_index_mean=np.mean(amoc_index,axis=1)
    data=[oht,aht,tas,amoc_index,phi]
    for i in range(3):
        reflatdata=data[i][:,:,j_lat[i]]
        data[i]=reflatdata
    
    #%% Plots
    titles=['OHT (in PW)','AHT (in PW)','TAS (in K)','AMOC Index (in Sv)', '$\phi$ in ($^\circ N$)']
       
    for i in range(len(data)):
        ax[i].set_title(titles[i],fontsize=fontsize-3)
        ensemble_mean=np.mean(data[i],axis=1)
        if m==3 and i==3: ensemble_mean=[None]*251
        ax[i].plot(years,ensemble_mean,alpha=1,color=colors[m],lw=2,label=modelname[m])     
        ax[i].tick_params(axis='both',labelsize=fontsize-3)
        if m ==3 and i==3:
            q=0  
        else:
            for j in range(len(data[i][0,:])):
                ax[i].plot(years,data[i][:,j],alpha=0.05,color=colors[m])           
    ax[0].legend(fontsize=fontsize-5)
ax[4].set_xlabel('Year',fontsize=fontsize-3)

#%%Observations
iel_obs_bootstrap=xr.open_dataset('../data/_archives/passive_microwave/iel_zm_yr_05deg_bil/iel_zm_yr_05deg_bil_nh_passive_microwave_NSIDC-0079.nc').iel_zm.to_numpy()
yr_iel_obs=np.linspace(1979,2022,44)
ax[0].plot(yr_oht_obs,oht_obs_yr[:,154],label='Obervations',color='black',lw=3)
ax[4].plot(yr_iel_obs,iel_obs_bootstrap,color='black',label='Observations',lw=3)
ax[0].legend(fontsize=fontsize-5)
plt.tight_layout()

#%%
directory='plots'
plots_path = os.path.join(directory,'time_series_combined.pdf')
plt.savefig(plots_path,dpi=1200)