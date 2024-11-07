# -*- coding: utf-8 -*-
"""
single mode time series combined plots
"""

import numpy as np
from Functions import *
from Funcs_importing import *
from Funcs_plots import *
import matplotlib.pyplot as plt
import xarray as xr
import os
#%%
modelname='CanESM5'
ref_lat_want=65
j_lat=np.zeros(4)
years=np.linspace(1850,2100,251)
if modelname=='CanESM5': modelfunc=import_CanESM5()
elif modelname=='UKESM1-0-LL': modelfunc=import_UKESM1_0()
elif modelname=='MPIESM1-2-HR': modelfunc=import_MPIESM1_2_HR()
elif modelname=='CNRM-CM6-1': modelfunc=import_CNRM_CM6_1()
#%% Importing and unit change
oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tas= modelfunc
omon,lat,lev,msft=import_Omon(0,'yr',modelname)
msft_sv=msft*(10**-9) #converting to svedrups
msft_sv_total=np.sum(msft_sv,axis=1) # depth intergrated
years=np.linspace(1850,2100,251)
avg_p=20 # years
#%% ref lat index
if modelname=='CNRM-CM6-1':
    j_lat[0] = np.argmin(abs(np.array(oht_hist.variables["ref_lat_n"]) - ref_lat_want))
else:
    j_lat[0] = np.argmin(abs(np.array(oht_hist.variables["ref_lat"]) - ref_lat_want))
j_lat[1] = np.argmin(abs(np.array(aht_hist.variables["ref_lat_n"]) - ref_lat_want))
j_lat[2] = np.argmin(abs(np.array(tas_hist.variables["ref_lat_n"]) - ref_lat_want))
j_lat[3] = np.argmin(abs(lat - ref_lat_want))
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
colors=['Blue','Red','Green','Orange','Purple']

fontsize=18
fig,ax=plt.subplots(2,3,sharex=True,figsize=(18,12))
ax=ax.ravel()
plt.suptitle('Time Series Plots for '+modelname+' at '+str(ref_lat_want)+'$^\circ N$',fontsize=fontsize+1)
for i in range(len(data)):
    ax[i].set_title(titles[i]+' vs Time',fontsize=fontsize-3)
    ensemble_mean=np.mean(data[i],axis=1)
    ax[i].plot(years,ensemble_mean,alpha=1,color='Black',lw=2,label='Ensemble mean')
    ax[i].plot(years,[None]*251,color=colors[i],label='Ensemble Members')
    ax[i].tick_params(axis='both',labelsize=fontsize-3)
    for j in range(len(data[i][0,:])):
        ax[i].plot(years,data[i][:,j],alpha=0.1,color=colors[i])
        if i>=3: ax[i].set_xlabel('Year',fontsize=fontsize-3)
    ax[i].legend(fontsize=fontsize-4)
plt.tight_layout()

#%%
# directory='plots'
# plots_path = os.path.join(directory,'time_series_combined_CanESM5.pdf')
# plt.savefig(plots_path,dpi=600)