# -*- coding: utf-8 -*-
"""
amoc index map- ref lat_vs time
"""
import xarray as xr
import numpy as np
from Functions import *
from Funcs_importing import *
from Funcs_plots import *
import matplotlib.pyplot as plt
import os
#%% import data
modelname=['CanESM5', 'UKESM1-0-LL','MPI-ESM1-2-HR']
lat=[None]*3
lev=[None]*3
msft=[None]*3
for i in range(3):
    if i==2:lat[i],lev[i],msft_temp=import_Omon(1,'yr',modelname[i])
    else:lat[i],lev[i],msft_temp=import_Omon(0,'yr',modelname[i])
    msft[i]=msft_temp*(10**-9) #converting to svedrups
#%%
amoc_index=[None]*3
for i in range(3):
    msft_mean=np.mean(msft[i],axis=3)
    amoc_index[i]=np.nanmax(msft_mean,axis=1)
#%%
fontsize=18
fig,ax=plt.subplots(3,1,sharex=True,figsize=(12,12),height_ratios = [1,1,1.4])
for i in range(3):
    extent=(lat[i][0],lat[i][-1],1850,2100)
    cmin=0
    plot=ax[i].contourf(amoc_index[i],extent=extent,levels=np.linspace(cmin,22,12),vmin=cmin,vmax=22,cmap='magma')
    # plot2=ax[i].contourf(amoc_index[i],extent=extent,levels=np.linspace(-1,0,50),vmin=-1,vmax=0,cmap='Greens')
    ax[i].set_title(modelname[i],fontsize=fontsize)
    ax[i].set_ylabel('Year',fontsize=fontsize-3)
    ax[i].tick_params(axis='both',labelsize=fontsize-4)
plt.colorbar(plot,ax=ax[2],location='bottom')
plt.suptitle('AMOC index maps',fontsize=fontsize)
ax[2].set_xlabel('Latitude ($^\circ N$)',fontsize=fontsize-3)

plt.tight_layout()

#%%
directory='plots'
plots_path = os.path.join(directory,'amoc_lat_time.pdf')
plt.savefig(plots_path,dpi=1200)