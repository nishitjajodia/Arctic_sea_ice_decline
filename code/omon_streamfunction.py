# -*- coding: utf-8 -*-
"""
streamfucntion variation plots
"""
import numpy as np
from Funcs_importing import *
import matplotlib.pyplot as plt
import os
#%%
model_name='CanESM5'
if model_name=='CanESM5':
    basins_name=['Atlantic','Pacific','Global']
else:
    basins_name=['Atlantic','Global','Pacific']
lat,lev,msft=import_Omon(0,'yr',model_name)
msft_sv=msft*(10**-9) #converting to svedrups
# msft_sv_total=np.sum(msft_sv,axis=1) # depth intergrated
years=np.linspace(1850,2100,251).astype(int)
msft_mm=np.mean(msft_sv,axis=3)
#%%

fig,ax=plt.subplots(1,2,figsize=(10,7),sharey=True)
#psi vs depth
year=2000
j_year=year-1850
ref_lat_range=np.arange(0,100,20)

for i in range(len(ref_lat_range)):
    j_lat=np.argmin(abs( lat- ref_lat_range[i]))
    ax[0].plot(lev,msft_mm[j_year,:,j_lat],label=ref_lat_range[i])
ax[0].grid()
ax[0].legend()
ax[0].set_xlabel('Depth (in m)')
ax[0].set_ylabel('Stream fucntion. $\psi$ (in SV)')
# ax[0].set_ylim(-3,13)
#%% psi vs lat
depths=np.array([0,100,500,1000,3000,5000])
for i in range(len(depths)):
    j_lev=np.argmin(abs( lev- depths[i]))
    ax[1].plot(lat,msft_mm[j_year,j_lev,:],label=depths[i])
ax[1].grid()
ax[1].legend()
ax[1].set_xlabel('Latitude (in $\degree N$)')

plt.tight_layout()
plt.show()