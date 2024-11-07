# -*- coding: utf-8 -*-
"""
Time series plots: 4 models, ensemble members abd ensemble means
- iel_zm vs time plots
- oht_zm vs time plots
"""

import numpy as np
from Functions import *
from Funcs_importing import *
from Funcs_plots import *
import matplotlib.pyplot as plt
import xarray as xr
import os

can_oht_hist, can_oht_ssp370, can_aht_hist, can_aht_ssp370, can_phi_hist, can_phi_ssp370, can_tas_hist, can_tas_ssp370,can_oht,can_aht,can_phi,can_tas=import_CanESM5()
ukesm_oht_hist, ukesm_oht_ssp370, ukesm_aht_hist, ukesm_aht_ssp370, ukesm_phi_hist, ukesm_phi_ssp370, ukesm_tas_hist, ukesm_tas_ssp370,ukesm_oht,ukesm_aht,ukesm_phi,ukesm_tas=import_UKESM1_0()
mpiesmHR_oht_hist, mpiesmHR_oht_ssp370, mpiesmHR_aht_hist, mpiesmHR_aht_ssp370, mpiesmHR_phi_hist, mpiesmHR_phi_ssp370, mpiesmHR_tas_hist, mpiesmHR_tas_ssp370,mpiesmHR_oht,mpiesmHR_aht,mpiesmHR_phi,mpiesmHR_tas=import_MPIESM1_2_HR()
cnrmcm16_oht_hist, cnrmcm16_oht_ssp370, cnrmcm16_aht_hist, cnrmcm16_aht_ssp370, cnrmcm16_phi_hist, cnrmcm16_phi_ssp370, cnrmcm16_tas_hist, cnrmcm16_tas_ssp370,cnrmcm16_oht,cnrmcm16_aht,cnrmcm16_phi,cnrmcm16_tas=import_CNRM_CM6_1()
oht_obs_mon,oht_obs_yr,yr_oht_obs,lat_oht_obs=import_oht_obs()
#%%
modelnames=['CanESM5','UKESM1-0-LL','MPI-ESM1-2-HR','CNRM-CN6-1']
colours=['orange','blue','green','red']
years=np.linspace(1850,2100,251)
fontsize=18
#%% iel vs time
phi_data=[can_phi,ukesm_phi,mpiesmHR_phi,cnrmcm16_phi]
phi_mean=np.zeros([251,4])
for i in range(4):
    phi_mean[:,i]=np.mean(phi_data[i],axis=1)
# iel_obs_nasa=xr.open_dataset('../data/_archives/passive_microwave/iel_zm_yr_05deg_bil/iel_zm_yr_05deg_bil_nh_passive_microwave_NSIDC-0051.nc').iel_zm.to_numpy()
iel_obs_bootstrap=xr.open_dataset('../data/_archives/passive_microwave/iel_zm_yr_05deg_bil/iel_zm_yr_05deg_bil_nh_passive_microwave_NSIDC-0079.nc').iel_zm.to_numpy()
yr_iel_obs=np.linspace(1979,2022,44)
# plt.plot(yr_iel_obs,iel_obs_nasa)

fig,ax=plt.subplots(2,1,sharex=True,figsize=(8,10))

#########
plinex=2015*np.ones(10)
pliney=np.linspace(65,90,10)
ax[0].plot(yr_iel_obs,iel_obs_bootstrap,color='black',label='Observations',lw=3)
ax[0].plot(plinex,pliney,color='black',linestyle='--')
ax[0].text(1970,86,'Historical',fontsize=12)
ax[0].text(2020,86,'ssp370 (future projection)',fontsize=12)
for i in range(4):
    ax[0].plot(years,phi_mean[:,i],alpha=0.8,color=colours[i],label=modelnames[i])
    
    for j in range(len(phi_data[i][0,:])):
        ax[0].plot(years,phi_data[i][:,j],alpha=0.1,lw=0.5,color=colours[i])
ax[0].legend(fontsize=12,loc='upper left')
ax[0].set_xlim(1850,2050)
ax[0].set_ylim(67,88)
ax[0].set_title('$\phi$ vs Time',fontsize=fontsize)
# ax[0].set_xlabel('Year',fontsize=fontsize-2)
ax[0].set_ylabel('$\phi (\degree N)$',fontsize=fontsize-2)

# directory='plots'
# plots_path = os.path.join(directory,'observations_vs_model_iel.pdf')
# plt.savefig(plots_path,dpi=600)

#%% oht vs time plots- at ref latitude 65N

ref_lat_want=65
ref_lats=[can_oht_hist.ref_lat,ukesm_oht_hist.ref_lat,mpiesmHR_oht_hist.ref_lat,cnrmcm16_oht_hist.ref_lat_n]
# j_ref_lat=np.zeros(4)
oht_data=[can_oht,ukesm_oht,mpiesmHR_oht,cnrmcm16_oht]
oht_mean=np.zeros([251,4])
for i in range(4):
    j_ref_lat=int(np.argmin(abs(ref_lats[i].to_numpy() - ref_lat_want)))
    tempdata=oht_data[i].copy()
    oht_data[i]=tempdata[...,j_ref_lat]
    oht_mean[:,i]=np.mean(oht_data[i],axis=1)   # ensemble mean


ax[1].plot(yr_oht_obs,oht_obs_yr[:,154],label='Obervations',color='black',lw=3)
plinex=2015*np.ones(10)
pliney=np.linspace(0.1,0.45,10)
ax[1].plot(plinex,pliney,color='black',linestyle='--')
# ax[1].text(2005,0.38,'Historical',fontsize=12)
# ax[1].text(2020,0.38,'ssp370 (future projection)',fontsize=12)
for i in range(4):
    ax[1].plot(years,oht_mean[:,i],alpha=0.8,color=colours[i],label=modelnames[i])
    
    for j in range(len(oht_data[i][0,:])):
        ax[1].plot(years,oht_data[i][:,j],alpha=0.1,color=colours[i],lw=0.5)

ax[1].set_xlim(1850,2100)
ax[1].set_ylim(0.18,0.37)
ax[1].set_title('OHT vs Time',fontsize=fontsize)
ax[1].set_xlabel('Year',fontsize=fontsize-2)
ax[1].set_ylabel('OHT (PW)',fontsize=fontsize-2)
plt.tight_layout()

directory='plots'
plots_path = os.path.join(directory,'observations_vs_model_phi_OHT.pdf')
plt.savefig(plots_path,dpi=1200)

#%% iel + oht plot
# for i in range(4):
#     fontsize=17
#     fig2,ax1=plt.subplots(figsize=(8,6))
#     plt.title('OHT and $\phi$ from '+modelnames[i],fontsize=fontsize)
#     ax1.plot(years,oht_mean[:,i],color='blue',label='OHT')
#     ax1.plot(years,[None]*251,color='red',label='$\phi$')
#     ax2=ax1.twinx()
#     ax2.plot(years,phi_mean[:,i],color='red',label='$\phi$')
#     ax1.legend(loc='upper left',fontsize=fontsize-3)
#     ax1.set_ylabel('Poleward Ocean Heat Transport (PW)',fontsize=fontsize-2)
#     ax2.set_ylabel('Ice Edge Latitude,$\phi$ ($^\circ N$)',fontsize=fontsize-2)
#     ax1.set_xlabel('Year',fontsize=fontsize-2)
#     ax1.tick_params(axis='both',labelsize=fontsize-3)
#     ax2.tick_params(axis='both',labelsize=fontsize-3)
#     plt.tight_layout()
    
#     # directory='plots'
#     # plots_path = os.path.join(directory,'OHT_iel_timeseries_'+modelnames[i]+'.pdf')
#     # plt.savefig(plots_path,dpi=600)
