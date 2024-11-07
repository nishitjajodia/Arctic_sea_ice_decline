"""
variable wise correlation heat maps_cimbined for models
"""
#%% Importing
import numpy as np
from Functions import *
from Funcs_importing import *
from Funcs_plots import *
import matplotlib.pyplot as plt
from scipy import stats
import os

can_oht_hist, can_oht_ssp370, can_aht_hist, can_aht_ssp370, can_phi_hist, can_phi_ssp370, can_tas_hist, can_tas_ssp370,can_oht,can_aht,can_phi,can_tas=import_CanESM5()
ukesm_oht_hist, ukesm_oht_ssp370, ukesm_aht_hist, ukesm_aht_ssp370, ukesm_phi_hist, ukesm_phi_ssp370, ukesm_tas_hist, ukesm_tas_ssp370,ukesm_oht,ukesm_aht,ukesm_phi,ukesm_tas=import_UKESM1_0()
mpiesmHR_oht_hist, mpiesmHR_oht_ssp370, mpiesmHR_aht_hist, mpiesmHR_aht_ssp370, mpiesmHR_phi_hist, mpiesmHR_phi_ssp370, mpiesmHR_tas_hist, mpiesmHR_tas_ssp370,mpiesmHR_oht,mpiesmHR_aht,mpiesmHR_phi,mpiesmHR_tas=import_MPIESM1_2_HR()
cnrmcm16_oht_hist, cnrmcm16_oht_ssp370, cnrmcm16_aht_hist, cnrmcm16_aht_ssp370, cnrmcm16_phi_hist, cnrmcm16_phi_ssp370, cnrmcm16_tas_hist, cnrmcm16_tas_ssp370,cnrmcm16_oht,cnrmcm16_aht,cnrmcm16_phi,cnrmcm16_tas=import_CNRM_CM6_1()
oht_obs_mon,oht_obs_yr,yr_oht_obs,lat_oht_obs=import_oht_obs()

#%% Model name and colors
modelnames=['CanESM5','UKESM1-0-LL','MPI-ESM1-2-HR','CNRM-CM6-1']
colours=['orange','blue','green','red']
years=np.linspace(1850,2100,251)
avg_p=20
#%% datasets
phi_data=[can_phi,ukesm_phi,mpiesmHR_phi,cnrmcm16_phi]
oht_data=[can_oht,ukesm_oht,mpiesmHR_oht,cnrmcm16_oht]
aht_data=[can_aht,ukesm_aht,mpiesmHR_aht,cnrmcm16_aht]
tas_data=[can_tas,ukesm_tas,mpiesmHR_tas,cnrmcm16_tas]
#%% ref lats
ref_lat_range=np.linspace(50,87,15)
pearsonr_corr=np.zeros([len(ref_lat_range),41,4]) #lats x times x models
p_values=np.zeros([len(ref_lat_range),41,4])      #lats x times x models
j_lat=np.zeros([len(ref_lat_range),4])            #lats x models
extent=[-avg_p,avg_p, ref_lat_range[0],ref_lat_range[-1]]
#%% OHT- phi
for j in range(len(ref_lat_range)):
    j_lat[j,3] = np.argmin(abs(np.array(cnrmcm16_oht_hist.variables["ref_lat_n"]) - ref_lat_range[j]))
    j_lat[j,0] = np.argmin(abs(np.array(can_oht_hist.variables["ref_lat"]) - ref_lat_range[j]))
    j_lat[j,1] = np.argmin(abs(np.array(ukesm_oht_hist.variables["ref_lat"]) - ref_lat_range[j]))
    j_lat[j,2] = np.argmin(abs(np.array(mpiesmHR_oht_hist.variables["ref_lat"]) - ref_lat_range[j]))\

j_lat=j_lat.astype(int)
fig,ax=plt.subplots(2,2,sharex=True,sharey=True,figsize=(16,12))
plot=[None]*4
ax=ax.ravel()
fontsize=18
plt.suptitle('Correlation maps for $r[\Delta OHT,\Delta \phi]$',fontsize=fontsize+2,y=0.93)
for model in range(4):
    
    delta_phi_fixed=np.average(phi_data[model] [151:151+avg_p,:],axis=0) -np.average(phi_data[model][150-avg_p:150,:],axis=0)
    for j in range(len(ref_lat_range)):
        for i in range(-avg_p,avg_p+1):
            delta_oht=np.average(oht_data[model] [151+i:151+avg_p+i,:,j_lat[j,model ]],axis=0) -np.average(oht_data[model][150-avg_p+i:150+i,:,j_lat[j,model ]],axis=0)
            pearsonr_corr[j,i+avg_p,model],p_values[j,i+avg_p,model]=stats.pearsonr(delta_oht,delta_phi_fixed)
    
    plot[model]=ax[model].contourf(pearsonr_corr[:,:,model],extent=extent,levels=np.linspace(-1, 1,25), cmap='RdBu_r',vmin=-1, vmax=1)
    ax[model].contourf(p_values[:,:,model],extent=extent,levels=np.linspace(0.1, 1,2),cmap='gray',alpha=0.2,hatches='.')
    ax[model].contour(p_values[:,:,model],extent=extent,levels=np.linspace(0.1, 1,2),colors='black')
    ax[model].set_title(modelnames[model],fontsize=fontsize)
    ax[model].tick_params(axis='both',labelsize=fontsize-4)
    if model%2==0: ax[model].set_ylabel('Reference Latitude ($^\circ N$)',fontsize=fontsize)
    if model>=2: ax[model].set_xlabel('Lead time of $\Delta \phi$ (Years) ',fontsize=fontsize)
plt.subplots_adjust(hspace=0.1, wspace=0.05)
cbar=fig.colorbar(plot[1], ax=ax, orientation='horizontal', fraction=0.046, pad = 0.08)
cbar.ax.tick_params(labelsize=fontsize-4)
#%%
# directory='plots'
# plots_path = os.path.join(directory,'pearsonr_pvalue_modelcombined_oht_phi.pdf')
# plt.savefig(plots_path,dpi=1200)