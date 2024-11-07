"""
MPI-EMS1-2-HR
"""
import xarray as xr
import numpy as np
from Functions import *
from Funcs_importing import *
from Funcs_plots import *
import matplotlib.pyplot as plt
import os

oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tas= import_MPIESM1_2_HR()
years=np.linspace(1850,2100,251)
#%%
r_hist,r_fut,init_hist,init_fut,p_hist,p_fut,f_hist,f_fut= irpf_indices(phi_hist,oht_hist,aht_hist,tas_hist,phi_ssp370,oht_ssp370,aht_ssp370,tas_ssp370)
bnd=oht_hist.bnd.to_numpy()
member=oht_hist.member.to_numpy()

#%% Calculating correlations
max_corr,max_corr_lag=np.zeros([10,4]),np.zeros([10,4])
lat_vary_toggle=np.array([0,0,1])  # oht,aht,tas | 0=off, 1=on
avg_p=20 # years
fig1, ax = plt.subplots(5,2,figsize=(20, 30))
ax = ax.ravel()
for i in range(10):
    ref_lat_want=65 + lat_vary_toggle*(-8 + 2*i)  # 1x3 matrix | oht, aht, tas | 55 to 74
    corr_o_phi, corr_a_phi, corr_o_a, corr_t_phi= lead_lag_corr(avg_p,ref_lat_want,oht_hist,aht_hist,tas_hist,phi,oht,aht,tas)
      
    max_corr[i,:],max_corr_lag[i,:]=Calc_max_corr(corr_o_phi,corr_a_phi,corr_o_a,corr_t_phi,avg_p)
#%%
    LL_plots(corr_o_phi, corr_a_phi, corr_o_a, corr_t_phi, avg_p,ref_lat_want,ax[i])
plt.tight_layout()

#%% Maxx corr vs ref latitude
corr_ref_lat_plots(57,73,max_corr,max_corr_lag,'MPIESM12HR_corr_refcary_o_a.png',savefig=0)
#%%
phi_mean=np.mean(phi,axis=1)

fig2,ax=plt.subplots(figsize=(10,7))
ax.plot(years,phi_mean,alpha=1,color='black',label='Ensemble Mean')

for i in range(10):
    if i==0:
        ax.plot(years,phi[:,i],alpha=0.2,color='blue',label='Ensemble Members')
    else:
        ax.plot(years,phi[:,i],alpha=0.2,color='blue')
    
plt.grid()
# plt.tight_layout()
ax.set_xlim(1850,2100)
ax.set_ylim(65,90)
ax.set_title('$\phi$ vs Time (year)- MPIESM1_2_HR')
ax.legend()
directory='plots'
plots_path = os.path.join(directory,'iel_time_plot_MPIESM1_2_HR')
plt.savefig(plots_path)