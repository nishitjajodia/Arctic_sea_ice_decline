"""
Can-ESM5
"""
import xarray as xr
import numpy as np
from Functions import *
from Funcs_importing import *
from Funcs_plots import *
import matplotlib.pyplot as plt
import os

oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tas= import_CanESM5()
lat,lev,msft=import_Omon(0,'yr','CanESM5')
msft_sv=msft*(10**-9) #converting to svedrups
msft_sv_total=np.sum(msft_sv,axis=1) # depth intergrated
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
    corr_o_phi, corr_a_phi, corr_o_a, corr_t_phi= lead_lag_corr(avg_p,ref_lat_want,oht_hist,aht_hist,tas_hist,phi,oht,aht,tas,'CanESM5')
      
    max_corr[i,:],max_corr_lag[i,:]=Calc_max_corr(corr_o_phi,corr_a_phi,corr_o_a,corr_t_phi,avg_p)
#%%
    LL_plots(corr_o_phi, corr_a_phi, corr_o_a, corr_t_phi, avg_p,ref_lat_want,ax[i])
# plt.tight_layout()

#%% Maxx corr vs ref latitude
corr_ref_lat_plots(57,73,max_corr,max_corr_lag,'CanESM5_corr_refcary_o_a.png',savefig=0)

#%%
amoc_index_lat=65
j_lat_psi = np.argmin(abs(lat - amoc_index_lat))
amoc_index=np.nanmax(msft_sv[:,:,j_lat_psi,:],axis=1)
amoc_index_mean=np.mean(amoc_index,axis=1)
#%%
phi_mean=np.mean(phi,axis=1)

fig2,ax2=plt.subplots(figsize=(10,7))
ax2.plot(years,phi_mean,alpha=1,color='blue',label=' Ice Edge Latitude')
ax3=ax2.twinx()
ax3.plot(years,amoc_index_mean,alpha=1,color='red',label='AMOC Index')
for i in range(25):
    ax3.plot(years,amoc_index[:,i],alpha=0.1,color='red')
    if i==0:
        ax2.plot(years,phi[:,i],alpha=0.2,color='blue')
    else:
        ax2.plot(years,phi[:,i],alpha=0.2,color='blue')
ax2.plot([],[],color='red',label='AMOC Index')
ax2.grid()
plt.tight_layout()
ax2.set_xlim(1850,2100)
ax2.set_ylim(65,90)
ax2.set_ylabel('Ice Edge Latitude ($\degree N$)',fontsize=15)
ax3.set_ylabel('AMOC Index (Sv)',fontsize=15)
ax2.set_xlabel('Years',fontsize=15)
# ax.set_title('$\phi$ vs Time (year)- CanESM5')
ax2.legend(fontsize=15)
# ax3.legend()
# directory='plots'
# plots_path = os.path.join(directory,'iel_time_plot_CanESM5')
# plt.savefig(plots_path)