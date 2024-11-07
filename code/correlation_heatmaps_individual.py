"""
LL- Correlation Heatmaps individualy
"""
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from Funcs_importing import *
from Functions import *
import os
#%% modename
modelname='CanESM5'
if modelname=='CanESM5': modelfunc=import_CanESM5()
elif modelname=='UKESM1-0-LL': modelfunc=import_UKESM1_0()
elif modelname=='MPIESM1-2-HR': modelfunc=import_MPIESM1_2_HR()
elif modelname=='CNRM-CM6-1': modelfunc=import_CNRM_CM6_1()
#%% Importing and unit change
oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tas= modelfunc
lat,lev,msft=import_Omon(0,'yr',modelname)
msft_sv=msft*(10**-9) #converting to svedrups
msft_sv_total=np.sum(msft_sv,axis=1) # depth intergrated
years=np.linspace(1850,2100,251)
avg_p=20 # years
#%% AMOC index
amoc_index_lat=65
j_lat_psi = np.argmin(abs(lat - amoc_index_lat))
amoc_index=np.nanmax(msft_sv[:,:,j_lat_psi,:],axis=1)
amoc_index_mean=np.mean(amoc_index,axis=1)
#%%
data=[oht,aht,tas,amoc_index,phi]
#%% Initialize corr arrays
ref_lat_range=np.linspace(55,80,15)
pearsonr_corr=np.zeros([len(ref_lat_range),41])
p_values=np.zeros([len(ref_lat_range),41])
j_lat=np.zeros([len(ref_lat_range),4])
for j in range(len(ref_lat_range)):
    if modelname=='CNRM-CM6-1':
        j_lat[j,0] = np.argmin(abs(np.array(oht_hist.variables["ref_lat_n"]) - ref_lat_range[j]))
    else:
        j_lat[j,0] = np.argmin(abs(np.array(oht_hist.variables["ref_lat"]) - ref_lat_range[j]))
    j_lat[j,1] = np.argmin(abs(np.array(aht_hist.variables["ref_lat_n"]) - ref_lat_range[j]))
    j_lat[j,2] = np.argmin(abs(np.array(tas_hist.variables["ref_lat_n"]) - ref_lat_range[j]))
j_lat=j_lat.astype(int)
#%% Choice of variable
varchoice=np.array([0,4])
#%% Calculating values
# delta_var2_fixed=np.average(data[varchoice[1]] [151:151+avg_p,:,j_lat[j,varchoice[1] ]],axis=0) -np.average(data[varchoice[1]][150-avg_p:150,:,j_lat[j,varchoice[1] ]],axis=0)
delta_var2_fixed=np.average(data[varchoice[1]] [151:151+avg_p,:],axis=0) -np.average(data[varchoice[1]][150-avg_p:150,:],axis=0)
xext=np.zeros(len(ref_lat_range))
for j in range(15):
    for i in range(-avg_p,avg_p+1):
        delta_var1=np.average(data[varchoice[0]] [151+i:151+avg_p+i,:,j_lat[j,varchoice[0] ]],axis=0) -np.average(data[varchoice[0]][150-avg_p+i:150+i,:,j_lat[j,varchoice[0] ]],axis=0)
        pearsonr_corr[j,i+avg_p],p_values[j,i+avg_p]=stats.pearsonr(delta_var1,delta_var2_fixed)
        xext[j]=-20 + np.argmax(pearsonr_corr[j,:])
#%% Plots
fig, ax = plt.subplots(2,sharex=True,figsize=(10,10))
plot=[None]*2
cbar=[None]*2
fontsize=18
extent=[-avg_p,avg_p, ref_lat_range[0],ref_lat_range[-1]]

ax[1].set_xlabel('Lead/Lag (Years)',fontsize=fontsize)
plot[0]=ax[0].contourf(pearsonr_corr,extent=extent,levels=np.linspace(-1, 1,50), cmap='RdBu_r',vmin=-1, vmax=1)
ax[0].plot(xext,ref_lat_range,color='black',alpha=1,linestyle=':',linewidth=3)
ax[0].set_title('(a) Correlation Map for $r[\Delta OHT, \Delta \phi]$',fontsize=fontsize-1)
ax[1].set_title('(b) p-value Map for $r[\Delta OHT, \Delta \phi]$',fontsize=fontsize-1)
plot[1]=ax[1].contourf(p_values,extent=extent,levels=np.linspace(0,0.05,5),cmap='Greens_r')
ax[1].contour(p_values,extent=extent,levels=np.linspace(0,0.05,2),colors='black',linestyles='--')
for i in range(2):
    ax[i].set_ylabel('Reference Latitude ($\degree N$)',fontsize=fontsize)
    ax[i].tick_params(axis='both',labelsize=fontsize-3)  # Change fontsize for x-axis
    # ax[i].yticks(fontsize=fontsize-3)  # Change fontsize for y-axis
    cbar[i]=fig.colorbar(plot[i],ax=ax[i],location='right')
    cbar[i].ax.tick_params(labelsize=fontsize-4)
# ax[0].grid()
plt.tight_layout()
#%%
# directory='plots'
# plots_path = os.path.join(directory,'pearsonr_pvalue_CanESM5_oht_phi.pdf')
# plt.savefig(plots_path,dpi=1200)
