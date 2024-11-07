"""
LL-correlation plots
"""
import xarray as xr
import numpy as np
from Functions import *
from Funcs_importing import *
from Funcs_plots import *
import matplotlib.pyplot as plt
from scipy import stats
import os
#%% choose model
modelname='CanESM5'

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
#%% ref lat indices
ref_lat_want=65 #degree N
if modelname=='CNRM-CM6-1':
    j_lat_oht = np.argmin(abs(np.array(oht_hist.variables["ref_lat_n"]) - ref_lat_want))
else:
    j_lat_oht = np.argmin(abs(np.array(oht_hist.variables["ref_lat"]) - ref_lat_want))
j_lat_aht = np.argmin(abs(np.array(aht_hist.variables["ref_lat_n"]) - ref_lat_want))
j_lat_tas = np.argmin(abs(np.array(tas_hist.variables["ref_lat_n"]) - ref_lat_want))

#%%amoc index
amoc_index_lat=65
j_lat_psi = np.argmin(abs(lat - amoc_index_lat))
amoc_index=np.nanmax(msft_sv[:,:,j_lat_psi,:],axis=1)
amoc_index_mean=np.mean(amoc_index,axis=1)
#%% Delta values- Fixed at 1980-2000 and 2001-2021
delta_phi_fixed= np.average(phi[151:151+avg_p,:],axis=0) -np.average(phi[150-avg_p:150,:],axis=0)
delta_aht_fixed=np.average(aht[151:151+avg_p,:,j_lat_aht],axis=0) -np.average(aht[150-avg_p:150,:,j_lat_aht],axis=0)
delta_oht_fixed=np.average(oht[151:151+avg_p,:,j_lat_oht],axis=0) -np.average(oht[150-avg_p:150,:,j_lat_oht],axis=0)
delta_tas_fixed=np.average(tas[151:151+avg_p,:,j_lat_tas],axis=0) -np.average(tas[150-avg_p:150,:,j_lat_tas],axis=0)
delta_amocindex_fixed= np.average(amoc_index[151:151+avg_p,:],axis=0) -np.average(amoc_index[150-avg_p:150,:],axis=0)
#%% Intialize corr arrays
corr_o_phi=np.zeros([41,2])  # corr b/w delta OHT and delta phi
corr_a_phi=np.zeros([41,2])  # corr b/w delta AHT and delta phi
corr_o_a=np.zeros([41,2])    # corr b/w delta OHT and delta AHT
corr_t_phi=np.zeros([41,2])  # corr b/w delta T and delta phi
corr_amoc_phi=np.zeros([41,2]) # corr b/w delta amoc index and delta phi
corr_amoc_oht=np.zeros([41,2]) # corr b/w delta amoc index and delta OHT
#%% Calculation of 
for i in range(-avg_p,avg_p+1):
    delta_phi= np.average(phi[151+i:151+avg_p+i,:],axis=0) -np.average(phi[150-avg_p+i:150+i,:],axis=0)
    delta_oht= np.average(oht[151+i:151+avg_p+i,:,j_lat_oht],axis=0) -np.average(oht[150-avg_p+i:150+i,:,j_lat_oht],axis=0)
    delta_aht= np.average(aht[151+i:151+avg_p+i,:,j_lat_aht],axis=0) -np.average(aht[150-avg_p+i:150+i,:,j_lat_aht],axis=0)
    delta_tas= np.average(tas[151+i:151+avg_p+i,:,j_lat_tas],axis=0) -np.average(tas[150-avg_p+i:150+i,:,j_lat_aht],axis=0)
    delta_amocindex= np.average(amoc_index[151+i:151+avg_p+i,:],axis=0) -np.average(amoc_index[150-avg_p+i:150+i,:],axis=0)
    
    corr_o_phi[i+avg_p,0]=stats.pearsonr(delta_oht,delta_phi_fixed)[0]
    corr_o_phi[i+avg_p,1]=stats.pearsonr(delta_oht_fixed,delta_phi)[0]
    
#%% Plots- both ways corr coeff 
# fig, ax = plt.subplots(figsize=(10,7))
# x=np.linspace(-avg_p, avg_p,(2*avg_p)+1)
# ax.plot(x,corr_o_phi[:,0],color='blue')
# ax.plot(x,corr_o_phi[:,1],color='red')
# ax.set_ylim(-1,1)
# ax.legend()
# ax.spines['bottom'].set_position('zero')
# ax.spines['right'].set_color('none')
# ax.spines['top'].set_color('none')
# ax.spines['left'].set_position(('outward', 10))
# ax.grid()












