# -*- coding: utf-8 -*-
"""
AMOC index
- calculate amoc index for chosen ref lat
- calculate correlation coeff of delta amoc index with oht and aht
- plot corr_coeff vs lead/lag time
- choose model- options: UKESM1-0-LL and CanESM5
"""

import xarray as xr
import numpy as np
from Functions import *
from Funcs_importing import *
from Funcs_plots import *
import matplotlib.pyplot as plt
import os
import scipy.stats as scistat
#%% Data import and format
model_name='CanESM5'
# model_name='UKESM1-0-LL'

if model_name=='CanESM5':
    members=25
    basins_name=['Atlantic','Pacific','Global']
    oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tascan_tas=import_CanESM5()
elif model_name=='UKESM1-0-LL':
    members=13
    basins_name=['Atlantic','Global','Pacific']
    oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tas= import_UKESM1_0()

omon,lat,lev,msft=import_Omon(0,'yr',model_name)
msft_sv=msft*(10**-9) #converting to svedrups
msft_sv_total=np.sum(msft_sv,axis=1) # depth intergrated
years=np.linspace(1850,2100,251)

#%%
avg_p=20 # years
delta_phi= np.average(phi[151:151+avg_p,:],axis=0) -np.average(phi[150-avg_p:150,:],axis=0)

#%%  Calculate AMOC index and correlation cooeffs for different lead/lag periods
amoc_index_lat=40
ref_lat_want=40
j_lat_psi = np.argmin(abs(lat - amoc_index_lat))
fig1, ax = plt.subplots()
amoc_index=np.nanmax(msft_sv[:,:,j_lat_psi,:],axis=1)  # calculate amoc index
# initialise arrays for correlation coeffs
corr_amoc_phi=np.zeros(41)
cov_amoc_phi=np.zeros(41)
corr_o_phi=np.zeros(41)
corr_amoc_oht=np.zeros(41)
p1=np.zeros(41)
p2=np.zeros(41)
# index of array with re_lat_req
j_lat_oht = np.argmin(abs(np.array(oht_hist.variables["ref_lat"]) - ref_lat_want))

# calc delta oht for fixed peried 1980-2000
delta_oht_fixed=np.average(oht[151:151+avg_p,:,j_lat_oht],axis=0) -np.average(oht[150-avg_p:150,:,j_lat_oht],axis=0)

# calculate correlation coeffecients
for i in range(-avg_p,avg_p+1):
    delta_amoc_index= np.average(amoc_index[151+i:151+i+avg_p,:],axis=0) -np.average(amoc_index[150+i-avg_p:150+i,:],axis=0)
    delta_oht= np.average(oht[151+i:151+avg_p+i,:,j_lat_oht],axis=0) -np.average(oht[150-avg_p+i:150+i,:,j_lat_oht],axis=0)
    corr_amoc_phi[i+avg_p]=np.corrcoef(delta_amoc_index,delta_phi,rowvar=False)[1,0]
    corr_o_phi[i+avg_p]=np.corrcoef(delta_oht,delta_phi,rowvar=False)[1,0]
    corr_amoc_oht[i+avg_p]=np.corrcoef(delta_oht_fixed,delta_amoc_index,rowvar=False)[1,0]
    p1[i+avg_p]=scistat.pearsonr(delta_amoc_index,delta_phi)[1]
    p2[i+avg_p]=scistat.pearsonr(delta_amoc_index,delta_oht)[1]
    # cov_amoc_phi[i+avg_p]=np.cov(delta_amoc_index,delta_phi,rowvar=False)[1,0]
    
#%% plots: correlation coeff vs lead/lag time
x=np.linspace(-avg_p, avg_p,(2*avg_p)+1)
ax.plot(x,corr_amoc_phi,label='r[$\Delta AMOC, \Delta \phi$]')
ax.plot(x,corr_amoc_oht,label='r[$\Delta AMOC, \Delta OHT$]')
# ax.plot(x,cov_amoc_phi,label='cov')

ax.set_ylim(-1,1)
ax.set_xlim(-20,20)
ax.grid()
plt.legend()
ax.spines['bottom'].set_position('zero')

ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
plt.show()
# ax.spines['left'].set_position(('outward', 10))
# ax.xaxis.set_ticks_position('bottom')
# ax.yaxis.set_ticks_position('left')

#%% stat tests plots
plt.plot(x,p1,color='blue')
plt.plot(x,p2,color='orange')
plt.grid()






















