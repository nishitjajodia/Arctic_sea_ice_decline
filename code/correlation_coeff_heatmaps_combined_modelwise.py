# -*- coding: utf-8 -*-
"""
Correlation heatmaps
Plot 2D heatmaps of correlation coeffs 
"""

import numpy as np
from Functions import *
from Funcs_importing import *
from Funcs_plots import *

can_oht_hist, can_oht_ssp370, can_aht_hist, can_aht_ssp370, can_phi_hist, can_phi_ssp370, can_tas_hist, can_tas_ssp370,can_oht,can_aht,can_phi,can_tas=import_CanESM5()
ukesm_oht_hist, ukesm_oht_ssp370, ukesm_aht_hist, ukesm_aht_ssp370, ukesm_phi_hist, ukesm_phi_ssp370, ukesm_tas_hist, ukesm_tas_ssp370,ukesm_oht,ukesm_aht,ukesm_phi,ukesm_tas=import_UKESM1_0()
mpiesmHR_oht_hist, mpiesmHR_oht_ssp370, mpiesmHR_aht_hist, mpiesmHR_aht_ssp370, mpiesmHR_phi_hist, mpiesmHR_phi_ssp370, mpiesmHR_tas_hist, mpiesmHR_tas_ssp370,mpiesmHR_oht,mpiesmHR_aht,mpiesmHR_phi,mpiesmHR_tas=import_MPIESM1_2_HR()
cnrmcm16_oht_hist, cnrmcm16_oht_ssp370, cnrmcm16_aht_hist, cnrmcm16_aht_ssp370, cnrmcm16_phi_hist, cnrmcm16_phi_ssp370, cnrmcm16_tas_hist, cnrmcm16_tas_ssp370,cnrmcm16_oht,cnrmcm16_aht,cnrmcm16_phi,cnrmcm16_tas=import_CNRM_CM6_1()
#%%
ref_lat_range=np.linspace(55,80,15)
avg_p=20  # years
lat_vary_toggle=np.array([1,1,1])  # oht,aht,tas | 0=off, 1=on
corr_o_phi=np.zeros([len(ref_lat_range),41])
corr_a_phi=np.zeros([len(ref_lat_range),41])
corr_o_a=np.zeros([len(ref_lat_range),41])
corr_t_phi=np.zeros([len(ref_lat_range),41])

savefigs=0
#%% CanESM5
for i in range(len(ref_lat_range)):
    ref_lat_want=lat_vary_toggle*ref_lat_range[i]
    corr_o_phi[i,:], corr_a_phi[i,:], corr_o_a[i,:], corr_t_phi[i,:]= lead_lag_corr(avg_p,ref_lat_want,can_oht_hist,can_aht_hist,can_tas_hist,can_phi,can_oht,can_aht,can_tas,'CanESM5')

relation_heatmap(corr_o_phi, corr_a_phi, corr_o_a, corr_t_phi,ref_lat_range,avg_p,'Correlation Coefficient heatmaps_CanESM5',savefigs)

#%% UKESM
for i in range(len(ref_lat_range)):
    ref_lat_want=lat_vary_toggle*ref_lat_range[i]
    corr_o_phi[i,:], corr_a_phi[i,:], corr_o_a[i,:], corr_t_phi[i,:]= lead_lag_corr(avg_p,ref_lat_want,ukesm_oht_hist,ukesm_aht_hist,ukesm_tas_hist,ukesm_phi,ukesm_oht,ukesm_aht,ukesm_tas,'UKESM1-0-LL')

relation_heatmap(corr_o_phi, corr_a_phi, corr_o_a, corr_t_phi,ref_lat_range,avg_p,'Correlation Coefficient heatmaps_UKESM10LL',savefigs)

#%% MPIESM1-2-HR
for i in range(len(ref_lat_range)):
    ref_lat_want=lat_vary_toggle*ref_lat_range[i]
    corr_o_phi[i,:], corr_a_phi[i,:], corr_o_a[i,:], corr_t_phi[i,:]= lead_lag_corr(avg_p,ref_lat_want,mpiesmHR_oht_hist,mpiesmHR_aht_hist,mpiesmHR_tas_hist,mpiesmHR_phi,mpiesmHR_oht,mpiesmHR_aht,mpiesmHR_tas,'MPI-ESM1-2-HR')

relation_heatmap(corr_o_phi, corr_a_phi, corr_o_a, corr_t_phi,ref_lat_range,avg_p,'Correlation Coefficient heatmaps_MPI-ESM1-2-HR',savefigs)

#%%CNRM-CM6-1

for i in range(len(ref_lat_range)):
    ref_lat_want=lat_vary_toggle*ref_lat_range[i]
    corr_o_phi[i,:], corr_a_phi[i,:], corr_o_a[i,:], corr_t_phi[i,:]= lead_lag_corr(avg_p,ref_lat_want,cnrmcm16_oht_hist,cnrmcm16_aht_hist,cnrmcm16_tas_hist,cnrmcm16_phi,cnrmcm16_oht,cnrmcm16_aht,cnrmcm16_tas,'CNRM-CM6-1')

relation_heatmap(corr_o_phi, corr_a_phi, corr_o_a, corr_t_phi,ref_lat_range,avg_p,'Correlation Coefficient heatmaps_CNRM-CM6-1',savefigs)

