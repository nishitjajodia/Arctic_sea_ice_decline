# -*- coding: utf-8 -*-
"""
Covariance heatmaps
Plot 2D heatmaps of correlation coeffs 
"""

import numpy as np
from functions import *
from ImportingFuncs import *

can_oht_hist, can_oht_ssp370, can_aht_hist, can_aht_ssp370, can_phi_hist, can_phi_ssp370, can_tas_hist, can_tas_ssp370,can_oht,can_aht,can_phi,can_tas=import_CanESM5()
ukesm_oht_hist, ukesm_oht_ssp370, ukesm_aht_hist, ukesm_aht_ssp370, ukesm_phi_hist, ukesm_phi_ssp370, ukesm_tas_hist, ukesm_tas_ssp370,ukesm_oht,ukesm_aht,ukesm_phi,ukesm_tas=import_UKESM0_1()
mpiesmHR_oht_hist, mpiesmHR_oht_ssp370, mpiesmHR_aht_hist, mpiesmHR_aht_ssp370, mpiesmHR_phi_hist, mpiesmHR_phi_ssp370, mpiesmHR_tas_hist, mpiesmHR_tas_ssp370,mpiesmHR_oht,mpiesmHR_aht,mpiesmHR_phi,mpiesmHR_tas=import_MPIESM1_2_HR()
#%%
ref_lat_range=np.linspace(55,80,15)
avg_p=20  # years
lat_vary_toggle=np.array([1,1,1])  # oht,aht,tas | 0=off, 1=on
cov_o_phi=np.zeros([len(ref_lat_range),41])
cov_a_phi=np.zeros([len(ref_lat_range),41])
cov_o_a=np.zeros([len(ref_lat_range),41])
cov_t_phi=np.zeros([len(ref_lat_range),41])

savefigs=0
#%% CanESM5
for i in range(len(ref_lat_range)):
    ref_lat_want=lat_vary_toggle*ref_lat_range[i]
    cov_o_phi[i,:], cov_a_phi[i,:], cov_o_a[i,:], cov_t_phi[i,:]= lead_lag_cov(avg_p,ref_lat_want,can_oht_hist,can_aht_hist,can_tas_hist,can_phi,can_oht,can_aht,can_tas)

relation_heatmap(cov_o_phi, cov_a_phi, cov_o_a, cov_t_phi,ref_lat_range,avg_p,'Covariance Coefficient heatmaps_CanESM5',savefigs,fixlim=0)
#%% UKESM
for i in range(len(ref_lat_range)):
    ref_lat_want=lat_vary_toggle*ref_lat_range[i]
    cov_o_phi[i,:], cov_a_phi[i,:], cov_o_a[i,:], cov_t_phi[i,:]= lead_lag_cov(avg_p,ref_lat_want,ukesm_oht_hist,ukesm_aht_hist,ukesm_tas_hist,ukesm_phi,ukesm_oht,ukesm_aht,ukesm_tas)

relation_heatmap(cov_o_phi, cov_a_phi, cov_o_a, cov_t_phi,ref_lat_range,avg_p,'Covariance Coefficient heatmaps_UKESM01LL',savefigs,fixlim=0)

#%% MPIESM1-2-HR
for i in range(len(ref_lat_range)):
    ref_lat_want=lat_vary_toggle*ref_lat_range[i]
    cov_o_phi[i,:], cov_a_phi[i,:], cov_o_a[i,:], cov_t_phi[i,:]= lead_lag_cov(avg_p,ref_lat_want,mpiesmHR_oht_hist,mpiesmHR_aht_hist,mpiesmHR_tas_hist,mpiesmHR_phi,mpiesmHR_oht,mpiesmHR_aht,mpiesmHR_tas)

relation_heatmap(cov_o_phi, cov_a_phi, cov_o_a, cov_t_phi,ref_lat_range,avg_p,'Covariance Coefficient heatmaps_MPI-ESM1-2-HR',savefigs,fixlim=0)
