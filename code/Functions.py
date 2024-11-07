"""
functions- General fucntions for calculating things
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
from scipy import stats

def irpf_indices(phi_hist,oht_hist,aht_hist,tas_hist,phi_ssp370,oht_ssp370,aht_ssp370,tas_ssp370):
    r_hist= np.array([phi_hist.realisation_id.to_numpy(),oht_hist.realisation_id.to_numpy(),aht_hist.realisation_id.to_numpy(),tas_hist.realisation_id.to_numpy()])
    r_fut= np.array([phi_ssp370.realisation_id.to_numpy(),oht_ssp370.realisation_id.to_numpy(),aht_ssp370.realisation_id.to_numpy(),tas_ssp370.realisation_id.to_numpy()])
    init_hist= np.array([phi_hist.initialisation_id.to_numpy(),oht_hist.initialisation_id.to_numpy(),aht_hist.initialisation_id.to_numpy(),tas_hist.initialisation_id.to_numpy()])
    init_fut= np.array([phi_ssp370.initialisation_id.to_numpy(),oht_ssp370.initialisation_id.to_numpy(),aht_ssp370.initialisation_id.to_numpy(),tas_ssp370.initialisation_id.to_numpy()])
    p_hist= np.array([phi_hist.physics_id.to_numpy(),oht_hist.physics_id.to_numpy(),aht_hist.physics_id.to_numpy(),tas_hist.physics_id.to_numpy()])
    p_fut= np.array([phi_ssp370.physics_id.to_numpy(),oht_ssp370.physics_id.to_numpy(),aht_ssp370.physics_id.to_numpy(),tas_ssp370.physics_id.to_numpy()])
    f_hist= np.array([phi_hist.forcing_id.to_numpy(),oht_hist.forcing_id.to_numpy(),aht_hist.forcing_id.to_numpy(),tas_hist.forcing_id.to_numpy()])
    f_fut= np.array([phi_ssp370.forcing_id.to_numpy(),oht_ssp370.forcing_id.to_numpy(),aht_ssp370.forcing_id.to_numpy(),tas_ssp370.forcing_id.to_numpy()])
    
    return r_hist,r_fut,init_hist,init_fut,p_hist,p_fut,f_hist,f_fut


def lead_lag_corr(avg_p,ref_lat_want,oht_hist,aht_hist,tas_hist,phi,oht,aht,tas,modelname):
    
    #indices at which ref lat occurs
    if modelname=='CNRM-CM6-1':
        j_lat_oht = np.argmin(abs(np.array(oht_hist.variables["ref_lat_n"]) - ref_lat_want[0]))
    else:
        j_lat_oht = np.argmin(abs(np.array(oht_hist.variables["ref_lat"]) - ref_lat_want[0]))
    j_lat_aht = np.argmin(abs(np.array(aht_hist.variables["ref_lat_n"]) - ref_lat_want[1]))
    j_lat_tas = np.argmin(abs(np.array(tas_hist.variables["ref_lat_n"]) - ref_lat_want[2]))

    delta_phi= np.average(phi[151:151+avg_p,:],axis=0) -np.average(phi[150-avg_p:150,:],axis=0)
    delta_aht_fixed=np.average(aht[151:151+avg_p,:,j_lat_aht],axis=0) -np.average(aht[150-avg_p:150,:,j_lat_aht],axis=0)
    corr_o_phi=np.zeros(41)  # corr b/w delta OHT and delta phi
    corr_a_phi=np.zeros(41)  # corr b/w delta AHT and delta phi
    corr_o_a=np.zeros(41)    # corr b/w delta OHT and delta AHT
    corr_t_phi=np.zeros(41)  # corr b/w delta T and delta phi

    for i in range(-avg_p,avg_p+1):  
        delta_oht= np.average(oht[151+i:151+avg_p+i,:,j_lat_oht],axis=0) -np.average(oht[150-avg_p+i:150+i,:,j_lat_oht],axis=0)
        delta_aht= np.average(aht[151+i:151+avg_p+i,:,j_lat_aht],axis=0) -np.average(aht[150-avg_p+i:150+i,:,j_lat_aht],axis=0)
        delta_tas= np.average(tas[151+i:151+avg_p+i,:,j_lat_tas],axis=0) -np.average(tas[150-avg_p+i:150+i,:,j_lat_aht],axis=0)
        corr_o_phi[i+avg_p]=np.corrcoef(delta_oht,delta_phi,rowvar=False)[1,0]
        corr_a_phi[i+avg_p]=np.corrcoef(delta_aht,delta_phi,rowvar=False)[1,0]
        corr_o_a[i+avg_p]=np.corrcoef(delta_oht,delta_aht_fixed,rowvar=False)[1,0]
        corr_t_phi[i+avg_p]=np.corrcoef(delta_tas,delta_phi,rowvar=False)[1,0]

    return corr_o_phi, corr_a_phi, corr_o_a, corr_t_phi

def lead_lag_cov(avg_p,ref_lat_want,oht_hist,aht_hist,tas_hist,phi,oht,aht,tas):
    
    #indices at which ref lat occurs
    j_lat_oht = np.argmin(abs(np.array(oht_hist.variables["ref_lat"]) - ref_lat_want[0]))
    j_lat_aht = np.argmin(abs(np.array(aht_hist.variables["ref_lat_n"]) - ref_lat_want[1]))
    j_lat_tas = np.argmin(abs(np.array(tas_hist.variables["ref_lat_n"]) - ref_lat_want[2]))

    delta_phi= np.average(phi[151:151+avg_p,:],axis=0) -np.average(phi[150-avg_p:150,:],axis=0)
    delta_aht_fixed=np.average(aht[151:151+avg_p,:,j_lat_aht],axis=0) -np.average(aht[150-avg_p:150,:,j_lat_aht],axis=0)
   
    cov_o_phi=np.zeros(41)  # cov b/w delta OHT and delta phi
    cov_a_phi=np.zeros(41)  # cov b/w delta AHT and delta phi
    cov_o_a=np.zeros(41)    # cov b/w delta OHT and delta AHT
    cov_t_phi=np.zeros(41)  # cov b/w delta T and delta phi
    for i in range(-avg_p,avg_p+1):  
        delta_oht= np.average(oht[151+i:151+avg_p+i,:,j_lat_oht],axis=0) -np.average(oht[150-avg_p+i:150+i,:,j_lat_oht],axis=0)
        delta_aht= np.average(aht[151+i:151+avg_p+i,:,j_lat_aht],axis=0) -np.average(aht[150-avg_p+i:150+i,:,j_lat_aht],axis=0)
        delta_tas= np.average(tas[151+i:151+avg_p+i,:,j_lat_tas],axis=0) -np.average(tas[150-avg_p+i:150+i,:,j_lat_aht],axis=0)

        cov_o_phi[i+avg_p]=np.cov(delta_oht,delta_phi,rowvar=False)[1,0]
        cov_a_phi[i+avg_p]=np.cov(delta_aht,delta_phi,rowvar=False)[1,0]
        cov_o_a[i+avg_p]=np.cov(delta_oht,delta_aht_fixed,rowvar=False)[1,0]
        cov_t_phi[i+avg_p]=np.cov(delta_tas,delta_phi,rowvar=False)[1,0]
    
    return cov_o_phi, cov_a_phi, cov_o_a, cov_t_phi
  
def LLpearson(avg_p,ref_lat_j,data1,data2,modelname):
    ''' data2 is fixed in time'''
    

def Calc_max_corr(corr_o_phi,corr_a_phi,corr_o_a,corr_t_phi,avg_p):
    max_corr,max_corr_lag=np.zeros(4),np.zeros(4)
    max_corr[0]=np.max(corr_o_phi)
    max_corr[1]=np.min(corr_a_phi)
    max_corr[2]=np.min(corr_o_a)
    max_corr[3]=np.max(corr_t_phi)
    
    max_corr_lag[0]= -avg_p + np.argmax(corr_o_phi)
    max_corr_lag[1]= -avg_p + np.argmin(corr_a_phi)
    max_corr_lag[2]= -avg_p + np.argmin(corr_o_a)
    max_corr_lag[3]= -avg_p + np.argmax(corr_t_phi)
    
    return max_corr,max_corr_lag
    
    


    