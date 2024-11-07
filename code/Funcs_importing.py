# -*- coding: utf-8 -*-
"""
Functions for importing model data
"""

import numpy as np
import xarray as xr
import h5py
#%%
def import_oht_obs():
    data=h5py.File('../data/ECCOv4r5_rc2_analysis/MHT/MHT.jld2',"r")
    oht_obs_mon=data['single_stored_object'][:]
    yr_oht_obs= np.linspace(1992, 2019,28)
    lat_oht_obs= np.arange(-89.0, 90.0, 1.0)
    
    oht_obs_yr=np.zeros([28,179])
    for i in range(28):
        oht_obs_yr[i,:]= np.mean(oht_obs_mon[i*12:(i*12)+12,:],axis=0)
    
    return oht_obs_mon,oht_obs_yr,yr_oht_obs,lat_oht_obs

def import_iel_lomg(modelname):
    phi_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/iel_mon_1deg_bil/historical/iel_mon_1deg_bil_historical_'+modelname+'.nc')
    phi_ssp370 = xr.open_dataset('../data/_tmp_to_archive/cmip6/iel_mon_1deg_bil/ssp370/iel_mon_1deg_bil_ssp370_'+modelname+'.nc')
    
    if modelname=='UKESM1-0-LL':
        rows_to_del=[4,5,6]
    else:
        rows_to_del=[]
    phi=np.delete( np.concatenate((phi_hist.iel_n.to_numpy(), phi_ssp370.iel_n.to_numpy())),rows_to_del,axis=1)
    
    return phi_hist,phi_ssp370,phi

def import_CanESM5():
    oht_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/oht_from_hfbasin_yr_gn/historical/oht_from_hfbasin_yr_gn_historical_CanESM5.nc')
    oht_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/oht_from_hfbasin_yr_gn/ssp370/oht_from_hfbasin_yr_gn_ssp370_CanESM5.nc')
    phi_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/iel_zm_yr_1deg_bil/historical/iel_zm_yr_1deg_bil_historical_CanESM5.nc')
    phi_ssp370 = xr.open_dataset('../data/_tmp_to_archive/cmip6/iel_zm_yr_1deg_bil/ssp370/iel_zm_yr_1deg_bil_ssp370_CanESM5.nc')
    aht_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/aht_from_net_flux_yr_gn/historical/aht_from_net_flux_yr_gn_historical_CanESM5.nc')
    aht_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/aht_from_net_flux_yr_gn/ssp370/aht_from_net_flux_yr_gn_ssp370_CanESM5.nc')
    tas_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/tas_area_mean_yr_gn/historical/tas_area_mean_yr_gn_historical_CanESM5.nc')
    tas_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/tas_area_mean_yr_gn/ssp370/tas_area_mean_yr_gn_ssp370_CanESM5.nc')
    
    rows_to_del=[]
    oht=np.delete( np.concatenate((oht_hist.oht.to_numpy(), oht_ssp370.oht.to_numpy())),rows_to_del,axis=1)
    aht=np.delete( np.concatenate((aht_hist.aht_n.to_numpy(), aht_ssp370.aht_n.to_numpy())),rows_to_del,axis=1)
    phi=np.delete( np.concatenate((phi_hist.iel_zm_n.to_numpy(), phi_ssp370.iel_zm_n.to_numpy())),rows_to_del,axis=1)
    tas=np.delete( np.concatenate((tas_hist.tas_n.to_numpy(), tas_ssp370.tas_n.to_numpy())),rows_to_del,axis=1)
    
    return oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tas

def import_UKESM1_0():
   oht_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/oht_from_hfbasin_yr_gn/historical/oht_from_hfbasin_yr_gn_historical_UKESM1-0-LL.nc')
   oht_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/oht_from_hfbasin_yr_gn/ssp370/oht_from_hfbasin_yr_gn_ssp370_UKESM1-0-LL.nc')
   phi_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/iel_zm_yr_1deg_bil/historical/iel_zm_yr_1deg_bil_historical_UKESM1-0-LL.nc')
   phi_ssp370 = xr.open_dataset('../data/_tmp_to_archive/cmip6/iel_zm_yr_1deg_bil/ssp370/iel_zm_yr_1deg_bil_ssp370_UKESM1-0-LL.nc')
   aht_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/aht_from_net_flux_yr_gn/historical/aht_from_net_flux_yr_gn_historical_UKESM1-0-LL.nc')
   aht_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/aht_from_net_flux_yr_gn/ssp370/aht_from_net_flux_yr_gn_ssp370_UKESM1-0-LL.nc')
   tas_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/tas_area_mean_yr_gn/historical/tas_area_mean_yr_gn_historical_UKESM1-0-LL.nc')
   tas_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/tas_area_mean_yr_gn/ssp370/tas_area_mean_yr_gn_ssp370_UKESM1-0-LL.nc')
   
   rows_to_del=[4,5,6]
   oht=np.delete( np.concatenate((oht_hist.oht.to_numpy(), oht_ssp370.oht.to_numpy())),rows_to_del,axis=1)
   aht=np.delete( np.concatenate((aht_hist.aht_n.to_numpy(), aht_ssp370.aht_n.to_numpy())),rows_to_del,axis=1)
   phi=np.delete( np.concatenate((phi_hist.iel_zm_n.to_numpy(), phi_ssp370.iel_zm_n.to_numpy())),rows_to_del,axis=1)
   tas=np.delete( np.concatenate((tas_hist.tas_n.to_numpy(), tas_ssp370.tas_n.to_numpy())),rows_to_del,axis=1)
    
   return oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tas

def import_UKESM1_1():
    oht_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/oht_from_hfbasin_yr_gn/historical/oht_from_hfbasin_yr_gn_historical_UKESM1-1-LL.nc')
    oht_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/oht_from_hfbasin_yr_gn/ssp370/oht_from_hfbasin_yr_gn_ssp370_UKESM1-1-LL.nc')
    phi_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/iel_zm_yr_1deg_bil/historical/iel_zm_yr_1deg_bil_historical_UKESM1-1-LL.nc')
    phi_ssp370 = xr.open_dataset('../data/_tmp_to_archive/cmip6/iel_zm_yr_1deg_bil/ssp370/iel_zm_yr_1deg_bil_ssp370_UKESM1-1-LL.nc')
    aht_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/aht_from_net_flux_yr_gn/historical/aht_from_net_flux_yr_gn_historical_UKESM1-1-LL.nc')
    aht_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/aht_from_net_flux_yr_gn/ssp370/aht_from_net_flux_yr_gn_ssp370_UKESM1-1-LL.nc')
    tas_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/tas_area_mean_yr_gn/historical/tas_area_mean_yr_gn_historical_UKESM1-1-LL.nc')
    tas_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/tas_area_mean_yr_gn/ssp370/tas_area_mean_yr_gn_ssp370_UKESM1-1-LL.nc')
    
    rows_to_del=[]
    oht=np.delete( np.concatenate((oht_hist.oht.to_numpy(), oht_ssp370.oht.to_numpy())),rows_to_del,axis=1)
    aht=np.delete( np.concatenate((aht_hist.aht_n.to_numpy(), aht_ssp370.aht_n.to_numpy())),rows_to_del,axis=1)
    phi=np.delete( np.concatenate((phi_hist.iel_zm_n.to_numpy(), phi_ssp370.iel_zm_n.to_numpy())),rows_to_del,axis=1)
    tas=np.delete( np.concatenate((tas_hist.tas_n.to_numpy(), tas_ssp370.tas_n.to_numpy())),rows_to_del,axis=1)
    
    return oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tas

def import_MPIESM1_2_HR():
    oht_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/oht_from_hfbasin_yr_gn/historical/oht_from_hfbasin_yr_gn_historical_MPI-ESM1-2-HR.nc')
    oht_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/oht_from_hfbasin_yr_gn/ssp370/oht_from_hfbasin_yr_gn_ssp370_MPI-ESM1-2-HR.nc')
    phi_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/iel_zm_yr_05deg_bil/historical/iel_zm_yr_05deg_bil_historical_MPI-ESM1-2-HR.nc')
    phi_ssp370 = xr.open_dataset('../data/_tmp_to_archive/cmip6/iel_zm_yr_05deg_bil/ssp370/iel_zm_yr_05deg_bil_ssp370_MPI-ESM1-2-HR.nc')
    aht_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/aht_from_net_flux_yr_gn/historical/aht_from_net_flux_yr_gn_historical_MPI-ESM1-2-HR.nc')
    aht_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/aht_from_net_flux_yr_gn/ssp370/aht_from_net_flux_yr_gn_ssp370_MPI-ESM1-2-HR.nc')
    tas_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/tas_area_mean_yr_gn/historical/tas_area_mean_yr_gn_historical_MPI-ESM1-2-HR.nc')
    tas_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/tas_area_mean_yr_gn/ssp370/tas_area_mean_yr_gn_ssp370_MPI-ESM1-2-HR.nc')
    
    rows_to_del=[]
    oht=np.delete( np.concatenate((oht_hist.oht.to_numpy(), oht_ssp370.oht.to_numpy())),rows_to_del,axis=1)
    aht=np.delete( np.concatenate((aht_hist.aht_n.to_numpy(), aht_ssp370.aht_n.to_numpy())),rows_to_del,axis=1)
    phi=np.delete( np.concatenate((phi_hist.iel_zm_n.to_numpy(), phi_ssp370.iel_zm_n.to_numpy())),rows_to_del,axis=1)
    tas=np.delete( np.concatenate((tas_hist.tas_n.to_numpy(), tas_ssp370.tas_n.to_numpy())),rows_to_del,axis=1)
    
    return oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tas

def import_CNRM_CM6_1():
    oht_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/oht_from_hfx_hfy_yr_cc_approx/historical/oht_from_hfx_hfy_yr_cc_approx_historical_CNRM-CM6-1.nc')
    oht_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/oht_from_hfx_hfy_yr_cc_approx/ssp370/oht_from_hfx_hfy_yr_cc_approx_ssp370_CNRM-CM6-1.nc')
    phi_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/iel_zm_yr_1deg_bil/historical/iel_zm_yr_1deg_bil_historical_CNRM-CM6-1.nc')
    phi_ssp370 = xr.open_dataset('../data/_tmp_to_archive/cmip6/iel_zm_yr_1deg_bil/ssp370/iel_zm_yr_1deg_bil_ssp370_CNRM-CM6-1.nc')
    aht_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/aht_from_net_flux_yr_gn/historical/aht_from_net_flux_yr_gn_historical_CNRM-CM6-1.nc')
    aht_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/aht_from_net_flux_yr_gn/ssp370/aht_from_net_flux_yr_gn_ssp370_CNRM-CM6-1.nc')
    tas_hist = xr.open_dataset('../data/_tmp_to_archive/cmip6/tas_area_mean_yr_gn/historical/tas_area_mean_yr_gn_historical_CNRM-CM6-1.nc')
    tas_ssp370=xr.open_dataset('../data/_tmp_to_archive/cmip6/tas_area_mean_yr_gn/ssp370/tas_area_mean_yr_gn_ssp370_CNRM-CM6-1.nc')

    rows_to_del=[]
    oht=np.delete( np.concatenate((oht_hist.oht_n.to_numpy(), oht_ssp370.oht_n.to_numpy())),rows_to_del,axis=1)
    aht=np.delete( np.concatenate((aht_hist.aht_n.to_numpy(), aht_ssp370.aht_n.to_numpy())),rows_to_del,axis=1)
    phi=np.delete( np.concatenate((phi_hist.iel_zm_n.to_numpy(), phi_ssp370.iel_zm_n.to_numpy())),rows_to_del,axis=1)
    tas=np.delete( np.concatenate((tas_hist.tas_n.to_numpy(), tas_ssp370.tas_n.to_numpy())),rows_to_del,axis=1)
    
    return oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tas

def import_Omon(basin,freq,model_name):
    if model_name=='UKESM1-0-LL':
        msft_name='msftyz'
    else:
        msft_name='msftmz'
    omon=np.load('../data/' + model_name + '_Omon/'+msft_name+'_Omon_from_basin'+str(basin) +'_'+freq+'_' + model_name + '.npz')
    lat=omon['lat']
    lev=omon['lev']
    msft=omon[msft_name+'_'+freq]
    # r_id=omon['realization_id']
    # i_id=omon['initialization_id']
    # p_id=omon['physics_id']
    # f_id=omon['forcings_id']
    
    return lat,lev,msft

def import_velocity_em(modelname,scenario,component):
    velocity=np.load('../data/'+modelname+'_Omon/'+component+'o_Omon_'+modelname+'_'+scenario+'_em.npz')
    latitude=velocity['latitude']
    longitude=velocity['longitude']
    lev=velocity['lev']
    if scenario=='historical': vel=velocity['vo_hist_yr']
    elif scenario=='ssp370': vel=velocity['vo_ssp370_yr']
    
    return vel,latitude,longitude,lev