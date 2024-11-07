# -*- coding: utf-8 -*-
"""
plots fpr explaining correlation across ensembles
use UKESM1-0-LL, OHT-AHT as an example
"""
import numpy as np
import matplotlib.pyplot as plt
from Funcs_importing import *
from Functions import *

oht_hist, oht_ssp370, aht_hist, aht_ssp370, phi_hist, phi_ssp370, tas_hist, tas_ssp370,oht,aht,phi,tas= import_UKESM1_0()
years=np.linspace(1850,2100,251)

r_hist,r_fut,init_hist,init_fut,p_hist,p_fut,f_hist,f_fut= irpf_indices(phi_hist,oht_hist,aht_hist,tas_hist,phi_ssp370,oht_ssp370,aht_ssp370,tas_ssp370)
bnd=oht_hist.bnd.to_numpy()
member=oht_hist.member.to_numpy()
#%%
ref_lat_want=65
j_lat_oht = np.argmin(abs(np.array(oht_hist.variables["ref_lat"]) - ref_lat_want))
j_lat_aht = np.argmin(abs(np.array(aht_hist.variables["ref_lat_n"]) - ref_lat_want))
oht_mm=np.zeros([251,13])
aht_mm=np.zeros([251,13])
for i in range(9,251-10):
    oht_mm[i,:]=np.mean(oht[i-10:i+10,:,j_lat_oht])
    aht_mm[i,:]=np.mean(oht[i-10:i+10,:,j_lat_aht])

#%%
fig,ax=plt.subplots(2,1)
ax[0].plot(oht_mm[151,:])
ax[1].plot(aht_mm[151,:])
