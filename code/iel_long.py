"""
zonally varying iel data
"""
from Funcs_importing import *
import numpy as np
import matplotlib.pyplot as plt
#%%
modelname='CanESM5'
phi_hist_mon,phi_ssp370_mon,phi_mon=import_iel_lomg(modelname)

#%%
phi_yr=np.zeros([251,25,360])
years=np.linspace(1850,2100,251)
for i in range(251):
    phi_yr[i,...]= np.mean(phi_mon[i*12:(i*12)+12,...],axis=0)

#%%
phi_yr_mean=np.mean(phi_yr,axis=1)
phi_mon_mean=np.mean(phi_mon,axis=1)
#%%
x=np.linspace(-179.5,179.5,360)

for i in range(3012):
    plt.plot(x,phi_mon_mean[i],color='red')
    plt.ylim(0,90)
    plt.xlim(-179.5,179.5)
    plt.legend()
    plt.grid()
    plt.show()
