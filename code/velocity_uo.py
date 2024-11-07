# -*- coding: utf-8 -*-
"""
uo data analysis
"""

from Funcs_importing import *
import numpy as np
import matplotlib.pyplot as plt
import os

#%%
modelname='CanESM5'
component='u'
uo_hist_yr,latitude,longitude,lev=import_velocity(modelname,'historical',component)
uo_ssp370_yr,latitude,longitude,lev=import_velocity(modelname,'ssp370',component)
#%%
u_surface=uo_hist_yr[:,10,:,:]

#%%
for i in range(165):
    fig,ax=plt.subplots(figsize=(15,10))
    extent=[0,359,0,90]
    plot=plt.contourf(u_surface[i,...],extent=extent,levels=np.linspace(-0.8,0.8,20),vmin=-0.8,vmax=0.8,cmap='RdBu_r')
    plt.colorbar(plot,ax=ax,location='bottom')
    ax.set_facecolor('grey')
    plt.show()