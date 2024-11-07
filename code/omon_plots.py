# -*- coding: utf-8 -*-
"""
OMON_plots
-depth vs latitude heatmaps for streamfucntion
-depth vs latitude heatmaps for delta(psi)
-3 basins: atlantic, pacific, global
- choose model- options: UKESM1-0-LL and CanESM5 
"""

import numpy as np
from Funcs_importing import *
import matplotlib.pyplot as plt
import os
#%% Data import and format
model_name='MPI-ESM1-2-HR'
if model_name=='CanESM5':
    basins_name=['Atlantic','Pacific','Global']
else:
    basins_name=['Atlantic','Global','Pacific']
for j in range(3):
    basin=j
    lat,lev,msft=import_Omon(basin,'yr',model_name)
    msft_sv=msft*(10**-9) #converting to svedrups
    msft_sv_total=np.sum(msft_sv,axis=1) # depth intergrated
    years=np.linspace(1850,2100,251).astype(int)
    #%% lev vs lat heatmap
    # member=10
    fig, ax = plt.subplots(2,1,figsize=(10,10), height_ratios = [1, 1.4])
    yr_i=[0,-1]
    for i in range(2):
        data=np.mean(msft_sv[yr_i[i],:,:,:],axis=2) #member averaged
        # data=msft_sv[yr_i[i],:,:,member]
        ax[i].set_ylabel('Depth (m)')
        ax[i].set_title(model_name+' - $\Psi$ for '+basins_name[j]+' Basin, year='+ str(years[yr_i[i]]))
        extent=[lat[0],lat[-1],-lev[0],-lev[-1]]
        plot=ax[i].contourf(data,extent=extent, cmap='turbo',levels=50,extend='both',alpha=1)
        # plot2=ax[i].contour(data,extent=extent,levels=50,colors='black',alpha=0.2)
    plt.colorbar(plot,ax=ax[i],location='bottom')
    ax[i].set_xlabel('latitude ($\degree N$)')
    plt.tight_layout()
    plt.show()
    # %% delta Omon heatmaps- member averaged
    # avg_p=20
    
    # title='$\Delta\Psi$ (in Svedrups) Member averaged, for '+basins_name[j]+' Basin - '+model_name
    # delta_omon= np.average(msft_sv[151:151+avg_p,:,:,:],axis=0) -np.average(msft_sv[150-avg_p:150,:,:,:],axis=0)
    
    # data=np.mean(delta_omon,axis=2)
    # fig, ax = plt.subplots(figsize=(10,7))
    # ax.set_title(title)
    # ax.set_xlabel('latitude ($\degree N$)')
    # ax.set_ylabel('$Depth (m)$')
    # extent=[lat[0],lat[-1],-lev[0],-lev[-1]]
    # plot=ax.contourf(data,extent=extent,levels=30,cmap='turbo',vmin=-1.5,vmax=1.5)
    # plt.colorbar(plot,ax=ax,location='bottom')
    # plt.tight_layout()
    
    # # directory='plots'
    # # plotname='Delta_msftmz_Heatmap_depthVSlat_'+basins_name[j]+'_basin_MemberAveraged_'+model_name
    # # plots_path = os.path.join(directory,plotname)
    # # plt.savefig(plots_path)
    
    # plt.show()

#%%
