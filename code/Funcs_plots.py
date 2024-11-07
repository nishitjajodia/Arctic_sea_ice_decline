# -*- coding: utf-8 -*-
"""
Plotting functions- Functions for making plots
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
#%%
def corr_ref_lat_plots(minlat,maxlat,max_corr,max_corr_lag,figname,savefig=0):
    x=np.linspace(minlat,maxlat,10)
    fig2, ax1 = plt.subplots(2,1,figsize=(10, 15))

    plot_colour=['blue','red','purple','cyan']
    plot_legend=['r[$\Delta$OHT(t),$\Delta\phi$]','r[$\Delta$AHT(t),$\Delta\phi$]','r[$\Delta$OHT(t),$\Delta$AHT(t)','r[$\Delta$T(t),$\Delta\phi$']

    for i in range(4):
        ax1[0].plot(x,max_corr[:,i],color=plot_colour[i],label=plot_legend[i])
        ax1[1].plot(x,max_corr_lag[:,i],color=plot_colour[i],label=plot_legend[i])
        
    ax1[0].set_ylim(-1,1)
    ax1[1].set_ylim(-22,22)
    ax1[0].set_title('Maximum correlation vs Reference latitude ($\degree N$)')
    ax1[1].set_title('Lead/Lag at maxx correlation (years) vs Reference latitude($\degree N$)')
    for i in range(2):
        # ax1[i].legend()
        ax1[i].spines['bottom'].set_position('zero')
        ax1[i].spines['right'].set_color('none')
        ax1[i].spines['top'].set_color('none')
        ax1[i].spines['left'].set_position(('outward', 10))
        ax1[i].grid()
        ax1[i].set_xlabel('Reference latitude')
        
        if savefig==1:
            plt.savefig('CanESM5_corr_refcary_o_a.png')
    plt.tight_layout()
    plt.show()
    
    
def relation_heatmap(rel_o_phi, rel_a_phi, rel_o_a, rel_t_phi,ref_lat_range,avg_p,title, savefig=0,fixlim=1):
    data=np.array([rel_o_phi, rel_a_phi, rel_o_a, rel_t_phi])
    fig, ax = plt.subplots(2,2,sharey=True,sharex=True,figsize=(15,10))
    ax = ax.ravel()
    fontsize=15
    fig.suptitle(title,fontsize=fontsize+3)
    extent=[-avg_p,avg_p, ref_lat_range[0],ref_lat_range[-1]]
    labels=['r[$\Delta$OHT(t),$\Delta\phi$]','r[$\Delta$AHT(t),$\Delta\phi$]','r[$\Delta$OHT(t),$\Delta$AHT(t)]','r[$\Delta$T(t),$\Delta\phi$]']
    ax[0].set_ylabel('Reference Latitude ($\degree N$)',fontsize=fontsize)
    ax[2].set_ylabel('Reference Latitude ($\degree N$)',fontsize=fontsize)
    yext=np.linspace(55,80,15)
    for i in range(4):  
        ax[i].set_title(labels[i],fontsize=fontsize)
        if i==0 or i==3:
            xext=-20 + np.argmax(data[i],axis=1)
        elif i==1 or i==2:
            xext=-20 + np.argmin(data[i],axis=1)
        ax[i].plot(xext,yext,color='black',alpha=1,linestyle='-.',linewidth=2)
        if i>=2:
            ax[i].set_xlabel('Lead/Lag (Years)',fontsize=fontsize)
        if fixlim==1:
            plot=ax[i].contourf(data[i],extent=extent,levels=np.linspace(-1, 1, 50), cmap='RdBu_r',vmin=-1, vmax=1)
        else:
            plot=ax[i].contourf(data[i],extent=extent,levels=np.linspace(np.min(data[i]), np.max(data[i]), 50), cmap='turbo',vmin=np.min(data[i]), vmax=np.max(data[i]))
        # plt.colorbar(plot,ax=ax[i],location='bottom')

    plt.tight_layout()
    if savefig==1:
        directory='plots'
        plots_path = os.path.join(directory,title+'.pdf')
        plt.savefig(plots_path)
    plt.show()

def LL_plots(corr_o_phi, corr_a_phi, corr_o_a, corr_t_phi, avg_p,ref_lat_want,ax):
    x=np.linspace(-avg_p, avg_p,(2*avg_p)+1)

    ax.plot(x,corr_o_phi,color='blue',label='r[$\Delta$OHT(t),$\Delta\phi$]')
    ax.plot(x,corr_a_phi,color='red',label='r[$\Delta$AHT(t),$\Delta\phi$]')
    ax.plot(x,corr_o_a,color='purple',label='r[$\Delta$OHT(t),$\Delta$AHT(t)')
    ax.plot(x,corr_t_phi,color='cyan',label='r[$\Delta$T(t),$\Delta\phi$')
    ax.set_ylim(-1,1)
    ax.legend()
    ax.spines['bottom'].set_position('zero')

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.spines['left'].set_position(('outward', 10))
    # ax.xaxis.set_ticks_position('bottom')
    # ax.yaxis.set_ticks_position('left')
    ax.grid()
    
    ax.text(-18,0.9,'$\phi_{ref,OHT}=$'+"{0:.2g}".format(ref_lat_want[0])+'$\degree N$',ha='center')
    ax.text(-18,0.8,'$\phi_{ref,AHT}=$'+"{0:.2g}".format(ref_lat_want[1])+'$\degree N$',ha='center')
    ax.text(-18,0.7,'$\phi_{ref,TAS}=$'+"{0:.2g}".format(ref_lat_want[2])+'$\degree N$',ha='center')