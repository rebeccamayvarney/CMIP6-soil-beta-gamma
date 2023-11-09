#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 2022

@author: Rebecca Varney (r.varney@exeter.ac.uk)

Script plots maps for dCs each C4MIP simulation.
"""

#%%

# Analysis imports
import numpy as np
import numpy.ma as ma
import iris
import iris.coord_categorisation
import glob
import warnings
from iris.experimental.equalise_cubes import equalise_attributes

# My functions
from rmv_cmip_analysis import combine_netCDF_cmip6, open_netCDF, annual_average

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker


#%%

fig = plt.figure(1, figsize=(34, 52))
gs = gspec.GridSpec(10, 4, figure=fig, width_ratios=[1, 1, 1, 0.075], hspace=0.1, wspace=0.1)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':3,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':52,
    'xtick.labelsize':52,
    'ytick.labelsize':52,
    'font.size':52,
}
plt.rcParams.update(params)


#%%

# C4MIP simulations
C4MIP_simulation = ['1pctCO2', '1pctCO2-bgc', '1pctCO2-rad']

# CMIP6 models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
n_models = len(cmip6_models)


# Loop through each C4MIP simulation being considered
for c4mip_option in range(0, len(C4MIP_simulation)):
    c4mip = C4MIP_simulation[c4mip_option]
    
    if c4mip=='1pctCO2':
        column = 0
    
    if c4mip=='1pctCO2-bgc':
        column = 1
        
    if c4mip=='1pctCO2-rad':
        column = 2


    # for loop for each CMIP6 model
    for model_i in range(0, n_models):
            model = cmip6_models[model_i]
            
            print(model, c4mip)
            
            # Soil Carbon (cSoil)
            cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/cSoil_Emon_'+model+'_'+c4mip+'*', model)
            cSoil_cube = open_netCDF(cSoil_cube)
            n_lat = cSoil_cube.coord('latitude').points
            n_lon = cSoil_cube.coord('longitude').points
            cSoil_cube = annual_average(cSoil_cube)
            cSoil_data = cSoil_cube.data
            # Litter Carbon (cLitter)
            if model=='ACCESS-ESM1-5' or model=='BCC-CSM2-MR' or model=='CanESM5' or model=='CESM2' or model=='CNRM-ESM2-1' or model=='IPSL-CM6A-LR' or model=='MIROC-ES2L' or model=='MPI-ESM1-2-LR' or model=='NorESM2-LM' or model=='GFDL-ESM4':
                cLitter_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/cLitter_Lmon_'+model+'_'+c4mip+'*', model)
                cLitter_cube = open_netCDF(cLitter_cube)
                cLitter_cube = annual_average(cLitter_cube)
                cLitter_data = cLitter_cube.data
                Cs_data = cSoil_data + cLitter_data
            else:
                Cs_data = cSoil_data.copy()
            #
            deltaCs = np.mean(Cs_data[130:140], axis=0) - np.mean(Cs_data[0:5], axis=0)
            deltaCs = ma.masked_where(deltaCs==0, deltaCs)


            #%% Plotting
            row = model_i
            ax = fig.add_subplot(gs[row, column], projection=ccrs.PlateCarree())
            
            if model=='ACCESS-ESM1-5' and c4mip=='1pctCO2':
                ax.set_title('(a) 1pctCO2', y=1.2, fontweight='bold')
            if model=='ACCESS-ESM1-5' and c4mip=='1pctCO2-bgc':
                ax.set_title('(b) 1pctCO2-bgc', y=1.2, fontweight='bold')
            if model=='ACCESS-ESM1-5' and c4mip=='1pctCO2-rad':
                ax.set_title('(c) 1pctCO2-rad', y=1.2, fontweight='bold')
                
            
            gl = ax.gridlines(linestyle='solid', color='white', linewidth=0.5, alpha=0.5)
            if c4mip=='1pctCO2':
                gl.yformatter=LATITUDE_FORMATTER
                gl.ylocator = mticker.FixedLocator([-60, -30, 0, 30, 60, 90])
                gl.ylabels_left=True
            if model=='UKESM1-0-LL':
                gl.xformatter=LONGITUDE_FORMATTER
                gl.xlocator = mticker.FixedLocator([-90, 0, 90])
                gl.xlabels_bottom=True
            ax.coastlines()
            # set up the x and y coordination
            lat = n_lat
            lon = n_lon
            x, y = np.meshgrid(lon, lat)
            
            print(np.min(deltaCs), np.max(deltaCs))
            line = np.arange(-10, 10, 0.5)
            diff = plt.contourf(x, y, deltaCs, line, cmap='BrBG_r', extend='both', transform=ccrs.PlateCarree(central_longitude=0))
            ax.set_ylim(-70,90)
            
            if c4mip=='1pctCO2':
                ax.text(-1.2, 0.5, model, transform=ax.transAxes, va='top', fontweight='bold', fontsize=52)

ax=fig.add_subplot(gs[:,3])
ax=plt.gca()
fig.colorbar(diff, ax, orientation='vertical').set_label(r'$\Delta C_{s}$ (kg C m$^{-2}$)')
            

# Save figure
fig.savefig('figures/fig2', bbox_inches='tight')
plt.close()
            