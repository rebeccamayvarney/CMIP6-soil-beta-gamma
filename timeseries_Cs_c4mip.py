#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 2022

@author: Rebecca Varney (r.varney@exeter.ac.uk)

Script plots timeseries for dCs each C4MIP simulation.

"""

#%%

# Analysis imports
import numpy as np
import iris
import iris.coord_categorisation
import glob
import warnings
from iris.experimental.equalise_cubes import equalise_attributes

# My functions
from rmv_cmip_analysis import combine_netCDF_cmip6, combine_netCDF_model, open_netCDF, annual_average, area_average, global_total_percentage

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
from matplotlib.lines import Line2D


#%% Setting up the figure

fig_figure1 = plt.figure(1, figsize=(48,32))
gs = gspec.GridSpec(3, 3, figure=fig_figure1, hspace=0.4, wspace=0.4)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':7,
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

# CMIP6 models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
model_colors = ['peachpuff', '#fb8072', '#80b1d3', 'dodgerblue', 'darkcyan', 'darkgreen', 'olive', 'gold', 'orange', 'darkseagreen']
n_models = len(cmip6_models)

# C4MIP scenarios
C4MIP_simulation = ['1pctCO2', '1pctCO2-bgc', '1pctCO2-rad'] 

# dCO2 array
co2 = np.zeros((150))
for i in range(0, 150):
        if i == 0:
             co2[i] = 285
        else:
            co2[i] = co2[i-1]*1.01


# Loop through each C4MIP simulation being considered
for c4mip_option in range(0, len(C4MIP_simulation)):
    c4mip = C4MIP_simulation[c4mip_option]

    # for loop for each CMIP6 model
    for model_i in range(0, n_models):
            model = cmip6_models[model_i]
            
            print(model, c4mip)
            
            # land fraction
            landfraction = combine_netCDF_model('/home/rmv203/DATA/cmip6_data/sftlf_fx_'+model+'_historical*', model)
            
            # cSoil
            cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/cSoil_Emon_'+model+'_'+c4mip+'*', model)
            cSoil_cube = open_netCDF(cSoil_cube)
            cSoil_cube = annual_average(cSoil_cube)
            cSoil_cube = global_total_percentage(cSoil_cube, landfrac=landfraction, latlon_cons=None)
            cSoil_data = cSoil_cube.data
            # cLitter
            if model=='ACCESS-ESM1-5' or model=='BCC-CSM2-MR' or model=='CanESM5' or model=='CESM2' or model=='CNRM-ESM2-1' or model=='IPSL-CM6A-LR' or model=='MIROC-ES2L' or model=='MPI-ESM1-2-LR' or model=='NorESM2-LM' or model=='GFDL-ESM4':
                cLitter_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/cLitter_Lmon_'+model+'_'+c4mip+'*', model)
                cLitter_cube = open_netCDF(cLitter_cube)
                cLitter_cube = annual_average(cLitter_cube)
                cLitter_cube = global_total_percentage(cLitter_cube, landfrac=landfraction, latlon_cons=None)
                cLitter_data = cLitter_cube.data
                #
                Cs_data = cSoil_data + cLitter_data
            else:
                Cs_data = cSoil_data.copy()
            # rolling mean
            Cs_data_rm = np.zeros((len(Cs_data)-4))
            for i in range(2, len(Cs_data)-2):
                Cs_average = np.mean(Cs_data[i-2:i+2])
                Cs_data_rm[i-2] = Cs_average
            # dCs
            deltaCs = (Cs_data_rm - Cs_data_rm[0])/Cs_data_rm[0]


            # Time dimension
            if model == 'ACCESS-ESM1-5':
                time_dimension = cSoil_cube.coord('year').points
                time_dimension = time_dimension - 100
            
            
            # Temperature
            tas_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/tas_Amon_'+model+'_'+c4mip+'*', model)
            tas_cube = open_netCDF(tas_cube)
            tas_cube = annual_average(tas_cube)
            tas_cube = area_average(tas_cube, [0, 360, -90,  90])
            tas_data = tas_cube.data - 273.15
            # rolling mean
            tas_data_rm = np.zeros((len(tas_data)-4))
            for i in range(2, len(tas_data)-2):
                if i <= 5:
                    tas_average = np.mean(tas_data[i-2:i+13])
                    tas_data_rm[i-2] = tas_average
                elif i > 130:
                    tas_average = np.mean(tas_data[i-8:i+2])
                    tas_data_rm[i-2] = tas_average
                else:
                    tas_average = np.mean(tas_data[i-15:i+5])
                    tas_data_rm[i-2] = tas_average
            #dT
            deltaT = tas_data_rm - tas_data_rm[0]
            
            
            ## ajusting for ESM differences amongst CMIP6
            time_dimension = time_dimension[0:135]
            deltaCs = deltaCs[0:135]
            co2 = co2[0:135]
            deltaT = deltaT[0:135]
            
            
            #%% Plotting      
            
            if c4mip=='1pctCO2':
                x_axis = time_dimension
            elif c4mip=='1pctCO2-bgc':
                x_axis = co2
            elif c4mip=='1pctCO2-rad':
                x_axis = deltaT

            ax = fig_figure1.add_subplot(gs[0, c4mip_option])
        
            if model=='UKESM1-0-LL' and c4mip=='1pctCO2':
                ax.set_title('(a) 1pctCO2', y=1.15, fontweight='bold')
                #plt.ylim((-25,700))
            if model=='UKESM1-0-LL' and c4mip=='1pctCO2-bgc':
                ax.set_title('(b) 1pctCO2-bgc', y=1.15, fontweight='bold')
                #plt.ylim((-25,700))
            if model=='UKESM1-0-LL' and c4mip=='1pctCO2-rad':
                ax.set_title('(c) 1pctCO2-rad', y=1.15, fontweight='bold')   
                #plt.ylim((-325,25))
            
            plt.plot(x_axis, deltaCs, color=model_colors[model_i])
            
            if c4mip=='1pctCO2':
                plt.xlabel('Time (yr)')
                plt.ylabel(r'$\Delta C_{s}$ (PgC)')
                plt.xlim((0,135))
            elif c4mip=='1pctCO2-bgc':
                plt.xlabel(r'CO$_{2}$ (ppm)')
                plt.ylabel(r'$\Delta C_{s}^{BGC}$ (PgC)')
                plt.xlim((285,1000))
            elif c4mip=='1pctCO2-rad':
                plt.xlabel(r'$\Delta$T ($^{\circ}$C)')
                plt.ylabel(r'$\Delta C_{s}^{RAD}$ (PgC)')
                plt.xlim((0,3))

            if model=='UKESM1-0-LL' and c4mip=='1pctCO2-bgc':
                handels = []
                handels.extend([Line2D([0,0],[0,0], linewidth=20, color='peachpuff', label='ACCESS-ESM1-5')])
                handels.extend([Line2D([0,0],[0,0], linewidth=20, color='#fb8072', label='BCC-CSM2-MR')])
                handels.extend([Line2D([0,0],[0,0], linewidth=20, color='#80b1d3', label='CanESM5')])
                handels.extend([Line2D([0,0],[0,0], linewidth=20, color='dodgerblue', label='CESM2')])
                handels.extend([Line2D([0,0],[0,0], linewidth=20, color='darkcyan', label='GFDL-ESM4')])
                handels.extend([Line2D([0,0],[0,0], linewidth=20, color='darkgreen', label='IPSL-CM6A-LR')])
                handels.extend([Line2D([0,0],[0,0], linewidth=20, color='olive', label='MIROC-ES2L')])
                handels.extend([Line2D([0,0],[0,0], linewidth=20, color='gold', label='MPI-ESM1-2-LR')])
                handels.extend([Line2D([0,0],[0,0], linewidth=20, color='orange', label='NorESM2-LM')])
                handels.extend([Line2D([0,0],[0,0], linewidth=20, color='darkseagreen', label='UKESM1-0-LL')])
                labels = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
                leg1 = ax.legend(handels, labels, loc='lower center', ncol=5, borderaxespad=0.2, bbox_to_anchor=(0.5, -0.8), fontsize=52)
                plt.gca().add_artist(leg1)
            
            print(model, c4mip, deltaCs[130])

#%%
fig_figure1.savefig('figures/fig1_frac', bbox_inches='tight')
plt.close()