#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 2022

@author: Rebecca Varney (r.varney@exeter.ac.uk)

Script showing the calculation of soil carbon-climate feedback parameters (gamma_s).
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


#%%

# C4MIP simulations
C4MIP_simulation = ['1pctCO2-rad'] #'1pctCO2', , '1pctCO2-rad'

# CMIP6 models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
model_colors = ['peachpuff', '#fb8072', '#80b1d3', 'dodgerblue', 'darkcyan', 'darkgreen', 'olive', 'gold', 'orange', 'darkseagreen']
n_models = len(cmip6_models)


# Loop through each C4MIP simulation being considered
for c4mip_option in range(0, len(C4MIP_simulation)):
    c4mip = C4MIP_simulation[c4mip_option]
    
    cmip6_gamma = np.zeros((len(cmip6_models)))
    cmip6_gamma_NEW = np.zeros((len(cmip6_models)))
    cmip6_gammaT = np.zeros((len(cmip6_models)))
    cmip6_deltaT = np.zeros((len(cmip6_models)))
    
    # Set up subplot figure
    fig = plt.figure(1, figsize=(32,68))
    gs = gspec.GridSpec(5, 2, figure=fig, width_ratios=[1, 1], hspace=0.8, wspace=0.7)
    n = 10
    column_1 = 0
    row_1 = 0
    n_columns_1 = 3
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'
    mpl.rcParams['xtick.top'] = True 
    mpl.rcParams['ytick.right'] = True
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    params = {
        'lines.linewidth':3,
        'axes.facecolor':'white',
        'xtick.color':'k',
        'ytick.color':'k',
        'axes.labelsize':84,
        'xtick.labelsize':84,
        'ytick.labelsize':84,
        'font.size':84,
    }
    plt.rcParams.update(params)
            

    # for loop for each CMIP6 model
    for model_i, a, in zip(range(n_models), range(n)):
            model = cmip6_models[model_i]
            
            print(model, c4mip)
            
            # land fraction
            landfraction = combine_netCDF_model('/home/rmv203/DATA/cmip6_data/sftlf_fx_'+model+'_historical*', model)
            
            # Soil Carbon (cSoil)
            cSoil_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/cSoil_Emon_'+model+'_'+c4mip+'*', model)
            cSoil_cube = open_netCDF(cSoil_cube)
            cSoil_cube = annual_average(cSoil_cube)
            cSoil_cube = global_total_percentage(cSoil_cube, landfrac=landfraction, latlon_cons=None)
            cSoil_data = cSoil_cube.data
            
            # Veg Carbon (cVeg)
            cVeg_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/cVeg_Lmon_'+model+'_'+c4mip+'*', model)
            cVeg_cube = open_netCDF(cVeg_cube)
            cVeg_cube = annual_average(cVeg_cube)
            cVeg_cube = global_total_percentage(cVeg_cube, landfrac=landfraction, latlon_cons=None)
            cVeg_data = cVeg_cube.data
        
            # Litter Carbon (cLitter)
            if model=='ACCESS-ESM1-5' or model=='BCC-CSM2-MR' or model=='CanESM5' or model=='CESM2' or model=='CNRM-ESM2-1' or model=='IPSL-CM6A-LR' or model=='MIROC-ES2L' or model=='MPI-ESM1-2-LR' or model=='NorESM2-LM' or model=='GFDL-ESM4':
                cLitter_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/cLitter_Lmon_'+model+'_'+c4mip+'*', model)
                cLitter_cube = open_netCDF(cLitter_cube)
                cLitter_cube = annual_average(cLitter_cube)
                cLitter_cube = global_total_percentage(cLitter_cube, landfrac=landfraction, latlon_cons=None)
                cLitter_data = cLitter_cube.data
                Cs_data = cSoil_data + cLitter_data + cVeg_data
            else:
                Cs_data = cSoil_data + cVeg_data
            #
            Cs_data_trimmed = Cs_data[0:70]
            
            # Temperature (tas)
            tas_cube = combine_netCDF_cmip6('/home/rmv203/DATA/CMIP_'+c4mip+'/tas_Amon_'+model+'_'+c4mip+'*', model)
            tas_cube = open_netCDF(tas_cube)
            tas_cube = annual_average(tas_cube)
            tas_cube = area_average(tas_cube - 273.15, [0, 360, -90,  90])
            tas_data = tas_cube.data
            tas_data_trimmed = tas_data[0:70]
            #
            deltaT = (np.mean(tas_data_trimmed[65:70]) - np.mean(tas_data_trimmed[0:5]))

            
            #print(row_1, column_1)
            ax = fig.add_subplot(gs[row_1, column_1])
            plt.plot(tas_data[0:140], Cs_data[0:140], 'k.', markersize=50)
            
            gamma_4 = (np.mean(Cs_data[135:140]) - np.mean(Cs_data[0:5])) / (np.mean(tas_data[135:140]) - np.mean(tas_data[0:5]))
            y4 = gamma_4*tas_data[0:140] + (Cs_data[0] - gamma_4*tas_data[0])
            plt.plot(tas_data[0:140], y4, 'seagreen', linewidth=15)#, alpha=0.75)
            
            # gamma
            #gamma = (Cs_data_trimmed[-1] - Cs_data_trimmed[0]) / (np.mean(tas_data_trimmed[65:70]) - np.mean(tas_data_trimmed[0:5]))
            gamma_2 = (np.mean(Cs_data_trimmed[65:70]) - np.mean(Cs_data_trimmed[0:5])) / (np.mean(tas_data_trimmed[65:70]) - np.mean(tas_data_trimmed[0:5]))
            #print(gamma)
            # lin reg
            m_1 = np.polyfit(tas_data_trimmed[0:70], Cs_data_trimmed[0:70], 1)
            m = m_1[0]
            #print(m)
            m_function = np.poly1d(m_1)
            sorted_tas = np.sort(tas_data_trimmed)
            #plt.plot(sorted_tas, m_function(sorted_tas), 'k', linewidth=10)
            y2 = gamma_2*tas_data_trimmed + (Cs_data[0] - gamma_2*tas_data[0])
            plt.plot(tas_data_trimmed, y2, 'darkseagreen', linewidth=15)#, alpha=0.75)
            
            cmip6_gamma[model_i] = m
#            cmip6_gamma_NEW[model_i] = gamma
            cmip6_gammaT[model_i] = m*deltaT
            cmip6_deltaT[model_i] = deltaT
            
            
            plt.xlabel(r'T ($^{\circ}$C)')
            plt.ylabel(r'$C_{s}^{RAD}$ (PgC)')
            plt.title(model)    
            plt.ylim((np.mean(Cs_data[130:140])*0.95, Cs_data[0]+15))
            
            # increase row and column 
            row_1 += 1 
            if row_1==5:
                column_1 += 1
                row_1 = 0
                
            if model == 'ACCESS-ESM1-5':
                plt.suptitle(r'(b) $\gamma_{s}$', y=0.95, x=0.45, fontweight='bold')
       
# legend
handels = []
handels.extend([Line2D([0,0],[0,0], linewidth=40, color='k')])
handels.extend([Line2D([0,0],[0,0], linewidth=40, color='darkseagreen')])
handels.extend([Line2D([0,0],[0,0], linewidth=40, color='seagreen')])
labels = [r'$C_{s}^{RAD}$', r'2xCO$_{2}$ $\gamma_{s}$', r'4xCO$_{2}$ $\gamma_{s}$']
leg1 = ax.legend(handels, labels, loc='lower center', ncol=3, bbox_to_anchor=(-0.4, -1.25), fontsize=84)
plt.gca().add_artist(leg1)
 

#%%
fig.savefig('figures/fig3b', bbox_inches='tight')
plt.close()

mean = np.mean(cmip6_gamma_NEW.data)
std = np.std(cmip6_gamma_NEW.data)
print('CMIP6 ensemble', mean, std)


#%%
np.save('saved_data/cmip6_gamma_veg_2xCO2.npy', cmip6_gamma.data)
np.save('saved_data/cmip6_NEWgamma_veg_2xCO2.npy', cmip6_gamma_NEW.data)
np.save('saved_data/cmip6_deltaTrad_2xCO2.npy', cmip6_deltaT.data)