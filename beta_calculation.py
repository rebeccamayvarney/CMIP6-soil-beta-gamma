#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 2022

@author: Rebecca Varney (r.varney@exeter.ac.uk)

Script showing the calculation of soil carbon-concentration feedback parameters (beta_s).

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
from rmv_cmip_analysis import combine_netCDF_cmip6, combine_netCDF_model, open_netCDF, annual_average, global_total_percentage

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec as gspec
from matplotlib.lines import Line2D


#%%

# C4MIP simulations
C4MIP_simulation = ['1pctCO2-bgc'] #'1pctCO2', , '1pctCO2-rad'

# CMIP6 models
cmip6_models = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
model_colors = ['peachpuff', '#fb8072', '#80b1d3', 'dodgerblue', 'darkcyan', 'darkgreen', 'olive', 'gold', 'orange', 'darkseagreen']
n_models = len(cmip6_models)


# Loop through each C4MIP simulation being considered
for c4mip_option in range(0, len(C4MIP_simulation)):
    c4mip = C4MIP_simulation[c4mip_option]
    
    cmip6_beta = np.zeros((len(cmip6_models)))
    cmip6_beta_NEW = np.zeros((len(cmip6_models)))
    cmip6_betaCO2 = np.zeros((len(cmip6_models)))
    
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
    
    # CO2 array
    co2_plot = np.zeros((140))
    for i in range(0, 140):
            if i == 0:
                 co2_plot[i] = 285
            else:
                co2_plot[i] = co2_plot[i-1]*1.01
                
    # CO2 array
    co2 = np.zeros((70))
    for i in range(0, 70):
            if i == 0:
                 co2[i] = 285
            else:
                co2[i] = co2[i-1]*1.01
            

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
                Cs_data = cSoil_data + cLitter_data #+ cVeg_data
            else:
                Cs_data = cSoil_data #+ cVeg_data
                
            Cs_data_trimmed = Cs_data[0:70]
            
            delta_Cs_data_2xco2 = Cs_data[-1] - Cs_data[0]
            #print(delta_Cs_data_2xco2, Cs_data[-1], Cs_data[0])

            deltaCO2 = (co2[-1] - co2[0])
            
            
            # Plotting
            ax = fig.add_subplot(gs[row_1, column_1])
            plt.plot(co2_plot, Cs_data[0:140], 'k.', markersize=50)

            # beta
            beta_4 = (np.mean(Cs_data[135:140]) - np.mean(Cs_data[0:5])) / (co2_plot[-1] - co2_plot[0])
            y4 = beta_4*co2_plot + (Cs_data[0] - beta_4*co2_plot[0])
            plt.plot(co2_plot, y4, 'seagreen', linewidth=15)#, alpha=0.75)

            # beta
            beta_2 = (np.mean(Cs_data_trimmed[65:70]) - np.mean(Cs_data_trimmed[0:5])) / (co2[-1] - co2[0])
            # lin reg
            m_1 = np.polyfit(co2[0:70], Cs_data_trimmed[0:70], 1)
            m = m_1[0]
            #print(m)
            m_function = np.poly1d(m_1)
            sorted_co2 = np.sort(co2)
            #plt.plot(sorted_co2, m_function(sorted_co2), 'k', linewidth=10)
            y2 = beta_2*co2 + (Cs_data_trimmed[0] - beta_2*co2[0])
            plt.plot(co2, y2, 'darkseagreen', linewidth=15)#, alpha=0.75)
        
            
            cmip6_beta[model_i] = m
#            cmip6_beta_NEW[model_i] = beta
            cmip6_betaCO2[model_i] = m*deltaCO2
            
            plt.xlabel(r'CO$_{2}$ (ppm)')
            plt.ylabel(r'$C_{s}^{BGC}$ (PgC)')
            plt.xlim((280,1140))
            plt.ylim((Cs_data[0]-10,  np.mean(Cs_data[135:140])*1.075))
            plt.title(model)            
            
            # increase row and column 
            row_1 += 1 
            if row_1==5:
                column_1 += 1
                row_1 = 0
                
            if model == 'ACCESS-ESM1-5':
                plt.suptitle(r'(a) $\beta_{s}$', y=0.95, x=0.45, fontweight='bold')
                
                
# legend
handels = []
handels.extend([Line2D([0,0],[0,0], linewidth=40, color='k')])
handels.extend([Line2D([0,0],[0,0], linewidth=40, color='darkseagreen')])
handels.extend([Line2D([0,0],[0,0], linewidth=40, color='seagreen')])
labels = [r'$C_{s}^{BGC}$', r'2xCO$_{2}$ $\beta_{s}$', r'4xCO$_{2}$ $\beta_{s}$']
leg1 = ax.legend(handels, labels, loc='lower center', ncol=3, bbox_to_anchor=(-0.4, -1.25), fontsize=84)
plt.gca().add_artist(leg1)
                
            
#%%
fig.savefig('figures/fig3a', bbox_inches='tight')
plt.close()

mean = np.mean(cmip6_beta_NEW.data)
std = np.std(cmip6_beta_NEW.data)
print('ensemble', mean, std)

#%%
#np.save('saved_data/cmip6_beta_land_4xCO2.npy', cmip6_beta.data)
#np.save('saved_data/cmip6_NEWbeta_land_4xCO2.npy', cmip6_beta_NEW.data)
#np.save('saved_data/cmip6_deltaCO2_4xCO2.npy', deltaCO2)
