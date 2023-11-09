#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 11:27:28 2022

@author: rmv203
"""

#%%
#import matplotlib
#matplotlib.use('Agg')

# Analysis imports
import numpy as np
import numpy.ma as ma
import csv
import netCDF4
from netCDF4 import Dataset
import iris
import iris.coord_categorisation
import glob
import warnings
from iris.experimental.equalise_cubes import equalise_attributes

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm, rcParams, colors
from matplotlib import gridspec as gspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
import matplotlib.path as mpat
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D


#%%

cmip6_beta_2xCO2_soil = np.load('saved_data/cmip6_NEWbeta_soil_2xCO2.npy')
cmip6_beta_4xCO2_soil = np.load('saved_data/cmip6_NEWbeta_soil_4xCO2.npy')

cmip6_gamma_2xCO2_soil = np.load('saved_data/cmip6_NEWgamma_soil_2xCO2.npy')
cmip6_gamma_4xCO2_soil = np.load('saved_data/cmip6_NEWgamma_soil_4xCO2.npy')

cmip6_beta_2xCO2_land = np.load('saved_data/cmip6_NEWbeta_land_2xCO2.npy')
cmip6_beta_4xCO2_land = np.load('saved_data/cmip6_NEWbeta_land_4xCO2.npy')

cmip6_gamma_2xCO2_land = np.load('saved_data/cmip6_NEWgamma_land_2xCO2.npy')
cmip6_gamma_4xCO2_land = np.load('saved_data/cmip6_NEWgamma_land_4xCO2.npy')

print('beta')
beta_mean = np.mean((cmip6_beta_4xCO2_soil/cmip6_beta_4xCO2_land)*100, axis=0)
beta_std = np.std((cmip6_beta_4xCO2_soil/cmip6_beta_4xCO2_land)*100, axis=0)
print((cmip6_beta_4xCO2_soil/cmip6_beta_4xCO2_land)*100, beta_mean, beta_std)

print('gamma')
gamma_mean = np.mean((cmip6_gamma_4xCO2_soil/cmip6_gamma_4xCO2_land)*100, axis=0)
gamma_std = np.std((cmip6_gamma_4xCO2_soil/cmip6_gamma_4xCO2_land)*100, axis=0)
print((cmip6_gamma_4xCO2_soil/cmip6_gamma_4xCO2_land)*100, gamma_mean, gamma_std)



#%%

fig_figure1 = plt.figure(1, figsize=(50,44))
gs = gspec.GridSpec(2, 2, figure=fig_figure1, hspace=0.6, wspace=0.4)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':5,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':68,
    'xtick.labelsize':68,
    'ytick.labelsize':68,
    'font.size':68,
}
plt.rcParams.update(params)


column = 0
row = 0
ax = fig_figure1.add_subplot(gs[row, column])

label_list = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
x = np.arange(10)  # the label locations
width = 0.4  # the width of the bars

rects1 = ax.bar(x - 1/5, cmip6_beta_2xCO2_soil, width, color='sandybrown', label='Cs')
ax.bar(x + 1/5, cmip6_beta_2xCO2_land, width, color='g')

ax.set_ylabel(r'$\beta$ (PgC ppm$^{-1}$)')
ax.set_ylim((0, 1.75))
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_title(r'(a) 2xCO$_{2}$', y=1.15, fontweight='bold')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')



column = 1
row = 0
ax = fig_figure1.add_subplot(gs[row, column])


rects1 = ax.bar(x - 1/5, cmip6_beta_4xCO2_soil, width, color='sandybrown', label='Cs')
ax.bar(x + 1/5, cmip6_beta_4xCO2_land, width, color='g')

ax.set_ylabel(r'$\beta$ (PgC ppm$^{-1}$)')
ax.set_ylim((0, 1.75))
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_title(r'(b) 4xCO$_{2}$', y=1.15, fontweight='bold')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')



column = 0
row = 1
ax = fig_figure1.add_subplot(gs[row, column])

label_list = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
x = np.arange(10)  # the label locations
width = 0.4  # the width of the bars

rects1 = ax.bar(x - 1/5, cmip6_gamma_2xCO2_soil, width, color='sandybrown', label='Cs')
ax.bar(x + 1/5, cmip6_gamma_2xCO2_land, width, color='g')

ax.set_ylabel(r'$\gamma$ (PgC $^{\circ}$C$^{-1}$)')
ax.set_ylim((-110, 0))
ax.set_xticks(x)
ax.set_xticklabels(label_list)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')



column = 1
row = 1
ax = fig_figure1.add_subplot(gs[row, column])


rects1 = ax.bar(x - 1/5, cmip6_gamma_4xCO2_soil, width, color='sandybrown', label='Cs')
ax.bar(x + 1/5, cmip6_gamma_4xCO2_land, width, color='g')

ax.set_ylabel(r'$\gamma$ (PgC $^{\circ}$C$^{-1}$)')
ax.set_ylim((-110, 0))
ax.set_xticks(x)
ax.set_xticklabels(label_list)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')




# legend
handels = []
handels.extend([Line2D([0,0],[0,0], linewidth=20, color='sandybrown', label=r'$\Delta C_{s}$')])
handels.extend([Line2D([0,0],[0,0], linewidth=20, color='g', label=r'$\Delta C_{s, NPP}$')])
labels = [r'Soil ($\beta_{s}$, $\gamma_{s}$)', r'Land ($\beta_{L}$, $\gamma_{L}$)']
leg1 = ax.legend(handels, labels, loc='lower center', ncol=7, bbox_to_anchor=(-0.32, -0.85), fontsize=68)#title='Model colors',
plt.gca().add_artist(leg1)


#%%
fig_figure1.tight_layout()
fig_figure1.savefig('figures/fig5', bbox_inches='tight')
plt.close()

