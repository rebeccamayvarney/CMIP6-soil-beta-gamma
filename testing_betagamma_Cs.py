#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 10:12:56 2022

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

cmip6_deltaCs_1pctCO2_2xCO2 = np.load('saved_data/cmip6_deltaCsfull_2xCO2.npy')
cmip6_deltaCs_1pctCO2_4xCO2 = np.load('saved_data/cmip6_deltaCsfull_4xCO2.npy')

cmip6_betaCO2_2xCO2 = np.load('saved_data/cmip6_NEWbeta_soil_2xCO2.npy')
cmip6_betaCO2_4xCO2 = np.load('saved_data/cmip6_NEWbeta_soil_4xCO2.npy')

cmip6_gammaT_2xCO2 = np.load('saved_data/cmip6_NEWgamma_soil_2xCO2.npy')
cmip6_gammaT_4xCO2 = np.load('saved_data/cmip6_NEWgamma_soil_4xCO2.npy')

cmip6_deltaT_2xCO2 = np.load('saved_data/cmip6_deltaTrad_2xCO2.npy')
cmip6_deltaT_4xCO2 = np.load('saved_data/cmip6_deltaTrad_4xCO2.npy')

cmip6_2xCO2 = np.load('saved_data/cmip6_deltaCO2_2xCO2.npy')
cmip6_4xCO2 = np.load('saved_data/cmip6_deltaCO2_4xCO2.npy')


r_squared = np.corrcoef(cmip6_betaCO2_4xCO2, cmip6_gammaT_4xCO2)[0, 1]**2
print(r_squared)

#
deltaCs_estimated_2xCO2 = cmip6_betaCO2_2xCO2*cmip6_2xCO2 + cmip6_gammaT_2xCO2*cmip6_deltaT_2xCO2
deltaCs_estimated_4xCO2 = cmip6_betaCO2_4xCO2*cmip6_4xCO2 + cmip6_gammaT_4xCO2*cmip6_deltaT_4xCO2

#print(deltaCs_estimated_2xCO2)
#print(cmip6_deltaCs_1pctCO2_2xCO2)
print((cmip6_deltaCs_1pctCO2_4xCO2/deltaCs_estimated_4xCO2)*100)

#%%

fig_figure1 = plt.figure(1, figsize=(42,40))
gs = gspec.GridSpec(2, 2, figure=fig_figure1, hspace=0.5, wspace=0.3)
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
params = {
    'lines.linewidth':5,
    'axes.facecolor':'white',
    'xtick.color':'k',
    'ytick.color':'k',
    'axes.labelsize':58,
    'xtick.labelsize':58,
    'ytick.labelsize':58,
    'font.size':58,
}
plt.rcParams.update(params)


column = 0
row = 0
ax = fig_figure1.add_subplot(gs[row, column])

label_list = ['ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CanESM5', 'CESM2', 'GFDL-ESM4', 'IPSL-CM6A-LR', 'MIROC-ES2L', 'MPI-ESM1-2-LR', 'NorESM2-LM', 'UKESM1-0-LL']
x = np.arange(10)  # the label locations
width = 0.4  # the width of the bars

rects1 = ax.bar(x - 1/5, cmip6_deltaCs_1pctCO2_2xCO2, width, color='b', label='Cs')
ax.bar(x + 1/5, deltaCs_estimated_2xCO2, width, color='lightseagreen')

ax.set_ylabel(r'$\Delta C_{s}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((0, 195))
ax.set_title(r'(a) 2xCO$_{2}$', y=1.09, fontweight='bold')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')



column = 1
row = 0
ax = fig_figure1.add_subplot(gs[row, column])


rects1 = ax.bar(x - 1/5, cmip6_deltaCs_1pctCO2_4xCO2, width, color='b', label='Cs')
ax.bar(x + 1/5, deltaCs_estimated_4xCO2, width, color='lightseagreen')

ax.set_ylabel(r'$\Delta C_{s}$ (PgC)')
ax.set_xticks(x)
ax.set_xticklabels(label_list)
ax.set_ylim((-100, 550))
ax.set_title(r'(b) 4xCO$_{2}$', y=1.09, fontweight='bold')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')




# legend
handels = []
handels.extend([Line2D([0,0],[0,0], linewidth=20, color='b', label=r'1% CO$_{2}$')])
#handels.extend([Line2D([0,0],[0,0], linewidth=20, color='lightsteelblue', label=r'BGC + RAD')])
handels.extend([Line2D([0,0],[0,0], linewidth=20, color='lightseagreen', label=r'BGC + RAD')])
labels = [r'$\Delta C_{s}$ 1% CO$_{2}$', r'$\beta_{s} \Delta CO_{2}$ + $\gamma_{s} \Delta$T']
leg1 = ax.legend(handels, labels, loc='lower center', ncol=7, bbox_to_anchor=(-0.25, -0.75), fontsize=58)#title='Model colors',
plt.gca().add_artist(leg1)


#%%
fig_figure1.tight_layout()
fig_figure1.savefig('figures/fig4', bbox_inches='tight')
plt.close()

