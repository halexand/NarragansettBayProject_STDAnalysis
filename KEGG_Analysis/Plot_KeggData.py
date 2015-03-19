#!/usr/bin/env python

'''
Created on November 9, 2013

STD_Function.py

A generalized collection of 

@author: harrietalexander
'''

import sys
import glob
import os
import numpy as np
import pandas as pd
import re
import csv
import matplotlib.pyplot as plt
import pickle
import matplotlib.mlab as mlab
from matplotlib.colors import LogNorm, NoNorm
import matplotlib.colors as mcolors
import matplotlib as mpl
#from hcluster import pdist, linkage, dendrogram
import Colormap as cmap

c = mcolors.ColorConverter().to_rgb
mpl.rcParams['pdf.fonttype'] = 42
## USE PANDAS TO LOAD DATA FRAME VERSION OF TWO COUNT DATA
#Throt_data= pd.read_csv('Throt_KassCounts_PercentTotalKO_NoUnknown.txt', sep='\t', index_col=0)
#Skcos_data= pd.read_csv('Skcos_KassCounts_PercentTotalKO_NoUnknown.txt', sep='\t', index_col=0)
All_data= pd.read_csv('redo/All_NoUnknown_TPM.tab', sep='\t', index_col=0)
## Convert to TPM
# Throt_data=np.multiply(Throt_data, 1000000)
# Skcos_data=np.multiply(Skcos_data, 1000000)
All_data=np.multiply(All_data, 1000000)

##Set Max and min values
# m=.01
# MM=Skcos_data.max().max()

## Get Columns of interest
colList=['SS1','SS2','SS3', 'SS4', 'SS5','TS1','TS2','TS3', 'TS4', 'TS5']
colList_incubation=['SA','SB','SC', 'SD', 'SE','TA','TB','TC', 'TD', 'TE']
All_insitu=All_data
## PLOT SKCOS
# 
# col_labels=list(Skcos_data.index)
# row_labels=list(Skcos_data.columns.values)
# fig,ax=plt.subplots()
# heatmap = ax.pcolor(Skcos_data, cmap=plt.cm.jet, norm=LogNorm(vmin=m, vmax=MM))
# 
# ax.set_xticks(np.arange(Skcos_data.shape[1])+0.5, minor=False)
# ax.set_yticks(np.arange(Skcos_data.shape[0])+0.5, minor=False)
# ax.invert_yaxis()
# ax.xaxis.tick_top()
# 
# 
# ax.set_xticklabels(row_labels, minor=False)
# ax.set_yticklabels(col_labels, minor=False)
# 
# plt.colorbar(heatmap)
# 
# 
# ## PLOT THROT
# col_labels=list(Throt_data.index)
# row_labels=list(Throt_data.columns.values)
# fig2,ax2=plt.subplots()
# heatmap2 = ax2.pcolor(Throt_data, cmap=plt.cm.jet, norm=LogNorm(vmin=m, vmax=MM))
# 
# ax2.set_xticks(np.arange(Throt_data.shape[1])+0.5, minor=False)
# ax2.set_yticks(np.arange(Throt_data.shape[0])+0.5, minor=False)
# ax2.invert_yaxis()
# ax2.xaxis.tick_top()
# 
# 
# ax2.set_xticklabels(row_labels, minor=False)
# ax2.set_yticklabels(col_labels, minor=False)
# 
# plt.colorbar(heatmap2)


# 
# ## PLOT COMBIN
# col_labels=list(All_insitu.index)
# row_labels=list(All_insitu.columns.values)
# fig3,ax3=plt.subplots()
# heatmap3 = ax3.pcolor(All_insitu, cmap=plt.cm.jet, norm=LogNorm(vmin=m, vmax=All_insitu.max().max()))
# 
# ax3.set_xticks(np.arange(All_insitu.shape[1])+0.5, minor=False)
# ax3.set_yticks(np.arange(All_insitu.shape[0])+0.5, minor=False)
# ax3.invert_yaxis()
# ax3.xaxis.tick_top()
# 
# 
# ax3.set_xticklabels(row_labels, minor=False)
# ax3.set_yticklabels(col_labels, minor=False)
# 
# plt.colorbar(heatmap3)


## PLOT COMBIN TOTAL MAPPING
m=1e-5
# rvb = make_colormap([c('#32cd32'), c('#00ffff'), 0.25, c('violet'), c('blue'), 0.75, c('blue')])

Sum_All_insitu=All_insitu.sum()
All_insitu_Percent=All_insitu/Sum_All_insitu #calculate the percent of total kegg mapped data
sckos=['SS1','SS2','SS3', 'SS4', 'SS5']
All_insitu_Percent['mean']=All_insitu_Percent[sckos].mean(skipna=1, axis=1) #calculate mean value for each class

All_insitu_Percent=All_insitu_Percent.sort(columns='mean', ascending=False)#Sort by the mean value

All_insitu_Percent=All_insitu_Percent.drop('mean',1) #drop mean column

#Blue to cyan to yellow
rvb = cmap.make_colormap([c('#000000'), c('#00FFFF'), 0.5, c('#00FFFF'), c('#FFFF00')])

# #yellow to blue
# rvb = cmap.make_colormap([c('#ffffd9'), c('#41b6c4'), 0.5, c('#41b6c4'), c('#081d58')])

# #Purple to green
# rvb = cmap.make_colormap([c('#4d9221'), c('#f7f7f7'), 0.5, c('#f7f7f7'), c('#c51b7d')])

# #White to purple
# rvb = cmap.make_colormap([c('#fcfbfd'), c('#9e9ac8'), 0.5, c('#9e9ac8'), c('#3f007d')])

# #White to green
# rvb = cmap.make_colormap([c('#e5f5e0'), c('#74c476'), 0.5, c('#74c476'), c('#00441b')])

# #White to orange
# rvb = cmap.make_colormap([c('#fff5eb'), c('#fd8d3c'), 0.5, c('#fd8d3c'), c('#7f2704')])

# #White to black
# rvb = cmap.make_colormap([c('#ffffff'), c('#969696'), 0.5, c('#969696'), c('#000000')])
# 
# #White to black
# rvb = plt.cm.hot_r


col_labels=list(All_insitu_Percent.index)
row_labels=list(All_insitu_Percent.columns.values)
fig3,ax3=plt.subplots()
heatmap3 = ax3.pcolor(All_insitu_Percent, cmap=rvb, norm=LogNorm(vmin=m, vmax=All_insitu_Percent.max().max()))
#heatmap3 = ax3.pcolor(All_insitu_Percent, cmap=plt.cm.jet, vmin=0, vmax=.25)
ax3.set_xticks(np.arange(All_insitu_Percent.shape[1])+0.5, minor=False)
ax3.set_yticks(np.arange(All_insitu_Percent.shape[0])+0.5, minor=False)
ax3.invert_yaxis()
ax3.xaxis.tick_top()


ax3.set_xticklabels(row_labels, minor=False)
ax3.set_yticklabels(col_labels, minor=False)
plt.colorbar(heatmap3)

# fig4,ax4=plt.subplots()

# plt.scatter(All_insitu_Percent[colList[0:4]], All_insitu_Percent[colList[5:9]],norm=LogNorm(vmin=m, vmax=All_insitu_Percent.max().max()))
#y=pdist(All_insitu_Percent.T)
#z=linkage(y)
#fig4,ax4=plt.subplots()
#dendrogram(z, labels=All_insitu_Percent.columns)
plt.show()