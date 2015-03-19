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
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

if len(sys.argv)<3:
	print "Usage: python Plot_KeggData_pathway.py KO_List_metabolic.tab KO_CountsFile.tab"
	sys.exit(1)

Kegg_file=sys.argv[1]
KeggCounts=sys.argv[2]	

## USE PANDAS TO LOAD DATA FRAME VERSION OF TWO COUNT DATA
# Throt_data= pd.read_csv('Throt_KassCounts_PercentTotalKO_NoUnknown.txt', sep='\t', index_col=0)
# Skcos_data= pd.read_csv('Skcos_KassCounts_PercentTotalKO_NoUnknown.txt', sep='\t', index_col=0)
# All_data= pd.read_csv('All_KassCounts_NoUnknown.txt', sep='\t', index_col=0)
All_data=pd.read_csv(KeggCounts, sep='\t', index_col=0)
## LOAD LIST OF KO NUMBERS FOR HEATMAP


with open(Kegg_file) as f:
    KeggValues = f.read().splitlines()

KeggValues=[x.strip() for x in KeggValues]

## Convert to TPM
# Throt_data=np.multiply(Throt_data, 1000000)
# Skcos_data=np.multiply(Skcos_data, 1000000)
All_data=np.multiply(All_data, 1000000)
##Set Max and min values
## Get Columns of interest
colList_insitu=['SS1','SS2','SS3', 'SS4', 'SS5','TS1','TS2','TS3', 'TS4', 'TS5']
colList_incubation=['SA','SB','SC', 'SD', 'SE','TA','TB','TC', 'TD', 'TE']
#colList=['SS1','SS2','SS3', 'SS4', 'SS5','TS1','TS2','TS3', 'TS4', 'TS5','SA','SB','SC', 'SD', 'SE','TA','TB','TC', 'TD', 'TE']
colList=['S1', 'S2','S3','S4','S5','A','B','C','D','E']
All_insitu=All_data[colList]



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

Sum_All_insitu=All_insitu.sum()
All_insitu_Percent=All_insitu/Sum_All_insitu #calculate the percent of total kegg mapped data

All_insitu_Percent['mean']=All_insitu_Percent.mean(skipna=1, axis=1) #calculate mean value for each class

All_insitu_Percent=All_insitu_Percent.sort(column='mean', ascending=False)#Sort by the mean value

All_insitu_Percent=All_insitu_Percent.drop('mean',1) #drop mean column

df=All_insitu_Percent.T
col=df.columns
KeggValues2=[]
for k in KeggValues:
	if k in col:
		KeggValues2.append(k)
All_insitu_Percent=df[KeggValues2].T

maxV=0.01
minV=1e-6
col_labels=list(All_insitu_Percent.index)
row_labels=list(All_insitu_Percent.columns.values)
fig3,ax3=plt.subplots()
#heatmap3 = ax3.pcolor(All_insitu_Percent, cmap=plt.cm.jet, norm=LogNorm(vmin=m, vmax=All_insitu_Percent.max().max()))
heatmap3 = ax3.pcolor(All_insitu_Percent, cmap=plt.cm.jet, norm=LogNorm(vmin=minV, vmax=maxV))
print All_insitu_Percent
print All_insitu_Percent.max().max()
ax3.set_xticks(np.arange(All_insitu_Percent.shape[1])+0.5, minor=False)
ax3.set_yticks(np.arange(All_insitu_Percent.shape[0])+0.5, minor=False)
ax3.invert_yaxis()
ax3.xaxis.tick_top()

ax3.set_xticklabels(row_labels, minor=False)
ax3.set_yticklabels(col_labels, minor=False)
plt.colorbar(heatmap3)

# fig4,ax4=plt.subplots()

# plt.scatter(All_insitu_Percent[colList[0:4]], All_insitu_Percent[colList[5:9]],norm=LogNorm(vmin=m, vmax=All_insitu_Percent.max().max()))


plt.show()