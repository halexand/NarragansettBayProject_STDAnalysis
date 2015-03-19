#!/usr/bin/env python

'''
Created on November 9, 2013

Comparison script for Throt and Skcos count data; calls upon general and std functions 

@author: harrietalexander
'''
import math
import sys
import glob
import os
import numpy as np
import re
import csv
import matplotlib.pylab as plt
import matplotlib.cm as cm
import STD_Functions as STD
import General_STD_analysis_control as GSTD
import pickle

#Relies upon the Genearl and STD_functions python files. 
#Location of pickled files
PickleOutDir="/Users/harrietalexander/Dropbox/NB_141020/PickledData/"

#Throt_filenames from pickle
Throt_Pk_Count=['Throt_raw', 'Throt_tpm', 'Throt_SGNC', 'Throt_GeneList']

#Skcos Raw Data

#Skcos filenames from pickle
Skcos_Pk_Count=['Skcos_raw', 'Skcos_tpm', 'Skcos_SGNC', 'Skcos_GeneList']

#Get STD_all-- from pickle. 
[Throt_raw, Throt_tpm, Throt_SGNC, Throt_GeneList]=STD.unPickleAll(Throt_Pk_Count, PickleOutDir)

[Skcos_raw, Skcos_tpm, Skcos_SGNC, Skcos_GeneList]=STD.unPickleAll(Skcos_Pk_Count, PickleOutDir)
N=0
print Skcos_GeneList['N_Genes'][1]
#---- P-related genes *************************************

snames=['S1', 'S2', 'S3', 'S4', 'S5']
inames=['+N','-N','+P','-P', 'C']
x=range(10)
outDir='/Users/harrietalexander/Dropbox/NB_141020/Figures/SingleGene_Blast/StableGenes/'


N+=1
F=plt.figure(N,figsize=(10,10))
ax=F.add_subplot(111)
# Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
ax.set_axis_bgcolor("none")

# Set common labels
ax.set_xlabel('Sample',fontsize=14)
ax.set_ylabel('TPM', fontsize=14)

ax1=F.add_subplot(211)
ax2=F.add_subplot(212)
# ax3=F.add_subplot(223)
# ax4=F.add_subplot(224)
#ax1.set_title('T. rotula')
ax1.set_xticks(x[:5])
ax1.set_xticklabels(snames)
#ax2.set_title('S. costatum')
ax2.set_xticks(x[:5])
ax2.set_xticklabels(snames)
ax1.set_xlim([-.5, 4.5])
ax2.set_xlim([4.5, 9.5])

# ax3.set_xticklabels(inames)
# ax3.set_xticks(x[5:])
# ax4.set_xticks(x[5:])
# ax4.set_xticklabels(inames)
ThrotSG=Throt_GeneList['Stable_genes']
print ThrotSG
NRG='0'
for g in ThrotSG:
	ax1.plot(x[:5],Throt_tpm[g][:5], ls='--',color=NRG, marker='o',mec=NRG, ms=12,linewidth=3)
	ax2.plot(x[5:],Throt_tpm[g][5:],ls='None',color=NRG, marker='o',mec=NRG, ms=12,linewidth=3)
SkcosSG=Skcos_GeneList['Stable_genes']
print SkcosSG
for g in SkcosSG:
	ax1.plot(x[:5],Skcos_tpm[g][:5], ls='-',color=NRG, marker='s',mec=NRG, ms=12,linewidth=3)
	ax2.plot(x[5:],Skcos_tpm[g][5:],ls='None',color=NRG, marker='s',mec=NRG, ms=12,linewidth=3)

ax1.set_ylim(bottom=0)
ax2.set_ylim(bottom=0)
ax2.set_ylim(bottom=0)
ax2.set_ylim(bottom=0)

#plt.suptitle(ref,fontsize=16)

plt.savefig(outDir+'StableGenes_tpm.pdf')

plt.show()
