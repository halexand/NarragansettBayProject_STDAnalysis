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
Throt_Pk_STD=['Throt_PupP_STD', 'Throt_PupN_STD', 'Throt_PdnP_STD', 'Throt_PdnN_STD', 'Throt_NupP_STD', 'Throt_NupN_STD', 'Throt_NdnP_STD', 'Throt_NdnN_STD', 'Throt_All_P', 'Throt_All_N']

#Skcos Raw Data

#Skcos filenames from pickle
Skcos_Pk_Count=['Skcos_raw', 'Skcos_tpm', 'Skcos_SGNC', 'Skcos_GeneList']
Skcos_Pk_STD=['Skcos_PupP_STD', 'Skcos_PupN_STD', 'Skcos_PdnP_STD', 'Skcos_PdnN_STD', 'Skcos_NupP_STD', 'Skcos_NupN_STD', 'Skcos_NdnP_STD', 'Skcos_NdnN_STD', 'Skcos_All_P', 'Skcos_All_N']

#Get STD_all-- from pickle. 
[Throt_raw, Throt_tpm, Throt_SGNC, Throt_GeneList]=STD.unPickleAll(Throt_Pk_Count, PickleOutDir)



[Skcos_raw, Skcos_tpm, Skcos_SGNC, Skcos_GeneList]=STD.unPickleAll(Skcos_Pk_Count, PickleOutDir)
[Throt_PupP_STD, Throt_PupN_STD, Throt_PdnP_STD, Throt_PdnN_STD, Throt_NupP_STD, Throt_NupN_STD, Throt_NdnP_STD, Throt_NdnN_STD, Throt_All_P, Throt_All_N]=STD.unPickleAll(Throt_Pk_STD, PickleOutDir)
[Skcos_PupP_STD, Skcos_PupN_STD, Skcos_PdnP_STD, Skcos_PdnN_STD, Skcos_NupP_STD, Skcos_NupN_STD, Skcos_NdnP_STD, Skcos_NdnN_STD, Skcos_All_P, Skcos_All_N]=STD.unPickleAll(Skcos_Pk_STD, PickleOutDir)



#---- Si-related genes *************************************

Pfile='/Users/harrietalexander/Dropbox/NB_141020/Reference_Genes_THaps.tab'
Preader=csv.reader(open(Pfile),delimiter='\t')

Phash={}
for line in Preader:
	key = line[0].strip()
	item=line[1].strip()
	if key in Phash: 
		Phash[key].append(item)
	else:
		Phash[key]=[item]


snames=['S1', 'S2', 'S3', 'S4', 'S5']
inames=['+N','-N','+P','-P', 'C']
x=range(10)
outDir='/Users/harrietalexander/Dropbox/NB_141020/FinalSingleGeneFigures/'
N=1
Throt_hash=Throt_tpm
Skcos_hash=Skcos_tpm
for ref in Phash:
	N+=1
	F=plt.figure(N, figsize=(10,10))
	ax=F.add_subplot(111)
	# Turn off axis lines and ticks of the big subplot
	ax.spines['top'].set_color('none')
	ax.spines['bottom'].set_color('none')
	ax.spines['left'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
	ax.set_axis_bgcolor("none")
# 
# 	# Set common labels

	ax1=F.add_subplot(211)
 	ax2=F.add_subplot(212)


 	ax1.set_xlabel('Sample',fontsize=14)
 	ax1.set_ylabel('TPM', fontsize=14)
	
	for item in Phash[ref]:
		ax1.set_xlim([-.5, 4.5])
		ax1.set_xticks(x[:5])
		ax1.set_xticklabels(snames)
		ax2.set_xlim([4.5, 9.5])
		ax2.set_xticks(x[5:])
		ax2.set_xticklabels(inames)
# 		ax1.set_ylim([-50,50])
		##ax1.set_ylim([0,10])
		#ax2.set_ylim([0,10])
		NUp='#000099'
		NDn='#00FFFF'
		PUp='#FF0000'
		PDn='#FF9900'
		NRG='0.75'
		if re.search('Throt', item):	
			if any(item==s for s in Throt_GeneList['P_Genes'][0]):
			 	ax1.plot(x[:5],Throt_hash[item][:5], ls='--', color=PUp, marker='o', mec=PUp, ms=12,linewidth=3)
			 	ax2.plot(x[5:],Throt_hash[item][5:],ls='None',color=PUp, marker='o', mec=PUp, ms=12,linewidth=3)
			elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
			 	ax1.plot(x[:5],Throt_hash[item][:5],ls='--', color=PDn, marker='o', mec=PDn, ms=12,linewidth=3)
			 	ax2.plot(x[5:],Throt_hash[item][5:],ls='None',color=PDn, marker='o', mec=PDn, ms=12,linewidth=3)
			elif any(item==s for s in Throt_GeneList['N_Genes'][0]):
			 	ax1.plot(x[:5],Throt_hash[item][:5], ls='--',color=NUp, marker='o', mec=NUp, ms=12,linewidth=3)
			 	ax2.plot(x[5:],Throt_hash[item][5:],ls='None', color=NUp, marker='o', mec=NUp, ms=12,linewidth=3)
			elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
			 	ax1.plot(x[:5],Throt_hash[item][:5], ls='--', color=NDn, marker='o', mec=NDn, ms=12,linewidth=3)
			 	ax2.plot(x[5:],Throt_hash[item][5:],ls='None',color=NDn, marker='o', mec=NDn, ms=12,linewidth=3)
			else:
				ax1.plot(x[:5],Throt_hash[item][:5], ls='--',color=NRG, marker='o',mec=NRG, ms=12,linewidth=3)
			 	ax2.plot(x[5:],Throt_hash[item][5:],ls='None',color=NRG, marker='o',mec=NRG, ms=12,linewidth=3)
			 	

		elif re.search('Skcos', item):
			if any(item==s for s in Skcos_GeneList['P_Genes'][0]):
			 	ax1.plot(x[:5],Skcos_hash[item][:5], ls='-', color=PUp, marker='s', mec=PUp, ms=12,linewidth=3)
# 			 	ax2.plot(x[5:],Skcos_hash[item][5:],ls='None',color=PUp, marker='s', mec=PUp, ms=12,linewidth=3)
			elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
			 	ax1.plot(x[:5],Skcos_hash[item][:5],ls='-', color=PDn, marker='s', mec=PDn, ms=12,linewidth=3)
# 			 	ax2.plot(x[5:],Skcos_hash[item][5:],ls='None',color=PDn, marker='s', mec=PDn, ms=12,linewidth=3)
			elif any(item==s for s in Skcos_GeneList['N_Genes'][0]):
			 	ax1.plot(x[:5],Skcos_hash[item][:5], ls='-',color=NUp, marker='s', mec=NUp, ms=12,linewidth=3)
# 			 	ax2.plot(x[5:],Skcos_hash[item][5:],ls='None', color=NUp, marker='s', mec=NUp, ms=12,linewidth=3)
			elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
			 	ax1.plot(x[:5],Skcos_hash[item][:5], ls='-', color=NDn, marker='s', mec=NDn, ms=12,linewidth=3)
# 			 	ax2.plot(x[5:],Skcos_hash[item][5:],ls='None',color=NDn, marker='s', mec=NDn, ms=12,linewidth=3)
			else:
				ax1.plot(x[:5],Skcos_hash[item][:5], ls='-',color=NRG, marker='s',mec=NRG, ms=12,linewidth=3)
# 			 	ax2.plot(x[5:],Skcos_hash[item][5:],ls='None',color=NRG, marker='s',mec=NRG, ms=12,linewidth=3)
	plt.suptitle(ref,fontsize=16)
	plt.savefig(outDir+ref+'_SingleGenes_tpm.pdf')


plt.show()

	