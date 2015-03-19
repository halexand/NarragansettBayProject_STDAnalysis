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
[Throt_PupP_STD, Throt_PupN_STD, Throt_PdnP_STD, Throt_PdnN_STD, Throt_NupP_STD, Throt_NupN_STD, Throt_NdnP_STD, Throt_NdnN_STD, Throt_All_P, Throt_All_N]=STD.unPickleAll(Throt_Pk_STD, PickleOutDir)
[Skcos_PupP_STD, Skcos_PupN_STD, Skcos_PdnP_STD, Skcos_PdnN_STD, Skcos_NupP_STD, Skcos_NupN_STD, Skcos_NdnP_STD, Skcos_NdnN_STD, Skcos_All_P, Skcos_All_N]=STD.unPickleAll(Skcos_Pk_STD, PickleOutDir)
[Skcos_raw, Skcos_tpm, Skcos_SGNC, Skcos_GeneList]=STD.unPickleAll(Skcos_Pk_Count, PickleOutDir)

N=0

# #---- P-related genes *************************************
# 
# Pfile='/Users/harrietalexander/Dropbox/NB_141020/Thaps_P_uniq.tab'
# Preader=csv.reader(open(Pfile),delimiter='\t')
# Phash={}
# for line in Preader:
# 	key = line[0]
# 	item=line[1]
# 	if key in Phash: 
# 		Phash[key].append(line[1])
# 	else:
# 		Phash[key]=[item]
# 
# snames=['S1', 'S2', 'S3', 'S4', 'S5','C']
# x=range(6)
# outDir='/Users/harrietalexander/Dropbox/NB_141020/Figures/SingleGene_Blast/PGenes/'
# 
# for ref in Phash:
# 	N+=1
# 	F=plt.figure(N)
# 	ax=F.add_subplot(111)
# 	# Turn off axis lines and ticks of the big subplot
# 	ax.spines['top'].set_color('none')
# 	ax.spines['bottom'].set_color('none')
# 	ax.spines['left'].set_color('none')
# 	ax.spines['right'].set_color('none')
# 	ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
# 	ax.set_axis_bgcolor("none")
# 	
# 	# Set common labels
# 	ax.set_xlabel('Sample',fontsize=14)
# 	ax.set_ylabel('STD-N', fontsize=14)
# 	ax1=F.add_subplot(121)
# 	ax2=F.add_subplot(122)
# 	for item in Phash[ref]:
# 		ax1.set_title('T. rotula')
# 		ax1.set_xticks(x)
# 		ax1.set_xticklabels(snames)
# 		ax2.set_title('S. costatum')
# 		ax2.set_xticks(x)
# 		ax2.set_xticklabels(snames)
# 		
# 		if re.search('Throt', item):	
# 			if any(item==s for s in Throt_GeneList['P_Genes'][0]):
# 			 	ax1.plot(x,Throt_All_N[item], color='r')
# 			elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
# 			 	ax1.plot(x,Throt_All_N[item], '--', color='r')
# 			elif any(item==s for s in Throt_GeneList['N_Genes'][0]):
# 			 	ax1.plot(x,Throt_All_N[item], color='b')
# 			elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
# 			 	ax1.plot(x,Throt_All_N[item], '--', color='b')
# 			else:
# 				ax1.plot(x,Throt_All_N[item], color='0.75')
# 			 	
# 
# 		elif re.search('Skcos', item):
# 			if any(item==s for s in Skcos_GeneList['P_Genes'][0]):
# 			 	ax2.plot(x,Skcos_All_N[item], color='r')
# 			elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
# 			 	ax2.plot(x,Skcos_All_N[item], '--', color='r')
# 			elif any(item==s for s in Skcos_GeneList['N_Genes'][0]):
# 			 	ax2.plot(x,Skcos_All_N[item], color='b')
# 			elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
# 			 	ax2.plot(x,Skcos_All_N[item], '--', color='b')
# 			else:
# 				ax2.plot(x,Skcos_All_N[item], color='0.75')
# 	plt.suptitle(ref,fontsize=16)
# 
# 	plt.savefig(outDir+ref+'_SingleGenes_STDN.pdf')
# 
# #---- N-related genes *************************************
# 
Pfile='/Users/harrietalexander/Dropbox/NB_141020/Thaps_N_uniq.tab'
Preader=csv.reader(open(Pfile),delimiter='\t')
Phash={}
for line in Preader:
	key = line[0]
	item=line[1]
	if key in Phash: 
		Phash[key].append(line[1])
	else:
		Phash[key]=[item]


snames=['S1', 'S2', 'S3', 'S4', 'S5','C']
x=range(6)
outDir='/Users/harrietalexander/Dropbox/NB_141020/Figures/SingleGene_Blast/NGenes/'

for ref in Phash:
	N+=1
	F=plt.figure(N)
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
	ax.set_ylabel('STD-N', fontsize=14)
	ax1=F.add_subplot(121)
	ax2=F.add_subplot(122)
	for item in Phash[ref]:
		ax1.set_title('T. rotula')
		ax1.set_xticks(x)
		ax1.set_xticklabels(snames)
		ax2.set_title('S. costatum')
		ax2.set_xticks(x)
		ax2.set_xticklabels(snames)
		ax1.set_ylim([-10,10])
		if re.search('Throt', item):	
			if any(item==s for s in Throt_GeneList['N_Genes'][0]):
			 	ax1.plot(x,Throt_All_N[item], color='b')
			elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
			 	ax1.plot(x,Throt_All_N[item], '--', color='b')
			
			elif any(item==s for s in Throt_GeneList['P_Genes'][0]):
			 	ax1.plot(x,Throt_All_N[item], color='r')
			elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
			 	ax1.plot(x,Throt_All_N[item], '--', color='r')
			else:
				ax1.plot(x,Throt_All_N[item], color='0.75')
			 	

		elif re.search('Skcos', item):
			if any(item==s for s in Skcos_GeneList['N_Genes'][0]):
			 	ax2.plot(x,Skcos_All_N[item], color='b')
			elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
			 	ax2.plot(x,Skcos_All_N[item], '--', color='b')

			elif any(item==s for s in Skcos_GeneList['P_Genes'][0]):
			 	ax2.plot(x,Skcos_All_N[item], color='r')
			elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
			 	ax2.plot(x,Skcos_All_N[item], '--', color='r')
			else:
				ax2.plot(x,Skcos_All_N[item], color='0.75')
	plt.suptitle(ref,fontsize=16)

	plt.savefig(outDir+ref+'_SingleGenes_STDN.pdf')
# 

# 
# #---- Si-related genes *************************************
# 
Pfile='/Users/harrietalexander/Dropbox/NB_141020/Thaps_DeathGene.tab'
Preader=csv.reader(open(Pfile),delimiter='\t')
Phash={}
for line in Preader:
	key = line[0]
	item=line[1]
	if key in Phash: 
		Phash[key].append(line[1])
	else:
		Phash[key]=[item]



outDir='/Users/harrietalexander/Dropbox/NB_141020/Figures/SingleGene_Blast/SiGenes/'
snames=['S1', 'S2', 'S3', 'S4', 'S5','C']
x=range(6)

for ref in Phash:
	N+=1
	F=plt.figure(N)
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
	ax.set_ylabel('STD-N', fontsize=14)
	ax1=F.add_subplot(121)
	ax2=F.add_subplot(122)
	for item in Phash[ref]:
		ax1.set_title('T. rotula')
		ax1.set_xticks(x)
		ax1.set_xticklabels(snames)
		ax2.set_title('S. costatum')
		ax2.set_xticks(x)
		ax2.set_xticklabels(snames)
		ax1.set_ylim([-10,10])
		if re.search('Throt', item):	
			if any(item==s for s in Throt_GeneList['N_Genes'][0]):
			 	ax1.plot(x,Throt_All_N[item], color='b')
			elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
			 	ax1.plot(x,Throt_All_N[item], '--', color='b')
			elif any(item==s for s in Throt_GeneList['P_Genes'][0]):
			 	ax1.plot(x,Throt_All_N[item], color='r')
			elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
			 	ax1.plot(x,Throt_All_N[item], '--', color='r')
			else:
				ax1.plot(x,Throt_All_N[item], color='0.75')
			 	

		elif re.search('Skcos', item):
			if any(item==s for s in Skcos_GeneList['N_Genes'][0]):
			 	ax2.plot(x,Skcos_All_N[item], color='b')
			elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
			 	ax2.plot(x,Skcos_All_N[item], '--', color='b')

			elif any(item==s for s in Skcos_GeneList['P_Genes'][0]):
			 	ax2.plot(x,Skcos_All_N[item], color='r')
			elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
			 	ax2.plot(x,Skcos_All_N[item], '--', color='r')
			else:
				ax2.plot(x,Skcos_All_N[item], color='0.75')
	plt.suptitle(ref,fontsize=16)

	plt.savefig(outDir+ref+'_SingleGenes_STDN.pdf')

plt.show()

	