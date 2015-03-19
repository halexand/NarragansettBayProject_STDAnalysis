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

# for s in Skcos_GeneList['N_Genes'][0]:
# 	print s
# print len(Skcos_GeneList['N_Genes'][0])
# N=0
# 
# #---- P-related genes *************************************
# 
# Pfile='/Users/harrietalexander/Dropbox/NB_141020/Thaps_P_uniq.tab'
# Preader=csv.reader(open(Pfile),delimiter='\t')
# for key in Throt_GeneList:
# 	print key
# Phash={}
# for line in Preader:
# 	key = line[0]
# 	item=line[1]
# 	if key in Phash: 
# 		Phash[key].append(line[1])
# 	else:
# 		Phash[key]=[item]
# 
# snames=['S1', 'S2', 'S3', 'S4', 'S5']
# inames=['+N','-N','+P','-P', 'C']
# x=range(10)
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
# 	ax.set_ylabel('TPM', fontsize=14)
# 	
# 	ax1=F.add_subplot(221)
# 	ax2=F.add_subplot(222)
# 	ax3=F.add_subplot(223)
# 	ax4=F.add_subplot(224)
# 	for item in Phash[ref]:
# 		ax1.set_title('T. rotula')
# 		ax1.set_xticks(x[:5])
# 		ax1.set_xticklabels(snames)
# 		ax2.set_title('S. costatum')
# 		ax2.set_xticks(x[:5])
# 		ax2.set_xticklabels(snames)
# 		ax3.set_xticklabels(inames)
# 		ax3.set_xticks(x[5:])
# 		ax4.set_xticks(x[5:])
# 		ax4.set_xticklabels(inames)
# 		
# 		if re.search('Throt', item):	
# 			if any(item==s for s in Throt_GeneList['P_Genes'][0]):
# 			 	ax1.plot(x[:5],Throt_tpm[item][:5], color='r')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='r', marker='o', mec='r', ms=8)
# 			elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
# 			 	ax1.plot(x[:5],Throt_tpm[item][:5], '--', color='r')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='r', marker='s', mec='r', ms=8)
# 			elif any(item==s for s in Throt_GeneList['N_Genes'][0]):
# 			 	ax1.plot(x[:5],Throt_tpm[item][:5], color='b')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='b', marker='o', mec='b', ms=8)
# 			elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
# 			 	ax1.plot(x[:5],Throt_tpm[item][:5], '--', color='b')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='b', marker='s', mec='b', ms=8)
# 			else:
# 				ax1.plot(x[:5],Throt_tpm[item][:5], color='0.75')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],color='0.5', ls='None', marker='x', ms=8)
# 			 	
# 
# 		elif re.search('Skcos', item):
# 			if any(item==s for s in Skcos_GeneList['P_Genes'][0]):
# 			 	ax2.plot(x[:5],Skcos_tpm[item][:5], color='r')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='r', marker='o', mec='r', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
# 			 	ax2.plot(x[:5],Skcos_tpm[item][:5], '--', color='r')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='r', marker='s', mec='r', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['N_Genes'][0]):
# 			 	ax2.plot(x[:5],Skcos_tpm[item][:5], color='b')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='b', marker='o', mec='b', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
# 			 	ax2.plot(x[:5],Skcos_tpm[item][:5], '--', color='b')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='b', marker='s', mec='b', ms=8)
# 			else:
# 				ax2.plot(x[:5],Skcos_tpm[item][:5], color='0.75')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],color='0.5', ls='None', marker='x', ms=8)
# 	plt.suptitle(ref,fontsize=16)
# 
# 	plt.savefig(outDir+ref+'_SingleGenes_TPM.pdf')
# 
# #---- N-related genes *************************************
# 
# Pfile='/Users/harrietalexander/Dropbox/NB_141020/Thaps_N_uniq.tab'
# Preader=csv.reader(open(Pfile),delimiter='\t')
# for key in Throt_GeneList:
# 	print key
# Phash={}
# for line in Preader:
# 	key = line[0]
# 	item=line[1]
# 	if key in Phash: 
# 		Phash[key].append(line[1])
# 	else:
# 		Phash[key]=[item]
# 
# 
# snames=['S1', 'S2', 'S3', 'S4', 'S5']
# inames=['+N','-N','+P','-P', 'C']
# x=range(10)
# outDir='/Users/harrietalexander/Dropbox/NB_141020/Figures/SingleGene_Blast/NGenes/'
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
# 	ax.set_ylabel('TPM', fontsize=14)
# 	
# 	ax1=F.add_subplot(221)
# 	ax2=F.add_subplot(222)
# 	ax3=F.add_subplot(223)
# 	ax4=F.add_subplot(224)
# 	for item in Phash[ref]:
# 		ax1.set_title('T. rotula')
# 		ax1.set_xticks(x[:5])
# 		ax1.set_xticklabels(snames)
# 		ax2.set_title('S. costatum')
# 		ax2.set_xticks(x[:5])
# 		ax2.set_xticklabels(snames)
# 		ax3.set_xticklabels(inames)
# 		ax3.set_xticks(x[5:])
# 		ax4.set_xticks(x[5:])
# 		ax4.set_xticklabels(inames)
# 		
# 		if re.search('Throt', item):	
# 			
# 			if any(item==s for s in Throt_GeneList['N_Genes'][0]):
# 			 	ax1.plot(x[:5],Throt_tpm[item][:5], color='b')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='b', marker='o', mec='b', ms=8)
# 			elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
# 			 	ax1.plot(x[:5],Throt_tpm[item][:5], '--', color='b')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='b', marker='s', mec='b', ms=8)
# 			elif any(item==s for s in Throt_GeneList['P_Genes'][0]):
# 			 	ax1.plot(x[:5],Throt_tpm[item][:5], color='r')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='r', marker='o', mec='r', ms=8)
# 			elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
# 			 	ax1.plot(x[:5],Throt_tpm[item][:5], '--', color='r')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='r', marker='s', mec='r', ms=8)
# 			else:
# 				ax1.plot(x[:5],Throt_tpm[item][:5], color='0.75')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],color='0.5', ls='None', marker='x', ms=8)
# 			 	
# 
# 		elif re.search('Skcos', item):
# 			if any(item==s for s in Skcos_GeneList['N_Genes'][0]):
# 			 	ax2.plot(x[:5],Skcos_tpm[item][:5], color='b')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='b', marker='o', mec='b', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
# 			 	ax2.plot(x[:5],Skcos_tpm[item][:5], '--', color='b')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='b', marker='s', mec='b', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['P_Genes'][0]):
# 			 	ax2.plot(x[:5],Skcos_tpm[item][:5], color='r')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='r', marker='o', mec='r', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
# 			 	ax2.plot(x[:5],Skcos_tpm[item][:5], '--', color='r')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='r', marker='s', mec='r', ms=8)
# 
# 			else:
# 				ax2.plot(x[:5],Skcos_tpm[item][:5], color='0.75')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],color='0.5', ls='None', marker='x', ms=8)
# 	plt.suptitle(ref,fontsize=16)
# 	plt.savefig(outDir+ref+'_SingleGenes_TPM.pdf')
# 
# 
# 
# #---- Si-related genes *************************************
# 
# Pfile='/Users/harrietalexander/Dropbox/NB_141020/Thaps_Si_uniq.tab'
# Preader=csv.reader(open(Pfile),delimiter='\t')
# for key in Throt_GeneList:
# 	print key
# Phash={}
# for line in Preader:
# 	key = line[0]
# 	item=line[1]
# 	if key in Phash: 
# 		Phash[key].append(line[1])
# 	else:
# 		Phash[key]=[item]
# 
# 
# snames=['S1', 'S2', 'S3', 'S4', 'S5']
# inames=['+N','-N','+P','-P', 'C']
# x=range(10)
# outDir='/Users/harrietalexander/Dropbox/NB_141020/Figures/SingleGene_Blast/SiGenes/'
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
# 	ax.set_ylabel('TPM', fontsize=14)
# 	
# 	ax1=F.add_subplot(221)
# 	ax2=F.add_subplot(222)
# 	ax3=F.add_subplot(223)
# 	ax4=F.add_subplot(224)
# 	for item in Phash[ref]:
# 		ax1.set_title('T. rotula')
# 		ax1.set_xticks(x[:5])
# 		ax1.set_xticklabels(snames)
# 		ax2.set_title('S. costatum')
# 		ax2.set_xticks(x[:5])
# 		ax2.set_xticklabels(snames)
# 		ax3.set_xticklabels(inames)
# 		ax3.set_xticks(x[5:])
# 		ax4.set_xticks(x[5:])
# 		ax4.set_xticklabels(inames)
# 		
# 		if re.search('Throt', item):	
# 			
# 			if any(item==s for s in Throt_GeneList['N_Genes'][0]):
# 			 	ax1.plot(x[:5],Throt_tpm[item][:5], color='b')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='b', marker='o', mec='b', ms=8)
# 			elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
# 			 	ax1.plot(x[:5],Throt_tpm[item][:5], '--', color='b')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='b', marker='s', mec='b', ms=8)
# 			elif any(item==s for s in Throt_GeneList['P_Genes'][0]):
# 			 	ax1.plot(x[:5],Throt_tpm[item][:5], color='r')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='r', marker='o', mec='r', ms=8)
# 			elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
# 			 	ax1.plot(x[:5],Throt_tpm[item][:5], '--', color='r')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='r', marker='s', mec='r', ms=8)
# 			else:
# 				ax1.plot(x[:5],Throt_tpm[item][:5], color='0.75')
# 			 	ax3.plot(x[5:],Throt_tpm[item][5:],color='0.5', ls='None', marker='x', ms=8)
# 			 	
# 
# 		elif re.search('Skcos', item):
# 			if any(item==s for s in Skcos_GeneList['N_Genes'][0]):
# 			 	ax2.plot(x[:5],Skcos_tpm[item][:5], color='b')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='b', marker='o', mec='b', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
# 			 	ax2.plot(x[:5],Skcos_tpm[item][:5], '--', color='b')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='b', marker='s', mec='b', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['P_Genes'][0]):
# 			 	ax2.plot(x[:5],Skcos_tpm[item][:5], color='r')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='r', marker='o', mec='r', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
# 			 	ax2.plot(x[:5],Skcos_tpm[item][:5], '--', color='r')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='r', marker='s', mec='r', ms=8)
# 
# 			else:
# 				ax2.plot(x[:5],Skcos_tpm[item][:5], color='0.75')
# 			 	ax4.plot(x[5:],Skcos_tpm[item][5:],color='0.5', ls='None', marker='x', ms=8)	
# 	plt.suptitle(ref,fontsize=16)
# 	plt.savefig(outDir+ref+'_SingleGenes_TPM.pdf')
# 
# 
# plt.show()
# 
# 	#---- Si-related genes *************************************

Pfile='/Users/harrietalexander/Dropbox/NB_141020/Up_Skele_genes2.tab'
Preader=csv.reader(open(Pfile),delimiter='\t')
for key in Throt_GeneList:
	print key
Phash={}
for line in Preader:
	key = line[0]
	item=line[1]
	if key in Phash: 
		Phash[key].append(line[1])
	else:
		Phash[key]=[item]

snames=['S1', 'S2', 'S3', 'S4', 'S5']
inames=['+N','-N','+P','-P', 'C']
x=range(10)
outDir='/Users/harrietalexander/Dropbox/NB_141020/Figures/SingleGene_Blast/N-relatedSkcos/'
N=0
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
	ax.set_ylabel('TPM', fontsize=14)
	
	ax1=F.add_subplot(221)
	ax2=F.add_subplot(222)
	ax3=F.add_subplot(223)
	ax4=F.add_subplot(224)

	for item in Phash[ref]:
		ax1.set_title('T. rotula')
		ax1.set_xticks(x[:5])
		ax1.set_xticklabels(snames)
		ax2.set_title('S. costatum')
		ax2.set_xticks(x[:5])
		ax2.set_xticklabels(snames)
		ax3.set_xticklabels(inames)
		ax3.set_xticks(x[5:])
		ax4.set_xticks(x[5:])
		ax4.set_xticklabels(inames)
		ax1.set_ylim([0,1000])	
		ax2.set_ylim([0,1000])	
		ax3.set_ylim([0,1000])	
		ax4.set_ylim([0,1000])	
		if re.search('Throt', item):	
			
			if any(item==s for s in Throt_GeneList['N_Genes'][0]):
			 	ax1.plot(x[:5],Throt_tpm[item][:5], color='b')
			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='b', marker='o', mec='b', ms=8)
			elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
			 	ax1.plot(x[:5],Throt_tpm[item][:5], '--', color='b')
			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='b', marker='s', mec='b', ms=8)
			elif any(item==s for s in Throt_GeneList['P_Genes'][0]):
			 	ax1.plot(x[:5],Throt_tpm[item][:5], color='r')
			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='r', marker='o', mec='r', ms=8)
			elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
			 	ax1.plot(x[:5],Throt_tpm[item][:5], '--', color='r')
			 	ax3.plot(x[5:],Throt_tpm[item][5:],ls='None',color='r', marker='s', mec='r', ms=8)
			else:
				ax1.plot(x[:5],Throt_tpm[item][:5], color='0.75')
			 	ax3.plot(x[5:],Throt_tpm[item][5:],color='0.5', ls='None', marker='x', ms=8)
			 	

		elif re.search('Skcos', item):
			if any(item==s for s in Skcos_GeneList['N_Genes'][0]):
			 	ax2.plot(x[:5],Skcos_tpm[item][:5], color='b')
			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='b', marker='o', mec='b', ms=8)
			elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
			 	ax2.plot(x[:5],Skcos_tpm[item][:5], '--', color='b')
			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='b', marker='s', mec='b', ms=8)
			elif any(item==s for s in Skcos_GeneList['P_Genes'][0]):
			 	ax2.plot(x[:5],Skcos_tpm[item][:5], color='r')
			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='r', marker='o', mec='r', ms=8)
			elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
			 	ax2.plot(x[:5],Skcos_tpm[item][:5], '--', color='r')
			 	ax4.plot(x[5:],Skcos_tpm[item][5:],ls='None',color='r', marker='s', mec='r', ms=8)

			else:
				ax2.plot(x[:5],Skcos_tpm[item][:5], color='0.75')
			 	ax4.plot(x[5:],Skcos_tpm[item][5:],color='0.5', ls='None', marker='x', ms=8)	
	plt.suptitle(ref,fontsize=16)
	plt.savefig(outDir+ref+'_SingleGenes_TPM.pdf')


plt.show()

	