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

# #Skcos filenames from pickle
# Skcos_Pk_Count=['Skcos_raw', 'Skcos_tpm', 'Skcos_SGNC', 'Skcos_GeneList']
# Skcos_Pk_STD=['Skcos_PupP_STD', 'Skcos_PupN_STD', 'Skcos_PdnP_STD', 'Skcos_PdnN_STD', 'Skcos_NupP_STD', 'Skcos_NupN_STD', 'Skcos_NdnP_STD', 'Skcos_NdnN_STD', 'Skcos_All_P', 'Skcos_All_N']

#Get STD_all-- from pickle. 
[Throt_raw, Throt_tpm, Throt_SGNC, Throt_GeneList]=STD.unPickleAll(Throt_Pk_Count, PickleOutDir)
[Throt_PupP_STD, Throt_PupN_STD, Throt_PdnP_STD, Throt_PdnN_STD, Throt_NupP_STD, Throt_NupN_STD, Throt_NdnP_STD, Throt_NdnN_STD, Throt_All_P, Throt_All_N]=STD.unPickleAll(Throt_Pk_STD, PickleOutDir)

# [Skcos_raw, Skcos_tpm, Skcos_SGNC, Skcos_GeneList]=STD.unPickleAll(Skcos_Pk_Count, PickleOutDir)

N=0

#---- P-related genes *************************************

Pfile='/Users/harrietalexander/Dropbox/NB_141020/Thaps_P_uniq.tab'
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


titles=['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Control']
#upregulated N figure
NU=plt.figure(1)
fNU1=NU.add_subplot(231, aspect = 'equal')
fNU2=NU.add_subplot(232, aspect = 'equal')
fNU3=NU.add_subplot(233, aspect = 'equal')
fNU4=NU.add_subplot(234, aspect = 'equal')
fNU5=NU.add_subplot(235, aspect = 'equal')
fNU6=NU.add_subplot(236, aspect = 'equal')

#downregulated N
ND=plt.figure(2)
fND1=ND.add_subplot(231, aspect = 'equal')
fND2=ND.add_subplot(232, aspect = 'equal')
fND3=ND.add_subplot(233, aspect = 'equal')
fND4=ND.add_subplot(234, aspect = 'equal')
fND5=ND.add_subplot(235, aspect = 'equal')
fND6=ND.add_subplot(236, aspect = 'equal')

NU.suptitle("N-genes")
ND.suptitle("P-genes")
#upregulatedP 	
PU=plt.figure(4)
fPU1=PU.add_subplot(231, aspect = 'equal')
fPU2=PU.add_subplot(232, aspect = 'equal')
fPU3=PU.add_subplot(233, aspect = 'equal')
fPU4=PU.add_subplot(234, aspect = 'equal')
fPU5=PU.add_subplot(235, aspect = 'equal')
fPU6=PU.add_subplot(236, aspect = 'equal')

PU.suptitle("Si-Genes")
allFigs=[NU, ND, PU]
NUfigs=[fNU1,fNU2,fNU3,fNU4,fNU5,fNU6]
NDfigs=[fND1,fND2,fND3,fND4,fND5,fND6]
PUfigs=[fPU1,fPU2,fPU3,fPU4,fPU5,fPU6]
allAxes=[NUfigs, NDfigs, PUfigs]

for t, tNU, tND, tPU, tPD in zip(titles, NUfigs, NDfigs, PUfigs):
	tNU.set_title(t)
	tPU.set_title(t)	
	tND.set_title(t)
	tPD.set_title(t)		

	if xylim:
		tNU.set_ylim((u,v))
		tNU.set_xlim((u,v))
		tPU.set_ylim((u,v))
		tPU.set_xlim((u,v))
		tND.set_ylim((u,v))
		tND.set_xlim((u,v))
		tPD.set_ylim((u,v))
		tPD.set_xlim((u,v))


for species, c, a, m in zip(STD_combo, color, trans, shape):

	Ns=[species[1], species[3], species[5], species[7]]
	# =PupN_STD, PdnN_STD, NupN_STD, NdnN_STD
	Ps=[species[0], species[2], species[4], species[6]]
	# =PupP_STD, PdnP_STD, NupP_STD, NdnP_STD
	loop=0
	z=0
	for lN, lP in zip(Ns, Ps):
		if loop==0:
			for key in lN:
				fPU1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
				fPU2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
				fPU3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
				fPU4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
				fPU5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
				fPU6.scatter(lN[key][5], lP[key][5], color=c, alpha=a, marker=m)
		elif loop==1:
			for key in lN:
				fNU1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
				fNU2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
				fNU3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
				fNU4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
				fNU5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
				fNU6.scatter(lN[key][5], lP[key][5], color=c, alpha=a, marker=m)

		else:
			for key in lN:
				fND1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
				fND2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
				fND3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
				fND4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
				fND5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
				fND6.scatter(lN[key][5], lP[key][5], color=c, alpha=a, marker=m)

		loop+=1
		
for fig in allAxes:
	for subplot in fig:
		subplot.axvline(x=0, color='black')
		subplot.axhline(y=0, color='black')
for fig in allFigs:
	fig.text(0.5, 0.05, 'STD-N', ha='center', va='center') #x axis
	fig.text(0.06, 0.5, 'STD-P', ha='center', va='center', rotation='vertical') #y axis	
if saveFig: 
	NU.savefig("NU.pdf")
	ND.savefig("ND.pdf")
	PU.savefig("PU.pdf")
	PD.savefig("PD.pdf")


 
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
# 	ax.set_ylabel('SGNC', fontsize=14)
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
# 			 	ax1.plot(x[:5],Throt_SGNC[item][:5], color='b')
# 			 	ax3.plot(x[5:],Throt_SGNC[item][5:],ls='None',color='b', marker='o', mec='b', ms=8)
# 			elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
# 			 	ax1.plot(x[:5],Throt_SGNC[item][:5], '--', color='b')
# 			 	ax3.plot(x[5:],Throt_SGNC[item][5:],ls='None',color='b', marker='s', mec='b', ms=8)
# 			elif any(item==s for s in Throt_GeneList['P_Genes'][0]):
# 			 	ax1.plot(x[:5],Throt_SGNC[item][:5], color='r')
# 			 	ax3.plot(x[5:],Throt_SGNC[item][5:],ls='None',color='r', marker='o', mec='r', ms=8)
# 			elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
# 			 	ax1.plot(x[:5],Throt_SGNC[item][:5], '--', color='r')
# 			 	ax3.plot(x[5:],Throt_SGNC[item][5:],ls='None',color='r', marker='s', mec='r', ms=8)
# 			else:
# 				ax1.plot(x[:5],Throt_SGNC[item][:5], color='0.75')
# 			 	ax3.plot(x[5:],Throt_SGNC[item][5:],color='0.5', ls='None', marker='x', ms=8)
# 			 	
# 
# 		elif re.search('Skcos', item):
# 			if any(item==s for s in Skcos_GeneList['N_Genes'][0]):
# 			 	ax2.plot(x[:5],Skcos_SGNC[item][:5], color='b')
# 			 	ax4.plot(x[5:],Skcos_SGNC[item][5:],ls='None',color='b', marker='o', mec='b', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
# 			 	ax2.plot(x[:5],Skcos_SGNC[item][:5], '--', color='b')
# 			 	ax4.plot(x[5:],Skcos_SGNC[item][5:],ls='None',color='b', marker='s', mec='b', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['P_Genes'][0]):
# 			 	ax2.plot(x[:5],Skcos_SGNC[item][:5], color='r')
# 			 	ax4.plot(x[5:],Skcos_SGNC[item][5:],ls='None',color='r', marker='o', mec='r', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
# 			 	ax2.plot(x[:5],Skcos_SGNC[item][:5], '--', color='r')
# 			 	ax4.plot(x[5:],Skcos_SGNC[item][5:],ls='None',color='r', marker='s', mec='r', ms=8)
# 
# 			else:
# 				ax2.plot(x[:5],Skcos_SGNC[item][:5], color='0.75')
# 			 	ax4.plot(x[5:],Skcos_SGNC[item][5:],color='0.5', ls='None', marker='x', ms=8)
# 	plt.suptitle(ref,fontsize=16)
# 	plt.savefig(outDir+ref+'_SingleGenes_SGNC.pdf')
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
# 	ax.set_ylabel('SGNC', fontsize=14)
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
# 			 	ax1.plot(x[:5],Throt_SGNC[item][:5], color='b')
# 			 	ax3.plot(x[5:],Throt_SGNC[item][5:],ls='None',color='b', marker='o', mec='b', ms=8)
# 			elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
# 			 	ax1.plot(x[:5],Throt_SGNC[item][:5], '--', color='b')
# 			 	ax3.plot(x[5:],Throt_SGNC[item][5:],ls='None',color='b', marker='s', mec='b', ms=8)
# 			elif any(item==s for s in Throt_GeneList['P_Genes'][0]):
# 			 	ax1.plot(x[:5],Throt_SGNC[item][:5], color='r')
# 			 	ax3.plot(x[5:],Throt_SGNC[item][5:],ls='None',color='r', marker='o', mec='r', ms=8)
# 			elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
# 			 	ax1.plot(x[:5],Throt_SGNC[item][:5], '--', color='r')
# 			 	ax3.plot(x[5:],Throt_SGNC[item][5:],ls='None',color='r', marker='s', mec='r', ms=8)
# 			else:
# 				ax1.plot(x[:5],Throt_SGNC[item][:5], color='0.75')
# 			 	ax3.plot(x[5:],Throt_SGNC[item][5:],color='0.5', ls='None', marker='x', ms=8)
# 			 	
# 
# 		elif re.search('Skcos', item):
# 			if any(item==s for s in Skcos_GeneList['N_Genes'][0]):
# 			 	ax2.plot(x[:5],Skcos_SGNC[item][:5], color='b')
# 			 	ax4.plot(x[5:],Skcos_SGNC[item][5:],ls='None',color='b', marker='o', mec='b', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
# 			 	ax2.plot(x[:5],Skcos_SGNC[item][:5], '--', color='b')
# 			 	ax4.plot(x[5:],Skcos_SGNC[item][5:],ls='None',color='b', marker='s', mec='b', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['P_Genes'][0]):
# 			 	ax2.plot(x[:5],Skcos_SGNC[item][:5], color='r')
# 			 	ax4.plot(x[5:],Skcos_SGNC[item][5:],ls='None',color='r', marker='o', mec='r', ms=8)
# 			elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
# 			 	ax2.plot(x[:5],Skcos_SGNC[item][:5], '--', color='r')
# 			 	ax4.plot(x[5:],Skcos_SGNC[item][5:],ls='None',color='r', marker='s', mec='r', ms=8)
# 
# 			else:
# 				ax2.plot(x[:5],Skcos_SGNC[item][:5], color='0.75')
# 			 	ax4.plot(x[5:],Skcos_SGNC[item][5:],color='0.5', ls='None', marker='x', ms=8)	
# 	plt.suptitle(ref,fontsize=16)
# 	plt.savefig(outDir+ref+'_SingleGenes_SGNC.pdf')
# 
# 
# plt.show()

	