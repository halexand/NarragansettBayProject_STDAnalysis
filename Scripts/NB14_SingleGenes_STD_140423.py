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
Skcos_Pk_Count=['Skcos_raw', 'Skcos_tpm', 'Skcos_SGNC', 'Skcos_GeneList']
Skcos_Pk_STD=['Skcos_PupP_STD', 'Skcos_PupN_STD', 'Skcos_PdnP_STD', 'Skcos_PdnN_STD', 'Skcos_NupP_STD', 'Skcos_NupN_STD', 'Skcos_NdnP_STD', 'Skcos_NdnN_STD', 'Skcos_All_P', 'Skcos_All_N']

#Get STD_all-- from pickle. 
[Throt_raw, Throt_tpm, Throt_SGNC, Throt_GeneList]=STD.unPickleAll(Throt_Pk_Count, PickleOutDir)
[Throt_PupP_STD, Throt_PupN_STD, Throt_PdnP_STD, Throt_PdnN_STD, Throt_NupP_STD, Throt_NupN_STD, Throt_NdnP_STD, Throt_NdnN_STD, Throt_All_P, Throt_All_N]=STD.unPickleAll(Throt_Pk_STD, PickleOutDir)
[Skcos_PupP_STD, Skcos_PupN_STD, Skcos_PdnP_STD, Skcos_PdnN_STD, Skcos_NupP_STD, Skcos_NupN_STD, Skcos_NdnP_STD, Skcos_NdnN_STD, Skcos_All_P, Skcos_All_N]=STD.unPickleAll(Skcos_Pk_STD, PickleOutDir)
[Skcos_raw, Skcos_tpm, Skcos_SGNC, Skcos_GeneList]=STD.unPickleAll(Skcos_Pk_Count, PickleOutDir)


# #---- P-related genes *************************************
# 
# ## Get P genes from BLAST
# Pfile='/Users/harrietalexander/Dropbox/NB_141020/Thaps_P_uniq.tab'
# Preader=csv.reader(open(Pfile),delimiter='\t')
# 
# Phash={}
# for line in Preader:
# 	key = line[0]
# 	item=line[1]
# 	if key in Phash: 
# 		Phash[key].append(line[1])
# 	else:
# 		Phash[key]=[item]
# P_blastGenes=[]
# for key in Phash:
# 	for item in Phash[key]:
# 		P_blastGenes.append(item)
# tmp=list(set(P_blastGenes))
# P_blastGenes=tmp
# 
# Throt_PB_N={}
# Throt_PB_P={}
# Skcos_PB_N={}
# Skcos_PB_P={}
# for item in P_blastGenes:
# 	if re.search('Throt', item):	
# 		Throt_PB_N[item]=Throt_All_N[item]
# 		Throt_PB_P[item]=Throt_All_P[item]
# 	else:
# 		Skcos_PB_N[item]=Skcos_All_N[item]
# 		Skcos_PB_P[item]=Skcos_All_P[item]
# 		
# 
# A=plt.figure(1)
# A1=A.add_subplot(121, aspect='equal')
# A1.set_title('T. rotula')
# 
# A2=A.add_subplot(122, aspect='equal')
# A2.set_title('S. costatum')
# 
# A1.set_xlim([-30,30])
# A1.set_ylim([-30,30])
# A2.set_xlim([-30,30])
# A2.set_ylim([-30,30])
# 
# ax=A.add_subplot(111)
# # Turn off axis lines and ticks of the big subplot
# ax.spines['top'].set_color('none')
# ax.spines['bottom'].set_color('none')
# ax.spines['left'].set_color('none')
# ax.spines['right'].set_color('none')
# ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
# ax.set_axis_bgcolor("none")
# 
# ax.set_xlabel('STD-N',fontsize=14)
# ax.set_ylabel('STD-P', fontsize=14)
# 
# for item in Throt_PB_N:
# 	Ns=Throt_PB_N[item]
# 	Ps=Throt_PB_P[item]
# 	if any(item==s for s in Throt_GeneList['N_Genes'][0]):
# 		A1.plot(Ns,Ps, color='b', marker='o', mec='b', ms=8)
# 	elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
# 		A1.plot(Ns,Ps, '--',color='b',  marker='o', mec='b', ms=8)
# 	elif any(item==s for s in Throt_GeneList['P_Genes'][0]):
# 		A1.plot(Ns,Ps, color='r', marker='o', mec='r', ms=8)	
# 	elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
# 		A1.plot(Ns,Ps, '--',color='r', marker='o', mec='r', ms=8)	
# 	else:
# 		A1.plot(Ns,Ps, ls='None',color='0.5', marker='x', ms=8)	
# for item in Skcos_PB_N:
# 	Ns=Skcos_PB_N[item]
# 	Ps=Skcos_PB_P[item]
# 	if any(item==s for s in Skcos_GeneList['N_Genes'][0]):
# 		A2.plot(Ns,Ps, color='b', marker='o', mec='b', ms=8)
# 	elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
# 		A2.plot(Ns,Ps, '--',color='b',  marker='o', mec='b', ms=8)
# 	elif any(item==s for s in Skcos_GeneList['P_Genes'][0]):
# 		A2.plot(Ns,Ps, color='r', marker='o', mec='r', ms=8)	
# 	elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
# 		A2.plot(Ns,Ps, '--',color='r', marker='o', mec='r', ms=8)	
# 	else:
# 		A2.plot(Ns,Ps, ls='None',color='0.5', marker='x', ms=8)	
# 
# for subplot in [A1,A2]:
# 	subplot.axvline(x=0, color='black')
# 	subplot.axhline(y=0, color='black')
# 
# A.suptitle('P-Genes',fontsize=16)
# 
# 
# 
# 
# 
# #---- N-related genes *************************************

Pfile='/Users/harrietalexander/Dropbox/NB_141020/FinalSingleGeneFigures/NovelN.tab'
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
P_blastGenes=[]
for key in Phash:
	for item in Phash[key]:
		P_blastGenes.append(item)
tmp=list(set(P_blastGenes))
P_blastGenes=tmp

Throt_PB_N={}
Throt_PB_P={}
Skcos_PB_N={}
Skcos_PB_P={}
for item in P_blastGenes:
	if re.search('Throt', item):	
		Throt_PB_N[item]=Throt_All_N[item]
		Throt_PB_P[item]=Throt_All_P[item]
	else:
		Skcos_PB_N[item]=Skcos_All_N[item]
		Skcos_PB_P[item]=Skcos_All_P[item]
		

A=plt.figure(2)
A1=A.add_subplot(111,aspect='equal')
A1.set_title('T. rotula')
# 
# A2=A.add_subplot(122,aspect='equal')
# A2.set_title('S. costatum')

A1.set_xlim([-30,30])
A1.set_ylim([-30,30])
# A2.set_xlim([-30,30])
# A2.set_ylim([-30,30])
# ax=A.add_subplot(111)
# # Turn off axis lines and ticks of the big subplot
# ax.spines['top'].set_color('none')
# ax.spines['bottom'].set_color('none')
# ax.spines['left'].set_color('none')
# ax.spines['right'].set_color('none')
# ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
# ax.set_axis_bgcolor("none")

A1.set_xlabel('STD-N',fontsize=14)
A1.set_ylabel('STD-P', fontsize=14)

NUp='#000099'
NDn='#00FFFF'
PUp='#FF0000'
PDn='#FF9900'

NRG='0.75'

for item in Throt_PB_N:
	Ns=Throt_PB_N[item]
	Ps=Throt_PB_P[item]
	if any(item==s for s in Throt_GeneList['P_Genes'][0]):
		A1.plot(Ns,Ps,ls='--', color=PUp, marker='o', mec=PUp, ms=12,linewidth=3)	
	elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
		A1.plot(Ns,Ps, ls='--', color=PDn, marker='o', mec=PDn, ms=12,linewidth=3)	
	elif any(item==s for s in Throt_GeneList['N_Genes'][0]):
		A1.plot(Ns,Ps, ls='--',color=NUp, marker='o', mec=NUp, ms=12,linewidth=3)
	elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
		A1.plot(Ns,Ps,ls='--', color=NDn, marker='o', mec=NDn, ms=12,linewidth=3)
	
	else:
		A1.plot(Ns,Ps, ls='--',color=NRG, marker='o',mec=NRG, ms=12,linewidth=3)	
for item in Skcos_PB_N:
	Ns=Skcos_PB_N[item]
	Ps=Skcos_PB_P[item]
	if any(item==s for s in Skcos_GeneList['P_Genes'][0]):
		A1.plot(Ns,Ps, ls='-', color=PUp, marker='s', mec=PUp, ms=12,linewidth=3)	
	elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
		A1.plot(Ns,Ps, ls='-', color=PDn, marker='s', mec=PDn, ms=12,linewidth=3)	
	elif any(item==s for s in Skcos_GeneList['N_Genes'][0]):
		A1.plot(Ns,Ps, ls='-',color=NUp, marker='s', mec=NUp, ms=12,linewidth=3)
	elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
		A1.plot(Ns,Ps, ls='-', color=NDn, marker='s', mec=NDn, ms=12,linewidth=3)
	else:
		A1.plot(Ns,Ps,  ls='-',color=NRG, marker='o',mec=NRG, ms=12,linewidth=3)	

for subplot in [A1]:
	subplot.axvline(x=0, color='black')
	subplot.axhline(y=0, color='black')

#A.suptitle('N-Genes',fontsize=16)

plt.show()
# 
# 
# 
# #---- NISP PLOT *************************************
# 
# Pfile='/Users/harrietalexander/Dropbox/NB_141020/FinalSingleGeneFigures/US2'
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
# P_blastGenes=[]
# for key in Phash:
# 	for item in Phash[key]:
# 		P_blastGenes.append(item)
# tmp=list(set(P_blastGenes))
# P_blastGenes=tmp
# 
# Throt_PB_N={}
# Throt_PB_P={}
# Skcos_PB_N={}
# Skcos_PB_P={}
# for item in P_blastGenes:
# 	if re.search('Throt', item):	
# 		Throt_PB_N[item]=Throt_All_N[item]
# 		Throt_PB_P[item]=Throt_All_P[item]
# 	else:
# 		Skcos_PB_N[item]=Skcos_All_N[item]
# 		Skcos_PB_P[item]=Skcos_All_P[item]
# 		
# 
# A=plt.figure(3, figsize=(10,10))
# A1=A.add_subplot(111)
# #A1.set_title('T. rotula')
# 
# #A2=A.add_subplot(122)# aspect='equal')
# #A2.set_title('S. costatum')
# 
# #A1.set_xlim([-50,50])
# #A1.set_ylim([-20,20])
# #A2.set_xlim([-50,10])
# #A2.set_ylim([-10,10])
# ax=A.add_subplot(111)
# Turn off axis lines and ticks of the big subplot
# ax.spines['top'].set_color('none')
# ax.spines['bottom'].set_color('none')
# ax.spines['left'].set_color('none')
# ax.spines['right'].set_color('none')
# ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
# ax.set_axis_bgcolor("none")
# 
# NUp='#000099'
# NDn='#00FFFF'
# PUp='#FF0000'
# PDn='#FF9900'
# 
# NRG='0.75'
# 
# A1.set_xlabel('STD-N',fontsize=14)
# A1.set_ylabel('STD-P', fontsize=14)
# names=['S1','S2','S3','S4','S5','C']
# for item in Throt_PB_N:
# 	Ns=Throt_PB_N[item]
# 	Ps=Throt_PB_P[item]
# 	if any(item==s for s in Throt_GeneList['P_Genes'][0]):
# 		A1.plot(Ns,Ps, ls='--', color=PUp, marker='o', mec=PUp, ms=12,linewidth=3)	
# 	elif any(item==s for s in Throt_GeneList['P_Genes'][1]):
# 		A1.plot(Ns,Ps, ls='--', color=PDn, marker='o', mec=PDn, ms=12,linewidth=3)	
# 		for l, ns,ps in zip(names,Ns,Ps):
# 			A1.text(ns,ps,l,fontsize=12)
# 	elif any(item==s for s in Throt_GeneList['N_Genes'][0]):
# 		A1.plot(Ns,Ps, ls='--',color=NUp, marker='o', mec=NUp, ms=12,linewidth=3)
# 	elif any(item==s for s in Throt_GeneList['N_Genes'][1]):
# 		A1.plot(Ns,Ps, ls='--', color=NDn, marker='o', mec=NDn, ms=12,linewidth=3)
# 	else:
# 		A1.plot(Ns,Ps, ls='--',color=NRG, marker='o',mec=NRG, ms=12,linewidth=3)	
# 
# for item in Skcos_PB_N:
# 	Ns=Skcos_PB_N[item]
# 	Ps=Skcos_PB_P[item]
# 	if any(item==s for s in Skcos_GeneList['P_Genes'][0]):
# 		A1.plot(Ns,Ps, ls='-', color=PUp, marker='s', mec=PUp, ms=12,linewidth=3)
# 	elif any(item==s for s in Skcos_GeneList['P_Genes'][1]):
# 		A1.plot(Ns,Ps, ls='-', color=PDn, marker='s', mec=PDn, ms=12,linewidth=3)	
# 		for l, ns,ps in zip(names,Ns,Ps):
# 			A1.text(ns,ps,l,fontsize=12)
# 	elif any(item==s for s in Skcos_GeneList['N_Genes'][0]):
# 		A1.plot(Ns,Ps, ls='-',color=NUp, marker='s', mec=NUp, ms=12,linewidth=3)
# 	elif any(item==s for s in Skcos_GeneList['N_Genes'][1]):
# 		A1.plot(Ns,Ps, ls='-', color=NDn, marker='s', mec=NDn, ms=12,linewidth=3)
# 	else:
# 		A1.plot(Ns,Ps, ls='-',color=NRG, marker='o',mec=NRG, ms=12,linewidth=3)	
# 
# for subplot in [A1]:  #,A2]:
# 	subplot.axvline(x=0, color='black')
# 	subplot.axhline(y=0, color='black')
# 	subplot.axvline(x=1, color='black',ls='--')
# 	subplot.axhline(y=1, color='black',ls='--')
# 
# #A.suptitle('Thaps_24435',fontsize=16)
# 
# 
# plt.show()
