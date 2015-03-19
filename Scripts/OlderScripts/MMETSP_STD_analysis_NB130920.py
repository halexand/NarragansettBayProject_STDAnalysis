#!/usr/bin/env python

'''
Created on November 9, 2013

GSO_STD_analysis

A generalized collection of 

@author: harrietalexander
'''

import sys
import glob
import os
import numpy as np
import re
import csv
import matplotlib.pylab as plt
import matplotlib.cm as cm
import STD_Functions as STD

##Load and modify raw count data

#Import all counts from tab delimited file (SeqID, Count, Count, Count, etc.)
CountFile="/Users/harrietalexander/Dropbox/NB_130920/SD-NB_MMETSP0013_count.tab"
Counts_raw=STD.importDict(CountFile)

#Calculate tpm
Counts_tpm=STD.TPMCalc(Counts_raw)

#Get list of tpm counts for plotting etc. 
TPM_list, Seq_keys=STD.makeLists(Counts_tpm, 10)


#Cutoff data based on histogram data
TPM_Pass_Hash, Seq_Pass_List, TPM_Fail_List, Seq_Fail_List=STD.cutDataDict(Counts_tpm, 2)

print "Total number of genes:", len(Seq_keys), "Genes at 2 tpm cut:", len(Seq_Pass_List)

# ##Create Histogram
# 
# #Each individual sample 
# 
# fig=plt.figure(1)
# titles=['S1','S2','S3', 'S4', 'S5', 'A', 'B', 'C', 'D', 'E']
# for x in xrange(10):
# 	ax=fig.add_subplot(2,5,x+1)
# 	ax.hist(TPM_list[x], 100, [0,5])
# # 	plt.ylim((0,18000))
# 	plt.title(titles[x])
# plt.show()
# 
# #Flat list (sum of all counts)
# TPM_list_flat=STD.flatList(TPM_list)
# fig=plt.figure(1)
# plt.hist(TPM_list_flat, 50, [0,4])
# plt.show()


##Get Stable Genes

#Import all ASC output data at 125 cutoff
PP_Dir="/Users/harrietalexander/Dropbox/NB_130920/MMETSP_ASC/PostProb/"
Tail_125="*_125.txt"


#Incubation and field Post-prob data in to dictionary
incuHashPP_125, fieldHashPP_125 = STD.importASC(PP_Dir,Tail_125)

#Get subhash of only genes that passed the pre-defined tpm cutoff
incuHashPP_125_cut=STD.subHash_fromKeyList(incuHashPP_125, Seq_Pass_List)
fieldHashPP_125_cut=STD.subHash_fromKeyList(fieldHashPP_125, Seq_Pass_List)

#Get the stable genes based on the cutoff value
StableGenes=STD.testPP(incuHashPP_125_cut, 0.1)

#Calculate SGNC
mean_StableGeneCount,stdev_StableGeneCount, SGNC_normalized=SGNC_Counts=STD.SGNCCalc(StableGenes, Counts_tpm)		

## Get Differential Regulated Genes

Pgenes, Ngenes=STD.difRegGenes(PP_Dir, Seq_Pass_List, 0.95)

## Calculate STD scores

#P up STD
PupP_STD, PupN_STD=STD.STD_Calc(Pgenes[0], SGNC_normalized, 5, 6, 7, 8, range(5))
#P down STD
PdnP_STD, PdnN_STD=STD.STD_Calc(Pgenes[1], SGNC_normalized, 5, 6, 7, 8, range(5))
#N up STD
NupP_STD, NupN_STD=STD.STD_Calc(Ngenes[0], SGNC_normalized, 5, 6, 7, 8, range(5))
#P down STD
NdnP_STD, NdnN_STD=STD.STD_Calc(Ngenes[1], SGNC_normalized, 5, 6, 7, 8, range(5))

## Plot STD scores 

#Plot with each in situ sample colored the same. 

colors=cm.rainbow(np.linspace(0,1,5))
fig=plt.figure(1)
x=0
titles=['P-up', 'P-down', 'N-up', 'N-down']
xylim=1 #Change to zero to plot on full scale (no scaling)
u=-10 #axis lower lim
v=10 #axis upper lim
for key in PupP_STD:
  	ax=fig.add_subplot(2,2,x+1)
	for P,N,c in zip(PupP_STD[key],PupN_STD[key], colors):
		ax.scatter(P, N, color=c)
plt.title(titles[x])
if xylim:
	plt.ylim((u,v))
	plt.xlim((u,v))
x+=1
for key in PdnP_STD:
  	ax=fig.add_subplot(2,2,x+1)
 	for P,N,c in zip(PdnP_STD[key],PdnN_STD[key], colors):
		ax.scatter(P, N, color=c)
plt.title(titles[x])
if xylim:
	plt.ylim((u,v))
	plt.xlim((u,v))
x+=1
for key in NupP_STD:
  	ax=fig.add_subplot(2,2,x+1)
 	for P,N,c in zip(NupP_STD[key],NupN_STD[key], colors):
		ax.scatter(P, N, color=c)
plt.title(titles[x])
if xylim:
	plt.ylim((u,v))
	plt.xlim((u,v))
x+=1
for key in NdnP_STD:
  	ax=fig.add_subplot(2,2,x+1)
 	for P,N,c in zip(NdnP_STD[key],NdnN_STD[key], colors):
		ax.scatter(P, N, color=c)		
plt.title(titles[x])
if xylim:
	plt.ylim((u,v))
	plt.xlim((u,v))
	
fig.text(0.5, 0.04, 'STD-P', ha='center', va='center')
fig.text(0.06, 0.5, 'STD-N', ha='center', va='center', rotation='vertical')



# ## Plot without color 
# fig=plt.figure(1)
# x=0
# for key in PupP_STD:
#  	ax=fig.add_subplot(2,2,x+1)
# 	Pdata=PupP_STD[key]
# 	Ndata=PupN_STD[key]
# 	ax.scatter(Pdata, Ndata)
# 	plt.title(titles[x])
# # 	plt.ylim((-10,10))
# # 	plt.xlim((-10,10))
# 
# x+=1
# for key in PdnP_STD:
#  	ax=fig.add_subplot(2,2,x+1)
# 	Pdata=PdnP_STD[key]
# 	Ndata=PdnN_STD[key]
# 	ax.scatter(Pdata, Ndata)
# 	plt.title(titles[x])
# # 	plt.ylim((-10,10))
# # 	plt.xlim((-10,10))
# # 
# x+=1
# for key in NupP_STD:
#  	ax=fig.add_subplot(2,2,x+1)
# 	Pdata=NupP_STD[key]
# 	Ndata=NupN_STD[key]
# 	ax.scatter(Pdata, Ndata)
# 	plt.title(titles[x])
# # 	plt.ylim((-10,10))
# # 	plt.xlim((-10,10))
# 
# x+=1
# for key in NdnP_STD:
#  	ax=fig.add_subplot(2,2,x+1)
# 	Pdata=NdnP_STD[key]
# 	Ndata=NdnN_STD[key]
# 	ax.scatter(Pdata, Ndata)
# 	plt.title(titles[x])
# # 	plt.ylim((-10,10))
# # 	plt.xlim((-10,10))
# fig.text(0.5, 0.04, 'STD-P', ha='center', va='center')
# fig.text(0.06, 0.5, 'STD-N', ha='center', va='center', rotation='vertical')

plt.show()
