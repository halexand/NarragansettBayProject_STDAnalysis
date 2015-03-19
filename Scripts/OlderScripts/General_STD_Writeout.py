#!/usr/bin/env python

'''
Created on November 9, 2013

General_STD_analysis
Functions to run general anlaysis on different species; call functions in differnet script

@author: harrietalexander
'''

import sys
import glob
import os
import numpy as np
import itertools
import re
import csv
import matplotlib.pylab as plt
import matplotlib.cm as cm
import matplotlib.backends.backend_pdf as PdfPages
import STD_Writeout_Functions as STD

def General_STD_analysis(CountFile, cc, PP_Dir, Tail_125, SeqName, RR):

## Input CountFile (tabdelimited raw counts), number of data collumns (counts not seqname), Post-prob directory, tail name for searching, and species name. 

	##Load and modify raw count data
	
	#Import all counts from tab delimited file (SeqID, Count, Count, Count, etc.)
	#CountFile="/Users/harrietalexander/Dropbox/NB_130920/SD-NB_GSO101_count.tab"
	Counts_raw=STD.importDict(CountFile)
	
	#Calculate tpm
	Counts_tpm=STD.TPMCalc(Counts_raw)
	
	
	#Get list of tpm counts for plotting etc. 
	TPM_list, Seq_keys=STD.makeLists(Counts_tpm, cc)
	
	
	#Cutoff data based on histogram data
	cutoff=2
	TPM_Pass_Hash, Seq_Pass_List, TPM_Fail_List, Seq_Fail_List=STD.cutDataDict(Counts_tpm, cutoff)
	
	print "Total number of genes:", len(Seq_keys), "Genes at ", cutoff, "tpm cut:", len(Seq_Pass_List)
	
	##Get Stable Genes
	
	#Import all ASC output data at 125 cutoff
	#PP_Dir="/Users/harrietalexander/Dropbox/NB_130920/GSO_ASC/PostProb/"
	#Tail_125="*_125.txt"
	
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
	PupP_STD, PupN_STD=STD.STD_Calc(Pgenes[0], SGNC_normalized, 5, 6, 7, 8, RR)
	#P down STD
	PdnP_STD, PdnN_STD=STD.STD_Calc(Pgenes[1], SGNC_normalized, 5, 6, 7, 8, RR)
	#N up STD
	NupP_STD, NupN_STD=STD.STD_Calc(Ngenes[0], SGNC_normalized, 5, 6, 7, 8, RR)
	#P down STD
	NdnP_STD, NdnN_STD=STD.STD_Calc(Ngenes[1], SGNC_normalized, 5, 6, 7, 8, RR)
	Stable_All=[StableGenes]
	print Stable_All
	STD_all=[PupP_STD, PupN_STD, PdnP_STD, PdnN_STD, NupP_STD, NupN_STD, NdnP_STD, NdnN_STD]
	## Plot STD scores 
	
	return STD_all

def STD_color_subplot(STD_all, xylim, u, v, fignum=None, SeqName=None, m=None, a=None):
	#List of lists with all STD outputs, xylim=yes or no limit on scale, lower limit, upper limit, figure number 
	#sequence name (species), m = shape for plotting, a = transparency for plotting
	#Plot with each in situ sample colored the same.
	
	colors=cm.rainbow(np.linspace(0,1,5))
	fig=plt.figure(fignum)
	x=0
	titles=['P-up', 'P-down', 'N-up', 'N-down']
# 	xylim=1
# 	u=-10
# 	v=10
	
	PupP_STD=STD_all[0]
	PupN_STD=STD_all[1]
	PdnP_STD=STD_all[2]
	PdnN_STD=STD_all[3]
	NupP_STD=STD_all[4]
	NupN_STD=STD_all[5]
	NdnP_STD=STD_all[6]
	NdnN_STD=STD_all[7]

	for key in PupP_STD:
		ax=fig.add_subplot(2,2,x+1)
		for P,N,c in zip(PupP_STD[key],PupN_STD[key], colors):
			ax.scatter(N, P, color=c, marker=m, alpha=a)
	plt.title(titles[x])
	if xylim:
		plt.ylim((u,v))
		plt.xlim((u,v))
	x+=1
	for key in PdnP_STD:
		ax=fig.add_subplot(2,2,x+1)
		for P,N,c in zip(PdnP_STD[key],PdnN_STD[key], colors):
			ax.scatter(N, P, color=c, marker=m, alpha=a)
	plt.title(titles[x])
	if xylim:
		plt.ylim((u,v))
		plt.xlim((u,v))
	x+=1
	for key in NupP_STD:
		ax=fig.add_subplot(2,2,x+1)
		for P,N,c in zip(NupP_STD[key],NupN_STD[key], colors):
			ax.scatter(N, P, color=c, marker=m, alpha=a)
	plt.title(titles[x])
	if xylim:
		plt.ylim((u,v))
		plt.xlim((u,v))
	x+=1
	for key in NdnP_STD:
		ax=fig.add_subplot(2,2,x+1)
		for P,N,c in zip(NdnP_STD[key],NdnN_STD[key], colors):
			ax.scatter(N, P, color=c, marker=m, alpha=a)		
	plt.title(titles[x])
	if xylim:
		plt.ylim((u,v))
		plt.xlim((u,v))
	fig.text(0.5, 0.95, SeqName, ha='center', va='center' )	
	fig.text(0.5, 0.05, 'STD-N', ha='center', va='center') #x axis
	fig.text(0.06, 0.5, 'STD-P', ha='center', va='center', rotation='vertical') #y axis
	return fig



def plot_by_day(STD_combo, xylim, u, v, SeqName=None, shape=None, trans=None, colorUp=None, colorDown=None): 
	#program will loop through multiple sets of STD scores for different organisms (combined with sequence name, shape marker, transparency, colors, etc.)
	#and will create a sample by sample/day by day plot of the NISP 
	#Breaks the figures into N and P 
	titles=['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5']
# 	fN1=plt.figure(1)
# 	fP1=plt.figure(2)
# 	fN2=plt.figure(3)
# 	fP2=plt.figure(4)
# 	fN3=plt.figure(5)
# 	fP3=plt.figure(6)
# 	fN4=plt.figure(7)
# 	fP4=plt.figure(8)
# 	fN5=plt.figure(9)
# 	fP5=plt.figure(10)
	N=plt.figure(2)
	fN1=N.add_subplot(151, aspect = 'equal')
	fN2=N.add_subplot(152, aspect = 'equal')
	fN3=N.add_subplot(153, aspect = 'equal')
 	fN4=N.add_subplot(154, aspect = 'equal')
 	fN5=N.add_subplot(155, aspect = 'equal')
 	N.suptitle("N-related Genes")
 	P=plt.figure(3)
	fP1=P.add_subplot(151, aspect = 'equal')
	fP2=P.add_subplot(152, aspect = 'equal')
	fP3=P.add_subplot(153, aspect = 'equal')
 	fP4=P.add_subplot(154, aspect = 'equal')
 	fP5=P.add_subplot(155, aspect = 'equal')
 	P.suptitle("P-related Genes")
 	Nfigs=[fN1,fN2,fN3,fN4,fN5]
 	Pfigs=[fP1,fP2,fP3,fP4,fP5]
 	for t, tN, tP in zip(titles, Nfigs, Pfigs):
		tN.set_title(t)
		tP.set_title(t)		
		if xylim:
			tN.set_ylim((u,v))
			tN.set_xlim((u,v))
			tP.set_ylim((u,v))
			tP.set_xlim((u,v))

	for species, cU, cD, a, m in zip(STD_combo, colorUp, colorDown, trans, shape):
# 		PupP_STD=STD_all[0]
# 		PupN_STD=STD_all[1]
# 		PdnP_STD=STD_all[2]
# 		PdnN_STD=STD_all[3]
# 		NupP_STD=STD_all[4]
# 		NupN_STD=STD_all[5]
# 		NdnP_STD=STD_all[6]
# 		NdnN_STD=STD_all[7]
		Ns=[species[1], species[3], species[5], species[7]]
		# =PupN_STD, PdnN_STD, NupN_STD, NdnN_STD
		Ps=[species[0], species[2], species[4], species[6]]
		# =PupP_STD, PdnP_STD, NupP_STD, NdnP_STD
		loop=0
		z=0
		for lN, lP in zip(Ns, Ps):
			if loop==0:
				for key in lN:
					fP1.scatter(lN[key][0], lP[key][0], color=cU, alpha=a, marker=m)
					fP2.scatter(lN[key][1], lP[key][1], color=cU, alpha=a, marker=m)
					fP3.scatter(lN[key][2], lP[key][2], color=cU, alpha=a, marker=m)
					fP4.scatter(lN[key][3], lP[key][3], color=cU, alpha=a, marker=m)
					fP5.scatter(lN[key][4], lP[key][4], color=cU, alpha=a, marker=m)
			elif loop==1:
				for key in lN:
					fP1.scatter(lN[key][0], lP[key][0], color=cD, alpha=a, marker=m)
					fP2.scatter(lN[key][1], lP[key][1], color=cD, alpha=a, marker=m)
					fP3.scatter(lN[key][2], lP[key][2], color=cD, alpha=a, marker=m)
					fP4.scatter(lN[key][3], lP[key][3], color=cD, alpha=a, marker=m)
					fP5.scatter(lN[key][4], lP[key][4], color=cD, alpha=a, marker=m)
			elif loop==2:
				for key in lN:
					fN1.scatter(lN[key][0], lP[key][0], color=cU, alpha=a, marker=m)
					fN2.scatter(lN[key][1], lP[key][1], color=cU, alpha=a, marker=m)
					fN3.scatter(lN[key][2], lP[key][2], color=cU, alpha=a, marker=m)
					fN4.scatter(lN[key][3], lP[key][3], color=cU, alpha=a, marker=m)
					fN5.scatter(lN[key][4], lP[key][4], color=cU, alpha=a, marker=m)
			else:
				for key in lN:
					fN1.scatter(lN[key][0], lP[key][0], color=cD, alpha=a, marker=m)
					fN2.scatter(lN[key][1], lP[key][1], color=cD, alpha=a, marker=m)
					fN3.scatter(lN[key][2], lP[key][2], color=cD, alpha=a, marker=m)
					fN4.scatter(lN[key][3], lP[key][3], color=cD, alpha=a, marker=m)
					fN5.scatter(lN[key][4], lP[key][4], color=cD, alpha=a, marker=m)

			loop+=1

def makeSubplot(r, c, figNum):
	U=plt.figure(figNum)
	R=range(r*c)
	fU=[]

	for x in R:
		z = int(str(r) + str(c) + str(x+1))
		t=U.add_subplot(z, aspect='equal')
		fU.append(t)
	allFig=[U,fU]
	return(allFig)
		
def plot_by_day2(STD_combo, xylim, u, v, saveFig, SeqName=None, shape=None, trans=None, color=None): 
	#program will loop through multiple sets of STD scores for different organisms (combined with sequence name, shape marker, transparency, colors, etc.)
	#and will create a sample by sample/day by day plot of the NISP 
	#Breaks the figures into N and P 
	titles=['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Control']
	#upregulated N figure
	NU=plt.figure(1)
	fNU1=NU.add_subplot(151, aspect = 'equal')
	fNU2=NU.add_subplot(152, aspect = 'equal')
	fNU3=NU.add_subplot(153, aspect = 'equal')
 	fNU4=NU.add_subplot(154, aspect = 'equal')
 	fNU5=NU.add_subplot(155, aspect = 'equal')
 	fNU6=NU
 	#downregulated N
 	ND=plt.figure(2)
	fND1=ND.add_subplot(151, aspect = 'equal')
	fND2=ND.add_subplot(152, aspect = 'equal')
	fND3=ND.add_subplot(153, aspect = 'equal')
 	fND4=ND.add_subplot(154, aspect = 'equal')
 	fND5=ND.add_subplot(155, aspect = 'equal')
 	NU.suptitle("UP-REGULATED N")
 	ND.suptitle("DOWN-REGULATED N")
	#upregulatedP 	
 	PU=plt.figure(4)
	fPU1=PU.add_subplot(151, aspect = 'equal')
	fPU2=PU.add_subplot(152, aspect = 'equal')
	fPU3=PU.add_subplot(153, aspect = 'equal')
 	fPU4=PU.add_subplot(154, aspect = 'equal')
 	fPU5=PU.add_subplot(155, aspect = 'equal')
 	#downregulatedP
	PD=plt.figure(3)
	fPD1=PD.add_subplot(151, aspect = 'equal')
	fPD2=PD.add_subplot(152, aspect = 'equal')
	fPD3=PD.add_subplot(153, aspect = 'equal')
 	fPD4=PD.add_subplot(154, aspect = 'equal')
 	fPD5=PD.add_subplot(155, aspect = 'equal')

 	PU.suptitle("UP-REGULATED P")
 	PD.suptitle("DOWN-REGULATED P")
 	allFigs=[NU, ND, PU, PD]
 	NUfigs=[fNU1,fNU2,fNU3,fNU4,fNU5]
 	NDfigs=[fND1,fND2,fND3,fND4,fND5]
	PUfigs=[fPU1,fPU2,fPU3,fPU4,fPU5]
	PDfigs=[fPD1,fPD2,fPD3,fPD4,fPD5]
	allAxes=[NUfigs, NDfigs, PUfigs, PDfigs]	
 	for t, tNU, tND, tPU, tPD in zip(titles, NUfigs, NDfigs, PUfigs, PDfigs):
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
			elif loop==1:
				for key in lN:
					fPD1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
					fPD2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
					fPD3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
					fPD4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
					fPD5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
			elif loop==2:
				for key in lN:
					fNU1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
					fNU2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
					fNU3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
					fNU4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
					fNU5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
			else:
				for key in lN:
					fND1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
					fND2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
					fND3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
					fND4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
					fND5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
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


def plot_by_day3(STD_combo, xylim, u, v, saveFig, SeqName=None, shape=None, trans=None, color=None, titles=None): 
	#program will loop through multiple sets of STD scores for different organisms (combined with sequence name, shape marker, transparency, colors, etc.)
	#and will create a sample by sample/day by day plot of the NISP 
	#Breaks the figures into N and P 
	titles=['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Control']
	#upregulated N figure
	allNU=makeSubplot(2, 3, 1)
	NU=allNU[0]
	NUfigs=allNU[1]
 	#downregulated N
	allND=makeSubplot(2, 3, 2)
	ND=allND[0]
	NDfigs=allND[1]
 	NU.suptitle("UP-REGULATED N")
 	ND.suptitle("DOWN-REGULATED N")
	#upregulatedP 	
	allPU=makeSubplot(2, 3, 3)
	PU=allPU[0]
	PUfigs=allPU[1]
 	#downregulatedP
	allPD=makeSubplot(2, 3, 4)
	PD=allPD[0]
	PDfigs=allPD[1]


 	PU.suptitle("UP-REGULATED P")
 	PD.suptitle("DOWN-REGULATED P")
 	allFigs=[NU, ND, PU, PD]
	allAxes=[NUfigs, NDfigs, PUfigs, PDfigs]	
 	for t, tNU, tND, tPU, tPD in zip(titles, NUfigs, NDfigs, PUfigs, PDfigs):
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
			elif loop==1:
				for key in lN:
					fPD1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
					fPD2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
					fPD3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
					fPD4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
					fPD5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
			elif loop==2:
				for key in lN:
					fNU1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
					fNU2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
					fNU3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
					fNU4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
					fNU5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
			else:
				for key in lN:
					fND1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
					fND2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
					fND3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
					fND4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
					fND5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
			loop+=1
			
	for fig in allAxes:
		for subplot in fig:
			subplot.axvline(x=0, color='black')
			subplot.axhline(y=0, color='black')
	for fig in allFigs:
		fig.text(0.5, 0.05, 'STD-N', ha='center', va='center') #x axis
		fig.text(0.06, 0.5, 'STD-P', ha='center', va='center', rotation='vertical') #y axis	
	if saveFig: 
		NU.savefig("NU_con.pdf")
		ND.savefig("ND_con.pdf")
		PU.savefig("PU_con.pdf")
		PD.savefig("PD_con.pdf")
		
def plot_by_day_Combined(STD_combo, xylim, u, v, saveFig, SeqName=None, shape=None, trans=None, color=None): 
	#program will loop through multiple sets of STD scores for different organisms (combined with sequence name, shape marker, transparency, colors, etc.)
	#and will create a sample by sample/day by day plot of the NISP 
	#Breaks the figures into N and P 
	titles=['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5']
	#upregulated N figure
	NU=plt.figure(1)
	fNU1=NU.add_subplot(151, aspect = 'equal')
	fNU2=NU.add_subplot(152, aspect = 'equal')
	fNU3=NU.add_subplot(153, aspect = 'equal')
 	fNU4=NU.add_subplot(154, aspect = 'equal')
 	fNU5=NU.add_subplot(155, aspect = 'equal')
 	#downregulated N
 	ND=plt.figure(2)
	fND1=ND.add_subplot(151, aspect = 'equal')
	fND2=ND.add_subplot(152, aspect = 'equal')
	fND3=ND.add_subplot(153, aspect = 'equal')
 	fND4=ND.add_subplot(154, aspect = 'equal')
 	fND5=ND.add_subplot(155, aspect = 'equal')
 	NU.suptitle("UP-REGULATED N")
 	ND.suptitle("DOWN-REGULATED N")
	#upregulatedP 	
 	PU=plt.figure(4)
	fPU1=PU.add_subplot(151, aspect = 'equal')
	fPU2=PU.add_subplot(152, aspect = 'equal')
	fPU3=PU.add_subplot(153, aspect = 'equal')
 	fPU4=PU.add_subplot(154, aspect = 'equal')
 	fPU5=PU.add_subplot(155, aspect = 'equal')
 	#downregulatedP
	PD=plt.figure(3)
	fPD1=PD.add_subplot(151, aspect = 'equal')
	fPD2=PD.add_subplot(152, aspect = 'equal')
	fPD3=PD.add_subplot(153, aspect = 'equal')
 	fPD4=PD.add_subplot(154, aspect = 'equal')
 	fPD5=PD.add_subplot(155, aspect = 'equal')

 	PU.suptitle("All N & P Regulated Genes")
 	PD.suptitle("DOWN-REGULATED P")
 	allFigs=[NU, ND, PU, PD]
 	NUfigs=[fNU1,fNU2,fNU3,fNU4,fNU5]
 	NDfigs=[fND1,fND2,fND3,fND4,fND5]
	PUfigs=[fPU1,fPU2,fPU3,fPU4,fPU5]
	PDfigs=[fPD1,fPD2,fPD3,fPD4,fPD5]
	allAxes=[NUfigs, NDfigs, PUfigs, PDfigs]	
 	for t, tNU, tND, tPU, tPD in zip(titles, NUfigs, NDfigs, PUfigs, PDfigs):
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
		Ps=[species[0], species[2], species[4], species[6]]
		lN=dict(Ns[0].items()+Ns[1].items()+Ns[2].items()+Ns[3].items())
		lP=dict(Ps[0].items()+Ps[1].items()+Ps[2].items()+Ps[3].items())
		for key in lN:
			fPU1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
			fPU2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
			fPU3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
			fPU4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
			fPU5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
# 		elif loop==1:
# 			for key in lN:
# 				fPD1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
# 				fPD2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
# 				fPD3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
# 				fPD4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
# 				fPD5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
# 		elif loop==2:
# 			for key in lN:
# 				fNU1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
# 				fNU2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
# 				fNU3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
# 				fNU4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
# 				fNU5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
# 		else:
# 			for key in lN:
# 				fND1.scatter(lN[key][0], lP[key][0], color=c, alpha=a, marker=m)
# 				fND2.scatter(lN[key][1], lP[key][1], color=c, alpha=a, marker=m)
# 				fND3.scatter(lN[key][2], lP[key][2], color=c, alpha=a, marker=m)
# 				fND4.scatter(lN[key][3], lP[key][3], color=c, alpha=a, marker=m)
# 				fND5.scatter(lN[key][4], lP[key][4], color=c, alpha=a, marker=m)
# 		loop+=1
			
	for fig in allAxes:
		for subplot in fig:
			subplot.axvline(x=0, color='black')
			subplot.axhline(y=0, color='black')
	for fig in allFigs:
		fig.text(0.5, 0.05, 'STD-N', ha='center', va='center') #x axis
		fig.text(0.06, 0.5, 'STD-P', ha='center', va='center', rotation='vertical') #y axis	

def plot_by_day_Combined2(STD_combo, xylim, u, v, saveFig, SeqName=None, shape=None, trans=None, color=None, titles=None): 
	#program will loop through multiple sets of STD scores for different organisms (combined with sequence name, shape marker, transparency, colors, etc.)
	#and will create a sample by sample/day by day plot of the NISP 
	#Breaks the figures into N and P 
	titles=['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Control']
	#upregulated N figure
	allNU=makeSubplot(2, 3, 1)
	NU=allNU[0]
	NUfigs=allNU[1]
 	#downregulated N
	allND=makeSubplot(2, 3, 2)
	ND=allND[0]
	NDfigs=allND[1]
 	NU.suptitle("UP-REGULATED N")
 	ND.suptitle("DOWN-REGULATED N")
	#upregulatedP 	
	allPU=makeSubplot(2, 3, 3)
	PU=allPU[0]
	PUfigs=allPU[1]
 	#downregulatedP
	allPD=makeSubplot(2, 3, 4)
	PD=allPD[0]
	PDfigs=allPD[1]


 	PU.suptitle("UP-REGULATED P")
 	PD.suptitle("DOWN-REGULATED P")
 	allFigs=[NU, ND, PU, PD]
	allAxes=[NUfigs, NDfigs, PUfigs, PDfigs]	
 	for t, tNU, tND, tPU, tPD in zip(titles, NUfigs, NDfigs, PUfigs, PDfigs):
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
		Ps=[species[0], species[2], species[4], species[6]]
		lN=dict(Ns[0].items()+Ns[1].items()+Ns[2].items()+Ns[3].items())
		lP=dict(Ps[0].items()+Ps[1].items()+Ps[2].items()+Ps[3].items())
		for key in lN:
			x=len(PUfigs)
			for f in range(x):
				PUfigs[f].scatter(lN[key][f], lP[key][f], color=c, alpha=a, marker=m)
			
	for fig in allAxes:
		for subplot in fig:
			subplot.axvline(x=0, color='black')
			subplot.axhline(y=0, color='black')
	for fig in allFigs:
		fig.text(0.5, 0.05, 'STD-N', ha='center', va='center') #x axis
		fig.text(0.06, 0.5, 'STD-P', ha='center', va='center', rotation='vertical') #y axis	



def MakeList(All_STD):	
	#GeneralFormat=['PupP_STD', 'PupN_STD', 'PdnP_STD', 'PdnN_STD', 'NupP_STD', 'NupN_STD', 'NdnP_STD', 'NdnN_STD']

	Ns=[All_STD[1], All_STD[3], All_STD[5], All_STD[7]]
	#Pup; Pdn; Nup; Ndn
	Ps=[All_STD[0], All_STD[2], All_STD[4], All_STD[6]]

	#Number of bins 
	Nnums=[[],[],[],[]]
	Pnums=[[],[],[],[]]
	count=0
	for N, P in zip(Ns, Ps): 
		st1=[]
		st2=[]
		st3=[]
		st4=[]
		st5=[]
		for key in N:
			item=N[key]
			st1.append(item[0])
			st2.append(item[1])
			st3.append(item[2])
			st4.append(item[3])
			st5.append(item[4])
		Ntmp=[st1, st2, st3, st4, st5]
		Nnums[count]=Ntmp
		st1=[]
		st2=[]
		st3=[]
		st4=[]
		st5=[]
		for key in P:
			item=P[key]
			st1.append(item[0])
			st2.append(item[1])
			st3.append(item[2])
			st4.append(item[3])
			st5.append(item[4])
		Ptmp=[st1, st2, st3, st4, st5]
		Pnums[count]=Ptmp
		count+=1
	Nums=[Nnums,Pnums]
	return Nums
	
def STD_Histogram(List, PlotName=None, SpeciesName=None,outputDir=None):
	#	Pup; Pdn; Nup; Ndn
	Pdif=[[],[],[],[],[]]
	Ndif=[[],[],[],[],[]]

	for x in range(5): 
		P1=List[0][x]
		P2=List[1][x]
		Ps=P1+P2
		N1=List[2][x]
		N2=List[3][x]
		Ns=N1+N2
		Pdif[x].append(Ps)
		Ndif[x].append(Ns)
	Ndif2=[]
	for x in Ndif:
		Ndif2.append(x[0])
	Pdif2=[]
	for x in Pdif:
		Pdif2.append(x[0])
	plt.figure()
	n, bins, patches = plt.hist(Pdif2, 100, normed=1, histtype='stepfilled', range=[-10,-1], label=['S1', 'S2', 'S3', 'S4','S5'])
	plt.title(SpeciesName+'P-regulated')
	plt.xlabel(PlotName)
	plt.legend()
	plt.savefig(outputDir+SpeciesName+PlotName+'_Preg.neg.hist.png')
	
	plt.figure()
	n, bins, patches = plt.hist(Ndif2, 100, normed=1, histtype='stepfilled', range=[-10,-1], label=['S1', 'S2', 'S3', 'S4','S5'])
	plt.title(SpeciesName+'N-regulated')
	plt.xlabel(PlotName)
	plt.legend()
	plt.savefig(outputDir+SpeciesName+PlotName+'_Nreg.neg.hist.png')

def STD_Histogram(List, PlotName=None, SpeciesName=None,outputDir=None):
	#	Pup; Pdn; Nup; Ndn
	Pdif=[[],[],[],[],[]]
	Ndif=[[],[],[],[],[]]

	for x in range(5): 
		P1=List[0][x]
		P2=List[1][x]
		Ps=P1+P2
		N1=List[2][x]
		N2=List[3][x]
		Ns=N1+N2
		Pdif[x].append(Ps)
		Ndif[x].append(Ns)
	Ndif2=[]
	for x in Ndif:
		Ndif2.append(x[0])
	Pdif2=[]
	for x in Pdif:
		Pdif2.append(x[0])
	plt.figure()
	n, bins, patches = plt.hist(Pdif2, 100, normed=1, histtype='stepfilled', range=[-10,-1], label=['S1', 'S2', 'S3', 'S4','S5'])
	plt.title(SpeciesName+'P-regulated')
	plt.xlabel(PlotName)
	plt.legend()
	plt.savefig(outputDir+SpeciesName+PlotName+'_Preg.neg.hist.png')
	
	plt.figure()
	n, bins, patches = plt.hist(Ndif2, 100, normed=1, histtype='stepfilled', range=[-10,-1], label=['S1', 'S2', 'S3', 'S4','S5'])
	plt.title(SpeciesName+'N-regulated')
	plt.xlabel(PlotName)
	plt.legend()
	plt.savefig(outputDir+SpeciesName+PlotName+'_Nreg.neg.hist.png')
	
def STD_Histogram_all(List, PlotName=None, SpeciesName=None,outputDir=None):
	#	Pup; Pdn; Nup; Ndn
	allList=[[],[],[],[],[]]

	for x in range(5): 
		P1=List[0][x]
		P2=List[1][x]
		N1=List[2][x]
		N2=List[3][x]
		all=N1+N2+P1+P2
		allList[x].append(all)
		
	AllDif2=[]
	for x in allList:
		AllDif2.append(x[0])
	plt.figure()
	n, bins, patches = plt.hist(AllDif2, 100, normed=1, histtype='stepfilled', range=[-10,-1], label=['S1', 'S2', 'S3', 'S4','S5'])
	plt.title(SpeciesName+'AllDifGenes')
	plt.xlabel(PlotName)
	plt.legend()
	plt.savefig(outputDir+SpeciesName+PlotName+'_All.neg.hist.png')
	
# def Ven_Comp_GeneLists(GeneList):
# 	#List of Quadrants; quadrant list of samples
# 	for Quad in GeneList: 
# 		lofSets=[]
# 		for S in Quad:
# 			s=set(S)
# 			lofSets.append(s)
# 		#Pairwise Comparison:
# 		L=len(lofSets)
# 		for x in range(L):
			
def graphAllGenes(All_STD):
	Ns=[All_STD[1], All_STD[3], All_STD[5], All_STD[7]]
	print Ns
	#Pup; Pdn; Nup; Ndn
	Ps=[All_STD[0], All_STD[2], All_STD[4], All_STD[6]]
	A=plt.figure()
	fA1=A.add_subplot(221, )#aspect = 'equal')
	fA2=A.add_subplot(222, )#aspect = 'equal')
	fA3=A.add_subplot(223, )#aspect = 'equal')
	fA4=A.add_subplot(224, )#aspect = 'equal')
	title=['P-up', 'P-dn', 'N-up', 'N-dn']
	print Ns[3]['Seq59467'], Ps[3]['Seq59467']
	figs=[fA1, fA2, fA3, fA4]
	for N, P, fig, t in zip(Ns, Ps, figs,title):
		colors=cm.rainbow(np.linspace(0,1,len(N)))
		for key,c in zip(N, colors): 
			n=N[key]
			p=P[key]
			fig.plot(n,p, color=c, marker='o', markersize=3)
			fig.set_xlim([-10,10])
			fig.set_ylim([-10,10])
		fig.set_title(t)
		fig.axvline(x=0, color='black')
		fig.axhline(y=0, color='black')
		
