#!/usr/bin/env python

'''
Created on November 9, 2013

Comparison script for Throt and Skcos count data; calls upon general and std functions 

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
import General_STD_analysis_control as GSTD
import pickle
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

## PRINT STD TO FILE

#Relies upon the Genearl and STD_functions python files. 
#Location of pickled files
PickleOutDir="/Users/harrietalexander/Dropbox/NB_141020/PickledData/"

#Throt_filenames from pickle
Throt_filename=['Throt_PupP_STD', 'Throt_PupN_STD', 'Throt_PdnP_STD', 'Throt_PdnN_STD', 'Throt_NupP_STD', 'Throt_NupN_STD', 'Throt_NdnP_STD', 'Throt_NdnN_STD']
 

#Skcos Raw Data

#Skcos filenames from pickle
Skcos_filename=['Skcos_PupP_STD', 'Skcos_PupN_STD', 'Skcos_PdnP_STD', 'Skcos_PdnN_STD', 'Skcos_NupP_STD', 'Skcos_NupN_STD', 'Skcos_NdnP_STD', 'Skcos_NdnN_STD']

#Get STD_all-- from pickle. 
Throt_STD_all=STD.unPickleAll(Throt_filename, PickleOutDir)
Skcos_STD_all=STD.unPickleAll(Skcos_filename, PickleOutDir)

SpeciesName=['Throt_', 'Skcos_']
TreatmentName=['DnMinusP_', 'UpMinusP_', 'DnMinusN_', 'UpMinusN_']
STDName=['STDN', 'STDP']
headerRow=['geneID', 'S1', 'S2', 'S3', 'S4', 'S5', 'Control']


allSpecies=[Throt_STD_all, Skcos_STD_all]
for species, SsName in zip(allSpecies, SpeciesName):
	STD_N_Dict={}
	STD_P_Dict={}
	Ns=[species[1], species[3], species[5], species[7]]
	# =PupN_STD, PdnN_STD, NupN_STD, NdnN_STD
	Ps=[species[0], species[2], species[4], species[6]]
	# =PupP_STD, PdnP_STD, NupP_STD, NdnP_STD
	for x in range(4):
		STD_N_Dict.update(Ns[x])
		STD_P_Dict.update(Ps[x])
		print STD_P_Dict
	AllNP=[Ns,Ps]
	for N, Sname in zip(AllNP, STDName):
		for S, Tname in zip(N, TreatmentName):
			fileName=SsName+Tname+Sname+'.tab'
			STD.Write_Dict_To_File(S, fileName,headerRow)
	STD.Write_Dict_To_File(STD_N_Dict, SsName+'_All_STDN.tab',headerRow)
	STD.Write_Dict_To_File(STD_P_Dict, SsName+'_All_STDP.tab',headerRow)

# ## PRINT SGNC TPM RAW TO FILE
# PickleOutDir="/Users/harrietalexander/Dropbox/NB_141020/PickledData/"
# 
# #Throt_filenames from pickle
# Throt_Pk_Count=['Throt_raw', 'Throt_tpm', 'Throt_SGNC', 'Throt_GeneList']
# 
# #Skcos Raw Data
# 
# #Skcos filenames from pickle
# Skcos_Pk_Count=['Skcos_raw', 'Skcos_tpm', 'Skcos_SGNC', 'Skcos_GeneList']
# 
# #Get STD_all-- from pickle. 
# [Throt_raw, Throt_tpm, Throt_SGNC, Throt_GeneList]=STD.unPickleAll(Throt_Pk_Count, PickleOutDir)
# 
# [Skcos_raw, Skcos_tpm, Skcos_SGNC, Skcos_GeneList]=STD.unPickleAll(Skcos_Pk_Count, PickleOutDir)
# 
# headerRow=['geneID', 'S1', 'S2', 'S3', 'S4', 'S5', '+N', '-N', '+P', '-P', 'Control']
# 
# STD.Write_Dict_To_File(Skcos_raw, 'Skcos_rawCounts.tab',headerRow)
# STD.Write_Dict_To_File(Throt_raw, 'Throt_rawCounts.tab',headerRow)
# STD.Write_Dict_To_File(Skcos_tpm, 'Skcos_tpm.tab',headerRow)
# STD.Write_Dict_To_File(Throt_tpm, 'Throt_tpm.tab',headerRow)
# STD.Write_Dict_To_File(Skcos_SGNC, 'Skcos_SGNC.tab',headerRow)
# STD.Write_Dict_To_File(Throt_SGNC, 'Throt_SGNC.tab',headerRow)

