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
import matplotlib.pyplot as pyplt
import matplotlib.cm as cm
import STD_Functions as STD
import General_STD_analysis_control as GSTD
import pickle
import pandas as pd
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


def color_dict(gradient):
    ''' Takes in a list of RGB sub-lists and returns dictionary of 
        colors in RGB and hex form for use in a graphing function 
        defined later on '''
    return {"hex":[RGB_to_hex(RGB) for RGB in gradient],
            "red":[RGB[0] for RGB in gradient],
            "green":[RGB[1] for RGB in gradient],
            "blue":[RGB[2] for RGB in gradient]}
            
def hex_to_RGB(hex):
    ''' "#FFFFFF" -> [255,255,255] '''
    # Pass 16 to the integer function for change of base
    return [int(hex[i:i+2], 16) for i in range(1,6,2)] 
    
def RGB_to_hex(RGB):
    ''' [255,255,255] -> "#FFFFFF" '''
    # Components need to be integers for hex to make sense
    RGB = [int(x) for x in RGB]
    return "#"+"".join(["0{0:x}".format(v) if v < 16 else 
                        "{0:x}".format(v) for v in RGB])
                        
def linear_gradient(start_hex, finish_hex="#FFFFFF", n=100):
    ''' returns a gradient list of (n) colors between two hex colors.
        start_hex and finish_hex should be the full six-digit color string, 
        inlcuding the number sign ("#FFFFFF") '''
    # Starting and ending colors in RGB form
    s = hex_to_RGB(start_hex)
    f = hex_to_RGB(finish_hex)
    # Initilize a list of the output colors with the starting color
    RGB_list = [s]
    # Calcuate a color at each evenly spaced value of t from 1 to n
    for t in range(1, n):
        # Interpolate RGB vector for color at the current value of t
        curr_vector = [ int(s[j] + (float(t)/(n-1))*(f[j]-s[j])) 
                        for j in range(3)]
        # Add it to our list of output colors
        RGB_list.append(curr_vector)

    return color_dict(RGB_list)


# wd='/Users/harrietalexander/Dropbox/NB_Paper/Supplemental_DataTables/STD/'
# os.chdir(wd)
# Skcos_STDN={}
# for file in glob.iglob('Skcos*STDN.tab'):
# 	name=file.replace('_STDN.tab','')
# 	Ar= pd.read_csv(file, sep='\t', index_col=0)
# 	Skcos_STDN[name]=Ar
# Skcos_STDP={}
# for file in glob.iglob('Skcos*STDP.tab'):
# 	name=file.replace('_STDP.tab','')
# 	Ar= pd.read_csv(file, sep='\t', index_col=0)
# 	Skcos_STDP[name]=Ar
# 
# Throt_STDN={}
# for file in glob.iglob('Throt*STDN.tab'):
# 	name=file.replace('_STDN.tab','')
# 	Ar= pd.read_csv(file, sep='\t', index_col=0)
# 	Throt_STDN[name]=Ar
# 
# Throt_STDP={}
# for file in glob.iglob('Throt*STDP.tab'):
# 	name=file.replace('_STDP.tab','')
# 	Ar= pd.read_csv(file, sep='\t', index_col=0)
# 	Throt_STDP[name]=Ar

wd='/Users/harrietalexander/Dropbox/NB_Paper/Supplemental_DataTables/STD/All/'
os.chdir
Throt_All_STDN= pd.read_csv(wd+'Throt__All_STDN.tab', sep='\t', index_col=0)
Throt_All_STDP= pd.read_csv(wd+'Throt__All_STDP.tab', sep='\t', index_col=0)
Skcos_All_STDN= pd.read_csv(wd+'Skcos__All_STDN.tab', sep='\t', index_col=0)
Skcos_All_STDP= pd.read_csv(wd+'Skcos__All_STDP.tab', sep='\t', index_col=0)

colList=['S1', 'S2', 'S3', 'S4', 'S5']
NU=plt.figure(1)
fNU1=NU.add_subplot(151, aspect = 'equal')
fNU2=NU.add_subplot(152, aspect = 'equal')
fNU3=NU.add_subplot(153, aspect = 'equal')
fNU4=NU.add_subplot(154, aspect = 'equal')
fNU5=NU.add_subplot(155, aspect = 'equal')
NUfigs=[fNU1,fNU2,fNU3,fNU4,fNU5]
shape_combo=["s", "o"]
trans_combo=[0.5, 0.5]
lims=[[-20,20],[-20,20],[-20,20],[-20,20],[-20,20]]
lightBlue=hex_to_RGB("#66FFFF")
darkBlue=hex_to_RGB("#000099")
lightOrange=hex_to_RGB("#FFDB4D")
darkOrange=hex_to_RGB("#FF6600")

colorN=["#000099", "#66FFFF"]
colorP=["#FF6600","#FFDB4D"]

Ncmap=linear_gradient(start_hex=colorN[0], finish_hex=colorN[1],n=10)
print Ncmap
# NDict=color_dict([lightBlue, darkBlue])
# Ncmap=mpl.colors.LinearSegmentedColormap('my_colormap',NDict,256)
# for fig in NUfigs:
# 		fig.axhline(y=1, color='#FF6600')
# 		fig.axvline(x=1, color='#000099')
# 		fig.axhspan(1,15, facecolor="#FFDB4D", alpha=0.15)
# 		fig.axvspan(1,15, facecolor="#66FFFF", alpha=0.15)


STDN=[Throt_All_STDN, Skcos_All_STDN]
STDP=[Throt_All_STDP, Skcos_All_STDP]
for SN, SP,m, cN, cP in zip(STDN, STDP, shape_combo, colorN, colorP):
	for x,fig,l in zip(colList,NUfigs, lims):
		Nquad=[]
		Pquad=[]
		Ncount=0
		Pcount=0
		dN=SN[x]
		dP=SP[x]
		print dN.shape[0]
		for N,P in zip(dN.T.iteritems(),dP.T.iteritems()) :
			if (N[1]<0.5) & (P[1]>0.5):
				Pquad.append(N[0])
				Pcount+=1
			elif (N[1]>0.5) & (P[1]<0.5):
				Nquad.append(N[0])
				Ncount+=1
		NN=dN.T[Nquad].median()
 		PN=dP.T[Nquad].median()
 		Ppercent=1. * Pcount/dN.shape[0]
		NP=dN.T[Pquad].median()
		PP=dP.T[Pquad].median()
 		Npercent=1. * Ncount/dN.shape[0]
# 		fig.scatter(NN,PN, color=Npercent, marker=m, alpha=1,s=Npercent*3000)
# 		fig.scatter(NP,PP, color=Ppercent, marker=m, alpha=1, s=Ppercent*3000)
		NCmap = plt.cm.get_cmap('jet', 11) 
		Nbounds=[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]
		#Nnorm = colors.BoundaryNorm(Nbounds, Ncmap.N)
				
		sc1=fig.scatter(1,-5, c=Npercent, marker=m,s=200, cmap=NCmap,vmin=0,vmax=1)
		sc2=fig.scatter(-5,1, c=Ppercent, marker=m, s=200, cmap=NCmap,vmin=0,vmax=1)
		fig.set_xlim(l)
		fig.set_ylim(l)
		fig.axvline(x=0, color='black')
		fig.axhline(y=0, color='black')




plt.colorbar(sc1)
plt.colorbar(sc2)
plt.show()
# print Skcos_All_STDN['S2']
# heatmap, xedges, yedges = np.histogram2d(Skcos_All_STDP['S2'],Skcos_All_STDN['S2'],bins=20, range=[[-40,40],[-40,40]])
# extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
# pyplt.imshow(heatmap, extent=extent, cmap=plt.cm.jet)
# pyplt.show()


##PLOT AVERAGES
# 
# STDN=[Throt_STDN,Skcos_STDN]
# STDP=[Throt_STDN,Skcos_STDP]
# Names=['Throt', 'Skcos']
# NU=plt.figure(1)
# fNU1=NU.add_subplot(161, aspect = 'equal')
# fNU2=NU.add_subplot(162, aspect = 'equal')
# fNU3=NU.add_subplot(163, aspect = 'equal')
# fNU4=NU.add_subplot(164, aspect = 'equal')
# fNU5=NU.add_subplot(165, aspect = 'equal')
# fNU6=NU.add_subplot(166, aspect = 'equal')
# NUfigs=[fNU1,fNU2,fNU3,fNU4,fNU5,fNU6]
# shape_combo=["^", "o"]
# trans_combo=[0.5, 0.5]
# colorN=["#66FFFF","#000099"]
# colorP=["#FFDB4D","#FF6600"]
# 
# PU=plt.figure(2)
# fPU1=PU.add_subplot(161, aspect = 'equal')
# fPU2=PU.add_subplot(162, aspect = 'equal')
# fPU3=PU.add_subplot(163, aspect = 'equal')
# fPU4=PU.add_subplot(164, aspect = 'equal')
# fPU5=PU.add_subplot(165, aspect = 'equal')
# fPU6=PU.add_subplot(166, aspect = 'equal')
# PUfigs=[fPU1,fPU2,fPU3,fPU4,fPU5,fPU6]
# shape_combo=["^", "o"]
# trans_combo=[0.5, 0.5]
# colorN=["#FFDB4D","#FF6600","#66FFFF","#000099"]
# 
# 
# # for SN, SP,s in zip(STDN, STDP, shape_combo):
# # 	for key,c in zip(SN, colorN):
# # 		print key
# # 		MeanP=SP[key].mean()
# # 		MeanN=SN[key].mean()
# # 		MedP=SP[key].median()
# # 		MedN=SN[key].median()
# # 		print MedP.max()
# # 		print MedP.min()
# # 		print MedN.max()
# # 		print MedN.min()
# # 		for i, fig in zip(MeanP.index, NUfigs):
# # 			fig.scatter(MeanP[i], MeanN[i], color=c, marker=s)
# # 		for i, fig in zip(MedP.index, PUfigs):
# # 			fig.scatter(MeanP[i], MeanN[i], color=c, marker=s)
# 
# for fig in NUfigs:
# 	fig.axvline(x=0, color='black')
# 	fig.axhline(y=0, color='black')
# 	fig.set_ylim([-40,40])
# 	fig.set_xlim([-40,40])
# 
# for fig in PUfigs:
# 	fig.axvline(x=0, color='black')
# 	fig.axhline(y=0, color='black')
# 	fig.set_ylim([-40,40])
# 	fig.set_xlim([-40,40])
	


#Relies upon the Genearl and STD_functions python files. 
#Location of pickled files
# PickleOutDir="/Users/harrietalexander/Dropbox/NB_141020/PickledData/"
# 
# #Throt_filenames from pickle
# Throt_filename=['Throt_PupP_STD', 'Throt_PupN_STD', 'Throt_PdnP_STD', 'Throt_PdnN_STD', 'Throt_NupP_STD', 'Throt_NupN_STD', 'Throt_NdnP_STD', 'Throt_NdnN_STD']
# 
# headerRow=['geneID', 'S1', 'S2', 'S3', 'S4', 'S5']
# #Skcos Raw Data
# 
# #Skcos filenames from pickle
# # Skcos_filename=['Skcos_PupP_STD', 'Skcos_PupN_STD', 'Skcos_PdnP_STD', 'Skcos_PdnN_STD', 'Skcos_NupP_STD', 'Skcos_NupN_STD', 'Skcos_NdnP_STD', 'Skcos_NdnN_STD']
# 
# #Get STD_all-- from pickle. 
# Throt_STD_all=STD.unPickleAll(Throt_filename, PickleOutDir)
# allSpecies=[Throt_STD_all]
# # Skcos_STD_all=STD.unPickleAll(Skcos_filename, PickleOutDir)
# 
# for species in allSpecies:
# 	Ns=[species[1], species[3], species[5], species[7]]
# 	# =PupN_STD, PdnN_STD, NupN_STD, NdnN_STD
# 	Ps=[species[0], species[2], species[4], species[6]]
# 	# =PupP_STD, PdnP_STD, NupP_STD, NdnP_STD
# 	for N, P in zip(Ns, Ps):
# 		Nlist=[]
# 		Ngenes=[]
# 		STD.Write_Dict_To_File(N,headerRow,'Test')
	
# 		for key in N:
# 			Nlist.append(N[key])
# 			Ngenes.append(key)
# 		Plist=[]
# 		Pgenes=[]
# 		for key in P:
# 			Plist.append(P[key])
# 			Pgenes.append(key)
# 		