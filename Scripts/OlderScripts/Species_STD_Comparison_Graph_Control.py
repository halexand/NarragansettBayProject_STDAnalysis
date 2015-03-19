#!/usr/bin/env python

'''
Created on November 9, 2013

Comparison script for GSO and MMETSP count data; calls upon general and std functions 

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

#Relies upon the Genearl and STD_functions python files. 
#Location of pickled files
PickleOutDir="/Users/harrietalexander/Dropbox/NB_141020/PickledData/"

#GSO_filenames from pickle
GSO_filename=['GSO_PupP_STD', 'GSO_PupN_STD', 'GSO_PdnP_STD', 'GSO_PdnN_STD', 'GSO_NupP_STD', 'GSO_NupN_STD', 'GSO_NdnP_STD', 'GSO_NdnN_STD']

#MMETSP Raw Data
MMETSP_CountFile="/Users/harrietalexander/Dropbox/NB_130920/SD-NB_MMETSP0013_count.tab"
MMETSP_PP_Dir="/Users/harrietalexander/Dropbox/NB_130920/MMETSP_ASC/PostProb/"
MMETSP_Tail_125="*_125.txt"

#MMETSP filenames from pickle
MMETSP_filename=['MMETSP_PupP_STD', 'MMETSP_PupN_STD', 'MMETSP_PdnP_STD', 'MMETSP_PdnN_STD', 'MMETSP_NupP_STD', 'MMETSP_NupN_STD', 'MMETSP_NdnP_STD', 'MMETSP_NdnN_STD']

#Get STD_all-- from pickle. 
GSO_STD_all=STD.unPickleAll(GSO_filename, PickleOutDir)

MMETSP_STD_all=STD.unPickleAll(MMETSP_filename, PickleOutDir)

#Scatter plots

STD_Combined=[MMETSP_STD_all, GSO_STD_all]
print STD_Combined
# #STD_Combined=[GSO_STD_all]
# shape_combo=["o", "D"]
# trans_combo=[0.5, 0.5]
# colorDown=["#FF9900","#66CCFF"]
# colorUp=["#CC6600","#3333FF"]
# SeqName=["T. rotula", "Skeletonema"]
# # #Colored plots with different colored dots for each day
# # # GSTD.STD_color_subplot(GSO_STD_all, 1, -10, 10, 1, "T. rotula", m="o", a=0.5)
# # # GSTD.STD_color_subplot(MMETSP_STD_all, 1, -10,10, 2, "Skeletonema", m="D", a=0.5)
# # plt.show()
# # 
# # #Plots with colored points for each different 
# # #GSTD.plot_by_day(STD_Combined, 1, -10, 10, SeqName=None, shape=shape_combo, trans=trans_combo, colorUp=colorUp, colorDown=colorDown) 
# # #GSTD.plot_by_day2(STD_Combined, 1, -10, 10, saveFig=0, SeqName=None, shape=shape_combo, trans=trans_combo, color=colorDown) 
# # 
# #GSTD.plot_by_day_Combined2(STD_Combined, 1, -10, 10, saveFig=0, SeqName=None, shape=shape_combo, trans=trans_combo, color=colorDown) 
# #Histogram
# GSO_Nums=GSTD.MakeList(GSO_STD_all)
# GSO_NList=GSO_Nums[0]
# #	Pup; Pdn; Nup; Ndn
# GSO_PList=GSO_Nums[1]
# 
# MMETSP_Nums=GSTD.MakeList(MMETSP_STD_all)
# MMETSP_NList=MMETSP_Nums[0]
# #	Pup; Pdn; Nup; Ndn
# MMETSP_PList=MMETSP_Nums[1]
# # PlotOutputDir='/Users/harrietalexander/Dropbox/NB_130920/STD_Histograms_Control/'
# # GSTD.STD_Histogram(GSO_NList, 'STD-N', 'GSO-', PlotOutputDir)
# # GSTD.STD_Histogram(GSO_PList, 'STD-P', 'GSO-', PlotOutputDir)
# # 
# # GSTD.STD_Histogram(MMETSP_NList, 'STD-N', 'MMETSP-', PlotOutputDir)
# # GSTD.STD_Histogram(MMETSP_PList, 'STD-P', 'MMETSP-', PlotOutputDir)
# # 
# # GSTD.STD_Histogram_all(GSO_NList, 'STD-N', 'GSO-', PlotOutputDir)
# # GSTD.STD_Histogram_all(GSO_PList, 'STD-P', 'GSO-', PlotOutputDir)
# # 
# # GSTD.STD_Histogram_all(MMETSP_NList, 'STD-N', 'MMETSP-', PlotOutputDir)
# # GSTD.STD_Histogram_all(MMETSP_PList, 'STD-P', 'MMETSP-', PlotOutputDir)
# # 
# #Plot quads look at how the variation shifts. 
# def divideList(list, y):
# 	newList=[]
# 	for sublist in list: 
# 		newSubList = [float(x)/float(y) for x in sublist]
# 		newList.append(newSubList)
# 	return newList
# 		
# 
# GSO_Quad=STD.count_Quadrants3(GSO_STD_all, 0.5, 0.5)
# GSO_geneList=GSO_Quad[0]
# GSO_quadnums=GSO_Quad[1]
# GSO_total=GSO_Quad[2]
# GSO_quadDiv=divideList(GSO_quadnums, GSO_total)
# # print len(set(GSO_geneList[3][0]+GSO_geneList[3][1]+GSO_geneList[3][2]+GSO_geneList[3][3]+GSO_geneList[3][4]))
# # print len(set(GSO_geneList[3][0]+GSO_geneList[3][2]+GSO_geneList[3][4]))
# # print len(set(GSO_geneList[3][0]))
# # print len(set(GSO_geneList[3][2]))
# # print len(set(GSO_geneList[3][4]))
# 
# MMETSP_Quad=STD.count_Quadrants3(MMETSP_STD_all, 0.5, 0.5)
# MMETSP_geneList=MMETSP_Quad[0]
# MMETSP_quadnums=MMETSP_Quad[1]
# MMETSP_total=MMETSP_Quad[2]
# MMETSP_quadDiv=divideList(MMETSP_quadnums, MMETSP_total)
# # print len(set(MMETSP_geneList[1][0]+MMETSP_geneList[1][1]+MMETSP_geneList[1][2]+MMETSP_geneList[1][3]+MMETSP_geneList[1][4]))
# # print len(set(MMETSP_geneList[1][1]+MMETSP_geneList[1][3]))
# # print len(set(MMETSP_geneList[1][1]))
# # print len(set(MMETSP_geneList[1][3]))
# # colors=['k', 'g', 'b', 'r']
# # for x in range(4):
# # 	c=colors[x]
# # 	plt.plot(GSO_quadDiv[x], color=c)
# # 	plt.plot(MMETSP_quadDiv[x], linestyle='dashed', color=c)
# # plt.legend(["GSO:Co-limited", "MMETSP:Co-limited", "GSO:P-limited", "MMETSP:P-limited", "GSO:Replete", "MMETSP:Replete","GSO:N-limited", "MMETSP:N-limited"])
# # plt.show()
# 
# #Plot with different range of values between 0 and 1 to see how it changes the plots. 
# # 
# #calculate over a range of cutoff values. Look at the variations in the output information 
# cutoffs=np.linspace(0.4, 0.6, 5)
# GSOCutoffs=[]
# MMETSPCutoffs=[]
# for c in cutoffs:
# 	GSO_Quad=STD.count_Quadrants3(GSO_STD_all, c, c)
# 	GSO_quadnums=GSO_Quad[1]
# 	GSO_total=GSO_Quad[2]
# 	GSO_quadDiv=divideList(GSO_quadnums, GSO_total)
# 	GSOCutoffs.append(GSO_quadDiv)
# 	MMETSP_Quad=STD.count_Quadrants3(MMETSP_STD_all, c, c)
# 	MMETSP_quadnums=MMETSP_Quad[1]
# 	MMETSP_total=MMETSP_Quad[2]
# 	MMETSP_quadDiv=divideList(MMETSP_quadnums, MMETSP_total)
# 	MMETSPCutoffs.append(MMETSP_quadDiv)
# # 
# # NU=plt.figure(1)
# # fNU1=NU.add_subplot(151)
# # fNU2=NU.add_subplot(152)
# # fNU3=NU.add_subplot(153)
# # fNU4=NU.add_subplot(154)
# # fNU5=NU.add_subplot(155)
# # NUfigs=[fNU1,fNU2,fNU3,fNU4,fNU5]
# # 
# # for i in range(6):
# # 	fig=NUfigs[i]
# # 	GSO_Plot=GSOCutoffs[i]
# # 	MMETSP_Plot=MMETSPCutoffs[i]
# # 	colors=['k', 'g', 'b', 'r']
# # 	stitle="Cutoff: "+str(cutoffs[i])
# # 	for x in range(4):
# # 		c=colors[x]
# # 		fig.plot([1,2,3,4,5,6], GSO_Plot[x], color=c)
# # 		fig.plot([1,2,3,4,5,6], MMETSP_Plot[x], linestyle='dashed', color=c)
# # 	fig.set_title(stitle)
# # 	fig.set_xticks([1, 2, 3, 4, 5,6])
# # 	fig.set_ylim([0, 0.8])
# # 	if i==0:
# # 		fig.set_ylabel('% Quadrent')
# # 	if i==2: 
# # 		fig.set_xlabel('Sample')
# 
# # fig.legend(["GSO:Co-limited", "MMETSP:Co-limited", "GSO:P-limited", "MMETSP:P-limited", "GSO:Replete", "MMETSP:Replete","GSO:N-limited", "MMETSP:N-limited"])
# # 
# # # Iterate through N and P different cutoffs and calculate the stdeves. 
# # 
# cutoffs=np.linspace(0.25, 0.75, 5)
# GSOCutoffs=[]
# MMETSPCutoffs=[]
# for n in cutoffs:
# 	for p in cutoffs:
# 		GSO_Quad=STD.count_Quadrants3(GSO_STD_all, n, p)
# 		GSO_quadnums=GSO_Quad[1]
# 		GSO_total=GSO_Quad[2]
# 		GSO_quadDiv=divideList(GSO_quadnums, GSO_total)
# 		GSOCutoffs.append(GSO_quadDiv)
# 		MMETSP_Quad=STD.count_Quadrants3(MMETSP_STD_all, n, p)
# 		MMETSP_quadnums=MMETSP_Quad[1]
# 		MMETSP_total=MMETSP_Quad[2]
# 		MMETSP_quadDiv=divideList(MMETSP_quadnums, MMETSP_total)
# 		MMETSPCutoffs.append(MMETSP_quadDiv)
# 
# def getMeanSTD(CutoffList):
# 	L1=[]
# 	L2=[]
# 	L3=[]
# 	L4=[]
# 	for list in CutoffList:
# 		L1.append(list[0])
# 		L2.append(list[1])
# 		L3.append(list[2])
# 		L4.append(list[3])
#  	La1=np.array(L1)
#  	La1Mean=np.mean(La1, axis=0)
#  	La1Std=np.std(La1, axis=0)
#  	La2=np.array(L2)
#  	La2Mean=np.mean(La2, axis=0)
#  	La2Std=np.std(La2, axis=0)
#  	La3=np.array(L3)
#  	La3Mean=np.mean(La3, axis=0)
#  	La3Std=np.std(La3, axis=0)
#  	La4=np.array(L4)
#  	La4Mean=np.mean(La4, axis=0)
#  	La4Std=np.std(La4, axis=0)
# 	Means=[La1Mean, La2Mean, La3Mean, La4Mean]
# 	Stds=[La1Std, La2Std, La3Std, La4Std]
# 	AllData=[Means, Stds]
# 	return AllData
# # 	
# # 	
# MMETSP_MeanSTD=getMeanSTD(MMETSPCutoffs)
# GSO_MeanSTD=getMeanSTD(GSOCutoffs)
# # 
# # 
# print len(MMETSPCutoffs[1][0])
# fig2=plt.figure(2)
# fig=fig2.add_subplot(111)
# colors=['k', 'g', 'b', 'r']
# for x in range(4):
# 	c=colors[x]
# 	fig.errorbar([1,2,3,4,5,6], GSO_MeanSTD[0][x], yerr=GSO_MeanSTD[1][x], color=c, marker='o', markersize=10, markeredgecolor=c)
# 	fig.errorbar([1,2,3,4,5,6], MMETSP_MeanSTD[0][x], yerr=MMETSP_MeanSTD[1][x], linestyle='dashed', color=c, marker='s', markersize=10, markeredgecolor=c)
# fig.set_xticks([1, 2, 3, 4, 5, 6])
# fig.set_xlim([0.5,6.5])
# fig.legend(["GSO:Co-limited", "MMETSP:Co-limited", "GSO:P-limited", "MMETSP:P-limited", "GSO:Replete", "MMETSP:Replete","GSO:N-limited", "MMETSP:N-limited"])
# fig.set_ylabel('% Quadrent')
# fig.set_xlabel('Sample')
# fig.set_title('Mean and Standard Deviation with cutoff between 0.25 and 0.75')
# # 
# # 
# N=6
# ind = np.arange(N)
# width=0.25
# U=plt.figure(4)
# fU1=U.add_subplot(221)
# fU2=U.add_subplot(222)
# fU3=U.add_subplot(223)
# fU4=U.add_subplot(224)
# Us=[fU1, fU2, fU3, fU4]
# titles=["Co-limited", "P-limited", "Replete", "N-limited"]
# for x in range(4):
# 	fig=Us[x]
# 	GSOU=fig.bar(ind, GSO_MeanSTD[0][x], width, yerr=GSO_MeanSTD[1][x], color='k', ecolor='k')
# 	MMETSPU=fig.bar(ind+width, MMETSP_MeanSTD[0][x], width, yerr=MMETSP_MeanSTD[1][x], color='w', ecolor='k')
# 	fig.set_ylabel('% Quadrant')
# 	fig.set_xlabel('Sample')
# 	fig.set_xticks(ind+width)
# 	fig.set_xticklabels( ('S1', 'S2', 'S3', 'S4', 'S5', 'Control') )
# 	fig.set_title(titles[x])	
# 	fig.set_ylim([0, 0.9])
# 	if x==1:
# 		fig.legend( (GSOU[0], MMETSPU[0]), ('T. rotula', 'Skeletonema') )
# 
# plt.show()
# 
# ## Graph all genes movement in space. 	
# 
# GSTD.graphAllGenes(GSO_STD_all)
# GSTD.graphAllGenes(MMETSP_STD_all)
# 
# 
# 
# 
# plt.show()
















