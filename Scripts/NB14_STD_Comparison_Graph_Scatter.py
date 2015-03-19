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

##  Scatter plots ********************************

STD_Combined=[Skcos_STD_all, Throt_STD_all]


#STD_Combined=[Throt_STD_all]
shape_combo=["o", "o"]
trans_combo=[0.5, 0.5]
colorN=["#66FFFF","#000099"]
colorP=["#FFDB4D","#FF6600"]

# GSTD.plot_by_day_Combined2(STD_Combined, 1, -10, 10, saveFig=1, SeqName=None, shape=shape_combo, trans=trans_combo, color=colorDown) 
# 
GSTD.plot_by_dayLine_combo(STD_Combined, 1,-15, 15, saveFig=1, SeqName=None, shape=shape_combo, trans=trans_combo, colorP=colorP, colorN=colorN)
plt.show()
# GSTD.plot_vs_Control(STD_Combined, 1, -10, 10, 0, 'P', SeqName=None, shape=shape_combo, trans=trans_combo, color=colorDown)



# #Histogram ************************
# 
# 
# Throt_Nums=GSTD.MakeList(Throt_STD_all)
# Throt_NList=Throt_Nums[0]
# #	Pup; Pdn; Nup; Ndn
# Throt_PList=Throt_Nums[1]
# 
# Skcos_Nums=GSTD.MakeList(Skcos_STD_all)
# Skcos_NList=Skcos_Nums[0]
# #	Pup; Pdn; Nup; Ndn
# Skcos_PList=Skcos_Nums[1]
# PlotOutputDir='/Users/harrietalexander/Dropbox/NB_141020/Figures/'
# GSTD.STD_Histogram(Throt_NList, 'STD-N', 'Throt-', PlotOutputDir)
# GSTD.STD_Histogram(Throt_PList, 'STD-P', 'Throt-', PlotOutputDir)
# 
# GSTD.STD_Histogram(Skcos_NList, 'STD-N', 'Skcos-', PlotOutputDir)
# GSTD.STD_Histogram(Skcos_PList, 'STD-P', 'Skcos-', PlotOutputDir)
# 
# GSTD.STD_Histogram_all(Throt_NList, 'STD-N', 'Throt-', PlotOutputDir)
# GSTD.STD_Histogram_all(Throt_PList, 'STD-P', 'Throt-', PlotOutputDir)
# 
# GSTD.STD_Histogram_all(Skcos_NList, 'STD-N', 'Skcos-', PlotOutputDir)
# GSTD.STD_Histogram_all(Skcos_PList, 'STD-P', 'Skcos-', PlotOutputDir)
# # 
# 


# #Quadrant Analysis ************************
# 
# #Plot quads look at how the variation shifts. 
# def divideList(list, y):
# 	newList=[]
# 	for sublist in list: 
# 		newSubList = [float(x)/float(y) for x in sublist]
# 		newList.append(newSubList)
# 	return newList
# 		
# 
# Throt_Quad=STD.count_Quadrants3(Throt_STD_all, 0.5, 0.5)
# Throt_geneList=Throt_Quad[0]
# Throt_quadnums=Throt_Quad[1]
# Throt_total=Throt_Quad[2]
# Throt_quadDiv=divideList(Throt_quadnums, Throt_total)
# # print len(set(Throt_geneList[3][0]+Throt_geneList[3][1]+Throt_geneList[3][2]+Throt_geneList[3][3]+Throt_geneList[3][4]))
# # print len(set(Throt_geneList[3][0]+Throt_geneList[3][2]+Throt_geneList[3][4]))
# # print len(set(Throt_geneList[3][0]))
# # print len(set(Throt_geneList[3][2]))
# # print len(set(Throt_geneList[3][4]))
# 
# Skcos_Quad=STD.count_Quadrants3(Skcos_STD_all, 0.5, 0.5)
# Skcos_geneList=Skcos_Quad[0]
# Skcos_quadnums=Skcos_Quad[1]
# Skcos_total=Skcos_Quad[2]
# Skcos_quadDiv=divideList(Skcos_quadnums, Skcos_total)
# # print len(set(Skcos_geneList[1][0]+Skcos_geneList[1][1]+Skcos_geneList[1][2]+Skcos_geneList[1][3]+Skcos_geneList[1][4]))
# # print len(set(Skcos_geneList[1][1]+Skcos_geneList[1][3]))
# # print len(set(Skcos_geneList[1][1]))
# # print len(set(Skcos_geneList[1][3]))
# # colors=['k', 'g', 'b', 'r']
# # for x in range(4):
# # 	c=colors[x]
# # 	plt.plot(Throt_quadDiv[x], color=c)
# # 	plt.plot(Skcos_quadDiv[x], linestyle='dashed', color=c)
# # plt.legend(["Throt:Co-limited", "Skcos:Co-limited", "Throt:P-limited", "Skcos:P-limited", "Throt:Replete", "Skcos:Replete","Throt:N-limited", "Skcos:N-limited"])
# # plt.show()
# 
# #Plot with different range of values between 0 and 1 to see how it changes the plots. 
# 
# #calculate over a range of cutoff values. Look at the variations in the output information 
# cutoffs=np.linspace(0.25, 0.75, 5)
# ThrotCutoffs=[]
# SkcosCutoffs=[]
# for c in cutoffs:
# 	Throt_Quad=STD.count_Quadrants3(Throt_STD_all, c, c)
# 	Throt_quadnums=Throt_Quad[1]
# 	Throt_total=Throt_Quad[2]
# 	Throt_quadDiv=divideList(Throt_quadnums, Throt_total)
# 	ThrotCutoffs.append(Throt_quadDiv)
# 	Skcos_Quad=STD.count_Quadrants3(Skcos_STD_all, c, c)
# 	Skcos_quadnums=Skcos_Quad[1]
# 	Skcos_total=Skcos_Quad[2]
# 	Skcos_quadDiv=divideList(Skcos_quadnums, Skcos_total)
# 	SkcosCutoffs.append(Skcos_quadDiv)
# 
# NU=plt.figure(1)
# fNU1=NU.add_subplot(151)
# fNU2=NU.add_subplot(152)
# fNU3=NU.add_subplot(153)
# fNU4=NU.add_subplot(154)
# fNU5=NU.add_subplot(155)
# NUfigs=[fNU1,fNU2,fNU3,fNU4,fNU5]
# 
# for i in range(6):
# 	fig=NUfigs[i]
# 	Throt_Plot=ThrotCutoffs[i]
# 	Skcos_Plot=SkcosCutoffs[i]
# 	colors=['k', 'g', 'b', 'r']
# 	stitle="Cutoff: "+str(cutoffs[i])
# 	for x in range(4):
# 		c=colors[x]
# 		fig.plot([1,2,3,4,5,6], Throt_Plot[x], color=c)
# 		fig.plot([1,2,3,4,5,6], Skcos_Plot[x], linestyle='dashed', color=c)
# 	fig.set_title(stitle)
# 	fig.set_xticks([1, 2, 3, 4, 5,6])
# 	fig.set_ylim([0, 0.8])
# 	if i==0:
# 		fig.set_ylabel('% Quadrent')
# 	if i==2: 
# 		fig.set_xlabel('Sample')
# 
# fig.legend(["Throt:Co-limited", "Skcos:Co-limited", "Throt:P-limited", "Skcos:P-limited", "Throt:Replete", "Skcos:Replete","Throt:N-limited", "Skcos:N-limited"])
# 
# Iterate through N and P different cutoffs and calculate the stdeves. 
# 
# cutoffs=np.linspace(0.25, 0.75, 5)
# print cutoffs
# ThrotCutoffs=[]
# SkcosCutoffs=[]
# for n in cutoffs:
# 	for p in cutoffs:
# 		Throt_Quad=STD.count_Quadrants3(Throt_STD_all, n, p)
# 		Throt_quadnums=Throt_Quad[1]
# 		Throt_total=Throt_Quad[2]
# 		Throt_quadDiv=divideList(Throt_quadnums, Throt_total)
# 		ThrotCutoffs.append(Throt_quadDiv)
# 		Skcos_Quad=STD.count_Quadrants3(Skcos_STD_all, n, p)
# 		Skcos_quadnums=Skcos_Quad[1]
# 		Skcos_total=Skcos_Quad[2]
# 		Skcos_quadDiv=divideList(Skcos_quadnums, Skcos_total)
# 		SkcosCutoffs.append(Skcos_quadDiv)
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
# 	AllData=[Means, Stds, [L1, L2, L3, L4]]
# 	return AllData
# 	
# 	
# Skcos_MeanSTD=getMeanSTD(SkcosCutoffs)
# Throt_MeanSTD=getMeanSTD(ThrotCutoffs)
# 
# fileName=['Throt_L1','Throt_L2', 'Throt_L3', 'Throt_L4']
# Sample=['S1','S2','S3','S4','S5','C']
# for x,f in zip(Throt_MeanSTD[2], fileName):
# 	test= file( f, "w" )
# 	test.write("\t".join(Sample))
# 	test.write( "\n" )
# 	for y in x: 
# 		test.write( "\t".join(str(z) for z in y) )
# 		test.write( "\n" )
# test.close()
# 
# 
# fileName=['Skcos_L1','Skcos_L2', 'Skcos_L3', 'Skcos_L4']
# for x,f in zip(Skcos_MeanSTD[2], fileName):
# 	test= file( f, "w" )
# 	test.write("\t".join(Sample))
# 	test.write( "\n" )
# 
# 	for y in x: 
# 		test.write( "\t".join(str(z) for z in y) )
# 		test.write( "\n" )
# test.close()
# 
# 
# #Getting data in the correct format for Tukey/ANOVA test
# fileName='Skcos_All'
# Quad=['L1','L2','L3','L4']
# Sample=['S1','S2','S3','S4','S5','C']
# test=file(fileName,'w')
# test.write('Quad')
# test.write('\t')
# test.write('Sample')
# test.write('\t')
# test.write('Percent')
# test.write('\n')
# for x,q in zip(Skcos_MeanSTD[2],Quad):
# 	for y in x: 
# 		for z,s in zip(y, Sample):
# 			ll=[q,s,z]
# 			test.write( "\t".join(str(zz) for zz in ll) )
# 			test.write( "\n" )
# test.close()
# 
# 	
# # 
# 
# #print len(SkcosCutoffs[1][0])
# fig2=plt.figure(2)
# fig=fig2.add_subplot(111)
# colors=['k', 'g', 'b', 'r']
# for x in range(4):
# 	c=colors[x]
# 	fig.errorbar([1,2,3,4,5,6], Throt_MeanSTD[0][x], yerr=Throt_MeanSTD[1][x], color=c, marker='o', markersize=10, markeredgecolor=c)
# 	fig.errorbar([1,2,3,4,5,6], Skcos_MeanSTD[0][x], yerr=Skcos_MeanSTD[1][x], linestyle='dashed', color=c, marker='s', markersize=10, markeredgecolor=c)
# fig.set_xticks([1, 2, 3, 4, 5, 6])
# fig.set_xlim([0.5,6.5])
# fig.legend(["Throt:Co-limited", "Skcos:Co-limited", "Throt:P-limited", "Skcos:P-limited", "Throt:Replete", "Skcos:Replete","Throt:N-limited", "Skcos:N-limited"])
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
# 	ThrotU=fig.bar(ind, Throt_MeanSTD[0][x], width, yerr=Throt_MeanSTD[1][x], color='k', ecolor='k')
# 	SkcosU=fig.bar(ind+width, Skcos_MeanSTD[0][x], width, yerr=Skcos_MeanSTD[1][x], color='w', ecolor='k')
# 	fig.set_ylabel('% Quadrant')
# 	fig.set_xlabel('Sample')
# 	fig.set_xticks(ind+width)
# 	fig.set_xticklabels( ('S1', 'S2', 'S3', 'S4', 'S5', 'Control') )
# 	fig.set_title(titles[x])	
# 	fig.set_ylim([0, 0.9])
# 	if x==1:
# 		fig.legend( (ThrotU[0], SkcosU[0]), ('T. rotula', 'Skeletonema') )
# 
# plt.show()
# 
# ## Graph all genes movement in space. 	
# 
# GSTD.graphAllGenes(Throt_STD_all)
# GSTD.graphAllGenes(Skcos_STD_all)
# 
# 
# 
# 
# plt.show()
















