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
import General_STD_analysis_131117 as GSTD
import pickle
#Relies upon the Genearl and STD_functions python files. 

def pickleIt(file_name,data, outdir):
	try:
		with open(outdir+file_name+".pickle", "wb") as output_file:
			pickle.dump(data, output_file,-1)
        	output_file.close()
	except Exception:
		print "Cannot open the file:",file_name
		
def unPickleAll(file_name, outdir):
	listHash=[]
	for name in file_name:
		dict={}
		dict=pickle.load(open(outdir+name+".pickle", "rb"))
		GSO_STD_all.append(dict)
	return listHash 

GSO_CountFile="/Users/harrietalexander/Dropbox/NB_130920/SD-NB_GSO101_count.tab"
GSO_PP_Dir="/Users/harrietalexander/Dropbox/NB_130920/GSO_ASC/PostProb/"
GSO_Tail_125="*_125.txt"
PickleOutDir="~/Dropbox/dyhrmanlab_mason/scripts/Harriet/Python_Scripts/STD_Plots/Pickle/"
# GSO_STD_all=GSTD.General_STD_analysis(GSO_CountFile, 10, GSO_PP_Dir, GSO_Tail_125, "T. rotula")

GSO_filename=['GSO_PupP_STD', 'GSO_PupN_STD', 'GSO_PdnP_STD', 'GSO_PdnN_STD', 'GSO_NupP_STD', 'GSO_NupN_STD', 'GSO_NdnP_STD', 'GSO_NdnN_STD']

for filename, data in zip(GSO_filename, GSO_STD_all):
	pickleIt(filename, data, PickleOutDir)

# 
# MMETSP_CountFile="/Users/harrietalexander/Dropbox/NB_130920/SD-NB_MMETSP0013_count.tab"
# MMETSP_PP_Dir="/Users/harrietalexander/Dropbox/NB_130920/MMETSP_ASC/PostProb/"
# MMETSP_Tail_125="*_125.txt"
# MMETSP_STD_all=GSTD.General_STD_analysis(MMETSP_CountFile, 10, MMETSP_PP_Dir, MMETSP_Tail_125, "Skeletonema")
# 
## STD_all=[PupP_STD, PupN_STD, PdnP_STD, PdnN_STD, NupP_STD, NupN_STD, NdnP_STD, NdnN_STD]

#STD_Combined=[MMETSP_STD_all, GSO_STD_all]
STD_Combined=[GSO_STD_all]
shape_combo=["o", "D"]
trans_combo=[0.5, 0.5]
colorDown=["#FF9900","#66CCFF"]
colorUp=["#CC6600","#3333FF"]
SeqName=["T. rotula", "Skeletonema"]
# GSTD.STD_color_subplot(GSO_STD_all, 1, -10, 10, 1, "T. rotula", m="o", a=0.5)
# # GSTD.STD_color_subplot(MMETSP_STD_all, 1, -10,10, 2, "Skeletonema", m="D", a=0.5)
# plt.show()
# 

#GSTD.plot_by_day(STD_Combined, 1, -10, 10, SeqName=None, shape=shape_combo, trans=trans_combo, colorUp=colorUp, colorDown=colorDown) 
#GSTD.plot_by_day2(STD_Combined, 1, -10, 10, saveFig=0, SeqName=None, shape=shape_combo, trans=trans_combo, color=colorDown) 

GSTD.plot_by_day_Combined(STD_Combined, 1, -10, 10, saveFig=0, SeqName=None, shape=shape_combo, trans=trans_combo, color=colorDown) 
plt.show()