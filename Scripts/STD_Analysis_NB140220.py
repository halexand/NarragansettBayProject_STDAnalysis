#!/usr/bin/env python

'''
Created on Feb 20, 2014

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

PickleOutDir="/Users/harrietalexander/Dropbox/DyhrmanLab_MASON/Scripts/Harriet/Python_Scripts/STD_Plots/Pickle/"

#GSO Raw Data
GSO_CountFile="/Users/harrietalexander/Dropbox/NB_141020/SD-NB-Throt.HTSeqCounts.tab"
GSO_PP_Dir="/Users/harrietalexander/Dropbox/NB_141020/ASC_Throt/PostProb/"
GSO_Tail_125="*_125.txt"




GSO_STD_all=GSTD.General_STD_analysis(GSO_CountFile, 10, GSO_PP_Dir, GSO_Tail_125, "T. rotula", RR=[0,1,2,3,4,9])

MMETSP_CountFile="/Users/harrietalexander/Dropbox/NB_130920/SD-NB_MMETSP0013_count.tab"
MMETSP_PP_Dir="/Users/harrietalexander/Dropbox/NB_130920/MMETSP_ASC/PostProb/"
MMETSP_Tail_125="*_125.txt"

MMETSP_STD_all=GSTD.General_STD_analysis(MMETSP_CountFile, 10, MMETSP_PP_Dir, MMETSP_Tail_125, "S. costatum", RR=[0,1,2,3,4,9])


STD_Combined=[GSO_STD_all, MMETSP_STD_all]
shape_combo=["o", "D"]
trans_combo=[0.5, 0.5]
colorDown=["#FF9900","#66CCFF"]
colorUp=["#CC6600","#3333FF"]
SeqName=["T. rotula", "Skeletonema"]
#Colored plots with different colored dots for each day
#GSTD.STD_color_subplot(GSO_STD_all, 1, -10, 10, 1, "T. rotula", m="o", a=0.5)
# GSTD.STD_color_subplot(MMETSP_STD_all, 1, -10,10, 2, "Skeletonema", m="D", a=0.5)
GSTD.plot_by_day(STD_Combined, 1, -10, 10, SeqName=None, shape=shape_combo, trans=trans_combo, colorUp=colorUp, colorDown=colorDown) 


plt.show()

#Plots with colored points for each different 
#GSTD.plot_by_day(STD_Combined, 1, -10, 10, SeqName=None, shape=shape_combo, trans=trans_combo, colorUp=colorUp, colorDown=colorDown) 
#GSTD.plot_by_day2(STD_Combined, 1, -10, 10, saveFig=0, SeqName=None, shape=shape_combo, trans=trans_combo, color=colorDown) 
