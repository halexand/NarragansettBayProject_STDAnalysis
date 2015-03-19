#!/usr/bin/env python

'''
Created on November 9, 2013

Analysis of NB_140220

Goal: 1) Create Pickled data containing the following: 

Dictionary: 
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
import STD_Functions as STD
import General_STD_analysis_control as GSTD

PickleOutDir="/Users/harrietalexander/Dropbox/NB_141020/PickledData/"
Throt_Count='/Users/harrietalexander/Dropbox/NB_141020/SD-NB-Throt.HTSeqCounts.tab'
Throt_Pk_STD=['Throt_PupP_STD', 'Throt_PupN_STD', 'Throt_PdnP_STD', 'Throt_PdnN_STD', 'Throt_NupP_STD', 'Throt_NupN_STD', 'Throt_NdnP_STD', 'Throt_NdnN_STD', 'Throt_All_P', 'Throt_All_N']
Throt_Pk_Count=['Throt_raw', 'Throt_tpm', 'Throt_SGNC', 'Throt_GeneList']
PP_Dir="/Users/harrietalexander/Dropbox/NB_141020/ASC_Throt/PostProb/"
Tail_125="*_125.txt"
SeqName="Throt"

GSTD.General_STD_analysis_Pickle(Throt_Count, Throt_Pk_Count, Throt_Pk_STD, PickleOutDir, 10, PP_Dir, Tail_125, SeqName, [0,1,2,3,4,9])

PickleOutDir="/Users/harrietalexander/Dropbox/NB_141020/PickledData/"
Skcos_Count='/Users/harrietalexander/Dropbox/NB_141020/SD-NB-Skcos.HTSeqCounts.tab'
Skcos_Pk_STD=['Skcos_PupP_STD', 'Skcos_PupN_STD', 'Skcos_PdnP_STD', 'Skcos_PdnN_STD', 'Skcos_NupP_STD', 'Skcos_NupN_STD', 'Skcos_NdnP_STD', 'Skcos_NdnN_STD', 'Skcos_All_P', 'Skcos_All_N']
Skcos_Pk_Count=['Skcos_raw', 'Skcos_tpm', 'Skcos_SGNC', 'Skcos_GeneList']
PP_Dir="/Users/harrietalexander/Dropbox/NB_141020/ASC_Skcos/PostProb/"
Tail_125="*_125.txt"
SeqName="Skcos"

GSTD.General_STD_analysis_Pickle(Skcos_Count, Skcos_Pk_Count, Skcos_Pk_STD, PickleOutDir, 10, PP_Dir, Tail_125, SeqName, [0,1,2,3,4,9])

