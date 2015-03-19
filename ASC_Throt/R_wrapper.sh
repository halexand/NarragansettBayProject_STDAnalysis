#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=16:dc,vmem=40gb,walltime=100:00:00
#PBS -M alexandh@indiana.edu
#PBS -m abe
#PBS -N NB_Throt_ASC
#PBS -j oe
#PBS -q batch




module load R/2.14.0 

wd=/N/dc2/projects/ldeo/140120/ASC_Throt
cd $wd
/N/soft/mason/R/2.14.0/bin/Rscript RunASC_Throt_140120.R
