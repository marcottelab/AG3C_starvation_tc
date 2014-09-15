# -*- coding: utf-8 -*-


# <codecell>

##author: John Houser
##RUN AG3C analysis pipeline

import numpy as np
import sys
import os

sys.path.append('./data_processing')
sys.path.append('./fit_time_course')
sys.path.append('./GO_analysis')


import format_data_deseq
import format_data
import format_data_wapex
import format_data_norm_to_length
import RNA_data_format 
#import cluster_kmeans
import time_course_fit
import classify_fits
import id_convert
import ChartReport
import FuncClust
import go_rise
import mean_of_genes_in_go_term_mRNA as mofg
#import DAVID



#the folder to store the results in
prefix='/media/HD1_/Documents/AG3C/Results/RNA_dseq/'
datatype='Test_pr'
normtype='apex'

#Run the data analysis

#format proteomics data
print 'formating the data by '+normtype

if normtype=='deseq':
	#need to run the "R" script that does the DEseq normalization
	#	os.system("~/R-3.0.3/bin/Rscript /media/HD1_/Documents/AG3C/data_processing/Rdeseq.r '"+prefix+"' 'rna_data' '.csv'")
	#os.wait()
	format_data_deseq.run(prefix,datatype)
elif normtype=='read_depth':
	format_data.run(prefix,datatype)
elif normtype=='norm_to_length': #normalize each gene by the length of the transcript
	format_data_norm_to_length.run(prefix,datatype)
elif normtype=='apex':
	format_data_wapex.run(prefix,datatype)

'''
#fit the time course (this takes ~30min)
print 'fitting the data now...this could take 30min (or more)'
time_course_fit.run(prefix,datatype)

#sort the time courses
print 'classify the fits'
classify_fits.run(prefix,datatype)


#find go enrichments
#first convert names into suitable format
#print 'converting the names into the Entrez format'

classes=['up-regulated','down-regulated','temp_up-regulated','temp_down-regulated','other','combined']

for c in classes:
	#print c
	id_convert.run(prefix,datatype,c)

#find go terms in individual groups
print 'finding the go terms in and individual group'

for i in classes:
	print prefix+datatype+'_'+i+'_conv.txt' 
	FuncClust.DAVIDenrich(listF=prefix+datatype+'_'+str(i)+'_conv.txt',idType = 'ENTREZ_GENE_ID',bgName='Escherichia coli',category = 'GOTERM_BP_FAT')

ChartReport.DAVIDenrich(listF=prefix+datatype+'_combined_conv.txt',idType = 'ENTREZ_GENE_ID',bgName='Escherichia coli',category = 'GOTERM_BP_FAT')


#sort the go terms based upon their rise time
print 'sort the go terms based upon their rise time'
go_rise.run(prefix,datatype)

#find the average time course of the responders inside a particular go term
#print 'find the average time course of the responders inside a particular go term'
#for c in classes:
#	mofg.run(True,prefix,datatype,c)
'''

