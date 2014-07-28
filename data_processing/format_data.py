# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

##author: John Houser
##format, average, and normalize AG3C data

import numpy as np
import pandas as pd
import csv
import pickle

def run(prefix,datatype):
	T=[3,4,5,6,8,24,48,24*7,24*14]

	#load the raw counts into pandas
	counts=pd.read_csv(prefix+datatype+'_data.csv',sep='\t',na_filter=True)
	#treat all zeros as missing (will be automatically handled by pandas)
	counts=counts[counts>0]	
	
	sig_change=[]
	counts=counts.set_index('gene_id')
	counts
	#filter out proteins with low counts
	for i in counts.index.values:
		if type(i)!=type(float('nan')): 

			if np.sum(counts.loc[i,:])>=10:
				sig_change.append(i)	
	
	counts=counts[counts.index.isin(sig_change)]
	
	counts=counts.reset_index()

	#normalize the counts
	for i in list(counts.columns.values):
		#print counts[i].div(2.0)
		if i!='gene_id':
			
			N=counts[i].sum(axis=1)
		
			counts[i]=counts[i].div(float(N))

	#average the counts
	av_counts=pd.DataFrame()
	st_counts=pd.DataFrame()
	av_counts['gene_id']=counts['gene_id']
	st_counts['gene_id']=counts['gene_id']

	for t in T:
		av_counts['t'+str(t)]=counts.filter(regex='t'+str(t)).mean(axis=1)
		st_counts['t'+str(t)]=counts.filter(regex='t'+str(t)).std(axis=1)
	av_counts=av_counts.set_index('gene_id')
	st_counts=st_counts.set_index('gene_id')

	#print av_counts.index.values
	#normalize to the max across the time course
	for i in av_counts.index.values:

		if type(i)!=type(float('nan')):
			N=av_counts.ix[i][av_counts.columns].max(axis=1)
		
			av_counts.ix[i]=av_counts.ix[i].div(N)
			st_counts.ix[i]=st_counts.ix[i].div(N)
	
	#filter out those proteins that are changing by less than 50%
	print len(av_counts)
	sig_change=[]
	for i in av_counts.index.values:
		if type(i)!=type(float('nan')):
			N=av_counts.ix[i][av_counts.columns].max(axis=1)/av_counts.ix[i][av_counts.columns].min(axis=1)
			if N>=1.5:
				sig_change.append(i)
			
	av_counts=av_counts[av_counts.index.isin(sig_change)]
	st_counts=st_counts[st_counts.index.isin(sig_change)]

	avout=[]
	ref_out=[]
	stdout=[]	
	for i in av_counts.index.values:
		avout.append(av_counts.ix[i])
		ref_out.append(i)
		stdout.append(st_counts.ix[i])	
	avout=np.nan_to_num(avout)
	stdout=np.nan_to_num(stdout)
	
	
	print len(avout),len(avout[0])
	#pickle.dump((avout,stdout,ref_out),open(prefix+datatype+'_data_format.p','wb'))
	

if __name__=='__main__':
	run('/media/HD1_/Documents/AG3C/Results/','prot')
