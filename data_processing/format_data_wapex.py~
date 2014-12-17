import numpy as np
import pandas as pd
import csv
import pickle

def run(prefix,datatype):
	if 'Test' in datatype:
		T=[3,4,5,6]
	else:
		T=[3,4,5,6,8,24,48,24*7,336]

	count_cutoff=10
	
	#load the raw counts into pandas
	counts=pd.read_csv(prefix+datatype+'_data.csv',sep='\t',na_filter=True)
	#counts=counts.set_index('gene_id')

	#treat all zeros as missing (will be automatically handled by pandas)
	counts=counts[counts>0]
	sig_counts=[]#pd.DataFrame()
	#sig_counts['gene_id']=counts['gene_id']
	#filter out proteins with low total counts	
	for i in counts['gene_id']:
		#print counts[counts['gene_id']==i].filter(regex='sample').sum(axis=1)

		row_counts=counts[counts['gene_id']==i].filter(regex='sample').sum(axis=1)
		
		if row_counts.values[0]>=count_cutoff:
			sig_counts.append(i)
	
	
	counts=counts.set_index('gene_id')
	counts=counts[counts.index.isin(sig_counts)]
	counts=counts.reset_index()

	#normalize the protein counts using APEX
	pkey,Oi_values=pickle.load(open('./Results/Oi_name_map.p'))
	namedict=pickle.load(open('./Results/accession_to_common.p'))

	
	OI=dict()

	for i in range(len(Oi_values)-1):
	    if ''!=Oi_values[i]:
		OI[pkey[i].lower()]=float(Oi_values[i+1])

	#print counts
	#Normalize based on APEX method
	#print OI
	for c in counts:
		norm=0
		
		for i in list(set(counts.index)):
			#print counts.ix[i][0]
                        if i in namedict:
				if namedict[i].lower() in OI:
					#not taking into account p as p=1 or close enough.	
					norm+=float(counts.ix[i][c]/(OI[namedict[i].lower()]+1e-50))
					
		for i in list(set(counts.index)):
			if i in namedict:
				if namedict[i].lower() in OI:
					
					counts.ix[i][c]=counts.ix[i][c]/(norm*OI[namedict[i].lower()])

		#N=OI[counts['gene_id'][c].lower()]
		#counts[c]=counts[c]/N
		
	
	#counts=counts.set_index('gene_id')
	#average the counts
	av_counts=pd.DataFrame()
	st_counts=pd.DataFrame()
	av_counts['gene_id']=counts['gene_id']
	st_counts['gene_id']=counts['gene_id']
	#av_counts=av_counts.set_index('gene_id')
	#st_counts=st_counts.set_index('gene_id')
	#counts.reset_index()
	#print av_counts
	#print counts
	for t in T:
		#av_counts['t'+str(t)]=['']*len(av_counts)
		#st_counts['t'+str(t)]=['']*len(av_counts)
		#for i in counts.filter(regex='t'+str(t)+'$').mean(axis=1).index:
		av_counts['t'+str(t)]=counts.filter(regex='t'+str(t)+'$').mean(axis=1)
		st_counts['t'+str(t)]=counts.filter(regex='t'+str(t)+'$').std(axis=1)
		if t==336:
			av_counts['t'+str(t)]=counts.filter(regex=str(36)+'$')

	sig_change=[]
	for i in counts.index:
		maximum=counts[counts.index==i].filter(regex='sample').max(axis=1)
		minimum=counts[counts.index==i].filter(regex='sample').min(axis=1)
		if np.array(maximum)/np.array(minimum)>=1.5:
			sig_change.append(i)
	
	av_counts=av_counts[av_counts.index.isin(sig_change)]
	st_counts=st_counts[st_counts.index.isin(sig_change)]
	
	av_counts=av_counts.set_index('gene_id')
	st_counts=st_counts.set_index('gene_id')

	avout=[]
	ref_out=[]
	stdout=[]	
	for i in av_counts.index:
		
		avout.append(av_counts.ix[i])
		ref_out.append(i)
		stdout.append(st_counts.ix[i])	
	#avout=np.nan_to_num(avout)
	#stdout=np.nan_to_num(stdout)
	
	#print avout
	#print ref_out
	pickle.dump((avout,stdout,ref_out),open(prefix+datatype+'_data_format_wapex.p','wb'))
	

if __name__=='__main__':
	run('/media/HD1_/Documents/AG3C/Results/RNA_dseq/','Test')


