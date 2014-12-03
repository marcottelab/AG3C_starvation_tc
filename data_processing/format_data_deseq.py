import numpy as np
import pandas as pd
import csv
import pickle

def run(prefix,datatype):
	if datatype=='Test':
		T=[3,4,5,6]
	else:
		T=[3,4,5,6,8,24,48,24*7,24*14]

	#load the raw counts into pandas
	counts=pd.read_csv(prefix+datatype+'_data.csv',sep='\t',na_filter=True)
	#treat all zeros as missing (will be automatically handled by pandas)
	counts=counts[counts>0]	
	
	
	#load the size factors normalize the size factors so that the average of the t3 (first) time point is=1
	for t in T:
		if t>3: 
			sf_temp=pd.read_csv(prefix+'size_factors/'+datatype+'_size_factors3_t'+str(t),',')
			sf_temp.columns=['ID','size_factors']	
			
			Norm=sf_temp[sf_temp['ID'].str.contains('t3$')].mean()
			#print Norm
			sf_temp['size_factors']=sf_temp['size_factors'].divide(Norm['size_factors'])
			#print sf_temp[sf_temp['ID'].str.contains('t3')]
			
			sf_temp=sf_temp[~sf_temp['ID'].str.contains('t3$')]
			
			if t==4:
				sf=sf_temp	
			else:
				sf=sf.append(sf_temp)

	#normalize the counts
	for i in sf['ID']:
		#print counts[i].div(2.0)	
		N=sf['size_factors'][sf['ID'].str.contains(i)]
		
		counts[i]=counts[i].div(float(N))
	#average the counts
	av_counts=pd.DataFrame()
	st_counts=pd.DataFrame()
	av_counts['gene_id']=counts['gene_id']
	st_counts['gene_id']=counts['gene_id']

	for t in T:
		av_counts['t'+str(t)]=counts.filter(regex='t'+str(t)+'$').mean(axis=1)
		st_counts['t'+str(t)]=counts.filter(regex='t'+str(t)+'$').std(axis=1)
	av_counts=av_counts.set_index('gene_id')
	st_counts=st_counts.set_index('gene_id')
	
	#normalize to the max across the time course
	for i in av_counts.index.values:
		N=av_counts.ix[i][av_counts.columns].max(axis=1)
		
		av_counts.ix[i]=av_counts.ix[i].div(N)
		st_counts.ix[i]=st_counts.ix[i].div(N)

	files=[]
	for i in range(len(T)-1):
		files.append(prefix+'size_factors/'+datatype+'_sig_gene_3_t'+str(T[i+1]))

	#filter out datapoints that aren't changing significantly
	sig_change=load_p_values(files,',',prefix+datatype+'_data.csv')
	
	av_counts=av_counts[av_counts.index.isin(sig_change)]
	st_counts=st_counts[st_counts.index.isin(sig_change)]
	
	
	avout=[]
	ref_out=[]
	stdout=[]	
	for i in av_counts.index.values:
		if 'ECB_t' not in i and 'ECB_r' not in i and 'RF' not in i:
			avout.append(av_counts.ix[i])
			ref_out.append(i)
			stdout.append(st_counts.ix[i])	
	avout=np.nan_to_num(avout)
	stdout=np.nan_to_num(stdout)
	
	#print avout
	#print ref_out
	pickle.dump((avout,stdout,ref_out),open(prefix+datatype+'_results/'+datatype+'_data_format.p','wb'))
	

def load_p_values(files,delim,ref_file):
        #import data from several different csv files. output to a list
        
	ref_dict=dict()
	with open(ref_file,'rb') as csvfile:
		read=csv.reader(csvfile,delimiter='\t',quotechar='"')
		q=1
		for row in read:
			if len(row)>=2:
				ref_dict[str(q)]=row[0]
				q+=1
	output=[]
	J=0
	for f in files:
		
		with open(f, 'rb') as csvfile:
		    read = csv.reader(csvfile, delimiter=delim, quotechar='"')
		    first=1
		    J+=1
		    for row in read:
			if first==1:
				row0=row
		
			#if row[0] not in output:
			else:		
				if row[1] in ref_dict:
					output.append(ref_dict[row[1]])
					
			if first==1:
				first=0
		   
        
	output=list(set(output))

        return output


if __name__=='__main__':
	run('/media/HD1_/Documents/AG3C/Results/RNA_dseq/','Test')


