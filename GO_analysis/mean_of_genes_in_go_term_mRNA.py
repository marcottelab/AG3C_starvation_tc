# -*- coding: utf-8 -*-# <nbformat>3.0</nbformat>

# <codecell>

from readcsv import read_csv
import numpy as np
import pickle
import pylab as plt
import sys

def run(plots,prefix,datatype,group):

	    

	filename=prefix+datatype+'_'+group+'_conv.txt.clustReport'
	delim='\t'
	output=read_csv(filename,delim)

	mapping=pickle.load(open(prefix+datatype+'_name_map.p','rb'))

	mapping_r=dict()
	for m in mapping:
		mapping_r[mapping[m]]=m
	#print mapping_r

	#put genes in a given annotation cluster into a matrix (rows are different clusters, columns are different proteins)
	genes=[]
	annotations=[]
	fdr=[]
	#print output[1]
	for i in range(len(output)):
	    #print output[i]
	    if len(output[i])==0 or i==0: #we are in a new annotation cluster
		if i!=0:
		    genes.append(list(set(genes_temp))) #only keep unique set of genes in an annotation cluster
		    annotations.append(annotation_temp)
		    fdr.append(fdr_temp)

		genes_temp=[]
		annotation_temp=[]
		fdr_temp=[]

	    #check if we are on a row with annotation output
	    else:
		    
		if 'GO' in output[i][0]:
		    genes_temp=genes_temp+output[i][5].split(', ')
		    annotation_temp.append(output[i][1])
		    fdr_temp.append(output[i][12])
	#print genes	    
	#get the preped data and average the protein time courses for that go-term
	(avout,stdout,av_reference)=pickle.load(open(prefix+datatype+'_data_format.p','rb'))
	av_ofav_ingo=[]
	for i in range(len(genes)):
	    a=list(set(map(str.lower,genes[i])))
	    temp_avs=[]
	    for j in range(len(a)):
		try:
			temp_avs.append(avout[list(map(str.lower,av_reference)).index(mapping_r[a[j]])])			
	    	except:
			None

	    av_ofav_ingo.append(np.average(temp_avs,axis=0))
	#print temp_avs
	#plot results
	
	#if plots==True:
	t=[3,4,5,6,8,12,24,24*7,14*24]

	av2=np.ndarray.transpose(np.array(av_ofav_ingo))
	#print len(av2),av2,len(t)
	if plots==True:	
		plt.figure()
		plt.semilogx(t,av2)

		#plt.legend(annotations,bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
		plt.xlabel('time(hr)')
		plt.ylabel('average relative expression')
		plt.savefig(prefix+group+'_m_of_go.pdf')
	
		plt.show()
	
	return av_ofav_ingo,annotations,fdr

if __name__=="__main__":
	print 'running...'
	plots=sys.argv[1]
	prefix=sys.argv[2]
	datatype=sys.argv[3]
	group=sys.argv[4]

	run(plots,prefix,datatype,group)

